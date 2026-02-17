

get_target_list <- function(filepath) {
	out <- 
		readxl::read_excel(filepath) |>
		dplyr::filter(Status != "3-Reserved (AMP)") |>
		tidyr::separate_rows(ENSGID, sep = ";")
	return(out)

}


prepare_ensembl_ids <- function(x) {

	# Convert from gene.symbol to ensembl.gene
	library(EnsDb.Hsapiens.v86)
	output <- 
		ensembldb::select(
			EnsDb.Hsapiens.v86, 
			keys = x, 
			keytype = "SYMBOL", 
			columns = columns(EnsDb.Hsapiens.v86)
		) |>
		dplyr::filter(grepl("ENSG", GENEID)) |>
		dplyr::filter(TXBIOTYPE == "protein_coding") |>
		dplyr::select(GENEID, SYMBOL, UNIPROTID) |>
		dplyr::distinct()

	# Test outupt
	if(all(x %in% output$SYMBOL)) {
		return(output)
	}else{
		logger::log_warn("Missing ensembl ID.")
		return(output)
	}
}	

merge_targets <- function(x, y) {

	# tar_load(current_targets); x = current_targets; tar_load(references_ensgid); y = references_ensgid
	out <- 
		dplyr::full_join(x, y, by = c("SYMBOL", "ENSGID" = "GENEID")) |>
		dplyr::select(SYMBOL, ENSGID) |>
		dplyr::distinct(.keep_all = TRUE)
	return(out)

}

get_hpa_filter <- function(x, y, p) {

	# tar_load(normal_ihc_data); x = normal_ihc_data; tar_load(param_ls); p = param_ls; tar_load(current_targets); y = current_targets
	out <- 
		x |>
		dplyr::filter(grepl(paste0(p$core_tissue, collapse = "|"), Tissue)) |>
		dplyr::filter(Gene %in% y$ENSGID)
	return(out)
}

hpa_download <- function(url_str){
	#' Download and load Human Protein Atlas data.
	#'
	#' This function downloads a data file from the Human Protein Atlas (HPA)
	#' and loads it into an R data frame. For details, see:
	#' https://www.proteinatlas.org/about/download
	#'
	#' @param url_str (string) The name of the file to download from the HPA.
	#'   Defaults to "rna_cancer_sample.tsv.gz".
	#' @param target_ids (data.frame) Output of get_ensembl_id()
	#'
	#' @return A data frame containing the loaded data, or NULL if the download
	#'   or loading fails.
	#'
	#' @examples
	#' \dontrun{
	#'   hpa_data <- hpa_download("rna_cancer_sample.tsv.gz")
	#'   if (!is.null(hpa_data)) {
	#'     print(head(hpa_data))
	#'   }
	#' }
	#'
	#' @export

	# Check input arguments
	checkmate::assertCharacter(
		url_str, 
		len = 1, 
		any.missing = FALSE, 
		null.ok = FALSE
	)	
		
	# Define paths
	url_id <- 
		stringr::str_glue("https://www.proteinatlas.org/download/tsv/{url_str}")
	dest_id <- 
		stringr::str_glue("../../input/hpa/{url_str}")
	
	# Download file if it doesn't exist
	if(file.exists(dest_id)) {
		logger::log_info("Download exists:", dest_id)
	}else{
		tryCatch({
		  download.file(url_id, destfile = dest_id)
		  logger::log_info("Downloaded data to:", dest_id)
		}, error = function(e) {
		  logger::log_error("Failed to download data:", e$message)
		  return(NULL) # Return NULL if download fails
		})		
	}
	
	# Load files
	# TODO filter to genes of interest on the fly to save memory 
	# TODO Save to rds and remove zip file, and check for rds
	tryCatch({
		input_hpa <- readr::read_tsv(dest_id)
	}, error = function(e) {
		logger::log_error("Failed to read data:", e$message)
		return(NULL) # Return NULL if reading fails
	})
	
	# Test output
	if (!exists("input_hpa") || is.null(input_hpa)) {
		logger::log_warn("No data to return.")
		return(NULL)
	}else{
		return(input_hpa)
	}


}


get_paxdb <- function(x, p) {

	# Homo sapiens data downloaded from https://pax-db.org/download

	# tar_load(current_targets); x = current_targets

	fof <- 
		list.files(
			path = "../../input/paxdb/9606",
			pattern = ".txt",
			full.names=TRUE
		)
	fof <- 	
		grep(
			pattern = tolower(paste0(p$core_tissue, collapse = "|")), 
			x = fof, 
			ignore.case = TRUE,
			value = TRUE
		)	
		
	# Translate these to uniprotids
	library(EnsDb.Hsapiens.v86)
	output <- 
		ensembldb::select(
			EnsDb.Hsapiens.v86, 
			keys = c("EGFR", "ERBB2", unique(x$SYMBOL)), 
			keytype = "SYMBOL", 
			columns = c("SYMBOL", "GENEID", "TXBIOTYPE", "PROTEINID")
		) |>
		dplyr::filter(grepl("ENSG", GENEID)) |>
		dplyr::filter(TXBIOTYPE == "protein_coding") |>
		dplyr::distinct()
		
	# Perform lookup	
	lambda <- function(f) {
		tmp <- 
			readr::read_tsv(f, comment = "#", col_names = FALSE, show_col_types = FALSE) 
		if (ncol(tmp) > 1) {
			tmp <- 
				tmp |>
				dplyr::rename_with(~"ID_COL", dplyr::contains("X2")) |>
				dplyr::rename_with(~"SYMBOL", dplyr::contains("X1")) |>
				dplyr::rename_with(~"abundance", dplyr::contains("X3")) |>
				dplyr::mutate(PROTEINID = gsub("9606.", "", ID_COL)) |>
				dplyr::filter(PROTEINID %in% output$PROTEINID) |>
				dplyr::select(SYMBOL, PROTEINID, abundance)
			return(tmp)		
		}else{
			return(NULL)
		}	
	}
	inputs <- lapply(X = fof, FUN = lambda)
	names(inputs) <- basename(fof)


	lambda <- function(f) {
		tmp <- grep("^\\s*#", readLines(f), value=TRUE)		
		return(grep("name: H.sapiens", tmp, value = TRUE))
	}
	metadata <- lapply(X = fof, FUN = lambda)	
	metadata <- gsub("#name: H.sapiens - ", "", unlist(metadata))
	metadata <- 
		data.frame(
			src = basename(fof),
			metadata = metadata
		) |>
		tibble::as_tibble()
	
	out <- 
		inputs |>
		dplyr::bind_rows(.id = "src") |>
		dplyr::left_join(metadata, by = "src") |>
		dplyr::relocate(src, .after = "metadata")
		
	return(out)
}



get_xls <- function(filepath, required_col = "Target gene", ...) {
	#' Read All Sheets from an Excel Workbook
	#'
	#' This function takes the file path of an Excel workbook and reads all sheets
	#' into a named list of tibbles (data frames).
	#'
	#' @param filepath A string specifying the path to the Excel file. The
	#'   function will stop if the file does not exist.
	#' @param required_col A string specifying a column name that must be present
	#'   in every sheet. If `NULL`, no column validation is performed.
	#'   Defaults to `"Target gene"`.
	#' @param ... Additional arguments to be passed to `readxl::read_excel()`.
	#'   This allows for customization of the reading process (e.g., `col_names`, `skip`).
	#'
	#' @return A tibble of HUGO gene symbols 
	#'
	#' @examples
	#' \dontrun{
	#'   # Create a dummy excel file for testing
	#'   library(writexl)
	#'   dummy_data1 <- data.frame(ID = 1:5, Value = rnorm(5))
	#'   dummy_data2 <- data.frame(Gene = c("BRCA1", "TP53"), FoldChange = c(1.5, -2.1))
	#'   sheets_to_write <- list("DataSet1" = dummy_data1, "GeneExpression" = dummy_data2)
	#'   write_xlsx(sheets_to_write, "dummy_workbook.xlsx")
	#'
	#'   # Use the function to read the workbook
	#'   all_sheets <- read_all_excel_sheets("dummy_workbook.xlsx")
	#'
	#'   # Access a specific sheet
	#'   print(all_sheets$DataSet1)
	#' }

	# --- Input Validation ---
	# 1. Check if the filepath is a single, non-missing string.
	# 2. Check if the file actually exists at the given path.
	assertions::assert_string(filepath)
	assertions::assert_file_exists(filepath)

	# --- Reading the Workbook ---
	sheet_names <- readxl::excel_sheets(filepath)
	all_targets <- 
		lapply(X = sheet_names, FUN = function(x) readxl::read_excel(filepath, sheet = x)) 
		
	# --- Column Validation ---
	if (!is.null(required_col)) {
		test <- 
			lapply(
				X = all_targets, 
				FUN = function(x) required_col %in% colnames(x)
			)
		if( !all(unlist(test)) ) {
        logger::log_error(
          glue::glue(
            "Validation Error: The required column {required_col} was not found in all sheets.") 
          )
		}
	}
		
	all_targets <- 
		lapply(X = all_targets, FUN = function(x){dplyr::select(x, all_of("Target gene"))}) |>
		dplyr::bind_rows() |>
		dplyr::select(`Target gene`) |>
		dplyr::distinct()

	# Return the HUGO gene symbols
	logger::log_info(glue::glue("Loaded {nrow(all_targets)} targets"))
	return(all_targets)

}


get_hpa_quant <- function(x, y, p) {

	# tar_load(hpa_quant_f); x = hpa_quant_f; tar_load(param_ls); p = param_ls; tar_load(current_targets); y = current_targets
	target_set <- 
		y |>
		tidyr::separate_rows(ENSGID, sep = ";") |>
		dplyr::distinct(ENSGID) |>
		dplyr::pull(ENSGID)

	out <-
		readr::read_tsv(x, comment = "#") |>
		dplyr::rename("ENSGID" = "Gene ID") |>
		dplyr::rename("SYMBOL" = "Gene Name") |>
		dplyr::select(SYMBOL, ENSGID, contains(p$core_tissue)) |>
		tidyr::pivot_longer(contains(p$core_tissue), names_to = "tissue", values_to = "abundance") |>
		dplyr::filter(ENSGID %in% target_set)
	return(out)
}

get_gtex_quant <- function(x, y, p) {

	# tar_load(gtex_quant_f); x = gtex_quant_f; tar_load(all_targets); y = all_targets; tar_load(param_ls); p = param_ls

	library(dplyr)

	out <- 
		readxl::read_excel(x, sheet = "F protein normalized abundance", skip = 2)
	meta <-
		readxl::read_excel(x, sheet = "A experimental design")

	# Reorder/annotate data
	annot <- 
		meta |>
		dplyr::select(`GTEx ID`, `Sample`) |>
		dplyr::distinct(.keep_all = TRUE)
	out1 <- 
		out |>
		dplyr::filter(gene.id %in% y$ENSGID) |>
		dplyr::select(!gene.id.full) |>
		dplyr::mutate_at(vars(contains("GTEX-")), ~ (as.numeric(.))) |>
		tidyr::pivot_longer(!gene.id, names_to = "GTEx ID", values_to = "abundance") |>
		dplyr::filter(!grepl("reference", `GTEx ID`)) |>
		tidyr::separate(col = `GTEx ID`, into = c("GTEx ID", "replicate_id"), sep = "\\.\\.\\.") |>
		dplyr::left_join(annot, by = "GTEx ID") |>
		dplyr::left_join(y, by = c("gene.id" = "ENSGID")) |>
		dplyr::rename("ENSGID" = "gene.id") |>
		dplyr::relocate("SYMBOL") |>
		dplyr::relocate("Sample", .after = "SYMBOL") |>
		dplyr::filter(grepl(paste0(p$core_tissue, collapse = "|"), Sample, ignore.case = TRUE))

	return(out1)

}

get_gtex_quant_fang <- function(y, p) {

	# tar_load(all_targets); y = all_targets; tar_load(param_ls); p = param_ls

	fof <- 
		list.files(
			path = "../../input/gtex/fang/",
			pattern = ".gz",
			full.names = TRUE
		)
	fof_ls <-
		lapply(X = fof, FUN = readr::read_tsv)
	names(fof_ls) <- 
		stringr::str_extract(basename(fof), "(?<=\\.)[A-Za-z]+(?=_)")
	out <-
		dplyr::bind_rows(fof_ls, .id = "tissue") |>
		dplyr::select(-chr, -start, -end) |>
		dplyr::mutate(ENSGID = tools::file_path_sans_ext(gene_id)) |>
		dplyr::filter(ENSGID %in% y$ENSGID) |>
		dplyr::filter(grepl(paste0(p$core_tissue, collapse = "|"), tissue, ignore.case = TRUE)) |>
		dplyr::select(ENSGID, tissue, contains("GTEX-")) |>
		tidyr::pivot_longer(!ENSGID:tissue, names_to = "GTEx ID", values_to = "abundance") |>
		dplyr::left_join(y, by = "ENSGID") |>
		dplyr::relocate(SYMBOL) |>
		dplyr::relocate(tissue, .after = "SYMBOL")

	return(tb = out)

}







