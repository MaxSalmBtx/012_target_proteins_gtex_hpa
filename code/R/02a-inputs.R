
## @knitr get_data
# This downlaods a subset of the HPA data

if(FALSE){
	logger::log_warn("DEBUGGING SCAFFOLDING STILL ON - will stop quarto render!!!")
	params <- list()
	params$target_id <- "CEACAM5"
}

negative_control <- c("GAPDH")
positive_control_adc <- c("ERBB2", "MET", "EGFR", "NECTIN4", "TACSTD2",  "CD47", "CD276", "FAP", "MSLN",  "DLL3", "CLDN18")
positive_control_radio <- c("FOLH1", "SSTR2")
portfolio <- c("CEACAM5", "MMP14", "EPHA2")
amplicon_targets <- c("HM13", "TM9SF4", "PROCR", "SLCO4A1", "SLC52A2", "SLC39A4", "NCSTN", "IGSF8", "F11R", "ADAM15", "MPZL1", "EFNA4", "LAMP1", "EPHB3", "LAMP3", "ABCC5", "SYPL1", "CLDN4", "SLC12A7")
ref_cols <- 
	c(
		"candidate" = BTX_COL$pink_out, 
		"positive_control_adc" = BTX_COL$blue_d_out, 
		"positive_control_radio" = BTX_COL$blue_l_out,
		"portfolio" = BTX_COL$green_out
	)

ensembl_id_ls <- list()
ensembl_id_ls$candidate <- get_ensembl_id(x = params$target_id) 
ensembl_id_ls$portfolio <- get_ensembl_id(x = portfolio) 
# ensembl_id_ls$negative_control <- get_ensembl_id(x = negative_control) # need a better one that does not change in tumours
ensembl_id_ls$positive_control_adc <- get_ensembl_id(x = positive_control_adc) 
ensembl_id_ls$positive_control_radio <- get_ensembl_id(x = positive_control_radio) 

ensembl_id <- 
	ensembl_id_ls |>
	dplyr::bind_rows(.id = "src") |>
	dplyr::distinct(.keep_all = TRUE) |>
	dplyr::filter(GENEBIOTYPE == "protein_coding")

# Identifier mappings
tmp <- 
	ensembl_id |>
		dplyr::filter(SYMBOL == params$target_id & src == "candidate") 
target_ensembl_ids <- 
	tmp |>
		dplyr::select(GENEID) |>
		dplyr::distinct()

ensembl_id <- 
	ensembl_id |>
		dplyr::filter(SYMBOL != params$target_id) |>
		dplyr::bind_rows(tmp) # In case candidate is in reference targets

ensembl_isoform_ids <-
	ensembl_id |>
	dplyr::filter(src == "candidate") |>
	dplyr::select(TXID) |>
	dplyr::distinct() |>
	unlist()

ensembl_prot_ids <-
	ensembl_id |>
	dplyr::filter(src == "candidate") |>
	dplyr::select(PROTEINID) |>
	dplyr::distinct() |>
	tidyr::drop_na() |>
	unlist() 
	

external_links <-
	list(
		genecards = glue::glue("https://www.genecards.org/cgi-bin/carddisp.pl?gene={params$target_id}"),
		ensembl = glue::glue("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g={params$target_id}"),
		ncbi = glue::glue("https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000142627{as.character(target_ensembl_ids)}"),
		expression_atlas = glue::glue("http://www.ebi.ac.uk/gxa/query?geneQuery={params$target_id}"),
		hpa = glue::glue("http://www.proteinatlas.org/{as.character(target_ensembl_ids)}-{params$target_id}/tissue")
	)
	
logger::log_info("Complete: get_data")
## ---- end

## @knitr get_uniprot
## UniProt data (all the data you could need, >250 entries)
up <- UniProt.ws(taxId=9606)
uniprot_data <- 
	UniProt.ws::select(
		up, 
		keys = params$target_id, 
		keytype = c("GeneCards", "Gene_Name")[1], # gives non-human outputs...
		columns = columns(up)
	)

protein_full_name <- gsub("\\(.*", "", uniprot_data[,"Protein.names"])

new_uniport_query <- 
	uniprot_data |>
		dplyr::select(any_of(c("Entry", "Reviewed", "Entry.Name", "Protein.names", "Gene.Names", "Transmembrane", "Signal.peptide", "Subcellular.location..CC.", "Domain..FT.", "Motif", "Protein.families", "Pathway")))
library(UniProtExtractR)
uniprot_extract <- 
	UniProtExtractR::uniprotextract(
		new_uniport_query, 
		map.up=NULL, 
		write.local=FALSE
	) # includes transmembrane count/location
logger::log_info("Complete: get_uniprot")
## ---- end

## @knitr get_opentargets

## OpenTargets data
open_targets_data <- get_open_targets_data(x = target_ensembl_ids)

europmc_plt <-
	open_targets_data$europePMCQuery_ls |>
		dplyr::group_by(`disease.name`) |>
		dplyr::summarise(n = sum(resourceScore)) |>
		dplyr::arrange(n) |>
		ggplot(aes(x = n, y = forcats::fct_reorder(`disease.name`, n))) +
		geom_col(fill = BTX_COL$blue_d_in, color = BTX_COL$blue_d_out) +
		theme(
			plot.background = element_rect(color = NA, fill = "grey97")
		) 
logger::log_info("Complete: get_opentargets")
## ---- end

## @knitr get_compartments
compartments_ls <- get_compartments(enspids = ensembl_prot_ids) 
compartments_k <-
	compartments_ls$`knowledge_full.tsv`
colnames(compartments_k) <- 
	c(
	"ensembl_id",
	"gene_name",
	"GO",
	"location",
	"source",
	"evidence",
	"confidence"
	)	
compartments_k <-
	compartments_k |>
	dplyr::select(-GO, - gene_name) |>
	dplyr::distinct(.keep_all = TRUE) |>
	dplyr::group_by(location) |>
	dplyr::arrange(desc(confidence))
logger::log_info("Complete: get_compartments")
## ---- end


## @knitr get_beacon

# Manual database dump (11/08/2025)
beacon_data_a <- readxl::read_excel("../../input/beacon/Treatments-11-Aug-2025-15-11-26-48793.xlsx", sheet = "Treatments")

# beacon_data |>
	# dplyr::select(`Drug Targets`, `Therapeutic Class`) |>
	# tidyr::separate_longer_delim(cols = c(`Drug Targets`), delim = ";")

top_30 <- 
	names(sort(table(unlist(strsplit(beacon_data_a[["Drug Targets"]], "; "))), decreasing = TRUE)[1:30])

## ---- end


## @knitr get_globaldata

# Manual retrieval from Wahl et al: https://pubmed.ncbi.nlm.nih.gov/34857619/
wahl_data <- 
	readxl::read_excel(
		"../../input/wahl_et_al.xlsx", 
		sheet = "Sheet1"
	)

read_globadata <- function(x) {
	tmp <- readxl::read_excel(x, col_names = FALSE, trim_ws = TRUE)
	idx <- which(tmp == "Highest Development Stage", arr.ind = TRUE)
	idx <- idx[, "row"]
	header <- tmp[idx, ]
	tmp <- tmp[(idx + 1):nrow(tmp), ]
	colnames(tmp) <- header
	out <- 
		tmp |>
		#dplyr::filter(grepl("Oncology", `Therapy Area`)) |>
		#dplyr::filter(grepl(";", `Target`)) |> # Filter to conjugates
		dplyr::select(any_of(c("Drug Id", "Target", "Drug Name", "Therapy Area", "Indication", "Molecule Type", "Development Stage", "Highest Development Stage")))
	return(out)
}

data_dir = "../../input/globaldata/"
fof <- list.files(path = data_dir, pattern = "\\.xlsx", all.files = TRUE, full.names = TRUE)
globaldata_tb <- 
	lapply(X = fof, FUN = read_globadata) |>
	dplyr::bind_rows()

idx <- 
	stringdist::amatch(protein_full_name, globaldata_tb$Target, maxDist = 10)
idx <- globaldata_tb$Target[idx] # fuzzy match to protein full name

globaldata_target_all_tb <- 
	globaldata_tb |>
		dplyr::filter(grepl("Oncology", `Therapy Area`)) |>
		dplyr::filter(grepl(idx, Target, ignore.case = TRUE))

# globaldata_target_tb <- read_globadata(x = "../../input/globaldata/08072025_DDR1.xlsx")

# Prepare for presentation
lvls <- 
	rev(c(
		"Marketed",
		"Phase III",
		"Phase II",
		"Phase I",
		"Preclinical",
		"Discovery",
		"Inactive"
	))

# globaldata_indications_tb <-
	# globaldata_target_all_tb |>
	# dplyr::group_by(`Drug Name`) |>
	# dplyr::summarise(Indication = gsub("\n", "; ", Indication)) |>
	# dplyr::summarise(indications = stringr::str_c(Indication, collapse = "; ")) |>
	# dplyr::mutate(Indications = deduplicate_multi_delimiter_strings(indications))	
# globaldata_stage_tb <-
	# globaldata_target_all_tb |>
		# dplyr::group_by(`Drug Name`) |>
		# dplyr::summarise(`Highest Development Stage` = unique(`Highest Development Stage`, collapse = "; ")) |>
		# dplyr::mutate(`Highest Development Stage` = factor(`Highest Development Stage`, levels = lvls))
# globaldata_modality_tb <-
	# globaldata_target_all_tb |>
		# dplyr::select(`Drug Name`, `Molecule Type`) 

# globaldata_gt <-
	# globaldata_indications_tb |>
	# dplyr::left_join(globaldata_stage_tb, by = "Drug Name") |>
	# dplyr::left_join(globaldata_modality_tb, by = "Drug Name") |>
	# dplyr::select(`Drug Name`, `Highest Development Stage`, `Indications`) |>
	# dplyr::distinct(.keep_all = TRUE) |>
	# dplyr::arrange(`Highest Development Stage`) 


globaldata_gt <-
  globaldata_target_all_tb |>
  dplyr::mutate(
    `Drug Name` = stringr::str_trim(`Drug Name`),
    `Highest Development Stage` = factor(`Highest Development Stage`, levels = lvls, ordered = TRUE)
  ) |>
  # Ensure the grouping key is not empty
  dplyr::filter(!is.na(`Drug Name`) & `Drug Name` != "") |>
  # Group by the clean drug name and create all summary columns in a single step.
  dplyr::group_by(`Drug Name`) |>
  dplyr::summarise(
    # A) Find the single HIGHEST development stage for each drug.
    # The max() function on a factor correctly uses the defined levels.
    `Highest Development Stage` = max(`Highest Development Stage`, na.rm = TRUE),
    # B) Consolidate and deduplicate indications into a single string.
    # This handles indications that are split by newlines or semicolons.
    Indications = {
      all_indications <- na.omit(Indication) |>
        stringr::str_split(pattern = ";\\s*|\\n") |> # Split by semicolon or newline
        unlist() |>
        stringr::str_trim()
      unique(all_indications) |> stringr::str_c(collapse = "; ")
    }
  ) |>
  dplyr::ungroup() |> # Best practice to ungroup after summarizing
  # 3. FINALIZE
  # Select the desired columns and arrange by the development stage.
  dplyr::select(`Drug Name`, `Highest Development Stage`, `Indications`) |>
  dplyr::arrange(`Highest Development Stage`)

logger::log_info("Complete: get_globaldata")
## ---- end


## @knitr get_globocan
## WHO Prevalence data
globocan_data <- get_globocan_data()
## ---- end




##############################################################################



## @knitr get_hpa


if (FALSE) {
	hpa_metadata_ls <- lapply(FUN = get_hpa_data, X = ensembl_id$GENEID)
	names(hpa_metadata_ls) <- ensembl_id$GENEID
	hpa_metadata_ls <-
		Filter(Negate(is.null), hpa_metadata_ls)	
}else{
	## Filter to target
	tmp_ids <- 
		unique(ensembl_id$GENEID[ensembl_id$SYMBOL %in% params$target_id])
	hpa_metadata_ls <- 
		lapply(
			FUN = get_hpa_data, 
			X = tmp_ids
		)	
	names(hpa_metadata_ls) <- tmp_ids
	hpa_metadata_ls <-
		Filter(Negate(is.null), hpa_metadata_ls)
}


coi <- c("Gene","Gene synonym","Ensembl","Gene description","Uniprot","Chromosome","Position","Disease involvement","RNA tissue specificity","RNA tissue distribution","RNA tissue specificity score","RNA tissue specific nTPM","RNA single cell type specificity","RNA single cell type distribution","RNA single cell type specificity score","RNA single cell type specific nTPM","RNA cancer specificity","RNA cancer distribution","RNA cancer specificity score","RNA cancer specific FPKM","RNA blood cell specificity","RNA blood cell distribution","RNA blood cell specificity score","RNA blood cell specific nTPM","RNA blood lineage specificity","RNA blood lineage distribution","RNA blood lineage specificity score","RNA blood lineage specific nTPM","RNA cell line specificity","RNA cell line distribution","RNA cell line specificity score","RNA cell line specific nTPM","RNA tissue cell type enrichment","Antibody","Reliability (IH)","Reliability (IF)","Subcellular location","Secretome location","Secretome function","CCD Protein","CCD Transcript","Blood concentration - Conc. blood IM [pg/L]","Blood concentration - Conc. blood MS [pg/L]","Blood expression cluster","Tissue expression cluster","Cell line expression cluster","Single cell expression cluster","Interactions","Subcellular main location","Subcellular additional location","Antibody RRID","Cancer prognostics - Bladder Urothelial Carcinoma (TCGA)","Cancer prognostics - Breast Invasive Carcinoma (TCGA)","Cancer prognostics - Breast Invasive Carcinoma (validation)","Cancer prognostics - Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (TCGA)","Cancer prognostics - Colon Adenocarcinoma (TCGA)","Cancer prognostics - Colon Adenocarcinoma (validation)","Cancer prognostics - Glioblastoma Multiforme (TCGA)","Cancer prognostics - Glioblastoma Multiforme (validation)","Cancer prognostics - Head and Neck Squamous Cell Carcinoma (TCGA)","Cancer prognostics - Kidney Chromophobe (TCGA)","Cancer prognostics - Kidney Renal Clear Cell Carcinoma (TCGA)","Cancer prognostics - Kidney Renal Clear Cell Carcinoma (validation)","Cancer prognostics - Kidney Renal Papillary Cell Carcinoma (TCGA)","Cancer prognostics - Liver Hepatocellular Carcinoma (TCGA)","Cancer prognostics - Liver Hepatocellular Carcinoma (validation)","Cancer prognostics - Lung Adenocarcinoma (TCGA)","Cancer prognostics - Lung Adenocarcinoma (validation)","Cancer prognostics - Lung Squamous Cell Carcinoma (TCGA)","Cancer prognostics - Lung Squamous Cell Carcinoma (validation)","Cancer prognostics - Ovary Serous Cystadenocarcinoma (TCGA)","Cancer prognostics - Ovary Serous Cystadenocarcinoma (validation)","Cancer prognostics - Pancreatic Adenocarcinoma (TCGA)","Cancer prognostics - Pancreatic Adenocarcinoma (validation)","Cancer prognostics - Prostate Adenocarcinoma (TCGA)","Cancer prognostics - Rectum Adenocarcinoma (TCGA)","Cancer prognostics - Rectum Adenocarcinoma (validation)","Cancer prognostics - Skin Cuteneous Melanoma (TCGA)","Cancer prognostics - Stomach Adenocarcinoma (TCGA)","Cancer prognostics - Testicular Germ Cell Tumor (TCGA)","Cancer prognostics - Thyroid Carcinoma (TCGA)","Cancer prognostics - Uterine Corpus Endometrial Carcinoma (TCGA)")


hpa_metadata_tb <- 
	hpa_metadata_ls |>
		dplyr::bind_rows() |>
		dplyr::distinct() |>
		dplyr::select(any_of(grep("Cancer prognostics", coi, value = TRUE, invert = TRUE)))

hpa_prognosis_tb <- 
	hpa_metadata_ls |>
		dplyr::bind_rows() |>
		dplyr::distinct() |>
		dplyr::select(any_of(grep("Cancer prognostics", coi, value = TRUE)))

hpa_single_cell_tb <- 
	hpa_metadata_ls |>
		dplyr::bind_rows() |>
		dplyr::distinct() |>
		dplyr::select(any_of(grep("single cell", coi, value = TRUE)))

hpa_rna_tissue_tb <- 
	hpa_metadata_ls |>
		dplyr::bind_rows() |>
		dplyr::distinct() |>
		dplyr::select(any_of(grep("RNA tissue", coi, value = TRUE)))

hpa_rna_cancer_tb <- 
	hpa_metadata_ls |>
		dplyr::bind_rows() |>
		dplyr::distinct() |>
		dplyr::select(any_of(grep("RNA cancer", coi, value = TRUE)))

hpa_rna_blood_tb <- 
	hpa_metadata_ls |>
		dplyr::bind_rows() |>
		dplyr::distinct() |>
		dplyr::select(any_of(grep("RNA blood cell", coi, value = TRUE)))

hpa_rna_cellline_tb <- 
	hpa_metadata_ls |>
		dplyr::bind_rows() |>
		dplyr::distinct() |>
		dplyr::select(any_of(grep("RNA cell line", coi, value = TRUE)))
logger::log_info("Complete: get_hpa")
## ---- end

## @knitr get_healthy_data
files <- 
	c(
		"normal_ihc_data.tsv.zip", 
		"normal_ihc_tissues.tsv.zip", 
		"rna_tissue_consensus.tsv.zip",
		"rna_tissue_consensus_tissues.tsv.zip",
		"rna_tissue_hpa.tsv.zip",
		"rna_tissue_hpa_tissues.tsv.zip",
		"transcript_rna_tissue.tsv.zip",
		"rna_tissue_hpa_samples.tsv.zip",
		"rna_tissue_gtex.tsv.zip",
		"rna_tissue_detail_gtex.tsv.zip",
		"rna_tissue_gtex_tissues.tsv.zip",
		"rna_tissue_detail_gtex_tissues.tsv.zip",
		"rna_tissue_detail_gtex_samples.tsv.zip",
		"transcript_rna_gtexretina.tsv.zip",
		"rna_tissue_fantom.tsv.zip",
		"rna_tissue_detail_fantom.tsv.zip",
		"rna_tissue_fantom_tissues.tsv.zip",
		"rna_tissue_detail_fantom_tissues.tsv.zip",
		"rna_tissue_detail_fantom_samples.tsv.zip"
	)
files <- 	
	lapply(
		X = files, 
		FUN = check_hpa_size,
		limit = 1 # Only retrieve files smaller than 1Gb
	)
files <- unlist(Filter(Negate(is.null), files))	
normal_ls <- 
	lapply(
		X = files, 
		FUN = hpa_download, 
		target_ids = ensembl_id
	)
names(normal_ls) <- gsub(".tsv.zip", "", files)
logger::log_info("Complete: get_healthy_data")
## ---- end

## @knitr prepare_healthy_data
tmp_gtex <-
	normal_ls$rna_tissue_gtex |>
		dplyr::filter(Gene %in% target_ensembl_ids$GENEID) 
tmp_hpa <-
	normal_ls$rna_tissue_hpa |>
		dplyr::filter(Gene %in% target_ensembl_ids$GENEID)
tmp_fantom <-
	normal_ls$rna_tissue_fantom |>
		dplyr::filter(Gene %in% target_ensembl_ids$GENEID)
tmp_ihc <-
	normal_ls$normal_ihc_data |>
		dplyr::filter(Gene %in% target_ensembl_ids$GENEID) |>
		dplyr::mutate(Tissue = tolower(Tissue)) |>
		dplyr::mutate(Level = factor(Level, levels = c("Not detected", "Low", "Medium", "High")))
		
rna_tissue_merged <- 
	dplyr::full_join(tmp_gtex, tmp_hpa, by = c("Gene", "Gene name", "Tissue"), suffix = c(".gtex", ".hpa")) |>
		dplyr::full_join(tmp_fantom, by = c("Gene", "Gene name", "Tissue")) |>
		dplyr::select(Gene, `Gene name`, `Tissue`, contains("nTPM"), contains("Scaled")) |>
	dplyr::full_join(tmp_ihc, by = c("Gene", "Gene name", "Tissue")) |>
	dplyr::left_join(normal_ls$rna_tissue_detail_gtex_tissues, by = "Tissue", relationship =
  "many-to-many") |>
	tidyr::replace_na(list(Organ = "NA")) |>
	tidyr::replace_na(list(`nTPM.gtex` = -10)) |>
	tidyr::replace_na(list(`nTPM.hpa` = -10)) |>
	dplyr::arrange(Level) 
logger::log_info("Complete: prepare_healthy_data")
## ---- end

## @knitr get_tumour_data
files <- 
	c(
		"rna_cancer_sample.tsv.gz", # >1Gb
		"cancer_cptac.tsv.zip",
		"cancer_cptac_cancers.tsv.zip",
		"cancer_data.tsv.zip",
		"cancer_cancers.tsv.zip",
		"cancer_prognostic_data.tsv.zip",
		"cancer_rna_tcga_cancers.tsv.zip"	
	)
files <- 	
	lapply(
		X = files, 
		FUN = check_hpa_size,
		limit = 1 # Only retrieve files smaller than 1Gb
	)
files <- unlist(Filter(Negate(is.null), files))
tumour_ls <- 
	lapply(
		X = files, 
		FUN = hpa_download, 
		target_ids = ensembl_id
	)
names(tumour_ls) <- gsub(".tsv.zip", "", files)	
logger::log_info("Complete: get_tumour_data")
## ---- end

## @knitr get_singlecell_data
files <- 
	c(
		"rna_single_cell_read_count.zip",
		"rna_single_cell_datasets.tsv.zip",
		"rna_single_cell_type.tsv.zip",
		"rna_single_cell_type.tsv.zip",
		"rna_single_cell_cluster.tsv.zip",
		"rna_single_cell_clusters.tsv.zip"
	)
files <- 	
	lapply(
		X = files, 
		FUN = check_hpa_size,
		limit = 1 # Only retrieve files smaller than 1Gb
	)
files <- unlist(Filter(Negate(is.null), files))	
# singlecell_ls <- 
	# lapply(
		# X = files, 
		# FUN = hpa_download, 
		# target_ids = ensembl_id
	# )
# names(singlecell_ls) <- gsub(".tsv.zip", "", files)
logger::log_info("Complete: get_singlecell_data")
## ---- end

## @knitr get_blood_data
files <- 
	c(
		"blood_pea_diseases.tsv.zip",
		"blood_pea_disease_de.tsv.zip",
		"blood_immunoassay_concentration.tsv.zip",
		"blood_ms_concentration.tsv.zip"
	)
files <- 	
	lapply(
		X = files, 
		FUN = check_hpa_size,
		limit = 1 # Only retrieve files smaller than 1Gb
	)
files <- unlist(Filter(Negate(is.null), files))	
blood_ls <- 
	lapply(
		X = files, 
		FUN = hpa_download, 
		target_ids = ensembl_id
	)
names(blood_ls) <- gsub(".tsv.zip", "", files)
logger::log_info("Complete: get_blood_data")
## ---- end


## @knitr get_cell_data
files <- 
	c(
		"rna_cell_line_cancer.tsv.zip",
		"rna_cell_line_cancer_cancers.tsv.zip",
		"rna_celline.tsv.zip",
		"rna_cell_line_tcga_comparison.tsv.zip",
		"cell_line_analysis_data.tsv.zip"
	)
files <- 	
	lapply(
		X = files, 
		FUN = check_hpa_size,
		limit = 1 # Only retrieve files smaller than 1Gb
	)
files <- unlist(Filter(Negate(is.null), files))	
celllines_ls <- 
	lapply(
		X = files, 
		FUN = hpa_download, 
		target_ids = ensembl_id
	)
names(celllines_ls) <- gsub(".tsv.zip", "", files)
logger::log_info("Complete: get_cell_data")
## ---- end


## @knitr get_tau_data
tau <- tau_download()
logger::log_info("Complete: get_tau_data")
## ---- end

## @knitr get_tcga_data
# N.B. Could download whole array directly; this allows subsetting
# N.B. Loads more data available (ICGC/PCAWG/GDC, single cell + modalities)
# TODO: add Head and Neck Cancer (Puram 2017), MAGIC (Medulloblastoma), Breast Cancer (Yau 2010), HMF (not in Xena), ...


# xena_tcga <- ucsc_tcga_download(ensembl_id_df = ensembl_id) # comprehensive! Esp. clinical data for BRCA.
xena_toil <- ucsc_toil_download(ensembl_id_df = ensembl_id) # for formal gtex vs tcga
xena_pancan <- ucsc_pancan_download(x = params$target_id) # batch corrected
xena_gdc_pancan <- ucsc_gdc_pancan_download(ensembl_id_df = ensembl_id) # largest dataset; for cohort-level analyses, e.g. SCNA ~ EXP or OS ~ EXP analysis but has no covariates (purity)
# xena_icgc <- ucsc_icgc_download(ensembl_id_df = ensembl_id) # ICGC
# xena_pcawg <- ucsc_pcawg_download(ensembl_id_df = ensembl_id) # PCAWG, purities
xena_met500 <- ucsc_met500_download(ensembl_id_df = ensembl_id) # 868 metastatic samples

## Pre-canned annotations from XenaShiny
tcga_purity <- UCSCXenaShiny::load_data("tcga_purity") # https://www.nature.com/articles/ncomms9971#Sec14
toil_info <- UCSCXenaShiny::load_data("toil_info")
tcga_tmb <- UCSCXenaShiny::load_data("tcga_tmb")
tcga_subtypes <- UCSCXenaShiny::load_data("tcga_subtypes")
tcga_gtex <- UCSCXenaShiny::load_data("tcga_gtex")
tcga_genome_instability <- UCSCXenaShiny::load_data("tcga_genome_instability")
pcawg_purity <- UCSCXenaShiny::load_data("pcawg_purity")



# vis_pcawg_dist(Gene = params$target_id)

# Visualize Correlation between Gene and Tumor Stemness
# UCSCXenaShiny::vis_gene_stemness_cor(Gene = params$target_id, cor_method = "spearman", data_type = "mRNA", Plot = "TRUE", opt_pancan = .opt_pancan )

# UCSCXenaShiny::vis_gene_immune_cor(Gene = params$target_id, cor_method = "spearman", data_type = "mRNA", Immune_sig_type = "Cibersort", Plot = "TRUE", opt_pancan = .opt_pancan )

# UCSCXenaShiny::vis_gene_drug_response_asso(Gene = params$target_id, x_axis_type =  "median.diff", output_form =  "ggplot2" )


# UCSCXenaShiny::vis_ccle_tpm(
	# Gene = params$target_id, data_type = c("mRNA", "protein", "cnv")[1], use_log = TRUE, opt_pancan = .opt_pancan
# )

## Survival analysis
# tcga_surv_get()

logger::log_info("Complete: get_tcga_data")
## ---- end


## @knitr toil_tumour_normal_plt

# List objects:
			# rsem_deseq_norm_count , 
			# rsem_norm_count , 
			# transcript_expected_count , 
			# rsem_isoform_tpm , 
			# clinical_data 

TcgaTargetGTEX_ge <- 
	xena_toil$rsem_deseq_norm_count |>
	tibble::as_tibble(rownames = "sample") # genes
TcgaTargetGTEX_te <- 
	xena_toil$transcript_expected_count |>
	tibble::as_tibble(rownames = "sample") # transcripts
TcgaTargetGTEX_phenotype <-
	xena_toil$clinical_data$TcgaTargetGTEX_phenotype.txt.gz |>
	janitor::clean_names() |>
	dplyr::mutate(primary_site = iconv(primary_site, from = "ISO-8859-1", to = "UTF-8")) |>
	dplyr::mutate(sample_type = iconv(sample_type, from = "ISO-8859-1", to = "UTF-8")) |>
	dplyr::filter(!grepl("blood", primary_site, ignore.case = TRUE)) |>
	dplyr::filter(!grepl("cell line", sample_type, ignore.case = TRUE)) |>
	tidyr::drop_na(primary_site)
TcgaTargetGTEX_phenotype$sample_type2 <-
	dplyr::recode(TcgaTargetGTEX_phenotype$sample_type,
		`Primary Tumor` = "primary", 
		`Additional - New Primary` = "primary", 
		`Primary Solid Tumor` = "primary", 
		`Solid Tissue Normal` = "normal",                         
		`Normal Tissue` = "normal",                         
		`Recurrent Tumor` = "recurrence",                            
		`Recurrent Solid Tumor` = "recurrence",                            
		`Recurrent Blood Derived Cancer - Bone Marrow` = "recurrence",     
		`Metastatic` = "metastatic",                                   
		`Additional Metastatic` = "metastatic",                                   
		`Control Analyte` = "control"
	) 
key_indication <- c("Lung", "Esophagus", "Stomach", "Colon", "Head and Neck region", "Pancreas", "Breast", "Bladder", "Ovary", "Prostate", "Liver", "Uterus", "Soft tissue,Bone")	

tcga_target_gtex_merged <- 	
	TcgaTargetGTEX_phenotype |>
	dplyr::left_join(TcgaTargetGTEX_ge, by = "sample") |>
	dplyr::mutate(priority = ifelse(primary_site %in% key_indication, "Core", "Ancillary")) |>
	dplyr::mutate(priority = factor(priority, levels = c("Core", "Ancillary"))) |>
	dplyr::left_join(TcgaTargetGTEX_te, by = "sample") 

# TcgaTargetGTEX_phenotype |>
	# dplyr::filter(study == "TARGET") |>
	# # dplyr::group_by(primary_disease_or_tissue) |>
	# dplyr::count() |>
	# View()
	
	# dplyr::select(detailed_category) |>
	# dplyr::distinct() |>
	# dplyr::arrange() |>
	# View()
	
# Ridgeline plot
# tcga_target_gtex_merged |>
	# dplyr::count(primary_site, detailed_category)

# TODO: embed legend (automatable?), in healthy tissue LHS
# TODO: add more quantiles to tumour plots
ridge_tumour_plt <- plot_expression_ridge(plt_tumour = TRUE)
ridge_healthy_plt <- plot_expression_ridge(plt_tumour = FALSE)

logger::log_info("Complete: toil_tumour_normal_plt")
## ---- end


## @knitr xena_met500
# xena_met500
# List objects:
	# ge (log2(fpkm+0.001))
	# clinical_data 

met500_ge <- 
	xena_met500$ge |>
	tibble::as_tibble(rownames = "sample") # genes
met500_phenotype <-
	xena_met500$clinical_data |>
	janitor::clean_names() |>
	dplyr::mutate(tissue = iconv(tissue, from = "ISO-8859-1", to = "UTF-8")) |>
	dplyr::filter(!grepl("blood", tissue, ignore.case = TRUE)) |>
	dplyr::filter(!grepl("cell line", tissue, ignore.case = TRUE)) |>
	tidyr::drop_na(tissue) |>
	dplyr::rename("sample" = "sample_id")
	
key_indication <- tolower(c("Lung", "Esophagus", "Stomach", "Colon", "Oral", "Pancreas", "Breast", "Bladder", "Ovary", "Prostate", "Liver", "Uterus", "Soft tissue,Bone"))	

met500_merged <- 	
	met500_phenotype |>
	dplyr::left_join(met500_ge, by = "sample") |>
	dplyr::mutate(priority = ifelse(tissue %in% key_indication, "Core", "Ancillary")) |>
	dplyr::mutate(priority = factor(priority, levels = c("Core", "Ancillary"))) 
	
ridge_met500_plt <- plot_expression_met500_ridge(x = met500_merged)
logger::log_info("Complete: xena_met500")
## ---- end


## @knitr xena_met500_enrichment

keep <-
	met500_merged |>
	dplyr::group_by(cohort) |>
	dplyr::count() |>
	dplyr::filter(n > 10)
to_plot <-
	met500_merged |>
	dplyr::filter(cohort %in% keep[["cohort"]])
X <- 
	to_plot |>
	dplyr::select(all_of(as.character(target_ensembl_ids))) |>
	dplyr::mutate("s2" := !!sym(as.character(target_ensembl_ids)) + 0.001) |>
	as.matrix()
rownames(X) <- paste0("g", 1:nrow(X))
tmp <-	
	to_plot |>
	dplyr::group_by(cohort) |>
	dplyr::group_data()	|>
	unnest(.rows) 
gs_ls <- 
	split(tmp$`.rows`, f = tmp$cohort) |>
	lapply(FUN = function(x){paste0("g", x)})
rel_enrichment <- gsva(ssgseaParam(X, gs_ls), verbose=FALSE)

met500_enrichment_plt <- 
	rel_enrichment |>
		as.data.frame() |>
		tibble::rownames_to_column() |>
		ggplot(aes(x = forcats::fct_reorder(rowname, s2), y = s2)) + 
		geom_col(color = BTX_COL$blue_d_out, fill = BTX_COL$blue_d_in) +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
		theme(plot.background = element_rect(color = NA, fill = "grey97")) 

logger::log_info("Complete: xena_met500_enrichment")

## ---- end


## @knitr ridge_emd

indications_tb <-
	tcga_target_gtex_merged |>
		dplyr::filter(study != "GTEX") |> 
		dplyr::count(primary_disease_or_tissue)
emd_t_ls <- emd_ls <- kidney_ls <- vector("list", nrow(indications_tb))
names(emd_t_ls) <- names(emd_ls) <- names(kidney_ls) <- indications_tb$primary_disease_or_tissue
for (ind in indications_tb$primary_disease_or_tissue) {
	emd_ls[[ind]] <- 
		get_emd(x = ind, ref = NULL, use_isoforms = FALSE) 
	emd_t_ls[[ind]] <-
		get_emd(x = ind, ref = NULL, use_isoforms = TRUE) 
	kidney_ls[[ind]] <-	
		get_emd(x = ind, ref = c("Kidney", "Salivary Gland"), use_isoforms = FALSE) 		
}

# Remove all values where expression is lower than in healthy (set to 0)
m <- prepare_emd_plot(x = emd_ls, gene_rename = TRUE)
m_kid <- prepare_emd_plot(x = kidney_ls, gene_rename = TRUE)

if (FALSE) {
	col_fun <- colorRamp2(range(m), c(BTX_COL$blue_d_in, BTX_COL$blue_d_out))
}else{
	col_fun <- colorRampPalette(brewer.pal(9, "Purples"))(255)
}

# Prepare heatmap annotations

# TODO: add prevalence from globocan to row annotation (fuzzy match non-trivial)
# Use ICD-10 codes (an NLP job), or those in Xena phenotype tables in download

# tcga_target_gtex_merged |> 
	# dplyr::select("detailed_category", "incidence", "mortality") |>
	# dplyr::distinct() |>
	# tidyr::drop_na(incidence)


emd_hm <- 
	ComplexHeatmap::Heatmap(
	  m$m, 
	  col = col_fun,
	  row_km = 6, # Cluster rows into 6 groups and add space between them
	  column_split = factor(colnames(m$m) == params$target_id),
	  name = "EMD (positive only)",
	  na_col = "black",
	  row_title = NULL,
	  column_title = NULL,
	  bottom_annotation = m$column_ha,
	  row_names_gp = gpar(fontsize = 10),
	  show_row_names = TRUE,
	  show_column_names = TRUE
	)

emd_kidney_hm <- 
	ComplexHeatmap::Heatmap(
	  m_kid$m, 
	  col = col_fun,
	  row_km = 6, # Cluster rows into 6 groups and add space between them
	  column_split = factor(colnames(m$m) == params$target_id),
	  name = "EMD (positive only)",
	  na_col = "black",
	  row_title = NULL,
	  column_title = NULL,
	  bottom_annotation = m_kid$column_ha,	  
	  row_names_gp = gpar(fontsize = 10),
	  show_row_names = TRUE,
	  show_column_names = TRUE
	)
logger::log_info("Complete: ridge_emd")

## ---- end


## @knitr toil_iso_tumour_normal_plt

emd_isoforms_m <- prepare_emd_plot(x = emd_t_ls) # From ridge_emd chunk
emd_isoforms_m <- emd_isoforms_m$m
target_isoform_ids <-
	ensembl_id |>
	dplyr::filter(src == "candidate") |>
	dplyr::select(TXID) |>
	dplyr::distinct() |>
	unlist()
idx <- 
	grepl(paste(target_isoform_ids, collapse = "|"), colnames(emd_isoforms_m))
emd_isoforms_m <- emd_isoforms_m[, idx]
colnames(emd_isoforms_m) <- tools::file_path_sans_ext(colnames(emd_isoforms_m))

transcript_tau <-
	tau$tau_transcript_V8 |>
	dplyr::mutate(transcript_id = tools::file_path_sans_ext(transcript_id)) |>
	dplyr::filter(transcript_id %in% colnames(emd_isoforms_m)) |>
	dplyr::select(transcript_id, tau)

transcript_tau <-
	transcript_tau[match(colnames(emd_isoforms_m), transcript_tau$transcript_id), ]

column_ha <- 
	ComplexHeatmap::HeatmapAnnotation(
		tau_healthy = anno_barplot(transcript_tau$tau)
	)

emd_isoforms_hm <- 
	ComplexHeatmap::Heatmap(
	  emd_isoforms_m, 
	  col = col_fun,
	  row_km = 6, # Cluster rows into 6 groups and add space between them
	  top_annotation = column_ha,
	  name = "EMD (positive only)",
	  na_col = "black",
	  row_title = "",
	  column_title = "",
	  row_names_gp = gpar(fontsize = 10),
	  show_row_names = TRUE,
	  show_column_names = TRUE
	)
logger::log_info("Complete: toil_iso_tumour_normal_plt")
## ---- end

## @knitr get_globocan_plt
filter_out <-
	c(
		"All cancers",
		"Non-Hodgkin lymphoma",
		"Hodgkin lymphoma",
		"Leukaemia",
		"Multiple myeloma",
		"Vulva",
		"Vagina",
		"Penis"
	)

to_plot <-
  globocan_data |>
	dplyr::filter(!(Label %in% filter_out)) |>
	dplyr::mutate(across(where(is.character), ~ stringr::str_replace_all(.x, "\xa0", " ")))

emd_tb <-
	emd_ls |>
	dplyr::bind_rows(.id = "TCGA") |> 
	dplyr::mutate(SYMBOL = tools::file_path_sans_ext(SYMBOL)) |> 
	dplyr::filter(SYMBOL %in% target_ensembl_ids)

emd_annot <- rep(NA, nrow(to_plot))	
for (i in 1:nrow(to_plot)) {
	bait <-
		gsub("; ", "|", to_plot[i, "TCGA"])
	emd_annot[i] <-	
		emd_tb |>
		dplyr::filter(grepl(bait, TCGA)) |>
		dplyr::summarise(n = max(emd)) |>
		unlist()
}
emd_annot[is.infinite(emd_annot)] <- NA

globocan_plt <- 
	to_plot |>
	dplyr::bind_cols(emd_annot) |>
	ggplot(aes(x = incidence, y = mortality, size = mortality, label = Label, fill = emd_annot, color = emd_annot)) +
    geom_point(shape = 21, alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "#ececec") +  
    coord_trans(x = "log10", y = "log10") +
	# scale_colour_gradient(name = "EMD", low = BTX_COL$blue_l_out, high = BTX_COL$blue_d_out) +
	# scale_fill_gradient(name = "EMD", low = BTX_COL$blue_l_in, high = BTX_COL$blue_d_in) +
	paletteer::scale_colour_paletteer_c(name = "EMD", "harrypotter::ravenclaw2", direction = -1) +
	paletteer::scale_fill_paletteer_c(name = "EMD", "harrypotter::ravenclaw2", direction = -1) +
	scale_x_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_continuous(breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
	ggrepel::geom_text_repel(size = 3.5) +
	theme(
		plot.background = element_rect(color = NA, fill = "grey97")
	) 
logger::log_info("Complete: get_globocan_plt")				
## ---- end


## @knitr plot_tau_data
get_refs <-
	tau$tau_gene_V8 |>
		dplyr::filter(gene_id %in% ensembl_id$GENEID)

get_refs <-
	dplyr::left_join(ensembl_id, get_refs, by = c("GENEID" = "gene_id")) |>
	na.omit() |>
	dplyr::distinct(GENEID, .keep_all = TRUE) |>
	dplyr::mutate(highlight = ifelse(src == "candidate", BTX_COL$pink_out, BTX_COL$blue_l_out)) |>
	dplyr::mutate(highlight2 = ref_cols[match(src, names(ref_cols))]) 
line_ref <- 
	get_refs |>
		dplyr::filter(SYMBOL == params$target_id)

tau_txt <- line_ref$tau

tau_plt <-
	tau$tau_gene_V8 |>
		dplyr::filter(!is.na(tau)) |>
		ggplot(aes(x = tau)) +
		geom_density(color = BTX_COL$gray_out) +
		geom_point(data = get_refs, aes(x = tau, y = 0, color = highlight2)) +
		geom_vline(data = line_ref, aes(xintercept = tau, color = highlight2)) +
		ggrepel::geom_text_repel(data = get_refs, aes(y = 0, x = tau, label = SYMBOL, color = highlight2), size = 3.5) +
		scale_color_identity() +
		# scale_color_manual(ref) +
		theme(
			plot.background = element_rect(color = NA, fill = "grey97")
		) +
		annotate(
			"richtext",
			x = c(0.15, 0.85),
			y = c(1.5, 3),
			label = c("'Housekeeping'<br>genes", "Tissue-specfic<br>genes"),
			fill = NA, label.size = 0, family = "lato", size = 2, vjust = 1,
		) 
logger::log_info("Complete: plot_tau_data")		
## ---- end 

## @knitr tau_ntpm
max_exp_gtex <- 
	normal_ls$rna_tissue_gtex |>
		dplyr::group_by(Gene) |> 
		dplyr::summarise(score = max(`nTPM`)) |>
		dplyr::filter(Gene %in% ensembl_id$GENEID) |>
		dplyr::rename(GENEID = Gene) |>
		dplyr::left_join(ensembl_id, by = "GENEID") |> 
		dplyr::select(GENEID, src, SYMBOL, score) |>
		dplyr::distinct(,.keep_all=TRUE) |>
		dplyr::mutate(highlight = ref_cols[match(src, names(ref_cols))]) |>
		tidyr::replace_na(list(SYMBOL = BTX_COL$gray_out))

to_plot <-
	tau$tau_gene_V8	|>
		dplyr::select(gene_id, tau) |>
		dplyr::rename(GENEID = gene_id) |>
		dplyr::filter(GENEID %in% ensembl_id$GENEID) |>
		dplyr::left_join(max_exp_gtex, by = "GENEID") |>
		dplyr::filter(SYMBOL != "GAPDH")

tau_plt2 <-
	to_plot |>
		ggplot(aes(x = tau, y = score, color = highlight, label = SYMBOL)) +
		geom_point() +
		ggrepel::geom_text_repel(size = 3.5) +
		scale_color_identity() +
		theme(
			plot.background = element_rect(color = NA, fill = "grey97")
		)
logger::log_info("Complete: tau_ntpm")				
## ---- end

## @knitr hm_healthy

to_plot <-
	rna_tissue_merged |>
	#tidyr::drop_na(Level) |> # This crashes if no IHC data available
	dplyr::ungroup() |>
	dplyr::mutate(grp = as.factor(Tissue)) |>
	dplyr::group_by(grp) |>
	dplyr::slice_max(order_by = Level, n = 1) |>
	dplyr::select(Organ, Tissue, Level, nTPM.gtex, nTPM.hpa) |>
	dplyr::distinct(.keep_all = TRUE) |>
	dplyr::left_join(wahl_data, by = "Tissue")


hm_healthy_plt <- 
	to_plot |>
	ggpubr::ggscatter(
		,
		x = "nTPM.gtex",
		y = "nTPM.hpa",
		color = "EBRT_limits_gy", 
		label = "Tissue", 
		repel = TRUE,
		facet.by = "Level",
		alpha = 0.75,
		ggtheme = theme_minimal()
	) + 
	# theme(legend.position = "none") +
	# scale_colour_viridis_c() + 
	paletteer::scale_colour_paletteer_c(name = "EBRT limits (Gy)", "harrypotter::ravenclaw2") +
	theme(
		plot.background = element_rect(color = NA, fill = "grey97"),
		panel.border = element_rect(color = BTX_COL$gray_out, fill = NA, size = 0.5)
	)

## GT table (to solve over-plotting issue)
hm_healthy_gt <- 
	to_plot |>
	dplyr::ungroup() |>
	dplyr::mutate(Level = as.character(Level)) |>
	tidyr::replace_na(list(Level = "Unknown")) |>
	dplyr::mutate(Level = stringr::str_c("IHC:", Level)) |> 
	dplyr::mutate(Level = factor(Level, levels = c("IHC:High", "IHC:Medium", "IHC:Low", "IHC:Not detected", "IHC:Unknown"))) |>
	dplyr::select(`Tissue`, `EBRT_limits_gy`, `Level`, `nTPM.gtex`, `nTPM.hpa`, ) |>
	dplyr::group_by(Level) |>
	gt() |>
	gtExtras::gt_plt_bar(
		column = nTPM.gtex,
		color = BTX_COL$blue_l_out,
		fill = BTX_COL$blue_l_in,
		keep_column = FALSE,
		scale_type = "number",
		text_color = "white"
	) |>
	gtExtras::gt_plt_bar(
		column = nTPM.hpa,
		color = BTX_COL$blue_d_in,
		fill = BTX_COL$blue_d_out,
		scale_type = "number",
		text_color = "white"
	) |>
	gt::data_color(
	  columns = `EBRT_limits_gy`,
	  target_columns = `Tissue`,
	  method = "numeric",
	  # domain = c(0, 1),
	  palette = "GnBu"
	) |>	
	gt::cols_align(columns = EBRT_limits_gy, align = "left") |>
	gt::cols_label(
		`EBRT_limits_gy` = "EBRT Limit (Gy)",
		`nTPM.hpa` = "HPA",
		`nTPM.gtex` = "GTEx"
	) |>
	gt::row_group_order(groups = c("IHC:High", "IHC:Medium", "IHC:Low", "IHC:Not detected", "IHC:Unknown")) |>
    gtExtras::gt_theme_guardian() |>
	gt::tab_source_note("Source: The Human Protein Atlas (v24.0); Wahl et al (2021)") 


logger::log_info("Complete: hm_healthy")				

## ---- end

## @knitr ihc_tb
# replicates Protein Expression Overview Barplot
# TODO: normal_ihc_tissues
tmp_ihc <-
	normal_ls$normal_ihc_data |>
		dplyr::filter(Gene %in% target_ensembl_ids$GENEID) |>
		dplyr::mutate(Tissue = tolower(Tissue)) |>
		dplyr::mutate(Level = factor(Level, levels = c("Not detected", "Low", "Medium", "High")))

switch_prot_exp_overview <- nrow(tmp_ihc) > 0
		
prot_exp_overview <- 
	tmp_ihc |>
		dplyr::group_by(Level) |> 
		dplyr::count(`Tissue`) |>
		ggplot(aes(x = `Tissue`, y = n)) +
		geom_col(aes(n, Tissue), fill = BTX_COL$blue_d_in, color = BTX_COL$blue_d_out, width = 0.6) +
		facet_wrap(~Level, scales = "free_y", drop = FALSE) +
		theme(
			plot.background = element_rect(color = NA, fill = "grey97")
		)		
## ---- end


## @knitr transcript_rna_tissue
tmp_isoform <-
	normal_ls$transcript_rna_tissue |>
		dplyr::filter(ensgid %in% unlist(target_ensembl_ids)) |>
		dplyr::select(ensgid, enstid, contains("TPM")) |>
		tidyr::pivot_longer(cols = contains("TPM"), names_to = "tissue", values_to = "TPM") |> 
		tidyr::separate(col = tissue, into = c(NA, "tissue", NA), sep = "\\.") |>
		dplyr::group_by(tissue, enstid) |>
		dplyr::summarise(median_TPM = median(TPM)) |>
		tidyr::pivot_wider(names_from = tissue, values_from = median_TPM) 
		
expression_filter <-		
	tmp_isoform |>
		dplyr::select(!enstid) |> 
		apply(FUN = function(x){mean(x) >= 0.5}, MARGIN = 1)
tmp_isoform <-		
	tmp_isoform[expression_filter, ] # filter transcripts
	
to_plot <- 
	tmp_isoform |>
		dplyr::select(!enstid) |>
		as.matrix() |>
		scale(center = TRUE, scale = TRUE)
rownames(to_plot) <- tmp_isoform$enstid

to_plot <-
	to_plot[, rowSums(is.na(t(to_plot))) < nrow(to_plot)] # filter out NaN cols

# Tissue clustering not actiated, to avoid kmeans errors 
switch_isoform_hm <- nrow(to_plot) > 1
isoform_hm <- 
	ComplexHeatmap::Heatmap(
	  t(to_plot), 
	  col = HEATMAP_COL,
	  name = "z-score (TPM)",
	  na_col = "black",
	  row_title = "",
	  column_title = "",
	  row_names_gp = gpar(fontsize = 6),
	  show_row_names = TRUE,
	  show_column_names = TRUE
	)
logger::log_info("Complete: transcript_rna_tissue")				
## ---- end


## @knitr retina_boxplot
tmp_retina <-
	normal_ls$transcript_rna_gtexretina |>
	dplyr::filter(ensgid %in% ensembl_id$GENEID) |>
	dplyr::select(ensgid, contains("TPM")) |>
	dplyr::distinct(ensgid, .keep_all = TRUE) |>
	dplyr::left_join(ensembl_id, by = c("ensgid" = "GENEID")) |>
	dplyr::distinct(ensgid, .keep_all = TRUE) |>	
	dplyr::filter(SYMBOL != "GAPDH") |> 
	dplyr::mutate(highlight2 = ref_cols[match(src, names(ref_cols))]) 

retina_plt <-
	tmp_retina |>
		tidyr::pivot_longer(cols = contains("TPM"), names_to = "sampleid", values_to = "TPM") |>
		ggplot(aes(x = SYMBOL, y = TPM, col = highlight2)) +
		geom_boxplot() +
		scale_color_identity() +
		coord_trans(y = "sqrt") + 
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
		theme(
			plot.background = element_rect(color = NA, fill = "grey97")
		)

logger::log_info("Complete: retina_boxplot")				
## ---- end


## @knitr tau_transcript_V8
get_refs <-
	tau$tau_transcript_V8 |>
		dplyr::filter(transcript_id %in% tmp_isoform$enstid) |>
		dplyr::left_join(tmp_isoform, by = c("transcript_id" = "enstid")) |>
		dplyr::select(transcript_id, tau) |>
		tidyr::drop_na(tau)

# annot_tbl <-
	# ensembl_id |>
	# dplyr::select(GENEID, SYMBOL) |>
	# dplyr::distinct(.keep_all = TRUE)

# get_refs <-
	# dplyr::left_join(annot_tbl, get_refs, by = c("GENEID" = "ensgid")) |>
	# na.omit() |>
	# dplyr::mutate(highlight = ifelse(SYMBOL == params$target_id, BTX_COL$pink_out, BTX_COL$blue_l_out))
tau_isoform_plt <-
	tau$tau_transcript_V8 |>
		dplyr::filter(!is.na(tau)) |>
		ggplot(aes(x = tau)) +
		geom_density(color = BTX_COL$gray_out) +
		geom_point(data = get_refs, aes(x = tau, y = 0), col = BTX_COL$pink_out) +
		ggrepel::geom_text_repel(data = get_refs, aes(y = 0, x = tau, label = transcript_id), size = 3.5, angle = 90) +
		scale_color_identity() 

logger::log_info("Complete: tau_transcript_V8")				
## ---- end


## @knitr plot_bld1
get_refs <-
	blood_ls$blood_ms_concentration  |>
		dplyr::filter(`ENSG ID` %in% ensembl_id$GENEID)
annot_tbl <-
	ensembl_id |>
	dplyr::select(GENEID, SYMBOL, src) |>
	dplyr::distinct(.keep_all = TRUE)
get_refs <-
	dplyr::left_join(annot_tbl, get_refs, by = c("GENEID" = "ENSG ID")) |>
	na.omit() |>
	dplyr::mutate(highlight = ifelse(SYMBOL == params$target_id, BTX_COL$pink_out, BTX_COL$blue_l_out)) |>
	dplyr::mutate(highlight2 = ref_cols[match(src, names(ref_cols))]) 

line_ref <- 
	get_refs |>
		dplyr::filter(SYMBOL == params$target_id)

bld_plt <-
	blood_ls$blood_ms_concentration |>
		dplyr::filter(!is.na(`Conc [pg/L]`)) |>
		ggplot(aes(x = `Conc [pg/L]`)) +
		geom_density(color = BTX_COL$gray_out) +
		geom_vline(data = line_ref, aes(xintercept = `Conc [pg/L]`, color = highlight2)) +
		geom_vline(aes(xintercept = quantile(`Conc [pg/L]`)[2]), lty = 2) + 
		geom_vline(aes(xintercept = quantile(`Conc [pg/L]`)[3]), lty = 2) + 
		geom_vline(aes(xintercept = quantile(`Conc [pg/L]`)[4]), lty = 2) + 
		scale_x_continuous(
			trans = scales::log10_trans(),
			breaks = scales::trans_breaks("log10", function(x) 10^x),
			labels = scales::trans_format("log10", scales::math_format(10^.x))
		) + 
		geom_point(data = get_refs, aes(y = 0, x = `Conc [pg/L]`, color = highlight2)) +
		ggrepel::geom_text_repel(data = get_refs, aes(y = 0, x = `Conc [pg/L]`, label = SYMBOL, color = highlight2), size = 3.5) +
		scale_color_identity() +
		theme(
			plot.background = element_rect(color = NA, fill = "grey97")
		)
logger::log_info("Complete: plot_bld1")				
		
## ---- end 



## @knitr plot_bld2
get_refs <-
	blood_ls$blood_immunoassay_concentration  |>
		dplyr::filter(`ENSG ID` %in% ensembl_id$GENEID)
annot_tbl <-
	ensembl_id |>
	dplyr::select(GENEID, SYMBOL) |>
	dplyr::distinct(.keep_all = TRUE)
get_refs <-
	dplyr::left_join(annot_tbl, get_refs, by = c("GENEID" = "ENSG ID")) |>
	na.omit() |>
	dplyr::mutate(highlight = ifelse(SYMBOL == params$target_id, BTX_COL$pink_out, BTX_COL$blue_l_out))

bld_plt2 <-
	blood_ls$blood_immunoassay_concentration |>
		dplyr::filter(!is.na(`conc [pg/L]`)) |>
		ggplot(aes(x = `conc [pg/L]`)) +
		geom_density(color = BTX_COL$gray_out) +
		geom_vline(aes(xintercept = quantile(`conc [pg/L]`)[2]), lty = 2) + 
		geom_vline(aes(xintercept = quantile(`conc [pg/L]`)[3]), lty = 2) + 
		geom_vline(aes(xintercept = quantile(`conc [pg/L]`)[4]), lty = 2) + 
		scale_x_continuous(
			trans = scales::log10_trans(),
			breaks = scales::trans_breaks("log10", function(x) 10^x),
			labels = scales::trans_format("log10", scales::math_format(10^.x))
		) + 
		geom_vline(data = get_refs, aes(xintercept = `conc [pg/L]`, color = highlight)) +
		ggrepel::geom_text_repel(data = get_refs, aes(y = 0, x = `conc [pg/L]`, label = SYMBOL, color = highlight), size = 3.5) +
		scale_color_identity() +
		theme(
			plot.background = element_rect(color = NA, fill = "grey97")
		) 
logger::log_info("Complete: plot_bld2")						
## ---- end 


## @knitr purity

tcga_target_gtex_merged_purity <- 
	tcga_purity |>
		dplyr::select(sample, cancer_type, CPE) |>
		dplyr::left_join(tcga_target_gtex_merged, by = "sample")

idx <- 
	ensembl_id |>
		dplyr::filter(src == "candidate") |>
		dplyr::distinct(GENEID) |>
		unlist() |>
		paste(collapse = "|")

idx <- 
	grep(idx, colnames(tcga_target_gtex_merged_purity), value = TRUE)

purity_plt <- 
	tcga_target_gtex_merged_purity |>
		tidyr::drop_na(CPE) |>
		dplyr::group_by(cancer_type) |>
		dplyr::summarize(corn = cor(CPE, !!!syms(idx), use = "pairwise.complete.obs")) |>
		dplyr::mutate(lbl = ifelse(corn < 0, BTX_COL$pink_out, BTX_COL$blue_l_out)) |>
		ggplot(aes(x = forcats::fct_reorder(cancer_type, corn), y = corn, fill = lbl)) + 
		geom_col() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
		theme(plot.background = element_rect(color = NA, fill = "grey97")) +
		scale_fill_identity()
logger::log_info("Complete: purity")				

#### ---- end



## @knitr compgen_plt

compgen_plt <- 
  open_targets_data$compGenomicsQuery_ls |> 
	dplyr::mutate(is_human = speciesName == "Human") |>
	dplyr::filter(speciesName %in% c("Human", "Chimpanzee", "Macaque", "Mouse", "Rabbit", "Guinea Pig", "Dog", "Pig", "Norway rat - BN/NHsdMcwi")) |>
	dplyr::mutate(lbl = ifelse(is_human, targetGeneSymbol, speciesName)) |>
	dplyr::mutate(pnl = ifelse(is_human, "Human", "Non-Human")) |>
	ggplot(aes(x = queryPercentageIdentity, y = targetPercentageIdentity, size = targetPercentageIdentity, label = lbl)) +
    geom_point(fill = BTX_COL$blue_d_out, color = BTX_COL$blue_d_out, shape = 21, alpha = 0.5) +
	ggrepel::geom_text_repel(size = 3.5) +
	facet_wrap(vars(pnl), scales = "free") +
	theme(
		#panel.spacing = unit(.05, "lines"),
		legend.position = "none",
		plot.background = element_rect(color = NA, fill = "grey97"),
        panel.border = element_rect(color = BTX_COL$gray_out, fill = NA, size = 0.5)
	)
logger::log_info("Complete: compgen_plt")				
#### ---- end


## @knitr depmap_plt
font_family <- "lato"

summary_stats_psite <- 
	open_targets_data$depMapQuery_ls |>
	dplyr::group_by(diseaseFromSource) |>
	dplyr::summarise(med = median(geneEffect, na.rm = TRUE)) |>
	dplyr::arrange(desc(med))
tissue_order <- summary_stats_psite$diseaseFromSource


depmap_plt <- 
  open_targets_data$depMapQuery_ls |> 
	dplyr::filter(!grepl("pleomorphic", diseaseFromSource, ignore.case = TRUE)) |>
	dplyr::filter(!grepl("lymphoma", diseaseFromSource, ignore.case = TRUE)) |>
	dplyr::arrange(desc(median(geneEffect, na.rm = TRUE))) |>
	dplyr::mutate(
		diseaseFromSource = factor(diseaseFromSource, levels = tissue_order)
	) |>	
	ggplot(aes(x = diseaseFromSource, y = geneEffect)) +
	ggdist::stat_halfeye(fill_type = "segments", alpha = 0.3) +
	ggdist::stat_interval() +
	stat_summary(geom = "point", fun = median) +
	# stat_summary(
		# fun.data = get_n, 
		# geom = "text", 
		# family = font_family, 
		# size = 3, color = "grey75") +
		geom_hline(yintercept = -1, col = "grey30", lty = "dashed") +
		annotate("text", x = 0.5, y = -1, label = "Strongly Selective",
			   family = font_family, size = 3, hjust = 0) +
		scale_x_discrete(labels = toupper) +
		scale_y_continuous(limits = c(-2, 2)) +
		scale_color_manual(values = MetBrewer::met.brewer("VanGogh3")) +
		coord_flip(clip = "off") +
		guides(col = "none") +
		  theme_minimal(base_family = font_family) +
		  theme(
			plot.background = element_rect(color = NA, fill = "grey97"),
			panel.grid = element_blank(),
			panel.grid.major.x = element_line(linewidth = 0.1, color = "grey75"),
			plot.title.position = "plot",
			plot.subtitle = ggtext::element_textbox_simple(
			  margin = margin(t = 4, b = 16), size = 10),
			plot.caption = ggtext::element_textbox_simple(
			  margin = margin(t = 12), size = 7
			),
			plot.caption.position = "plot",
			axis.text.y = element_text(hjust = 0, margin = margin(r = -10)),
			plot.margin = margin(4, 4, 4, 4)
	  ) +
	stat_n_text(family = font_family, size = 3, color = "grey75") 
logger::log_info("Complete: depmap_plt")				
#### ---- end


## @knitr pea_volcano
pea_volcano_plt <-
	blood_ls$blood_pea_disease_de |>
	dplyr::filter(`ENSG ID` %in% target_ensembl_ids$GENEID) |>
	dplyr::filter(`p-value adjusted` <= 0.05) |>
	dplyr::arrange(`p-value adjusted`) |>
	dplyr::select(-Gene, -`ENSG ID`) |>
	dplyr::mutate(log10p = -log10(`p-value adjusted`)) |>
	dplyr::mutate(lbl = ifelse(grepl("cancer|oma", Disease, ignore.case = TRUE), BTX_COL$pink_out, BTX_COL$blue_d_out)) |>
	ggplot(aes(x = logFC, y = log10p, size = log10p, label = Disease)) +
		geom_point(aes(fill = lbl, color = lbl), shape = 21, alpha = 0.5) +
		geom_vline(xintercept = 0, color = BTX_COL$pink_out, lty = 2) +		
		ggrepel::geom_text_repel(size = 3.5) +
		theme(
			plot.background = element_rect(color = NA, fill = "grey97")
	) +
	theme(legend.position = "none") +
	scale_color_identity() + 
	scale_fill_identity() 
logger::log_info("Complete: pea_volcano")				
#### ---- end


## @knitr getLinkedOmicsData
# Will need to manually download CPTAC sets from:


# data('linkedOmics.data')
# linkedOmics.data |>
    # dplyr::select(c("project","OMICS Dataset","Unit","Samples","Atrributes","url")) |>
	# dplyr::filter(grepl("Proteome", `OMICS Dataset`))
	
# TCGA_COAD_protein <- 
	# TCGAbiolinks::getLinkedOmicsData(
		# project = "TCGA-OV",
		# dataset = "RPPA (Gene level)"
	# )

#### ---- end


## @knitr cancer_cptac

cancer_cptac_plt <-
	tumour_ls$cancer_cptac |>
		dplyr::filter(`Gene` %in% target_ensembl_ids$GENEID) |>
		# dplyr::select(-Gene, -`ENSG ID`) |>
		dplyr::mutate(log10p = -log10(`p-value adjusted`)) |>
		ggplot(aes(x = logFC, y = log10p, size = log10p, label = Cancer)) +
		geom_point(fill = BTX_COL$blue_d_out, color = BTX_COL$blue_d_out, shape = 21, alpha = 0.5) +
		geom_vline(xintercept = 0, color = BTX_COL$pink_out, lty = 2) +		
		ggrepel::geom_text_repel(size = 3.5) +
		theme(
			plot.background = element_rect(color = NA, fill = "grey97")
	) +
	theme(legend.position = "none") 

switch_cancer_cptac <- 
		tumour_ls$cancer_cptac |>
		dplyr::filter(`Gene` %in% target_ensembl_ids$GENEID) |>
		nrow() > 0
logger::log_info("Complete: cancer_cptac")				

## ---- end

## @knitr ihc_cancer_tb
# tumour_ls$cancer_cancers |> View() # Organ buckets
# replicates Protein Expression Overview Barplot

lvls <-
	tumour_ls$cancer_data |>
		dplyr::filter(Gene %in% target_ensembl_ids$GENEID) |>
		dplyr::filter(Cancer != "lymphoma") |> 	
		dplyr::arrange(High) |>
		dplyr::select(Cancer) |>
		unlist()

tmp_ihc <-
	tumour_ls$cancer_data |>
		dplyr::filter(Gene %in% target_ensembl_ids$GENEID) |>
		dplyr::filter(Cancer != "lymphoma") |> 
		tidyr::pivot_longer(cols = `High`:`Not detected`, names_to = "Level", values_to = "count") |>
		dplyr::mutate(Level = factor(Level, levels = c("Not detected", "Low", "Medium", "High"))) |>
		dplyr::mutate(Cancer = factor(Cancer, levels = lvls)) 

switch_ihc_cancer_tb <- nrow(tmp_ihc) > 0
		
prot_exp_c_overview <- 
	tmp_ihc |>
		ggplot(aes(y = Cancer, x = count, fill = Level)) +
		geom_bar(position = "stack", stat = "identity") + 
		paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
		theme(
			plot.background = element_rect(color = NA, fill = "grey97")
		)
logger::log_info("Complete: ihc_cancer_tb")				

## ---- end




