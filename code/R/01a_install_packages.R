# author      = "Max Salm"
# copyright   = "MIT"
# license     = "NA"
# version     = "1.0.0"
# status      = "Development"
# Aim         = Simple script to install all necessary R packages

### Check web access first
if ( Sys.info()["sysname"] == "Linux" ) {
	cmd <- "wget -p https://www.google.com/"
	test_connection <- 
		system(cmd, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
	if (length(test_connection) == 0) {
		cat("Connection to internet verified.\n")
		unlink("www.google.com", recursive = TRUE)
	}else{
		stop("No connection to internet found.\n")
	}
}

### Helper functions


check_package <- function(
  package_names, 
  dest_dir, 
  fresh_install = FALSE, 
  repository = "CRAN"
){

	#' Check and install R packages from various repositories.
	#'
	#' This function checks if a set of R packages are installed in a specified
	#' local repository. If not, it installs them from CRAN, or BIOC.
	#' It also allows for a fresh install, removing previous versions if needed.
	#'
	#' @param package_names A character vector of package names to check and install.
	#' @param dest_dir The destination directory for the local repository.
	#' @param fresh_install A logical value indicating whether to remove previous
	#'   installations of the packages. Defaults to FALSE.
	#' @param repository A character string specifying the repository to use.
	#'   Must be one of "CRAN" or "BIOC". Defaults to "CRAN".
	#'
	#' @return None (invisibly). The function installs packages and loads them.
	#'
	#' @examples
	#' \dontrun{
	#' # Check and install packages from CRAN
	#' checkPackage(
	#'   package_names = c("ggplot2", "dplyr"),
	#'   dest_dir = "~/R/my_packages"
	#' )
	#'
	#' # Perform a fresh install from CRAN
	#' checkPackage(
	#'   package_names = c("MASS", "nlme"),
	#'   dest_dir = "~/R/my_packages",
	#'   fresh_install = TRUE,
	#'   repository = "CRAN"
	#' )
	#' }
	#'
	#' @import checkmate
	#' @export

	# Argument checks using checkmate
	checkmate::assertCharacter(package_names, min.len = 1, any.missing = FALSE)
	checkmate::assertDirectory(dest_dir, access = "w")
	checkmate::assertLogical(fresh_install, len = 1)
	checkmate::assertChoice(repository, choices = c("CRAN", "BIOC"))

	# Check if the local repository directory exists.
	# If it doesn't, create it and update the library paths.
    if (!dir.exists(dest_dir)){
        dir.create(dest_dir, recursive = TRUE, mode = "0777")
        .libPaths(dest_dir)  # Update library paths to include the new directory
    } 
        
	# Get a vector of installed packages in the local repository.
    installed_packages  <- installed.packages(lib.loc = dest_dir)[, "Package"]

	# If fresh_install is TRUE, remove previously installed packages.
    if (fresh_install) {
        cat("Removing previous installations.\n")
        packages_to_remove  <- 
			package_names[package_names %in% installed_packages]
        remove.packages(packages_to_remove, lib = dest_dir)
    }

	# Determine which packages need to be newly installed.
    new_packages <- package_names[!(package_names %in% installed_packages)]
    if (length(new_packages) > 0) {
        message("Installing:")
        message(paste0(new_packages, "\n"))
		# Install packages based on the specified repository.		
        if(repository == "CRAN"){
            install.packages(new_packages, 
                            dependencies = TRUE, 
                            repos = "https://cran.rstudio.com/", 
                            lib = dest_dir)      
        }else if(repository == "BIOC"){
		if (!require("BiocManager", quietly = TRUE))
			install.packages("BiocManager")
		BiocManager::install(version = "3.22") # compatible with R 4.4.0				
            BiocManager::install(
              pkgs = new_packages, 
              lib = dest_dir, 
              suppressUpdates=TRUE, 
              ask=FALSE
            )     
        }else{
            stop("checkPackage(): Select from CRAN/BIOC")
        }                
    }else{
        logger::log_info("All packages found.\n")
    }

	# Load all requested packages.   
	package_load_status <- 
		lapply(
			package_names, 
			function(pkg) {
				tryCatch(
					{
						library(pkg, character.only = TRUE, lib.loc = dest_dir)
						TRUE # Return TRUE if loading succeeds
					}, 
					error = function(e) {
						cat(paste("Failed to load ", pkg, ":", e$message, "\n"))
						FALSE # Return FALSE if loading fails
					}
				)
			}
		)

	if (any(!unlist(package_load_status))) {
		logger::log_warn("Some packages failed to load. Check output for details.")		
	}
}

#########################
### Base dependencies ###
#########################
if (Sys.getenv("BIOENV_IMAGE") == "") {
    .libPaths(c("/home/m.salm/R/personal_lib", .libPaths()))
}
if (Sys.getenv("USER") == "m.salm") {
  personal_repo <- "../../code/R/r_lib"
} else {
  personal_repo <- .libPaths()[1]
}

base_pkgs <- c("quarto", "RColorBrewer", "pander", "checkmate", "docstring", "logger", "testthat", "mockery", "readxl", "janitor", "stringdist", "targets", "tarchetypes", "assertions", "overlapping")
check_package(
	package_names = base_pkgs, 
    dest_dir = personal_repo, 
    fresh_install = FALSE, 
    repository = "CRAN"
)

####################
### Bioconductor ###
####################
bioc_pkgs <- c("EnsDb.Hsapiens.v86", "ComplexHeatmap", "curatedTCGAData", "TCGAbiolinks", "TCGAbiolinksGUI.data", "EMDomics", "UniProt.ws", "org.Hs.eg.db", "AnnotationDbi", "GSVA", "vulcan")
check_package(
	package_names = bioc_pkgs, 
    dest_dir = personal_repo, 
    fresh_install = FALSE, 
    repository = "BIOC"
)


############
### CRAN ###
############
mran_pkgs <- c("rmarkdown", "readr", "Hmisc", "pander", "R.utils", "ggrepel", "UCSCXenaShiny", "sparklyr", "sparklyr.nested", "httr", "httr2", "jsonlite", "DescTools", "corrr", "ggcorrplot", "ppcor", "RankAggreg", "TopKLists", "RobustRankAggreg", "mousetrap")
check_package(
	package_names = mran_pkgs, 
    dest_dir = personal_repo, 
    fresh_install = FALSE, 
    repository = "CRAN"
)

table_pkgs <- c("kableExtra", "pander", "gt", "gtExtras", "gtsummary")
check_package(
	package_names = table_pkgs, 
    dest_dir = personal_repo, 
    fresh_install = FALSE, 
    repository = "CRAN"
)
    
plot_pkgs <- 
	c("ggplot2", "viridis", "ggpubr",  "magick", "showtext", "gplots", "ggExtra", "ggdist", "MetBrewer", "EnvStats", "ggtext", "colorRamp2", "iheatmapr", "legendry", "tidyclust", "GGally", "gghighlight")
check_package(
	package_names = plot_pkgs, 
    dest_dir = personal_repo, 
    fresh_install = FALSE, 
    repository = "CRAN"
)


# remotes::install_github("ropensci/UCSCXenaTools") # Needed to get the most up-to-date version, bug fix
# remotes::install_github("jespermaag/gganatogram") # organ maps
# remotes::install_github("CT-Data-Haven/stylehaven") # Logo annotations
# remotes::install_github("alex-bio/UniProtExtractR") # Tidy access to UniProtKB
# remotes::install_github("hrbrmstr/ggchicklet") # stacked barplots
# remotes::install_github("steveneschrich/surfaceome") # surfaceome
# remotes::install_github("hughjonesd/ggmagnify")
