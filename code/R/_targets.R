# _targets.R file

## To execute this in an R session:
# tar_make(callr_function = NULL)

library(targets)
library(tarchetypes)
# Source functions from the 'function_lib' directory
tar_source(
	files = "function_lib",
	envir = targets::tar_option_get("envir"),
	change_directory = FALSE
)
# Set global package dependencies for targets
tar_option_set(
	packages = c("quarto", "logger", "readr", "readxl", "janitor", "dplyr", "tidyr", "checkmate", "ensembldb", "EnsDb.Hsapiens.v86", "org.Hs.eg.db", "AnnotationDbi", "UCSCXenaTools", "ggplot2", "ggrepel", "paletteer", "plotly", "pander", "iheatmapr", "ppcor", "corrr", "ggcorrplot", "gtExtras", "ComplexHeatmap"),
	memory = "transient",
	garbage_collection = TRUE,
	workspace_on_error = TRUE
)


list(
	core_data_ls,
	hpa_healthy_ihc_ls,
	protein_healthy_abundance_ls,
	hpa_gtex_healthy_quant_ls
)
