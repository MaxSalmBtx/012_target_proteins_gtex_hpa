
core_data_ls <- 
	list(
		tar_target(param_file, "target_params.json", format = "file"),
		tar_target(param_ls, jsonlite::read_json(param_file, simplifyVector = TRUE)),
		tar_target(references_ensgid, prepare_ensembl_ids(x = c(param_ls$references))),
		tar_target(current_targets_f, "../../input/Bicycle_Merged_Target_List_Reserved_targets_06_02_26_FINAL.xlsx", format = "file"),
		tar_target(current_targets, get_target_list(current_targets_f)),
		tar_target(all_targets, merge_targets(x = current_targets, y = references_ensgid))
	)

hpa_healthy_ihc_ls <-
	list(
		tar_target(normal_ihc_data, hpa_download(url_str = "normal_ihc_data.tsv.zip")),
		tar_target(hpa_filtered, get_hpa_filter(x = normal_ihc_data, y = all_targets, p = param_ls)),
		tar_target(report_ihc_healthy_t, class(hpa_filtered)[1] == "tibble"),		
		# Specific report for this output
	    tarchetypes::tar_quarto(
			report_ihc_healthy,
			path = "report_ihc_healthy.qmd",
			quiet = TRUE,
			execute_params = list(dummy_param = report_ihc_healthy_t)
		)		
	)

protein_healthy_abundance_ls <-
	list(
		tar_target(paxdb_tb, get_paxdb(x = all_targets, p = param_ls)),
		tar_target(report_paxdb_t, class(paxdb_tb)[1] == "tibble"),
	    tarchetypes::tar_quarto(
			report_paxdb_healthy,
			path = "report_paxdb_healthy.qmd",
			quiet = FALSE,
			execute_params = list(dummy_param = report_paxdb_t)
		)
	)



hpa_gtex_healthy_quant_ls <- 
	list(
		tar_target(hpa_quant_f, "../../input/hpa/E-PROT-29-query-results.tsv", format = "file"),
		tar_target(hpa_quant, get_hpa_quant(x = hpa_quant_f, y = all_targets, p = param_ls)),
		# GTEx section [1]
		tar_target(gtex_quant_f, "../../input/gtex/jiang/Table_S1_gene_info_at_protein_level.xlsx", format = "file"),
		tar_target(gtex_quant, get_gtex_quant(gtex_quant_f, y = all_targets, p = param_ls)),
		# GTEx section [2]
		tar_target(gtex_quant_fang, get_gtex_quant_fang(y = all_targets, p = param_ls)),
		# Specific report for this output
		tar_target(report_proteomics_hpa_gtex_t, class(hpa_quant)[1] == "tibble"),		
	    tarchetypes::tar_quarto(
			report_proteomics_hpa_gtex,
			path = "report_proteomics_hpa_gtex.qmd",
			quiet = FALSE,
			execute_params = list(dummy_param = report_proteomics_hpa_gtex_t)
		)
	)

