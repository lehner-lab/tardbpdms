
#' tardbpdms
#'
#' Main analysis script.
#'
#' @param startStage Start at a specified analysis stage (default:1)
#' @param stopStage Stop at a specified analysis stage (default:0 i.e. no stop condition)
#' @param base_dir Base directory for all output file (default:NB private CRG server path; change accordingly)
#' @param DMS2structure_path Path to DMS2structure repository (default:NB private CRG server path; change accordingly)
#' @param rerun_fitness re-run fitness analysis? (default:F)
#' @param rerun_epistasis re-run epistasis analysis? (default:F)
#' @param rerun_structure re-run structure analysis? (default:F)
#' @param numCores Number of available CPU cores (default:10)
#'
#' @return Nothing
#' @export
tardbpdms <- function(
  startStage=1,
  stopStage=0,
  base_dir = "/users/blehner/afaure/DMS/Results/tardbpdms_proj",
  DMS2structure_path = "/users/blehner/afaure/DMS/Code/DMS2structure",
  rerun_fitness = F,
  rerun_epistasis = F,
  rerun_structure = F,
  numCores = 10
  ){

	colour_scheme <- list(
		"shade 0" = list(
			"#F4270C",
			"#F4AD0C",
			"#1B38A6",
			"#09B636"),
		"shade 1" = list(
			"#FFB0A5",
			"#FFE4A5",
			"#9DACE3",
			"#97E9AD"),
		"shade 2" = list(
			"#FF6A56",
			"#FFCB56",
			"#4C63B7",
			"#43C766"),
		"shade 3" = list(
			"#A31300",
			"#A37200",
			"#0C226F",
			"#007A20"),
		"shade 4" = list(
			"#410800",
			"#412D00",
			"#020B2C",
			"#00300D"))

  #First and last analysis stages
  first_stage <- startStage
  last_stage <- stopStage

	#DiMSum counts to fitness
	stagenum <- 1
	#290 DMS library
	tardbpdms_dimsumcounts_to_fitness(
	  counts_path = file.path(base_dir, "misc", "DiMSum_counts", "290", "BB_TARDBP_290_2017-06-13_variant_data_merge.RData"),
	  outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_dimsumcounts_to_fitness_290", stagenum=stagenum, base_dir=base_dir),
	  generations_path = file.path(base_dir, "misc", "BB_TARDBP_290_2017-06-13_generations.txt"),
	  all_reps = c(1,3,4), #Exclude replicate 2
	  numCores = numCores,
	  execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum) & rerun_fitness))
	#332 DMS library
	tardbpdms_dimsumcounts_to_fitness(
	  counts_path = file.path(base_dir, "misc", "DiMSum_counts", "332", "BB_TARDBP_332_2018-04-06_Q25_variant_data_merge.RData"),
	  outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_dimsumcounts_to_fitness_332", stagenum=stagenum, base_dir=base_dir),
	  generations_path = file.path(base_dir, "misc", "BB_TARDBP_332_2018-04-06_Q25_generations.txt"),
	  all_reps = c(2,3,4), #Exclude replicate 1
	  numCores = numCores,
	  execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum) & rerun_fitness))

	#Quality control plots
	stagenum <- 2
	tardbpdms_quality_control(
		toxicity_list = list(
			"290" = file.path(base_dir, "001_tardbpdms_dimsumcounts_to_fitness_290"),
			"332" = file.path(base_dir, "001_tardbpdms_dimsumcounts_to_fitness_332")),
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_quality_control", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Combine toxicity estimates from 290 and 332 experiments (regions)
	stagenum <- 3
	toxicity_norm_dt <- tardbpdms_combine_toxicity(
		toxicity_list = list(
			"290" = file.path(base_dir, "002_tardbpdms_quality_control", "290"),
			"332" = file.path(base_dir, "002_tardbpdms_quality_control", "332")),
		growth_rate_file = file.path(base_dir, "misc", "single_mutants_growth_rates.txt"),
		disease_mut_file = file.path(base_dir, "misc", "tdp43_reported_mutations.txt"),
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_combine_toxicity", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Calculate single and double mutant effects from AA PCA (using single mutants)
	stagenum <- 4
	toxicity_aaprop_dt <- tardbpdms_aa_properties_mutant_effects(
		toxicity_dt = toxicity_norm_dt,
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_aa_properties_mutant_effects", stagenum=stagenum, base_dir=base_dir),
		aaprop_file = file.path(base_dir, "misc", "amino_acid_properties", "amino_acid_properties_annotated_supplementary.txt"),
		aaprop_file_selected = file.path(base_dir, "misc", "amino_acid_properties", "selected.amino_acid_properties.txt"),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Calculate single and double mutant effects from aggregation tool predictions (using single mutants)
	stagenum <- 5
	toxicity_aggtools_dt <- tardbpdms_agg_tools_mutant_effects(
		toxicity_dt = toxicity_aaprop_dt,
		aggtool_results_file = file.path(base_dir, "misc", "aggregation_tools_results.RData"),
		aggscale_file = file.path(base_dir, "misc", "amino_acid_properties", "aggregation_scales.txt"),
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_agg_tools_mutant_effects", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Single mutant heatmaps
	stagenum <- 6
	tardbpdms_single_mutant_heatmaps(
		toxicity_dt = toxicity_aggtools_dt,
		disease_mut_file = file.path(base_dir, "misc", "tdp43_reported_mutations.txt"),
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_single_mutant_heatmaps", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Human disease mutations
	stagenum <- 7
	tardbpdms_human_disease_mutations(
		toxicity_dt = toxicity_aggtools_dt,
		missense_AF_file = file.path(base_dir, "misc", "tdp43_gnomAD_r2.0.2_missense_allelefreqs.txt"),
		disease_mut_file = file.path(base_dir, "misc", "tdp43_reported_mutations.txt"),
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_human_disease_mutations", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Dot plots showing explained variance of models to predict variant toxicity.
	stagenum <- 8
	tardbpdms_toxicity_model_summary(
		toxicity_dt = toxicity_aggtools_dt,
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_toxicity_model_summary", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Violin plots showing (residual) toxicity versus number of introduced AAs of various properties
	stagenum <- 9
	tardbpdms_num_introduced_aa_violins(
		toxicity_dt = toxicity_aggtools_dt,
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_num_introduced_aa_violins", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Hydrophobicity of WT and toxicity hotspot line plot
	stagenum <- 10
	tardbpdms_wt_hydrophobicity(
		toxicity_dt = toxicity_aggtools_dt,
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_wt_hydrophobicity", stagenum=stagenum, base_dir=base_dir),
		aaprop_file = file.path(base_dir, "misc", "amino_acid_properties", "amino_acid_properties_annotated_supplementary.txt"),
		aaprop_file_selected = file.path(base_dir, "misc", "amino_acid_properties", "selected.amino_acid_properties.txt"),
		wtfasta_path = file.path(base_dir, "misc", "TARDBP.fa"),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Epistasis analysis
	stagenum <- 11
	tardbpdms_epistasis_analysis(
		toxicity_path = file.path(base_dir, "003_tardbpdms_combine_toxicity", "combine_toxicity_values.RData"),
		miscpath = file.path(base_dir, "misc"),
		DMS2structure_path = DMS2structure_path,
		numCores = numCores,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum) & rerun_epistasis))

	#Secondary structure predictions
	stagenum <- 12
	tardbpdms_secondary_structure_predictions(
		toxicity_dt = toxicity_aggtools_dt,
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_secondary_structure_predictions", stagenum=stagenum, base_dir=base_dir),
		miscpath = file.path(base_dir, "misc"),
		DMS2structure_path = DMS2structure_path,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)),
		rerun_structure = rerun_structure)

	#Guenther structure propensities
	stagenum <- 13
	tardbpdms_guenther_structure_propensities(
		result_dir = file.path(base_dir, "misc", "misc_guenther_structures"),
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_guenther_structure_propensities", stagenum=stagenum, base_dir=base_dir),
		miscpath = file.path(base_dir, "misc"),
		DMS2structure_path = DMS2structure_path,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)),
		rerun_structure = rerun_structure)

	#PWI heatmaps
	stagenum <- 14
	tardbpdms_PWI_heatmaps(
		toxicity_dt = toxicity_aggtools_dt,
		PWI_dir_list = list(
			"290" = file.path(base_dir, "misc", "misc_epistasis_analysis/doubles_cond_290/processed_data/"),
			"332" = file.path(base_dir, "misc", "misc_epistasis_analysis/doubles_cond_332/processed_data/")),
		outpath = tardbpdms__format_dir(dir_suffix="_tardbpdms_PWI_heatmaps", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))
}
