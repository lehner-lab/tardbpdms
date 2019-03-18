
#' tardbpdms_agg_tools_mutant_effects
#'
#' Calculate single and double mutant effects from aggregation tool predictions (using single mutants).
#'
#' @param toxicity_dt data.table with single and double mutant codes (required)
#' @param outpath output path for plots and saved objects (required)
#' @param aggtool_results_file path to aggregation tool results file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return A data.table with single and double mutant effects from aggregation tool predictions
#' @export
#' @import data.table
tardbpdms_agg_tools_mutant_effects <- function(
  toxicity_dt,
  outpath,
  aggtool_results_file,
  execute = TRUE
  ){

	#Return previous results if analysis not executed
	if(!execute){
		load(file.path(outpath, "agg_tools_mutant_effects.RData"))
		return(dms_dt_aggtool)
	}

	#Display status
	message(paste("\n\n*******", "running stage: tardbpdms_agg_tools_mutant_effects", "*******\n\n"))

	#Create output directory
	tardbpdms__create_dir(tardbpdms_dir = outpath)

	### Load data
	###########################

	load(aggtool_results_file)
	#Remove WT
	agg_df <- result_list_final_logfc[rownames(result_list_final_logfc)!="wt",]
	#Add single mutant code
	rownames(agg_df) <- gsub("mut_", "", rownames(agg_df))

	### Single mutant effects on aggregation tool predictions
	###########################

	singles_dt <- copy(toxicity_dt[Nmut_aa==1 & !STOP,])
	#Add aggregation tool predictions
	singles_dt <- cbind(singles_dt, as.data.table(agg_df[singles_dt[,mut_code],]))

	### Double mutant effects on aggregation tool predictions
	###########################

	doubles_dt <- copy(toxicity_dt[Nmut_aa==2 & !STOP,])
	#Add aggregation tool predictions (sum of effects of singles)
	doubles_dt <- cbind(doubles_dt, as.data.table(agg_df[doubles_dt[, mut_code1],] + agg_df[doubles_dt[, mut_code2],]))

	### Merge 
	###########################

	dms_dt_aggtool <- rbind(
		copy(toxicity_dt[!(Nmut_aa==2 & !STOP) & !(Nmut_aa==1 & !STOP)]),
		singles_dt,
		doubles_dt,
		fill = T)

	#RData object
	save(dms_dt_aggtool, file = file.path(outpath, "agg_tools_mutant_effects.RData"))
	
	#Return normalised toxicity data.table
	return(dms_dt_aggtool)
}

