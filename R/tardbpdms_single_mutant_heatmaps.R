
#' tardbpdms_single_mutant_heatmaps
#'
#' Combine toxicity estimates from two different DMS experiments (regions).
#'
#' @param toxicity_dt data.table with single mutant toxicity values (required)
#' @param outpath output path for plots and saved objects (required)
#' @param disease_mut_file table of human disease mutations and classifications (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
tardbpdms_single_mutant_heatmaps <- function(
  toxicity_dt,
  disease_mut_file,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

	#Return previous results if analysis not executed
	if(!execute){
		return()
	}

  #Display status
  message(paste("\n\n*******", "running stage: tardbpdms_single_mutant_heatmaps", "*******\n\n"))

	#Create output directory
	tardbpdms__create_dir(tardbpdms_dir = outpath)

	#Subset to single AA mutants only
	dms_dt <- copy(toxicity_dt[Nmut_aa==1])

	#Toxicity hotspot positions
	mean_toxicity <- dms_dt[STOP==F,mean(abs(toxicity))]
	Pos_abs_hotspot <- dms_dt[STOP==F,.(hotspot = mean(abs(toxicity))>mean_toxicity),by=Pos_abs][hotspot==T,Pos_abs]
	dms_dt[, hotspot := as.numeric(Pos_abs %in% Pos_abs_hotspot)]

	#Disease mutations
	dis_mut <- read.table(disease_mut_file, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
	fALS_muts <- rownames(dis_mut)[!is.na(dis_mut$fALS)]
	sALS_muts <- rownames(dis_mut)[!is.na(dis_mut$sALS)]

	### Heatmap of transformed data
	###########################

	#Toxicity
	tardbpdms__single_mutant_heatmap(
		input_df = as.data.frame(dms_dt), 
		variable_name="toxicity", 
		output_file = file.path(outpath, 'heatmap_toxicity_center_scale.pdf'),
		x_annotation="sequence", 
		xaxis_angle=0, 
		na_colour="white", 
		na_text="-", 
		colour_low=colour_scheme[["shade 0"]][[3]], 
		colour_high=colour_scheme[["shade 0"]][[1]])

	#Toxicity - cluster positions
	tardbpdms__single_mutant_heatmap(
		input_df = as.data.frame(dms_dt), 
		variable_name="toxicity", 
		output_file = file.path(outpath, 'heatmap_toxicity_center_scale_colcluster.pdf'),
		x_annotation="both", 
		xaxis_angle=90, 
		na_colour="white", 
		na_text="-", 
		colour_low=colour_scheme[["shade 0"]][[3]], 
		colour_high=colour_scheme[["shade 0"]][[1]],
		cluster="column",
		width=24,
		height=8)

	#Toxicity (ALS mutations)
	tardbpdms__single_mutant_heatmap(
		input_df = as.data.frame(dms_dt), 
		variable_name="toxicity", 
		mut_dict=list("F" = fALS_muts, "S" = sALS_muts),
		output_file = file.path(outpath, 'heatmap_toxicity_center_scale_hmuts.pdf'),
		x_annotation="sequence", 
		xaxis_angle=0, 
		na_colour="white", 
		na_text="-", 
		colour_low=colour_scheme[["shade 0"]][[3]], 
		colour_high=colour_scheme[["shade 0"]][[1]])

	#Residual toxicity (after controlling for hydrophobicity and position)
	dms_dt[STOP==F, toxicity_resid := lm(toxicity ~ .*., data = .SD)[['residuals']], .SDcols = c('toxicity', 'PC1 (Hydrophobicity)', 'hotspot')]
	# dms_dt[STOP==F, toxicity_resid := loess(toxicity ~ .SD[[2]] + .SD[[3]] + .SD[[2]]*.SD[[3]], data = .SD)[['residuals']], .SDcols = c('toxicity', 'PC1 (Hydrophobicity)', 'hotspot')]
	dms_dt[, toxicity_resid_sig := p.adjust(2*pnorm(-abs(toxicity_resid/sigma)), method = "bonferroni")<0.01]
	tardbpdms__single_mutant_heatmap(
		input_df = as.data.frame(dms_dt), 
		variable_name="toxicity_resid", 
		mut_dict=list("*" = dms_dt[toxicity_resid_sig==T,mut_code]),
		output_file = file.path(outpath, 'heatmap_resid_toxicity_center_scale.pdf'),
		x_annotation="sequence", 
		xaxis_angle=0, 
		na_colour="white", 
		na_text="-", 
		colour_low=colour_scheme[["shade 0"]][[3]], 
		colour_high=colour_scheme[["shade 0"]][[1]])

}

