
#' tardbpdms_secondary_structure_predictions
#'
#' Secondary structure predictions.
#'
#' @param toxicity_dt data.table with single and double mutant toxicity values (required)
#' @param outpath output path for plots and saved objects (required)
#' @param miscpath path to misc scripts and data directory (required)
#' @param DMS2structure_path Path to DMS2structure repository (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#' @param rerun_structure re-run structure analysis? (default:F)
#'
#' @return Nothing
#' @export
#' @import data.table
tardbpdms_secondary_structure_predictions <- function(
  toxicity_dt,
  outpath,
  miscpath,
  DMS2structure_path,
  colour_scheme,
  execute = TRUE,
  rerun_structure = F
  ){

	#Return previous results if analysis not executed
	if(!execute){
		return()
	}

  #Display status
  message(paste("\n\n*******", "running stage: tardbpdms_secondary_structure_predictions", "*******\n\n"))

	#Create output directory
	tardbpdms__create_dir(tardbpdms_dir = outpath)

	### Re-run structure analysis
	###########################

	if(rerun_structure){

	  #Display status
	  message(paste("\n\n*******", "re-running structure predictions (this might take a while)", "*******\n\n"))

	  #Create output directory
	  tardbpdms__create_dir(tardbpdms_dir = file.path(miscpath, "misc_secondary_structure"))

		#Run misc script on command-line
	  system(paste0(
	  	file.path(miscpath, "scripts", "tardbpdms__secondary_structure.R"),
	  	" -o ",
	  	file.path(miscpath, "misc_secondary_structure"),
	  	" -d ",
	  	DMS2structure_path,
	  	" -e ",
	  	file.path(miscpath, "misc_epistasis_analysis")))
	}

	### Setup
	###########################

	#Load secondary structure predictions
	load(file.path(miscpath, "misc_secondary_structure", "result_list.RData"))
	ss_dt <- rbind(
		result_list[["PWI_cond_290"]][["secondary_structure_score"]][["association_score"]][, region := "290"][, score_type := "association_score"],
		result_list[["PWI_cond_290"]][["secondary_structure_score"]][["posE_pcor"]][, region := "290"][, score_type := "posE_pcor"],
		result_list[["PWI_cond_290"]][["secondary_structure_score"]][["negE_pcor"]][, region := "290"][, score_type := "negE_pcor"],
		result_list[["PWI_cond_332"]][["secondary_structure_score"]][["association_score"]][, region := "332"][, score_type := "association_score"],
		result_list[["PWI_cond_332"]][["secondary_structure_score"]][["posE_pcor"]][, region := "332"][, score_type := "posE_pcor"],
		result_list[["PWI_cond_332"]][["secondary_structure_score"]][["negE_pcor"]][, region := "332"][, score_type := "negE_pcor"])
	ss_dt[, Pos_abs := as.numeric(region)+Pos1-1]

	#Single AA mutants only (no STOPs)
	singles_dt <- copy(toxicity_dt[Nmut_aa==1 & !STOP])
	#Double AA mutants only (no STOPs)
	doubles_dt <- copy(toxicity_dt[Nmut_aa==2 & !STOP])
	doubles_dt[, fitness := fitness_cond]

	#Fitness hotspot positions
	mean_fitness <- singles_dt[,mean(abs(fitness))]
	Pos_abs_hotspot <- singles_dt[Nmut_aa==1 & !STOP,.(hotspot = mean(abs(fitness))>mean_fitness),by=Pos_abs][hotspot==T,Pos_abs]

	#Mutant positions
	Pos_abs_all <- unique(singles_dt[,Pos_abs])

	#Structure coordinates
	NTD_coords <- c(1, 78)
	RRM1_coords <- c(101, 176)
	RRM2_coords <- c(191, 262)
	PRLD_coords <- c(274, 414)
	LARKS_coords <- c(312, 317)
	LARKS_coords_focus <- c(312, 314)
	helix_coords <- c(321, 330)
	helix_coords_focus <- c(326, 330)

	### Secondary structure significance - association score (both positive and negative epistasis)
	###########################

	#Line plot
	plot_dt <- copy(ss_dt)[score_type == "association_score",]
	plot_dt[, alpha_helix_pval := -log10(alpha_p_seed0)]
	plot_dt[, beta_strand_pval := -log10(beta_p_seed0)]
	plot_df <- reshape2::melt(as.data.frame(plot_dt[,.(Pos_abs, alpha_helix_pval, beta_strand_pval)]), id = 'Pos_abs')
	data_min <- min(plot_df[,'value']) + 0.2*min(plot_df[,'value'])
	data_max <- max(plot_df[,'value']) + 0.2*max(plot_df[,'value'])
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(Pos_abs, value, color=variable)) +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(Pos_abs_hotspot), xmax=max(Pos_abs_hotspot), ymin=data_min, ymax=data_max), fill = "lightgrey", linetype = 0, alpha = 0.1) +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(LARKS_coords), xmax=max(LARKS_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(helix_coords), xmax=max(helix_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	  ggplot2::geom_vline(xintercept = median(Pos_abs_all), linetype = 2) +
	  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2) +
	  ggplot2::geom_line() +
	  ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(4, 2)]))
  }
	ggplot2::ggsave(file=file.path(outpath, 'tdp43_ss_pred_association_score_all.pdf'), width=8, height=3)

	### Secondary structure significance - positive and negative association score
	###########################

	#Line plot
	plot_dt <- copy(ss_dt)[score_type %in% c("posE_pcor", "negE_pcor"),]
	plot_dt[, alpha_helix_pval := -log10(alpha_p_seed0)]
	plot_dt[, beta_strand_pval := -log10(beta_p_seed0)]
	plot_df <- reshape2::melt(as.data.frame(plot_dt[,.(Pos_abs, alpha_helix_pval, beta_strand_pval, score_type)]), id = c('Pos_abs', 'score_type'))
	data_min <- min(plot_df[,'value']) + 0.2*min(plot_df[,'value'])
	data_max <- max(plot_df[,'value']) + 0.2*max(plot_df[,'value'])
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(Pos_abs, value, color=variable)) +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(Pos_abs_hotspot), xmax=max(Pos_abs_hotspot), ymin=data_min, ymax=data_max), fill = "lightgrey", linetype = 0, alpha = 0.1) +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(LARKS_coords), xmax=max(LARKS_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(helix_coords), xmax=max(helix_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	  ggplot2::geom_vline(xintercept = median(Pos_abs_all), linetype = 2) +
	  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2) +
	  ggplot2::geom_line(ggplot2::aes(linetype = score_type)) +
	  ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(4, 2)]))
  }
	ggplot2::ggsave(file=file.path(outpath, 'tdp43_ss_pred_association_score_pos_neg.pdf'), width=8, height=3)

	#Dotplots
	#Adjust p-values
	plot_dt2 <- rbind(
		copy(plot_dt)[, ss_type := "beta_strand"][, value := beta_strand_pval],
		copy(plot_dt)[, ss_type := "alpha_helix"][, value := alpha_helix_pval])[, value_adjust := -log10(p.adjust(10^-value, method = "BH"))]
	#Subset to LARKS and helix positions
	plot_dt2[Pos_abs %in% seq(LARKS_coords[1], LARKS_coords[2]), feature := "LARKS"]
	plot_dt2[feature=="LARKS" & ss_type=="alpha_helix", feature_cat := "LARKS_alpha"]
	plot_dt2[feature=="LARKS" & ss_type=="beta_strand", feature_cat := "LARKS_beta"]
	plot_dt2[Pos_abs %in% seq(helix_coords[1], helix_coords[2]), feature := "helix"]
	plot_dt2[feature=="helix" & ss_type=="alpha_helix", feature_cat := "helix_alpha"]
	plot_dt2[feature=="helix" & ss_type=="beta_strand", feature_cat := "helix_beta"]
	plot_df <- as.data.frame(plot_dt2[feature=="LARKS" | feature=="helix",])
	plot_df[,"feature_cat"] <- factor(plot_df[,"feature_cat"], levels = c("LARKS_beta", "LARKS_alpha", "helix_beta", "helix_alpha"))
	plot_df[,"score_type"] <- factor(plot_df[,"score_type"], levels = c("posE_pcor", "negE_pcor"))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(feature_cat, value_adjust, fill=score_type)) +
	  ggplot2::geom_dotplot(binaxis = "y", stackdir = "center") +
	  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2) +
	  ggplot2::theme_classic() +
	  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(3,1)]))
  }
	suppressMessages(ggplot2::ggsave(file=file.path(outpath, 'tdp43_ss_pred_association_score_pos_neg_features_dotplot.pdf'), width=4, height=4))

	### Secondary structure significance - association score (both positive and negative epistasis) - deepcontact
	###########################

	# ss_dt <- rbind(
	# 	result_list[["PWI_cond_290_deepcontact"]][["secondary_structure_score"]][["association_score"]][, region := "290"][, score_type := "association_score"],
	# 	result_list[["PWI_cond_332_deepcontact"]][["secondary_structure_score"]][["association_score"]][, region := "332"][, score_type := "association_score"])
	# ss_dt[, Pos_abs := as.numeric(region)+Pos1-1]

	# #Line plot
	# plot_dt <- copy(ss_dt)[score_type == "association_score",]
	# plot_dt[, alpha_helix_pval := -log10(alpha_p_seed0)]
	# plot_dt[, beta_strand_pval := -log10(beta_p_seed0)]
	# plot_df <- reshape2::melt(as.data.frame(plot_dt[,.(Pos_abs, alpha_helix_pval, beta_strand_pval)]), id = 'Pos_abs')
	# data_min <- min(plot_df[,'value']) + 0.2*min(plot_df[,'value'])
	# data_max <- max(plot_df[,'value']) + 0.2*max(plot_df[,'value'])
	# d <- ggplot2::ggplot(plot_df, ggplot2::aes(Pos_abs, value, color=variable)) +
	#   ggplot2::geom_rect(ggplot2::aes(xmin=min(Pos_abs_hotspot), xmax=max(Pos_abs_hotspot), ymin=data_min, ymax=data_max), fill = "lightgrey", linetype = 0, alpha = 0.1) +
	#   ggplot2::geom_rect(ggplot2::aes(xmin=min(LARKS_coords), xmax=max(LARKS_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	#   ggplot2::geom_rect(ggplot2::aes(xmin=min(helix_coords), xmax=max(helix_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	#   ggplot2::geom_vline(xintercept = median(Pos_abs_all), linetype = 2) +
	#   ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2) +
	#   ggplot2::geom_line() +
	#   ggplot2::theme_classic()
 #  if(!is.null(colour_scheme)){
 #    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(4, 2)]))
 #  }
	# ggplot2::ggsave(file=file.path(outpath, 'tdp43_ss_pred_association_score_all_deepcontact.pdf'), width=8, height=3)


}








