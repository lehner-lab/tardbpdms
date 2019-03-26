
#' tardbpdms_agg_tools_mutant_effects
#'
#' Calculate single and double mutant effects from aggregation tool predictions (using single mutants).
#'
#' @param toxicity_dt data.table with single and double mutant codes (required)
#' @param outpath output path for plots and saved objects (required)
#' @param aggtool_results_file path to aggregation tool results file (required)
#' @param aggscale_file path to aggregation scale file (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return A data.table with single and double mutant effects from aggregation tool predictions
#' @export
#' @import data.table
tardbpdms_agg_tools_mutant_effects <- function(
  toxicity_dt,
  outpath,
  aggtool_results_file,
  aggscale_file,
  colour_scheme,
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

	# ### Aggregation scales correlation with single mutants in hotspot
	# ###########################

	# #Single AA mutants only (no STOPs)
	# singles_dt_aggscale <- copy(singles_dt)
	# #Fitness hotspot positions
	# mean_toxicity <- singles_dt_aggscale[STOP==F,mean(abs(toxicity))]
	# Pos_abs_hotspot <- singles_dt_aggscale[Nmut_aa==1 & !STOP,.(hotspot = mean(abs(toxicity))>mean_toxicity),by=Pos_abs][hotspot==T,Pos_abs]
	# singles_dt_aggscale[, hotspot := as.numeric(Pos_abs %in% Pos_abs_hotspot)*2]

	# #Single mutant AA scales
	# singles_dt_aggscale <- tardbpdms__aa_properties_singles_loadings(
	# 	input_dt = singles_dt_aggscale,
 #  	aa_properties_file = aggscale_file)

 #  #AA properties PCA
 #  exp_pca_list <- tardbpdms__aa_properties_pca(
 #  	aa_properties_file = aggscale_file, 
 #  	return_evidences = T)
 #  agg_evidences <- exp_pca_list[["evidences"]]

	# #Single mutants in hotspot aggregation scales
	# singles_dt_agg <- singles_dt_aggscale[hotspot == 2,.SD,,.SDcols = names(singles_dt_aggscale)[grep("toxicity$|_score$", names(singles_dt_aggscale))]]
	# names(singles_dt_agg) <- gsub("_score$", "", names(singles_dt_agg))

	# #Barplot of R-squared for all scales
	# plot_df <- data.frame(
	# 	r_squared = cor(singles_dt_agg)[1,-1]^2,
	# 	scale = paste0(unlist(names(singles_dt_agg)[-1]), ":", unlist(agg_evidences[names(singles_dt_agg)[-1]])))
	# plot_df[,"scale"] <- factor(plot_df[,"scale"], levels = plot_df[order(plot_df[,"r_squared"], decreasing = T),"scale"])
	# d <- ggplot2::ggplot(plot_df, ggplot2::aes(scale, r_squared)) +
	#   ggplot2::geom_bar(stat = "identity") +
	#   ggplot2::theme_bw() +
 #    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
	# ggplot2::ggsave(file=file.path(outpath, 'aggregation_scales_correlations.pdf'), width=5, height=10)

	# ### Scatterplot of Tartaglia et al. predicted aggregation rate vs toxicity for single mutants in hotspot
	# ###########################

	# #Required equation values
	# eq_vals_list <- list(
	# 	"beta" = 1-c(0.35, 0.34, 0.72, 0.35, 0.4, 0.34, 0.37, 0.3, 0.06, 0.11, 0.6, 0.25, 0.47, 0.24, 0.26, 0.13, 0.13, 0.1, 0.32, 0),
	# 	"ASAa" = c(89, 119, 48, 61, 44, 53, 102, 44, 74, 144, 47, 48, 67, 195, 117, 175, 117, 140, 137, 0),
	# 	"ASAp" = c(107, 48, 58, 77, 69, 91, 49, 36, 28, 43, 0, 87, 0, 27, 43, 0, 0, 0, 0, 0),
	# 	"D" = c(5.92, 6.86, 2.64, 4.93, 4.03, 3.98, 2.86, 2.05, 2.06, 2.02, 0, 1.65, 0, 1.12, 0.66, 0, 0, 0, 0, 0))
	# eq_vals <- as.data.frame(do.call("cbind", eq_vals_list))
	# rownames(eq_vals) <- unlist(strsplit("RKDENQHSTYGCAWMFVILP", ""))
	# #Aromatic residues
	# eq_vals[,"A"] <- as.numeric(rownames(eq_vals) %in% unlist(strsplit("HFWYV", "")))
	# #Charged residues
	# eq_vals[,"C"] <- as.numeric(rownames(eq_vals) %in% unlist(strsplit("RK", ""))) - as.numeric(rownames(eq_vals) %in% unlist(strsplit("DE", "")))

	# #Charged/dipole side-chains
	# p_sc <- unlist(strsplit("KRDEHNQSTYCWM", ""))
	# a_sc <- rownames(eq_vals)[!rownames(eq_vals) %in% p_sc]

	# #Calculate
	# singles_dt_aggscale[WT_AA %in% a_sc & Mut %in% a_sc, phi_h := eq_vals[Mut,"ASAa"]/eq_vals[WT_AA,"ASAa"]]
	# singles_dt_aggscale[WT_AA %in% p_sc & Mut %in% p_sc, phi_h := eq_vals[WT_AA,"ASAp"]/eq_vals[Mut,"ASAp"]]
	# singles_dt_aggscale[WT_AA %in% a_sc & Mut %in% p_sc, phi_h := 1/eq_vals[Mut,"D"]]
	# singles_dt_aggscale[WT_AA %in% p_sc & Mut %in% a_sc, phi_h := eq_vals[WT_AA,"D"]]
	# singles_dt_aggscale[, phi_b := eq_vals[Mut,"beta"]/eq_vals[WT_AA,"beta"]]
	# singles_dt_aggscale[, phi_a := exp(eq_vals[Mut,"A"]-eq_vals[WT_AA,"A"])]
	# singles_dt_aggscale[, phi_c := exp(-0.5*(eq_vals[Mut,"C"]-eq_vals[WT_AA,"C"]))]
	# singles_dt_aggscale[, agg_rate := log(phi_h*phi_b*phi_a*phi_c)]

	# #Scatter plots - singles and doubles in hotspot
	# plot_df <- as.data.frame(singles_dt_aggscale[hotspot == 2 & !is.infinite(agg_rate),.SD,,.SDcols = c("toxicity", "agg_rate")])
 #  plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
 #  temp_cor <- as.data.table(plot_df)[,.(cor = round(cor(.SD[[1]], .SD[[2]]), 2), n = length(.SD[[1]]))]
	# d <- ggplot2::ggplot(plot_df, ggplot2::aes(agg_rate, toxicity)) +
 #    ggplot2::stat_binhex(bins=50) +
	#   ggplot2::xlab("Aggregation rate (Tartaglia et al., 2004)") +
	#   ggplot2::ylab("Toxicity") +
	#   ggplot2::theme_classic() +
	#   ggplot2::geom_smooth(method = "lm", linetype = 2, se = F, color = "black") +
	#   # ggplot2::coord_cartesian(ylim=c(-0.4, 0.4)) +
	#   ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 2, y = 0.3, size = 2) +
	#   ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 20))
	# ggplot2::ggsave(file=file.path(outpath, 'toxicity_vs_aggregation_rate_scatter_hotspot.pdf'), width=4.5, height=3)

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

