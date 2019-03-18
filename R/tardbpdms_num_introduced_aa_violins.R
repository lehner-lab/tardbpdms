
#' tardbpdms_num_introduced_aa_violins
#'
#' Violin plots showing (residual) toxicity versus number of introduced AAs of various properties.
#'
#' @param toxicity_dt data.table with single and double mutant toxicity values and their mutant effects from AA PCA (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
tardbpdms_num_introduced_aa_violins <- function(
  toxicity_dt,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

	#Return previous results if analysis not executed
	if(!execute){
		return()
	}

  #Display status
  message(paste("\n\n*******", "running stage: tardbpdms_num_introduced_aa_violins", "*******\n\n"))

	#Create output directory
	tardbpdms__create_dir(tardbpdms_dir = outpath)

	### Setup
	###########################

	#Single AA mutants only (no STOPs)
	singles_dt <- copy(toxicity_dt[Nmut_aa==1 & !STOP])
	#Double AA mutants only (no STOPs)
	doubles_dt <- copy(toxicity_dt[Nmut_aa==2 & !STOP])
	doubles_dt[, toxicity := toxicity_cond]

	#Type of introduced AAs
	AA_classes <- list(
		hydrophobic = unlist(strsplit("ACILMFWYV", "")),
		aromatic = unlist(strsplit("HFWYV", "")),
		charged = unlist(strsplit("RDEK", "")))

	#Classify singles and doubles according to number of introducted AAs
	suppressWarnings(for(i in names(AA_classes)){
	  singles_dt[, (i) := 0,,.SDcols=i]
	  singles_dt[WT_AA %in% AA_classes[[i]] & !Mut %in% AA_classes[[i]], (i) := -1,,.SDcols=i]
	  singles_dt[!WT_AA %in% AA_classes[[i]] & Mut %in% AA_classes[[i]], (i) := 1,,.SDcols=i]
	  i1 <- paste0(i, 1)
	  doubles_dt[, (i1) := 0,,.SDcols=i1]
	  doubles_dt[WT_AA1 %in% AA_classes[[i]] & !Mut1 %in% AA_classes[[i]], (i1) := -1,,.SDcols=i1]
	  doubles_dt[!WT_AA1 %in% AA_classes[[i]] & Mut1 %in% AA_classes[[i]], (i1) := 1,,.SDcols=i1]
	  i2 <- paste0(i, 2)
	  doubles_dt[, (i2) := 0,,.SDcols=i2]
	  doubles_dt[WT_AA2 %in% AA_classes[[i]] & !Mut2 %in% AA_classes[[i]], (i2) := -1,,.SDcols=i2]
	  doubles_dt[!WT_AA2 %in% AA_classes[[i]] & Mut2 %in% AA_classes[[i]], (i2) := 1,,.SDcols=i2]
	  doubles_dt[, (i) := .SD[[1]]+.SD[[2]],,.SDcols=c(i1, i2)]
	})

	#Type of propensity
	propensity_classes <- list(
		helix_propensity = "AGADIR",
		aggregation_propensity = "Zyggregator")

	#Classify singles and doubles according to propensity change
	suppressWarnings(for(i in names(propensity_classes)){
		j <- propensity_classes[[i]]
		#Singles - direction
	  singles_dt[, (i) := "no change",,.SDcols=i]
	  singles_dt[j>0, (i) := "increase",,.SDcols=c(i, j)]
	  singles_dt[j<0, (i) := "decrease",,.SDcols=c(i, j)]
	  #Doubles - direction
	  doubles_dt[, (i) := "no change",,.SDcols=i]
	  doubles_dt[eval(as.name(j))>0, (i) := "increase",,.SDcols=c(i, j)]
	  doubles_dt[eval(as.name(j))<0, (i) := "decrease",,.SDcols=c(i, j)]
	  #Doubles - combination
	  k <- paste0(i, "_combination")
	  singles_df <- as.data.frame(singles_dt)
	  rownames(singles_df) <- singles_dt[,mut_code]
	  doubles_dt1 <- as.data.table(singles_df[doubles_dt[,mut_code1],])
	  doubles_dt2 <- as.data.table(singles_df[doubles_dt[,mut_code2],])
	  doubles_dt[, (k) := "both unchanged",,.SDcols=k]
	  doubles_dt[doubles_dt1[,.SD[[1]],.SDcols=j]>0 & doubles_dt2[,.SD[[1]],.SDcols=j]>0, (k) := "both increase",,.SDcols=c(i, j, k)]
	  doubles_dt[doubles_dt1[,.SD[[1]],.SDcols=j]<0 & doubles_dt2[,.SD[[1]],.SDcols=j]<0, (k) := "both decrease",,.SDcols=c(i, j, k)]
	  doubles_dt[(doubles_dt1[,.SD[[1]],.SDcols=j]*doubles_dt2[,.SD[[1]],.SDcols=j])<0, (k) := "increase/decrease",,.SDcols=c(i, j, k)]
	})

	#Toxicity hotspot positions
	mean_toxicity <- singles_dt[STOP==F,mean(abs(toxicity))]
	Pos_abs_hotspot <- singles_dt[Nmut_aa==1 & !STOP,.(hotspot = mean(abs(toxicity))>mean_toxicity),by=Pos_abs][hotspot==T,Pos_abs]
	singles_dt[, hotspot := as.numeric(Pos_abs %in% Pos_abs_hotspot)*2]
	doubles_dt[, hotspot := as.numeric(Pos_abs1 %in% Pos_abs_hotspot) + as.numeric(Pos_abs2 %in% Pos_abs_hotspot)]

	#Residual toxicity controlling for hydrophobicity and position (and interactions between these)
	tox_dt <- rbind(singles_dt, doubles_dt, fill = T)
	tox_dt[, summary(lm(toxicity ~ .*., data = .SD))[["r.squared"]], .SDcols = c('toxicity', 'PC1 (Hydrophobicity)', 'hotspot')]
	tox_dt[, toxicity_resid := lm(toxicity ~ .*., data = .SD)[['residuals']], .SDcols = c('toxicity', 'PC1 (Hydrophobicity)', 'hotspot')]
	tox_dt[hotspot==max(hotspot), toxicity_resid_hotspot := lm(toxicity ~ ., data = .SD)[['residuals']], .SDcols = c('toxicity', 'PC1 (Hydrophobicity)')]
	
	#Plot colours
	temp_col <- c(
		colour_scheme[["shade 0"]][[3]], 
		colour_scheme[["shade 1"]][[3]],
		"white",
		colour_scheme[["shade 1"]][[1]],
		colour_scheme[["shade 0"]][[1]])

	### Violin plots - hydrophobic, aromatic, charged
	###########################

	#Violin plots - toxicity before and after controlling for hydrophobicity
	feat_dt <- tox_dt[,.(toxicity, toxicity_resid, hydrophobic, aromatic, charged)]
	plot_dt <- rbind(
		feat_dt[,.(toxicity, hydrophobic, aromatic, charged, control="none")],
		feat_dt[,.(toxicity = toxicity_resid, hydrophobic, aromatic, charged, control="hydrophobicity")])
	plot_df <- reshape2::melt(as.data.frame(plot_dt), id = c("toxicity", "control"))
	plot_df[,"value"] <- factor(plot_df[,"value"], levels = -2:2)
	plot_df[,"control"] <- factor(plot_df[,"control"], levels = c("none", "hydrophobicity"))
	names(temp_col) <- 2:-2
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(value, toxicity, fill = value)) +
	  ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
	  ggplot2::xlab("Number of introduced amino acids") +
	  ggplot2::ylab("Toxicity") +
	  ggplot2::theme_bw() +
	  ggplot2::facet_grid(control ~ variable) +
	  ggplot2::scale_fill_manual(values = temp_col) +
	  ggplot2::geom_hline(yintercept = 0, linetype = 2)
	ggplot2::ggsave(file=file.path(outpath, 'toxicity_vs_num_hydrophobicity_aromaticity_charge_before_after_control.pdf'), width=8, height=5.5)

	# #Predictive power of aromaticity after controlling for hydrophobicity
	# tox_dt[,cor(toxicity, aromatic)^2]
	# tox_dt[,cor(toxicity_resid, aromatic)^2]

	#Violin plots - toxicity before and after controlling for hydrophobicity
	feat_dt <- tox_dt[hotspot==max(hotspot),.(toxicity, toxicity_resid_hotspot, hydrophobic, aromatic, charged)]
	plot_dt <- rbind(
		feat_dt[,.(toxicity, hydrophobic, aromatic, charged, control="none")],
		feat_dt[,.(toxicity = toxicity_resid_hotspot, hydrophobic, aromatic, charged, control="hydrophobicity")])
	plot_df <- reshape2::melt(as.data.frame(plot_dt), id = c("toxicity", "control"))
	plot_df[,"value"] <- factor(plot_df[,"value"], levels = -2:2)
	plot_df[,"control"] <- factor(plot_df[,"control"], levels = c("none", "hydrophobicity"))
	names(temp_col) <- 2:-2
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(value, toxicity, fill = value)) +
	  ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
	  ggplot2::xlab("Number of introduced amino acids") +
	  ggplot2::ylab("Toxicity") +
	  ggplot2::theme_bw() +
	  ggplot2::facet_grid(control ~ variable) +
	  ggplot2::scale_fill_manual(values = temp_col) +
	  ggplot2::geom_hline(yintercept = 0, linetype = 2)
	ggplot2::ggsave(file=file.path(outpath, 'toxicity_vs_num_hydrophobicity_aromaticity_charge_before_after_control_hotspot.pdf'), width=8, height=5.5)

	### Violin plots - helix and aggregation propensity combinations
	###########################

	#Violin plots - toxicity before and after controlling for hydrophobicity
	feat_dt <- tox_dt[Nmut_aa==2,.(toxicity, toxicity_resid, hotspot, helix_propensity_combination, aggregation_propensity_combination)]
	plot_dt <- rbind(
		feat_dt[,.(toxicity, hotspot, helix_propensity_combination, aggregation_propensity_combination, control="none")],
		feat_dt[,.(toxicity = toxicity_resid, hotspot, helix_propensity_combination, aggregation_propensity_combination, control="hydrophobicity")])
	plot_df <- reshape2::melt(as.data.frame(plot_dt), id = c("toxicity", "hotspot", "control"))
	#Remove underscores
	plot_df[,"variable"] <- gsub("_", " ", plot_df[,"variable"])
	#Remove doubles where both unchanged
	plot_df <- plot_df[plot_df[,"value"]!="both unchanged",]
	plot_df[,"value"][plot_df[,"hotspot"]==max(plot_df[,"hotspot"]) & plot_df[,"value"]!="increase/decrease"] <- paste0(plot_df[,"value"][plot_df[,"hotspot"]==max(plot_df[,"hotspot"]) & plot_df[,"value"]!="increase/decrease"], ' (hotspot)')
	plot_df[,"value"] <- factor(plot_df[,"value"], levels = c("both decrease (hotspot)", "both decrease", "increase/decrease", "both increase", "both increase (hotspot)"))
	names(temp_col) <- rev(c("both decrease (hotspot)", "both decrease", "increase/decrease", "both increase", "both increase (hotspot)"))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(value, toxicity, fill = value)) +
	  ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
	  ggplot2::xlab("Predicted effect of introduced amino acids") +
	  ggplot2::ylab("Toxicity") +
	  ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
	  ggplot2::facet_grid(control ~ variable) +
	  ggplot2::scale_fill_manual(values = temp_col) +
	  ggplot2::geom_hline(yintercept = 0, linetype = 2)
	ggplot2::ggsave(file=file.path(outpath, 'toxicity_vs_predicted_helix_aggregation_propensity_combination_before_after_control.pdf'), width=6, height=5)

	### Scatter plots - hydrophobicity
	###########################

	#Scatter plots - all
	feat_dt <- tox_dt[,.SD,,.SDcols = c("Nmut_aa", "toxicity", "PC1 (Hydrophobicity)", "hotspot")]
	plot_df <- as.data.frame(feat_dt)
	colnames(plot_df)[3] <- "hydrophobicity"
  plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
  temp_cor <- plyr::ddply(plot_df, c("hotspot", "Nmut_aa"), plyr::summarize, cor = round(cor(hydrophobicity, toxicity), 2), n = length(hydrophobicity))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(hydrophobicity, toxicity)) +
    ggplot2::stat_binhex(bins=50) +
	  ggplot2::xlab("Hydrophobicity") +
	  ggplot2::ylab("Toxicity") +
	  ggplot2::theme_bw() +
	  ggplot2::facet_grid(hotspot ~ Nmut_aa) +
	  ggplot2::geom_smooth(method = "lm", linetype = 2, se = F, color = "black") +
	  # ggplot2::coord_cartesian(ylim=c(-0.4, 0.4)) +
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 30, y = -0.3, size = 2) +
	  ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 20))
	ggplot2::ggsave(file=file.path(outpath, 'toxicity_vs_hydrophobicity_scatter_all.pdf'), width=5, height=5)

	#Scatter plots - singles and doubles in hotspot
	feat_dt <- tox_dt[hotspot==max(hotspot),.SD,,.SDcols = c("toxicity", "PC1 (Hydrophobicity)")]
	plot_df <- as.data.frame(feat_dt)
	colnames(plot_df)[2] <- "hydrophobicity"
  plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
  temp_cor <- feat_dt[,.(cor = round(cor(.SD[[1]], .SD[[2]]), 2), n = length(.SD[[1]]))]
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(hydrophobicity, toxicity)) +
    ggplot2::stat_binhex(bins=50) +
	  ggplot2::xlab("Hydrophobicity") +
	  ggplot2::ylab("Toxicity") +
	  ggplot2::theme_classic() +
	  ggplot2::geom_smooth(method = "lm", linetype = 2, se = F, color = "black") +
	  # ggplot2::coord_cartesian(ylim=c(-0.4, 0.4)) +
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 30, y = -0.3, size = 2) +
	  ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 20))
	ggplot2::ggsave(file=file.path(outpath, 'toxicity_vs_hydrophobicity_scatter_hotspot.pdf'), width=4.5, height=3)

	### Scatter plots - aggregation propensity
	###########################

	#Scatter plots - all
	feat_dt <- tox_dt[,.SD,,.SDcols = c("Nmut_aa", "toxicity", "Zyggregator", "hotspot")]
	plot_df <- as.data.frame(feat_dt)
  plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
  temp_cor <- plyr::ddply(plot_df, c("hotspot", "Nmut_aa"), plyr::summarize, cor = round(cor(Zyggregator, toxicity), 2), n = length(Zyggregator))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(Zyggregator, toxicity)) +
    ggplot2::stat_binhex(bins=50) +
	  ggplot2::xlab("Zyggregator") +
	  ggplot2::ylab("Toxicity") +
	  ggplot2::theme_bw() +
	  ggplot2::facet_grid(hotspot ~ Nmut_aa) +
	  ggplot2::geom_smooth(method = "lm", linetype = 2, se = F, color = "black") +
	  # ggplot2::coord_cartesian(ylim=c(-0.4, 0.4)) +
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 0.25, y = -0.3, size = 2) +
	  ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 20))
	ggplot2::ggsave(file=file.path(outpath, 'toxicity_vs_zyggregator_scatter_all.pdf'), width=5, height=5)

	#Scatter plots - singles and doubles in hotspot
	feat_dt <- tox_dt[,.(toxicity, toxicity_resid, Zyggregator)]
	plot_dt <- rbind(
		feat_dt[,.(toxicity, Zyggregator, control="none")],
		feat_dt[,.(toxicity = toxicity_resid, Zyggregator, control="hydrophobicity")])
	plot_df <- reshape2::melt(as.data.frame(plot_dt), id = c("toxicity", "Zyggregator", "control"))
	plot_df[,"control"] <- factor(plot_df[,"control"], levels = c("none", "hydrophobicity"))
  plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
  temp_cor <- plyr::ddply(plot_df, c("control"), plyr::summarize, cor = round(cor(Zyggregator, toxicity), 2), n = length(Zyggregator))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(Zyggregator, toxicity)) +
    ggplot2::stat_binhex(bins=50) +
	  ggplot2::xlab("Zyggregator") +
	  ggplot2::ylab("Toxicity") +
	  ggplot2::theme_classic() +
	  ggplot2::facet_grid(control ~ .) +
	  ggplot2::geom_smooth(method = "lm", linetype = 2, se = F, color = "black") +
	  # ggplot2::coord_cartesian(ylim=c(-0.4, 0.4)) +
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 0.25, y = -0.3, size = 2) +
	  ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 20))
	ggplot2::ggsave(file=file.path(outpath, 'toxicity_vs_zyggregator_scatter_hotspot_before_after_control.pdf'), width=5, height=5.5)

	#Scatter plots - singles and doubles in hotspot
	feat_dt <- tox_dt[hotspot==max(hotspot),.(toxicity, toxicity_resid_hotspot, Zyggregator)]
	plot_dt <- rbind(
		feat_dt[,.(toxicity, Zyggregator, control="none")],
		feat_dt[,.(toxicity = toxicity_resid_hotspot, Zyggregator, control="hydrophobicity")])
	plot_df <- reshape2::melt(as.data.frame(plot_dt), id = c("toxicity", "Zyggregator", "control"))
	plot_df[,"control"] <- factor(plot_df[,"control"], levels = c("none", "hydrophobicity"))
  plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
  temp_cor <- plyr::ddply(plot_df, c("control"), plyr::summarize, cor = round(cor(Zyggregator, toxicity), 2), n = length(Zyggregator))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(Zyggregator, toxicity)) +
    ggplot2::stat_binhex(bins=50) +
	  ggplot2::xlab("Zyggregator") +
	  ggplot2::ylab("Toxicity") +
	  ggplot2::theme_classic() +
	  ggplot2::facet_grid(control ~ .) +
	  ggplot2::geom_smooth(method = "lm", linetype = 2, se = F, color = "black") +
	  # ggplot2::coord_cartesian(ylim=c(-0.4, 0.4)) +
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 0.25, y = -0.3, size = 2) +
	  ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 20))
	ggplot2::ggsave(file=file.path(outpath, 'toxicity_vs_zyggregator_scatter_hotspot_before_after_control_hotspot.pdf'), width=5, height=5.5)

}

