
#' tardbpdms_wt_hydrophobicity
#'
#' Hydrophobicity of WT and toxicity hotspot line plot.
#'
#' @param toxicity_dt data.table with single and double mutant toxicity values and their mutant effects from AA PCA (required)
#' @param outpath output path for plots and saved objects (required)
#' @param aaprop_file path to amino acid properties file (required)
#' @param aaprop_file_selected path to file with selected subset of identifiers
#' @param wtfasta_path path to WT FASTA sequence file (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
tardbpdms_wt_hydrophobicity <- function(
  toxicity_dt,
  outpath,
  aaprop_file,
  aaprop_file_selected,
  wtfasta_path,
  colour_scheme,
  execute = TRUE
  ){

	#Return previous results if analysis not executed
	if(!execute){
		return()
	}

  #Display status
  message(paste("\n\n*******", "running stage: tardbpdms_wt_hydrophobicity", "*******\n\n"))

	#Create output directory
	tardbpdms__create_dir(tardbpdms_dir = outpath)

	### Setup
	###########################

	#Single AA mutants only (no STOPs)
	singles_dt <- copy(toxicity_dt[Nmut_aa==1 & !STOP])

	#Fitness hotspot positions
	mean_toxicity <- singles_dt[,mean(abs(toxicity))]
	Pos_abs_hotspot <- singles_dt[Nmut_aa==1 & !STOP,.(hotspot = mean(abs(toxicity))>mean_toxicity),by=Pos_abs][hotspot==T,Pos_abs]

	#Mutant positions
	Pos_abs_all <- unique(singles_dt[,Pos_abs])

	#Structure coordinates
	NTD_coords <- c(1, 78)
	RRM1_coords <- c(101, 176)
	RRM2_coords <- c(191, 262)
	PRLD_coords <- c(274, 414)

	#Wildtype sequence
	wt_seq <- paste0(system(paste0("cat ", wtfasta_path), inter = T)[-1], collapse = "")

	### Amino acid properties data
	###########################

	#PCA
  exp_pca <- tardbpdms__aa_properties_pca(aa_properties_file = aaprop_file, selected_identifiers = unlist(fread(aaprop_file_selected, header = F)))
	pc_mat <- exp_pca$x

	#Add names and change signs of top PCs
	top_PCs <- 5
	top_PC_signs <- c(-1, -1, 1, -1, -1)
	top_PC_names <- c("Hydrophobicity", "Helix propensity", "Commonness", "Linker propensity", "Beta-sheet propensity")
	colnames(pc_mat)[grep("^PC", colnames(pc_mat))][1:5] <- paste0(colnames(pc_mat)[grep("^PC", colnames(pc_mat))][1:5], " (", top_PC_names, ")")
	pc_mat <- scale(pc_mat, scale=c(top_PC_signs, rep(1, dim(pc_mat)[1]-top_PCs)), center=F)

	### WT sequence in PC space
	###########################

	#WT loadings - all
	wt_pc <- as.data.frame(pc_mat[unlist(strsplit(wt_seq, "")),])
	rownames(wt_pc) <- 1:dim(wt_pc)[1]
	#Plot
	plot_df <- reshape2::melt(cbind(wt_pc, data.frame(position = 1:dim(wt_pc)[1])), id = "position")
	temp_min <- tapply(plot_df[,"value"], plot_df[,"variable"], min)
	temp_max <- tapply(plot_df[,"value"], plot_df[,"variable"], max)
	plot_df[,"data_min"] <- temp_min[plot_df[,"variable"]]+0.2*temp_min[plot_df[,"variable"]]
	plot_df[,"data_max"] <- temp_max[plot_df[,"variable"]]+0.2*temp_max[plot_df[,"variable"]]
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(position, value, color=variable)) +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(NTD_coords), xmax=max(NTD_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 2, color = "black") +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(RRM1_coords), xmax=max(RRM1_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 2, color = "black") +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(RRM2_coords), xmax=max(RRM2_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 2, color = "black") +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(Pos_abs_hotspot), xmax=max(Pos_abs_hotspot), ymin=data_min, ymax=data_max), fill = "lightgrey", linetype = 0, alpha = 0.1) +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(Pos_abs_all), xmax=max(Pos_abs_all), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(PRLD_coords), xmax=max(PRLD_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 2, color = "black") +
	  ggplot2::geom_vline(xintercept = median(Pos_abs_all)) +
	  ggplot2::geom_line() +
	  ggplot2::geom_smooth(span = 0.05, ggplot2::aes(fill = variable), method = 'loess', formula = 'y ~ x') +
	  ggplot2::theme_classic() + ggplot2::facet_grid(variable~., scales = "free_y")
	ggplot2::ggsave(file=file.path(outpath, 'tdp43_wt_pc_loadings.pdf'), width=14, height=40)

	#WT loadings - hydrophobicity
	singles_df <- singles_dt[,.(position = Pos_abs, value = abs(toxicity)*100, variable = "toxicity")]
	plot_df <- plot_df[grepl("PC1 ", plot_df[,"variable"]),]
	hydro_min <- min(plot_df[,"value"])+0.2*min(plot_df[,"value"])
	hydro_max <- max(plot_df[,"value"])+0.2*max(plot_df[,"value"])
	#Test correlation between wt hydrophobicity and mean absolute toxicity per position
	plot_df_test <- plot_df[plot_df[,"position"] %in% Pos_abs_all & grepl("PC1 ", plot_df[,"variable"]),]
	test_dt <- singles_dt[,.(mean_tox = mean(abs(toxicity))),by=Pos_abs]
	test_dt[, wt_hydrophobicity := plot_df_test[as.character(test_dt[,Pos_abs]),"value"]]
	tox_cor <- test_dt[,.(R = cor(wt_hydrophobicity, mean_tox), n = .N)]
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(position, value, color=variable)) +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(NTD_coords), xmax=max(NTD_coords), ymin=hydro_min, ymax=hydro_max), fill = NA, linetype = 2, color = "black") +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(RRM1_coords), xmax=max(RRM1_coords), ymin=hydro_min, ymax=hydro_max), fill = NA, linetype = 2, color = "black") +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(RRM2_coords), xmax=max(RRM2_coords), ymin=hydro_min, ymax=hydro_max), fill = NA, linetype = 2, color = "black") +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(Pos_abs_hotspot), xmax=max(Pos_abs_hotspot), ymin=hydro_min, ymax=hydro_max), fill = "lightgrey", linetype = 0, alpha = 0.1) +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(Pos_abs_all), xmax=max(Pos_abs_all), ymin=hydro_min, ymax=hydro_max), fill = NA, linetype = 1, color = "black") +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(PRLD_coords), xmax=max(PRLD_coords), ymin=hydro_min, ymax=hydro_max), fill = NA, linetype = 2, color = "black") +
	  ggplot2::geom_vline(xintercept = median(Pos_abs_all)) +
	  # ggplot2::geom_line(alpha = 0.5) +
	  ggplot2::geom_smooth(span = 0.05, ggplot2::aes(fill = variable), method = 'loess', formula = 'y ~ x') +
	  ggplot2::geom_smooth(data = singles_df, span = 0.15, ggplot2::aes(fill = variable), method = 'loess', formula = 'y ~ x') +
    ggplot2::annotate("text", label = paste0("R = ", round(tox_cor[,"R"], 2), " (", tox_cor[,"n"], ")") , x = 350, y = 15) +
	  ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = c(colour_scheme[["shade 0"]][[1]], "black"))
    d <- d + ggplot2::scale_colour_manual(values = c(colour_scheme[["shade 0"]][[1]], "black"))
  }
	ggplot2::ggsave(file=file.path(outpath, 'tdp43_wt_hydrophobicity_loadings.pdf'), width=8, height=3)

	save(singles_dt, singles_df, plot_df, file = file.path(outpath, 'tardbpdms_wt_hydrophobicity.RData'))
}

