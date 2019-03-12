
#' tardbpdms_quality_control
#'
#' Quality control plots including scatterplots comparing toxicity estimates between replicates.
#'
#' @param toxicity_list named list of folder paths with toxicity estimates (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
tardbpdms_quality_control <- function(
  toxicity_list,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

	#Return previous results if analysis not executed
	if(!execute){
		return()
	}

	#Create output directory
	tardbpdms__create_dir(tardbpdms_dir = outpath)
	
	### Load data
	###########################

	#DMS variant data (for individual replicates)
	dms_dt1 <- tardbpdms__load_toxicity(toxicity_list[[1]], object=T, read_filter=F)
	dms_dt2 <- tardbpdms__load_toxicity(toxicity_list[[2]], object=T, read_filter=F)

	#Rename column names consistently
	names(dms_dt1)[grep('toxicity[1-4]$', names(dms_dt1))] <- paste0("toxicity", 1:length(grep('toxicity[1-4]$', names(dms_dt1))))
	names(dms_dt2)[grep('toxicity[1-4]$', names(dms_dt2))] <- paste0("toxicity", 1:length(grep('toxicity[1-4]$', names(dms_dt2))))
	names(dms_dt1)[grep('toxicity[1-4]_cond$', names(dms_dt1))] <- paste0("toxicity", 1:length(grep('toxicity[1-4]_cond$', names(dms_dt1))), "_cond")
	names(dms_dt2)[grep('toxicity[1-4]_cond$', names(dms_dt2))] <- paste0("toxicity", 1:length(grep('toxicity[1-4]_cond$', names(dms_dt2))), "_cond")

	#Combine regions
	dms_dt1[, region := names(toxicity_list)[1]]
	dms_dt2[, region := names(toxicity_list)[2]]
	dms_dt <- rbind(
		dms_dt1[,.SD,,.SDcols = c("Nmut_aa", "region", "toxicity", "toxicity_cond", "toxicity_uncorr", "STOP", "mean_count", "is.reads0",
			names(dms_dt1)[grep('toxicity[1-4]$', names(dms_dt1))], names(dms_dt1)[grep('toxicity[1-4]_cond$', names(dms_dt1))])], 
		dms_dt2[,.SD,,.SDcols = c("Nmut_aa", "region", "toxicity", "toxicity_cond", "toxicity_uncorr", "STOP", "mean_count", "is.reads0",
			names(dms_dt2)[grep('toxicity[1-4]$', names(dms_dt2))], names(dms_dt2)[grep('toxicity[1-4]_cond$', names(dms_dt2))])])

	### Mean counts vs uncorrected toxicity
	###########################

	#Set uncorrected toxicity to toxicity for singles
	dms_dt[Nmut_aa==1, toxicity_uncorr := toxicity]
	#Set conditional toxicity to toxicity for singles
	dms_dt[Nmut_aa==1, toxicity_cond := toxicity]
	plot_df <- reshape2::melt(as.data.frame(dms_dt[,.(mean_count, toxicity_uncorr, toxicity_cond, region)]), id = c("mean_count", "region"))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(mean_count, value)) +
		ggplot2::scale_x_log10() +
    ggplot2::stat_binhex(bins=30) +
	  ggplot2::xlab("Mean input read count") +
	  ggplot2::ylab("Toxicity") +
	  ggplot2::geom_hline(yintercept = 0, linetype=2) +
	  ggplot2::geom_vline(xintercept = 10, linetype=2, colour = colour_scheme[["shade 0"]][[1]]) +
	  ggplot2::theme_bw() +
	  ggplot2::facet_grid(variable ~ region)
	suppressWarnings(suppressMessages(ggplot2::ggsave(file=file.path(outpath, 'toxicity_vs_mean_input_read_count.pdf'), width=5, height=4)))

	### Scatterplots comparing toxicity estimates between replicates
	###########################

	#Singles, doubles and STOPs mutants with defined toxicity
	dms_dt <- dms_dt[is.reads0==T & (Nmut_aa==1 | Nmut_aa==2) & (!is.na(toxicity) | !is.na(toxicity_cond))]
	# dms_dt <- dms_dt[(Nmut_aa==1 | Nmut_aa==2) & (!is.na(toxicity) | !is.na(toxicity_cond))]
	dms_dt[Nmut_aa==2, toxicity1 := toxicity1_cond]
	dms_dt[Nmut_aa==2, toxicity2 := toxicity2_cond]
	dms_dt[Nmut_aa==2, toxicity3 := toxicity3_cond]
	# dms_dt[Nmut_aa==2, toxicity4 := toxicity4_cond]

	#Convert infinite toxicity values to NA
	dms_dt[is.infinite(toxicity1), toxicity1 := NA]
	dms_dt[is.infinite(toxicity2), toxicity2 := NA]
	dms_dt[is.infinite(toxicity3), toxicity3 := NA]
	# dms_dt[is.infinite(toxicity4), toxicity4 := NA]

	#Number of variants with defined toxicity values in each replicate
	print(paste0("All AA variants (single, double AA substitutions): ", dms_dt[,.N], " (", paste0(unlist(dms_dt[,.N,Nmut_aa][,2]), collapse = ", "), ")"))
	print(paste0("Rep1 AA variants (single, double AA substitutions): ", dms_dt[!is.na(toxicity1),.N], " (",paste0(unlist(dms_dt[!is.na(toxicity1),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	print(paste0("Rep2 AA variants (single, double AA substitutions): ", dms_dt[!is.na(toxicity2),.N], " (",paste0(unlist(dms_dt[!is.na(toxicity2),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	print(paste0("Rep3 AA variants (single, double AA substitutions): ", dms_dt[!is.na(toxicity3),.N], " (",paste0(unlist(dms_dt[!is.na(toxicity3),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	# dms_dt[!is.na(toxicity4),.N]

	#Plot replicate toxicity correlation
	num_rep = 3
	# num_rep = 4
	tardbpdms__ggpairs_binhex(
		input_dt = dms_dt[Nmut_aa %in% c(1,2) & region=="290",.SD,.SDcols = paste0('toxicity', c(1:num_rep))], 
		# input_dt_upper = dms_dt[Nmut_aa==1,.SD,.SDcols = paste0('toxicity', c(1:num_rep))],
		output_file = file.path(outpath, "replicate_scatter_binhex_290.pdf"),
		xlab = "Toxicity",
		ylab = "Toxicity",
		width = 5,
		height = 5,
		label_size = 2,
  	plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]]),
  	colour_limits = c(0, 100))
	tardbpdms__ggpairs_binhex(
		input_dt = dms_dt[Nmut_aa %in% c(1,2) & region=="332",.SD,.SDcols = paste0('toxicity', c(1:num_rep))], 
		# input_dt_upper = dms_dt[Nmut_aa==1,.SD,.SDcols = paste0('toxicity', c(1:num_rep))],
		output_file = file.path(outpath, "replicate_scatter_binhex_332.pdf"),
		xlab = "Toxicity",
		ylab = "Toxicity",
		width = 5,
		height = 5,
		label_size = 2,
  	plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]]),
  	colour_limits = c(0, 100))

	#Plot replicate 1 and 2 correlation
	plot_dt <- copy(dms_dt)[Nmut_aa %in% c(1,2) & !is.na(toxicity1) & !is.na(toxicity2)]
  plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
  temp_cor <- plyr::ddply(plot_dt, c("region"), plyr::summarize, cor = round(cor(toxicity1, toxicity2, use = "pairwise.complete.obs"), 2), n = length(toxicity1))
  d <- ggplot2::ggplot(plot_dt, ggplot2::aes(toxicity1, toxicity2)) +
  	ggplot2::stat_binhex(bins=50) +
  	ggplot2::geom_abline(linetype = 2) + 
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 0, y = 0, size = 2) +
    ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 100)) +
    ggplot2::ylab("Toxicity (rep2)") +
    ggplot2::xlab("Toxicity (rep1)") +
    ggplot2::theme_classic() + 
    ggplot2::facet_grid(region ~.)
  ggplot2::ggsave(file=file.path(outpath, "replicate1_replicate2_scatter.pdf"), width=5, height=5)

	### Scale toxicity by median standard deviation of AA mutants in each library
	###########################

	tardbpdms__scale_toxicity(
	  toxicity_path=toxicity_list[[1]],
	  outpath=file.path(outpath, names(toxicity_list)[1]))

	tardbpdms__scale_toxicity(
	  toxicity_path=toxicity_list[[2]],
	  outpath=file.path(outpath, names(toxicity_list)[2]))

	### Load data
	###########################

	#DMS variant data (for individual replicates)
	dms_dt1 <- tardbpdms__load_toxicity(file.path(outpath, names(toxicity_list)[1]), object=T, read_filter=F)
	dms_dt2 <- tardbpdms__load_toxicity(file.path(outpath, names(toxicity_list)[2]), object=T, read_filter=F)

	#Rename column names consistently
	names(dms_dt1)[grep('toxicity[1-4]$', names(dms_dt1))] <- paste0("toxicity", 1:length(grep('toxicity[1-4]$', names(dms_dt1))))
	names(dms_dt2)[grep('toxicity[1-4]$', names(dms_dt2))] <- paste0("toxicity", 1:length(grep('toxicity[1-4]$', names(dms_dt2))))
	names(dms_dt1)[grep('toxicity[1-4]_cond$', names(dms_dt1))] <- paste0("toxicity", 1:length(grep('toxicity[1-4]_cond$', names(dms_dt1))), "_cond")
	names(dms_dt2)[grep('toxicity[1-4]_cond$', names(dms_dt2))] <- paste0("toxicity", 1:length(grep('toxicity[1-4]_cond$', names(dms_dt2))), "_cond")

	#Combine regions
	dms_dt1[, region := names(toxicity_list)[1]]
	dms_dt2[, region := names(toxicity_list)[2]]
	dms_dt <- rbind(
		dms_dt1[,.SD,,.SDcols = c("Nmut_aa", "region", "toxicity", "toxicity_cond", "toxicity_uncorr", "STOP", "mean_count", "is.reads0",
			names(dms_dt1)[grep('toxicity[1-4]$', names(dms_dt1))], names(dms_dt1)[grep('toxicity[1-4]_cond$', names(dms_dt1))])], 
		dms_dt2[,.SD,,.SDcols = c("Nmut_aa", "region", "toxicity", "toxicity_cond", "toxicity_uncorr", "STOP", "mean_count", "is.reads0",
			names(dms_dt2)[grep('toxicity[1-4]$', names(dms_dt2))], names(dms_dt2)[grep('toxicity[1-4]_cond$', names(dms_dt2))])])

	### Mean counts vs uncorrected toxicity
	###########################

	#Set uncorrected toxicity to toxicity for singles
	dms_dt[Nmut_aa==1, toxicity_uncorr := toxicity]
	#Set conditional toxicity to toxicity for singles
	dms_dt[Nmut_aa==1, toxicity_cond := toxicity]
	plot_df <- reshape2::melt(as.data.frame(dms_dt[,.(mean_count, toxicity_uncorr, toxicity_cond, region)]), id = c("mean_count", "region"))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(mean_count, value)) +
		ggplot2::scale_x_log10() +
    ggplot2::stat_binhex(bins=30) +
	  ggplot2::xlab("Mean input read count") +
	  ggplot2::ylab("Toxicity") +
	  ggplot2::geom_hline(yintercept = 0, linetype=2) +
	  ggplot2::geom_vline(xintercept = 10, linetype=2, colour = colour_scheme[["shade 0"]][[1]]) +
	  ggplot2::theme_bw() +
	  ggplot2::facet_grid(variable ~ region)
	suppressWarnings(suppressMessages(ggplot2::ggsave(file=file.path(outpath, 'toxicity_vs_mean_input_read_count_scale.pdf'), width=5, height=4)))

	### Scatterplots comparing toxicity estimates between replicates
	###########################

	#Singles, doubles and STOPs mutants with defined toxicity
	dms_dt <- dms_dt[is.reads0==T & (Nmut_aa==1 | Nmut_aa==2) & (!is.na(toxicity) | !is.na(toxicity_cond))]
	# dms_dt <- dms_dt[(Nmut_aa==1 | Nmut_aa==2) & (!is.na(toxicity) | !is.na(toxicity_cond))]
	dms_dt[Nmut_aa==2, toxicity1 := toxicity1_cond]
	dms_dt[Nmut_aa==2, toxicity2 := toxicity2_cond]
	dms_dt[Nmut_aa==2, toxicity3 := toxicity3_cond]
	# dms_dt[Nmut_aa==2, toxicity4 := toxicity4_cond]

	#Convert infinite toxicity values to NA
	dms_dt[is.infinite(toxicity1), toxicity1 := NA]
	dms_dt[is.infinite(toxicity2), toxicity2 := NA]
	dms_dt[is.infinite(toxicity3), toxicity3 := NA]
	# dms_dt[is.infinite(toxicity4), toxicity4 := NA]

	#Number of variants with defined toxicity values in each replicate
	print(paste0("All AA variants (single, double AA substitutions): ", dms_dt[,.N], " (", paste0(unlist(dms_dt[,.N,Nmut_aa][,2]), collapse = ", "), ")"))
	print(paste0("Rep1 AA variants (single, double AA substitutions): ", dms_dt[!is.na(toxicity1),.N], " (",paste0(unlist(dms_dt[!is.na(toxicity1),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	print(paste0("Rep2 AA variants (single, double AA substitutions): ", dms_dt[!is.na(toxicity2),.N], " (",paste0(unlist(dms_dt[!is.na(toxicity2),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	print(paste0("Rep3 AA variants (single, double AA substitutions): ", dms_dt[!is.na(toxicity3),.N], " (",paste0(unlist(dms_dt[!is.na(toxicity3),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	# dms_dt[!is.na(toxicity4),.N]

	#Plot replicate toxicity correlation
	num_rep = 3
	# num_rep = 4
	tardbpdms__ggpairs_binhex(
		input_dt = dms_dt[Nmut_aa %in% c(1,2) & region=="290",.SD,.SDcols = paste0('toxicity', c(1:num_rep))], 
		# input_dt_upper = dms_dt[Nmut_aa==1,.SD,.SDcols = paste0('toxicity', c(1:num_rep))],
		output_file = file.path(outpath, "replicate_scatter_binhex_290_scale.pdf"),
		xlab = "Toxicity",
		ylab = "Toxicity",
		width = 5,
		height = 5,
		label_size = 2,
  	plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]]),
  	colour_limits = c(0, 100))
	tardbpdms__ggpairs_binhex(
		input_dt = dms_dt[Nmut_aa %in% c(1,2) & region=="332",.SD,.SDcols = paste0('toxicity', c(1:num_rep))], 
		# input_dt_upper = dms_dt[Nmut_aa==1,.SD,.SDcols = paste0('toxicity', c(1:num_rep))],
		output_file = file.path(outpath, "replicate_scatter_binhex_332_scale.pdf"),
		xlab = "Toxicity",
		ylab = "Toxicity",
		width = 5,
		height = 5,
		label_size = 2,
  	plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]]),
  	colour_limits = c(0, 100))

	#Plot replicate 1 and 2 correlation
	plot_dt <- copy(dms_dt)[Nmut_aa %in% c(1,2) & !is.na(toxicity1) & !is.na(toxicity2)]
  plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
  temp_cor <- plyr::ddply(plot_dt, c("region"), plyr::summarize, cor = round(cor(toxicity1, toxicity2, use = "pairwise.complete.obs"), 2), n = length(toxicity1))
  d <- ggplot2::ggplot(plot_dt, ggplot2::aes(toxicity1, toxicity2)) +
  	ggplot2::stat_binhex(bins=50) +
  	ggplot2::geom_abline(linetype = 2) + 
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 0, y = 0, size = 2) +
    ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 100)) +
    ggplot2::ylab("Toxicity (rep2)") +
    ggplot2::xlab("Toxicity (rep1)") +
    ggplot2::theme_classic() + 
    ggplot2::facet_grid(region ~.)
  ggplot2::ggsave(file=file.path(outpath, "replicate1_replicate2_scatter_scale.pdf"), width=5, height=5)

}

