
#' tardbpdms__scale_toxicity
#'
#' Scale toxicity estimates between replicates.
#'
#' @param toxicity_path folder path to toxicity estimates (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
tardbpdms__scale_toxicity <- function(
  toxicity_path,
  outpath
  ){

	#Create output directory
	tardbpdms__create_dir(tardbpdms_dir = outpath)
	
	### Load data
	###########################

	#DMS variant data (for individual replicates)
	dms_dt <- tardbpdms__load_toxicity(toxicity_path, object=T, read_filter=F)

	### Scale fitness by maximum standard deviation of AA mutants in each library
	###########################

	#Standard deviation of AA mutants in each replicate
	dms_dt_sd <- list()
	num_reps <- as.numeric(gsub("fitness", "", names(dms_dt)[grep("fitness[1-4]$", names(dms_dt))]))
	for(j in num_reps){
		mut_aa_fitness <- unlist(dms_dt[Nmut_aa!=0,.SD,,.SDcols = c(paste0("fitness", j), paste0("fitness", j, "_cond"))])
		mut_aa_fitness[is.infinite(mut_aa_fitness)] <- NA
		dms_dt_sd[[as.character(j)]] <- sd(mut_aa_fitness, na.rm = T)
	}

	#Scale fitness and error
	for(j in num_reps){
		for(p in c("fitness", "sigma")){
			for(s in c("", "_cond", "_uncorr")){
				this_column <- paste0(p, j, s)
				dms_dt[, (this_column) := .SD[[1]]*max(unlist(dms_dt_sd))/dms_dt_sd[[as.character(j)]],,.SDcols = this_column]
			}
		}
	}
	
	### Estimate replicate error
	###########################

	#Set uncorrected fitness and sigma for variants that are not doubles (fitness[1-4]_uncorr and sigma[1-4]_uncorr to fitness[1-4] and sigma[1-4])
	for(i in paste0("fitness", num_reps, "_uncorr")){
		j <- gsub('_uncorr','',i)
		dms_dt[Nmut_aa!=2,(i) := .SD[[1]],,.SDcols=j]
		i <- gsub('fitness','sigma',i)
		j <- gsub('fitness','sigma',j)
		dms_dt[Nmut_aa!=2,(i) := .SD[[1]],,.SDcols=j]
	}

	#Replicate error
	num_reps_str <- paste0(num_reps, collapse = "")

	#Calculate
	dms_dt[,var_fitness := rowSums((rowMeans(.SD[,1:nchar(num_reps_str)],na.rm=T) - .SD[,1:nchar(num_reps_str)])^2/(.SD[,(nchar(num_reps_str)+1):(2*nchar(num_reps_str))]^2),na.rm=T) / 
	              rowSums(1/(.SD[,(nchar(num_reps_str)+1):(2*nchar(num_reps_str))]^2 ),na.rm=T),
	            ,.SDcols = c(grep(names(dms_dt),pattern=paste0("fitness[", num_reps_str, "]_uncorr")),grep(names(dms_dt),pattern=paste0("sigma[", num_reps_str, "]_uncorr")))]
	dms_dt[,avg_sigma := rowMeans(.SD[,(nchar(num_reps_str)+1):(2*nchar(num_reps_str))],na.rm=T),
	            ,.SDcols = c(grep(names(dms_dt),pattern=paste0("fitness[", num_reps_str, "]_uncorr")),grep(names(dms_dt),pattern=paste0("sigma[", num_reps_str, "]_uncorr")))]
	dms_dt[,isNA := rowSums(is.na(.SD)),,.SDcols = grep(names(dms_dt),pattern=paste0("fitness[", num_reps_str, "]_uncorr"))]
	dms_dt[isNA==0 & var_fitness != Inf & avg_sigma != Inf,avg_sigma_fit := loess(var_fitness ~ avg_sigma,span=0.75)$fitted]
	replicate_error <- dms_dt[,min(sqrt(avg_sigma_fit),na.rm=T)]

	#Plot average sigma versus fitness replicate error
	d <- ggplot2::ggplot(dms_dt[isNA==0 & !is.infinite(avg_sigma)],ggplot2::aes(avg_sigma)) +
	  ggplot2::geom_point(ggplot2::aes(y=sqrt(var_fitness))) +
	  ggplot2::geom_line(ggplot2::aes(y=sqrt(avg_sigma_fit)),color="red") +
	  ggplot2::geom_line(ggplot2::aes(y=sqrt(avg_sigma^2 + replicate_error^2)),color="darkgreen") +
	  ggplot2::geom_abline(color="yellow",linetype=2) + 
	  ggplot2::scale_x_log10() + 
	  ggplot2::scale_y_log10()
	suppressWarnings(suppressMessages(ggplot2::ggsave(file.path(outpath, "replicateerror_vs_avgsigma.pdf"), width = 7, height = 5)))

	#print(replicate_error)

	### Merge replicate fitness values
	###########################

	#### singles
	fitness_rx = dms_dt[Nmut_aa!=2,.SD,.SDcols = grep(paste0("fitness[", num_reps_str, "]$"),colnames(dms_dt))]
	sigma_rx = sqrt(dms_dt[Nmut_aa!=2,.SD,.SDcols = grep(paste0("sigma[", num_reps_str, "]$"),colnames(dms_dt))]^2 + 
	                  matrix(replicate_error^2,nrow = dim(fitness_rx)[1],ncol = dim(fitness_rx)[2]))
	dms_dt[Nmut_aa!=2,fitness := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
	dms_dt[Nmut_aa!=2,sigma := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
	d <- ggplot2::ggplot(dms_dt[Nmut_aa==1],ggplot2::aes(fitness,sigma)) + 
	  ggplot2::geom_hex() + 
	  ggplot2::scale_y_log10() + 
	  ggplot2::coord_cartesian(ylim=c(0.01,1))
	ggplot2::ggsave(file.path(outpath, "sigma_vs_fitness_singles.pdf"), width = 5, height = 5)

	#### doubles
	#uncorrected fitness
	fitness_rx = dms_dt[Nmut_aa==2,.SD,.SDcols = grep(paste0("fitness[", num_reps_str, "]_uncorr"),colnames(dms_dt))]
	sigma_rx = sqrt(dms_dt[Nmut_aa==2,.SD,.SDcols = grep(paste0("sigma[", num_reps_str, "]_uncorr"),colnames(dms_dt))]^2 + 
	                  matrix(replicate_error^2,nrow = dim(fitness_rx)[1],ncol = dim(fitness_rx)[2]))
	dms_dt[Nmut_aa==2,fitness_uncorr := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
	dms_dt[Nmut_aa==2,sigma_uncorr := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
	d <- ggplot2::ggplot(dms_dt[Nmut_aa==2 & !is.na(fitness_uncorr),],ggplot2::aes(fitness_uncorr,sigma_uncorr)) + 
	  ggplot2::geom_hex() + 
	  ggplot2::scale_y_log10() + 
	  ggplot2::coord_cartesian(ylim=c(0.01,1))
	ggplot2::ggsave(file.path(outpath, "sigma_vs_fitness_doubles_uncorr.pdf"), width = 5, height = 5)

	#### doubles
	#conditioned fitness
	fitness_rx = dms_dt[Nmut_aa==2,.SD,.SDcols = grep(paste0("fitness[", num_reps_str, "]_cond"),colnames(dms_dt))]
	sigma_rx = sqrt(dms_dt[Nmut_aa==2,.SD,.SDcols = grep(paste0("sigma[", num_reps_str, "]_cond"),colnames(dms_dt))]^2 + 
	                  matrix(replicate_error^2,nrow = dim(fitness_rx)[1],ncol = dim(fitness_rx)[2]))
	dms_dt[Nmut_aa==2,fitness_cond := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
	dms_dt[Nmut_aa==2,sigma_cond := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
	d <- ggplot2::ggplot(dms_dt[Nmut_aa==2,],ggplot2::aes(fitness_cond,sigma_cond)) + 
	  ggplot2::geom_hex() + 
	  ggplot2::scale_y_log10() + 
	  ggplot2::coord_cartesian(ylim=c(0.01,1))
	suppressWarnings(suppressMessages(ggplot2::ggsave(file.path(outpath, "sigma_vs_fitness_doubles_cond.pdf"), width = 5, height = 5)))

	### Output replicate data files
	###########################

	#Save objects
	#Merge replicate fitness values separately for each library
	singles_silent <- copy(dms_dt)[Nmut_aa!=2,.SD,,.SDcols = names(dms_dt)[!grepl("toxicity", names(dms_dt))]]
	doubles <- copy(dms_dt)[Nmut_aa==2,.SD,,.SDcols = names(dms_dt)[!grepl("toxicity", names(dms_dt))]]
	save(singles_silent, doubles, file = file.path(outpath, "DMS_processed_data.RData"))

	##### finalize data.tables
	silent = singles_silent[Nmut_aa==0,.(Pos,WT_AA,Mut,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,fitness,sigma)]
	singles = singles_silent[Nmut_aa==1,.(Pos,WT_AA,Mut,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,fitness,sigma)]

	#for doubles #add single mutant fitness/sigma values to double mutant table
	doubles[,fitness1 := singles[Pos == Pos1 & Mut == Mut1,fitness],.(Pos1,Mut1)]
	doubles[,sigma1 := singles[Pos == Pos1 & Mut == Mut1,sigma],.(Pos1,Mut1)]
	doubles[,fitness2 := singles[Pos == Pos2 & Mut == Mut2,fitness],.(Pos2,Mut2)]
	doubles[,sigma2 := singles[Pos == Pos2 & Mut == Mut2,sigma],.(Pos2,Mut2)]

	doubles = doubles[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,
	                     fitness1,sigma1,fitness2,sigma2,
	                     fitness_uncorr,sigma_uncorr,
	                     fitness_cond,sigma_cond)]

	#Exclude variants with STOP codons from downstream fitness analyses
	silent[,is.fitness := TRUE]
	singles[,is.fitness := !STOP]
	doubles[,is.fitness := !STOP]

	#write data to files
	write.table(x = silent, file = file.path(outpath, "DMS_silent.txt"),
	            quote = F,row.names = F, col.names = T)
	write.table(x = singles, file = file.path(outpath, "DMS_singles.txt"),
	            quote = F,row.names = F, col.names = T)
	write.table(x = doubles, file = file.path(outpath, "DMS_doubles.txt"),
	            quote = F,row.names = F, col.names = T)

}

