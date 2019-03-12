
#' tardbpdms_combine_toxicity
#'
#' Combine toxicity estimates from two different DMS experiments (regions).
#'
#' @param toxicity_list named list of folder paths with toxicity estimates (required)
#' @param growth_rate_file path to growth rate file (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return A data.table with normalised toxicity data
#' @export
#' @import data.table
tardbpdms_combine_toxicity <- function(
  toxicity_list,
  growth_rate_file,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

	#Return previous results if analysis not executed
	if(!execute){
		load(file.path(outpath, "combine_toxicity_values.RData"))
		return(dms_dt_center_scale)
	}

	#Create output directory
	tardbpdms__create_dir(tardbpdms_dir = outpath)
	
	### Load data
	###########################

	#Growth rate data
	gr_df <- read.table(growth_rate_file, header = F, skip = 1, stringsAsFactors = F, sep = "\t", quote = "", row.names = 1)[,1:4]
	#Mean growth rate table
	mean_gr_df <- data.frame(growth_rate = apply(log(gr_df), 1, mean, na.rm = T))
	mean_gr_df[is.na(mean_gr_df)] <- 0
	sd_gr_df <- data.frame(growth_rate = apply(log(gr_df), 1, sd, na.rm = T))
	sd_gr_df[is.na(sd_gr_df)] <- 0

	#DMS variant data 
	dms_dt1 <- tardbpdms__load_toxicity(toxicity_list[[1]])
	dms_dt2 <- tardbpdms__load_toxicity(toxicity_list[[2]])

	#Combine regions
	dms_dt1[, region := names(toxicity_list)[1]]
	dms_dt2[, region := names(toxicity_list)[2]]
	dms_dt <- rbind(dms_dt1, dms_dt2)
	#Remove missing fitness values
	dms_dt <- dms_dt[!is.na(fitness) | (!is.na(fitness_cond) & !is.na(fitness1) & !is.na(fitness2)),]
	#Absolute position (singles)
	dms_dt[, Pos_abs := Pos+as.numeric(region)-1]
	#Absolute position (doubles)
	dms_dt[, Pos_abs1 := Pos1+as.numeric(region)-1]
	dms_dt[, Pos_abs2 := Pos2+as.numeric(region)-1]
	#Mutation code (singles)
	dms_dt[, mut_code := paste0(WT_AA, Pos_abs, Mut)]
	#Mutation code (doubles)
	dms_dt[, mut_code1 := paste0(WT_AA1, Pos_abs1, Mut1)]
	dms_dt[, mut_code2 := paste0(WT_AA2, Pos_abs2, Mut2)]

	#Add growth rate data
	dms_dt[, growth_rate := mean_gr_df[mut_code,]]
	dms_dt[, growth_rate_low := mean_gr_df[mut_code,]-1.96*sd_gr_df[mut_code,]]
	dms_dt[, growth_rate_high := mean_gr_df[mut_code,]+1.96*sd_gr_df[mut_code,]]

	#Annotate number and type of mutations
	dms_dt[, mut_type := paste0("Nmut_aa: ", Nmut_aa, ", Nmut_codons: ", Nmut_codons)]
	dms_dt[STOP==T, mut_type := "STOP"]
	dms_dt[Nmut_codons >= 3 & Nmut_aa==0, mut_type := paste0("Nmut_aa: 0, Nmut_codons: >=3")]
	dms_dt[, mut_type := factor(mut_type, levels = unique(mut_type[order(Nmut_aa+10*as.numeric(STOP), Nmut_codons, decreasing = F)]))]

	### 1. Untransformed data
	###########################

	#Plot distributions split by number and type of mutation
	tardbpdms__plot_toxicity_distributions(
		input_dt = dms_dt, 
		output_file = file.path(outpath, 'single_double_syn_mutant_toxicity_hist.pdf'),
		colour_scheme = c("black", "darkgrey"))

	#Plot correlation between toxicity estimates and growth rate
	tardbpdms__plot_toxicity_vs_growthrate(
		input_dt = dms_dt, 
		output_file = file.path(outpath, 'single_mutant_toxicity_growthrate.pdf'),
		colour_scheme = c("black", "darkgrey"))
	#Plot correlation between toxicity estimates and growth rate (separately for each region)
	for(i in 1:length(toxicity_list)){
		tardbpdms__plot_toxicity_vs_growthrate(
			input_dt = dms_dt[region==names(toxicity_list)[i],], 
			output_file = file.path(outpath, paste0('single_mutant_toxicity_growthrate_', names(toxicity_list)[i], '.pdf')),
		colour_scheme = c("black", "darkgrey"))
	}

	### 2. Center fitness scores using weighted mean of single codon silent mutations
	###########################

	#Center fitness scores using weighted mean of single codon silent mutations
	# silent_codon_center <- dms_dt[mut_type=="Nmut_aa: 0, Nmut_codons: 1",.(mean=median(fitness, na.rm = T)),region]
	silent_codon_center <- dms_dt[mut_type=="Nmut_aa: 0, Nmut_codons: 1",.(mean=weighted.mean(fitness, sigma^-2, na.rm = T)),region]
	dms_dt_center <- copy(dms_dt)
	for(i in names(dms_dt_center)[grep("^fitness", names(dms_dt_center))]){
	  dms_dt_center[region==names(toxicity_list)[1],(i) := .SD - silent_codon_center[region==names(toxicity_list)[1],mean],,.SDcols=i]
	  dms_dt_center[region==names(toxicity_list)[2],(i) := .SD - silent_codon_center[region==names(toxicity_list)[2],mean],,.SDcols=i]
    j <- gsub('fitness','toxicity',i)
    dms_dt_center[, (j) := -.SD[[1]],,.SDcols = i]
	}

	#Plot distributions split by number and type of mutation
	tardbpdms__plot_toxicity_distributions(
		input_dt = dms_dt_center, 
		output_file = file.path(outpath, 'single_double_syn_mutant_toxicity_hist_center.pdf'),
		colour_scheme = c("black", "darkgrey"))

	#Plot correlation between toxicity estimates and growth rate
	tardbpdms__plot_toxicity_vs_growthrate(
		input_dt = dms_dt_center, 
		output_file = file.path(outpath, 'single_mutant_toxicity_growthrate_center.pdf'),
		colour_scheme = c("black", "darkgrey"))

	### 3. Scale fitness scores using ratio between weighted mean of single STOP codon fitness
	###########################

	#Scale fitness scores using ratio between weighted mean of single STOP codon fitness
	stop_center <- dms_dt_center[mut_type=="STOP" & Nmut_aa==1,.(mean=weighted.mean(fitness, sigma^-2, na.rm = T)),region]
	dms_dt_center_scale <- copy(dms_dt_center)
	for(i in names(dms_dt_center_scale)[grep("^fitness", names(dms_dt_center_scale))]){
	  dms_dt_center_scale[region==names(toxicity_list)[2], (i) := .SD * stop_center[region==names(toxicity_list)[1],mean]/stop_center[region==names(toxicity_list)[2],mean],,.SDcols=i]
	  j <- gsub("fitness", "sigma", i)
	  dms_dt_center_scale[region==names(toxicity_list)[2], (j) := .SD * stop_center[region==names(toxicity_list)[1],mean]/stop_center[region==names(toxicity_list)[2],mean],,.SDcols=j]
    j <- gsub('fitness','toxicity',i)
    dms_dt_center_scale[, (j) := -.SD[[1]],,.SDcols = i]
	}

	#Plot distributions split by number and type of mutation
	tardbpdms__plot_toxicity_distributions(
		input_dt = dms_dt_center_scale, 
		output_file = file.path(outpath, 'single_double_syn_mutant_toxicity_hist_center_scale.pdf'),
		colour_scheme = c("black", "darkgrey"))

	#Plot correlation between toxicity estimates and growth rate
	tardbpdms__plot_toxicity_vs_growthrate(
		input_dt = dms_dt_center_scale, 
		output_file = file.path(outpath, 'single_mutant_toxicity_growthrate_center_scale.pdf'),
		colour_scheme = c("black", "darkgrey"))

	### Toxicity w.r.t. position (singles)
	###########################

	#Repeat for before normalisation and after centering and scaling
	dms_dt_list <- list(
		"before_norm" = dms_dt,
		"center" = dms_dt_center,
		"center_scale" = dms_dt_center_scale
		)
	for(i in names(dms_dt_list)){
		tardbpdms__plot_toxicity_vs_position(
			input_dt = dms_dt_list[[i]][STOP==T & Nmut_aa==1,], 
			use_sigma=T, 
			output_file=file.path(outpath, paste0('toxicity_vs_position_STOP_', i, '.pdf')), 
			ylab="Relative fitness (log)", 
			label_offset = 0.01)
		tardbpdms__plot_toxicity_vs_position(
			input_dt = dms_dt_list[[i]][STOP==F & Nmut_aa==1,], 
			use_sigma=F, 
			output_file=file.path(outpath, paste0('toxicity_vs_position_single_AA_mutants_', i, '.pdf')), 
			ylab="Relative fitness (log)", 
			label_offset = 0.01)
		tardbpdms__plot_toxicity_vs_position(
			input_dt = copy(dms_dt_list[[i]])[, toxicity := abs(toxicity)][STOP==F & Nmut_aa==1,], 
			use_sigma=F, 
			output_file=file.path(outpath, paste0('toxicityabs_vs_position_single_AA_mutants_', i, '.pdf')), 
			ylab="Absolute relative fitness (log)", 
			span = 0.25, 
			label_offset = 0.01,
			plot_hotspot = T)
	}

	### Bounds on fitness
	###########################

	#Lower bound on fitness
	lower_bound_F <- list()
	for(i in 1:length(toxicity_list)){
		doubles_surpassing_upper_bound_F <- dms_dt_center_scale[is.fitness & Nmut_aa==2 & region==names(toxicity_list)[i] & fitness_cond>stop_center[region==names(toxicity_list)[i],mean],.N,]
		doubles_positive_F <- dms_dt_center_scale[is.fitness & Nmut_aa==2 & region==names(toxicity_list)[i] & fitness_cond>0,.N,]
		lower_bound_quantile <- doubles_surpassing_upper_bound_F / doubles_positive_F
		lower_bound_F[[names(toxicity_list)[i]]] <- quantile(dms_dt_center_scale[is.fitness & Nmut_aa==2 & region==names(toxicity_list)[i] & fitness_cond<0,fitness_cond,], lower_bound_quantile)
	}

	#Plot centered and scaled distributions split by number and type of mutation with bounds
	temp_cols <- as.list(c("black", "darkgrey", "black"))
	names(temp_cols) <- c(names(toxicity_list)[1], names(toxicity_list)[2], "STOP")
	dms_dt_center_scale_boundaries <- list(
		-lower_bound_F[[1]], 
		-lower_bound_F[[2]],
		"STOP" = -stop_center[region==names(toxicity_list)[1],mean])
	names(dms_dt_center_scale_boundaries)[1:2] <- names(toxicity_list)
	tardbpdms__plot_toxicity_distributions(
		input_dt = dms_dt_center_scale, 
		output_file = file.path(outpath, "single_double_syn_mutant_toxicity_hist_center_scale_bounds.pdf"), 
		boundary_list=dms_dt_center_scale_boundaries, 
		boundary_colour_list=temp_cols,
		colour_scheme = c("black", "darkgrey"))

	#Plot centered and scaled distributions split by number and type of mutation with bounds
	temp_cols <- as.list(c("black", "darkgrey", "black"))
	names(temp_cols) <- c(names(toxicity_list)[1], names(toxicity_list)[2], "STOP")
	dms_dt_center_scale_boundaries <- list(
		-lower_bound_F[[1]], 
		-lower_bound_F[[2]],
		"STOP" = -stop_center[region==names(toxicity_list)[1],mean])
	names(dms_dt_center_scale_boundaries)[1:2] <- names(toxicity_list)
	tardbpdms__plot_toxicity_distributions(
		input_dt = dms_dt_center_scale, 
		output_file = file.path(outpath, "single_double_syn_mutant_toxicity_hist_center_scale_upperbound.pdf"), 
		boundary_list=dms_dt_center_scale_boundaries[3], 
		boundary_colour_list=temp_cols[3],
		colour_scheme = c("black", "darkgrey"))

	#Plot centered and scaled distributions split by number and type of mutation with bounds (x axis limits according to 99% of single and double mutant data)
	xlim_focus <- quantile(unlist(dms_dt_center_scale[Nmut_aa==1 | Nmut_aa==2,.(fitness, fitness_cond)]), probs = c(0.005, 0.995), na.rm = T)
	xlim_focus[1] <- -(xlim_focus[1]-0.15)
	xlim_focus[2] <- -(xlim_focus[2]+0.15)
	tardbpdms__plot_toxicity_distributions(
		input_dt = dms_dt_center_scale[Nmut_aa %in% c(1, 2) & !STOP,], 
		output_file = file.path(outpath, 'single_double_syn_mutant_toxicity_hist_center_scale_bounds_xlimfocus.pdf'),
		boundary_list=dms_dt_center_scale_boundaries, 
		boundary_colour_list=temp_cols,
		xlim = xlim_focus,
		colour_scheme = c("black", "darkgrey"),
		width = 4,
		height = 3)

	#Plot centered and scaled distributions split by number and type of mutation with STOP bound (x axis limits according to 99% of single and double mutant data)
	xlim_focus <- quantile(unlist(dms_dt_center_scale[Nmut_aa==1 | Nmut_aa==2,.(fitness, fitness_cond)]), probs = c(0.005, 0.995), na.rm = T)
	xlim_focus[1] <- -(xlim_focus[1]-0.15)
	xlim_focus[2] <- -(xlim_focus[2]+0.15)
	tardbpdms__plot_toxicity_distributions(
		input_dt = dms_dt_center_scale[Nmut_aa %in% c(1, 2) & !STOP,], 
		output_file = file.path(outpath, 'single_double_syn_mutant_toxicity_hist_center_scale_upperbound_xlimfocus.pdf'),
		boundary_list=dms_dt_center_scale_boundaries[3], 
		boundary_colour_list=temp_cols[3],
		xlim = xlim_focus,
		colour_scheme = c("black", "darkgrey"),
		width = 4,
		height = 3)

	#Plot centered and scaled distributions (for singles and doubles combined) with STOP bound (x axis limits according to 99% of single and double mutant data)
	tardbpdms__plot_toxicity_distributions(
		input_dt = dms_dt_center_scale[Nmut_aa %in% c(1, 2) & !STOP,], 
		output_file = file.path(outpath, 'single_double_syn_mutant_toxicity_hist_center_scale_upperbound_xlimfocus_combined.pdf'),
		boundary_list=dms_dt_center_scale_boundaries[3], 
		boundary_colour_list=temp_cols[3],
		xlim = xlim_focus,
		colour_scheme = c("black", "darkgrey"),
		width = 4,
		height = 3,
		facet_by_mut_type = F)

	### Mutation statistics
	###########################

	#Variant classification statistics
	plot_dt <- dms_dt_center_scale
	plot_dt[region == names(toxicity_list)[1],.N,mut_type]
	plot_dt[region == names(toxicity_list)[2],.N,mut_type]

	#Significantly more/less toxic variants
	tox_dt <- copy(dms_dt_center_scale)
	tox_dt[Nmut_aa==1, toxicity_z := toxicity/sigma]
	tox_dt[Nmut_aa==2, toxicity_z := toxicity_cond/sigma_cond]
	tox_dt[, toxicity_sig := p.adjust(2*pnorm(-abs(toxicity_z)), method = "BH")<0.05]
	#Singles
	print(paste0("Significantly more toxic single AA variants: ", tox_dt[Nmut_aa==1 & toxicity_sig & toxicity>0,.N,]))
	print(paste0("Significantly less toxic single AA variants: ", tox_dt[Nmut_aa==1 & toxicity_sig & toxicity<0,.N,]))
	#Doubles
	print(paste0("Significantly more toxic double AA variants: ", tox_dt[Nmut_aa==2 & toxicity_sig & toxicity_cond>0,.N,]))
	print(paste0("Significantly less toxic double AA variants: ", tox_dt[Nmut_aa==2 & toxicity_sig & toxicity_cond<0,.N,]))
	#Both
	print(paste0("Significantly more toxic AA variants: ", tox_dt[((Nmut_aa==1 & toxicity>0) | (Nmut_aa==2 & toxicity_cond>0)) & toxicity_sig,.N,]))
	print(paste0("Significantly less toxic AA variants: ", tox_dt[((Nmut_aa==1 & toxicity<0) | (Nmut_aa==2 & toxicity_cond<0)) & toxicity_sig,.N,]))

	### Save combined fitness
	###########################

	#Write all single AA mutant tables to output files
	fwrite(dms_dt[Nmut_aa==1 & !STOP,], file = file.path(outpath, "combine_toxicity_values_singles_raw.txt"), sep = "\t")
	fwrite(dms_dt_center[Nmut_aa==1 & !STOP,], file = file.path(outpath, "combine_toxicity_values_singles_center.txt"), sep = "\t")
	fwrite(dms_dt_center_scale[Nmut_aa==1 & !STOP,], file = file.path(outpath, "combine_toxicity_values_singles_center_scale.txt"), sep = "\t")
	
	#RData object
	save(dms_dt, dms_dt_center, dms_dt_center_scale, dms_dt_center_scale_boundaries, file = file.path(outpath, "combine_toxicity_values.RData"))

	#Supplementary data file
	fwrite(dms_dt_center_scale, file = file.path(outpath, "supplementary_table_normalised_toxicity_estimates_all.txt"), sep = "\t")

	#Return normalised toxicity data.table
	return(dms_dt_center_scale)
}

