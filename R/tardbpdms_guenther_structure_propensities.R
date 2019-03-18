
#' tardbpdms_guenther_structure_propensities
#'
#' Guenther structure predictions.
#'
#' @param result_dir directory with Guenther structure match results (required)
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
tardbpdms_guenther_structure_propensities <- function(
  result_dir,
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
  message(paste("\n\n*******", "running stage: tardbpdms_guenther_structure_propensities", "*******\n\n"))

	#Create output directory
	tardbpdms__create_dir(tardbpdms_dir = outpath)

	### Re-run structure analysis
	###########################

	if(rerun_structure){

	  #Display status
	  message(paste("\n\n*******", "re-running structure propensity calculations (this might take a while)", "*******\n\n"))

	  #Create output directory
	  tardbpdms__create_dir(tardbpdms_dir = file.path(miscpath, "misc_guenther_structures"))

		#Run misc script on command-line
	  system(paste0(
	  	file.path(miscpath, "scripts", "tardbpdms__guenther_structures.R"),
	  	" -o ",
	  	file.path(miscpath, "misc_guenther_structures"),
	  	" -d ",
	  	DMS2structure_path,
	  	" -e ",
	  	file.path(miscpath, "misc_epistasis_analysis"),
	  	" -p ",
	  	file.path(miscpath, "pdb")))
	}

	### Setup
	###########################

	#List of structures
	pdb_list <- list("5wkd" = list("aa_seq" = "GNNQGSN", "idx_start" = 300), 
	  "5whn" = list("aa_seq" = "NFGAFS", "idx_start" = 312), 
	  "5whp" = list("aa_seq" = "NFGTFS", "idx_start" = 312),
	  "5wkb" = list("aa_seq" = "NFGEFS", "idx_start" = 312),  
	  # "6cf4" = list("aa_seq" = "NFGTFS", "idx_start" = 312),  
	  "6cew" = list("aa_seq" = "AMMAAA", "idx_start" = 321), 
	  "6cb9" = list("aa_seq" = "AALQSS", "idx_start" = 328), 
	  "6cfh" = list("aa_seq" = "SWGMMGMLASQ", "idx_start" = 333), 
	  "5wia" = list("aa_seq" = "GNNSYS", "idx_start" = 370))

	kernel_scores_kernal_width <- list()
	for(i in names(pdb_list)){
	  kernel_scores_kernal_width[[i]] <- fread(file.path(result_dir, "pdb", i, "processed_data", "rand_strategy_kernal_width_kernel_structure_propensity.txt"))
	}

	### Kernel smooth scores - randomisation strategy: kernal width
	###########################

	#Distinct segments
	pdb_segs <- c("5wkd", "5whn", "6cew", "6cb9", "6cfh", "5wia")
	#LARKS mutants
	pdb_LARKS <- c("5whn", "5whp", "5wkb")
	pdb_sets <- list("segs" = pdb_segs, "LARKS" = pdb_LARKS)
	for(pdb_set_name in names(pdb_sets)){
	  pdb_names <- pdb_sets[[pdb_set_name]]
	  for(score_name in c("association_score_norm")){
	    tardbpdms__guenther_structure_propensity_plot(
	    	kernel_scores=kernel_scores_kernal_width, 
	    	outpath=file.path(outpath, paste0("kernel_propensity_kernal_width_", pdb_set_name, "_", score_name, ".pdf")),
	    	colour_scheme=colour_scheme,
	    	pdb_names=pdb_names, 
	    	score_type=score_name, 
	    	pdb_list=pdb_list, 
	    	position_offset=290)
	   }
	}

}








