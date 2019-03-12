
#' tardbpdms__load_toxicity
#'
#' Load variant toxicity estimates.
#'
#' @param dms_dir path to variant toxicity directory (required)
#' @param object load from R data object (default:F)
#' @param read_filter filter reads based on is.reads0 (default:T)
#'
#' @return a variant toxicity data.table
#' @export
tardbpdms__load_toxicity <- function(
  dms_dir,
  object=F,
  read_filter=T
  ){
  #Load variant toxicity data
  if(object){
  	load(file.path(dms_dir, "DMS_processed_data.RData"))
  	silent = singles_silent[Nmut_aa==0]
  	singles = singles_silent[Nmut_aa==1] 
  }else{
    silent <- data.table::fread(file.path(dms_dir, "DMS_silent.txt"))
	  singles <- data.table::fread(file.path(dms_dir, "DMS_singles.txt"))
	  doubles <- data.table::fread(file.path(dms_dir, "DMS_doubles.txt"))
  }
  #Toxicity - silent
  for(i in names(silent)[grep("^fitness", names(silent))]){
    j <- gsub('fitness','toxicity',i)
    silent[, (j) := -.SD[[1]],,.SDcols = i]
  }
  #Toxicity - singles
  for(i in names(singles)[grep("^fitness", names(singles))]){
    j <- gsub('fitness','toxicity',i)
    singles[, (j) := -.SD[[1]],,.SDcols = i]
  }
  #Toxicity - doubles
  for(i in names(doubles)[grep("^fitness", names(doubles))]){
    j <- gsub('fitness','toxicity',i)
    doubles[, (j) := -.SD[[1]],,.SDcols = i]
  }
  #Combine
  if(read_filter){
    return(rbind(singles, silent, doubles, fill = T)[is.reads0==T,])
  }else{
    return(rbind(singles, silent, doubles, fill = T))
  }
}
