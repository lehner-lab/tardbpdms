
#' tardbpdms__aa_properties
#'
#' Create results folder for analysis script plots and saved objects.
#'
#' @param aa_properties_file path to amino acid properties file (required)
#' @param selected_identifiers path to file with selected subset of identifiers
#'
#' @return results from PCA
#' @export
tardbpdms__aa_properties <- function(
  aa_properties_file,
  selected_identifiers = NULL){
  #Load amino acid properties
  aa_df <- read.table(file=aa_properties_file, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  #Subset to selected amino acid properties if supplied
  if(!is.null(selected_identifiers)){
    aa_df <- aa_df[rownames(aa_df) %in% selected_identifiers | grepl("BBSUPPL", rownames(aa_df)),]
  }
  #Evidences list
  aa_evidences <- as.list(aa_df[,1])
  names(aa_evidences) <- rownames(aa_df)
  #Reformat and add AA identities
  aa_df <- aa_df[,-1]
  #Remove properties with NAs
  aa_evidences <- aa_evidences[apply(is.na(aa_df), 1, sum)==0]
  aa_df <- aa_df[apply(is.na(aa_df), 1, sum)==0,]
  #Normalise
  aa_mat <- scale(t(aa_df))
  return(aa_mat)
}

