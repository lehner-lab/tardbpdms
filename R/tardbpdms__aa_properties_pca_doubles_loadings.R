
#' tardbpdms__aa_properties_pca_doubles_loadings
#'
#' Create results folder for analysis script plots and saved objects.
#'
#' @param input_dt data.table with double mutant codes (required)
#' @param aa_properties_file path to amino acid properties file (required)
#' @param selected_identifiers path to file with selected subset of identifiers
#'
#' @return results from PCA
#' @export
tardbpdms__aa_properties_pca_doubles_loadings <- function(
  input_dt,
  aa_properties_file,
  selected_identifiers = NULL
  ){
  #AA properties PCA
  exp_pca <- tardbpdms__aa_properties_pca(aa_properties_file = aa_properties_file, selected_identifiers = selected_identifiers)
  #Mut codes
  mut_code1 <- input_dt[,paste0(WT_AA1, Mut1)]
  mut_code2 <- input_dt[,paste0(WT_AA2, Mut2)]
  #List of mutation effects on AA properties
  mut_list <- list()
  for(i in unique(c(mut_code1, mut_code2))){
    mut_list[[i]] <- exp_pca$x[substr(i, 2, 2),] - exp_pca$x[substr(i, 1, 1),]
  }
  #Double mutant effects on AA properties (sum)
  doubles_daaprop <- data.table::as.data.table(do.call("rbind", mut_list[mut_code1]) + do.call("rbind", mut_list[mut_code2]))
  names(doubles_daaprop) <- paste0(names(doubles_daaprop), "_score_sum")
  #Append to input data table
  return(cbind(input_dt, doubles_daaprop))
}
