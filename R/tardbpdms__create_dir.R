
#' tardbpdms__create_dir
#'
#' Create results folder for analysis script plots and saved objects.
#'
#' @param tardbpdms_dir directory path string (required)
#' @param execute whether or not the system command will be executed (required)
#' @param message message string (optional, default: NULL i.e. no message displayed)
#' @param overwrite_dir delete directory if already exists (optional, default: TRUE)
#'
#' @return Nothing
#' @export
tardbpdms__create_dir <- function(
  tardbpdms_dir, 
  execute = TRUE, 
  message = NULL, 
  overwrite_dir = TRUE){
  if(!execute){
    return()
  }
  if(!is.null(message)){
    message(paste("\n\n\n*******", message, "*******\n\n\n"))
  }
  if(dir.exists(tardbpdms_dir) & overwrite_dir){
    unlink(tardbpdms_dir, recursive = TRUE)
  }
  dir.create(tardbpdms_dir, showWarnings = FALSE)
  return(tardbpdms_dir)
}
