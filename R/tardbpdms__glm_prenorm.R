
#' tardbpdms__glm_prenorm
#'
#' Normalise data.frame columns.
#'
#' @param x data.frame to be normalised (required)
#'
#' @return A normalised data.frame
#' @export
tardbpdms__glm_prenorm <- function(
  x){
  for(i in 1:dim(x)[2]){
    if(length(unique(na.omit(x[,i])))>2 & !is.factor(x[,i])){
      x[,i]<-scale(x[,i])
    }
  }
  x
}
