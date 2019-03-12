
#' tardbpdms__plot_toxicity_distributions
#'
#' Plot toxicity distributions (and optional vertical lines highlighting specific toxicity values).
#'
#' @param input_dt a data.table with columns: STOP, toxicity_cond, region, mut_type (required)
#' @param output_file path to output plot (required)
#' @param boundary_list a list of toxicity values to highlight with vertical lines on plot
#' @param boundary_colour_list a list of colours corresponding the boundary_list vertical lines
#' @param xlim x axis limits
#' @param colour_scheme colour scheme list
#' @param width plot width
#' @param height plot height
#' @param facet_by_mut_type facet by mut_type (default:T)
#'
#' @return Nothing
#' @export
tardbpdms__plot_toxicity_distributions <- function(
  input_dt, 
  output_file, 
  boundary_list=NULL, 
  boundary_colour_list=NULL,
  xlim=NULL,
  colour_scheme=NULL,
  width=6,
  height=13,
  facet_by_mut_type=T
  ){
  plot_dt <- copy(input_dt)[STOP==FALSE & !is.na(toxicity_cond), toxicity := toxicity_cond]
  d <- ggplot2::ggplot(plot_dt, ggplot2::aes(toxicity, ..density.., color = region)) +
    ggplot2::geom_density(data = plot_dt[!is.na(toxicity),]) +
    # ggplot2::geom_density(data = input_dt[STOP==FALSE & !is.na(toxicity_cond),], ggplot2::aes(toxicity_cond, ..density.., color = region)) +
    ggplot2::xlab("Toxicity")
  if(!is.null(boundary_list)){
    for(i in names(boundary_list)){
      d <- d + ggplot2::geom_vline(xintercept = boundary_list[[i]], linetype=2, colour = boundary_colour_list[[i]])
    }
  }
  if(!is.null(xlim)){
    d <- d + ggplot2::coord_cartesian(xlim = xlim)
  }
  d <- d + ggplot2::theme_classic()
  if(facet_by_mut_type){
    d <- d + lemon::facet_rep_grid(mut_type~., scales = "free_y", repeat.tick.labels = 'bottom')
  }
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme))
  }
  ggplot2::ggsave(file=output_file, width=width, height=height)
}
