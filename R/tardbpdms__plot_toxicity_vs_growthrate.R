
#' tardbpdms__plot_toxicity_vs_growthrate
#'
#' Plot toxicity versus growth rate.
#'
#' @param input_dt a data.table with columns: growth_rate, growth_rate_low, growth_rate_high, toxicity, sigma, region (required)
#' @param output_file path to output plot (required)
#' @param colour_scheme colour scheme list
#' @param width plot width
#' @param height plot height
#'
#' @return Nothing
#' @export
tardbpdms__plot_toxicity_vs_growthrate <- function(
  input_dt, 
  output_file,
  colour_scheme=NULL,
  width=4,
  height=3
  ){
  x_axis_lims <- max(abs(unlist(input_dt[!is.na(growth_rate),.(toxicity-1.96*sigma, toxicity+1.96*sigma)])))
  y_axis_lims <- max(abs(unlist(input_dt[!is.na(growth_rate),.(growth_rate_low, growth_rate_high)])))
  temp_cor <- input_dt[!is.na(growth_rate),cor(toxicity, growth_rate)]
  d <- ggplot2::ggplot(input_dt[!is.na(growth_rate),], ggplot2::aes(toxicity, growth_rate, color = region)) +
    ggplot2::geom_abline(slope = -1, linetype = 2) + 
    # ggplot2::geom_hline(yintercept = 0, linetype = 3) + 
    # ggplot2::geom_vline(xintercept = 0, linetype = 3) +
    ggplot2::geom_smooth(method = lm, se = FALSE, size = 1) + 
    ggplot2::geom_linerange(ggplot2::aes(ymin=growth_rate_low, ymax=growth_rate_high), alpha = 0.5) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin=toxicity-1.96*sigma, xmax=toxicity+1.96*sigma), alpha = 0.5) +
    ggplot2::geom_point(size = 1) +
    ggplot2::coord_cartesian(xlim = c(-x_axis_lims, x_axis_lims), ylim = c(-y_axis_lims, y_axis_lims)) +
    ggplot2::ylab("Relative growth rate (log)") +
    ggplot2::xlab("Toxicity") +
    ggplot2::annotate("text", label = paste0("R = ", round(temp_cor, 2)) , x = x_axis_lims-0.2, y = -y_axis_lims, size = 4) + 
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme))
  }
  ggplot2::ggsave(file=output_file, width=width, height=height)
}
