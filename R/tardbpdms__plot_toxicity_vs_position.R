
#' tardbpdms__plot_toxicity_vs_position
#'
#' Plot toxicity versus position.
#'
#' @param input_dt a data.table with columns: toxicity, sigma, region, Pos (required)
#' @param output_file path to output plot (required)
#' @param ylab y-axis title (required)
#' @param use_sigma whether to use error in 'sigma' column (default:F i.e. error estimated from toxicity distribution)  
#' @param span loess smoother span
#' @param label_offset annotation label offset
#' @param plot_rect whether to plot
#' @param plot_hotspot whether to plot hotspot
#'
#' @return Nothing
#' @export
#' @import data.table
tardbpdms__plot_toxicity_vs_position <- function(
  input_dt, 
  output_file, 
  ylab, 
  use_sigma = F, 
  span = 0.75, 
  label_offset = 0.05, 
  plot_rect = F,
  plot_hotspot = F
  ){
  plot_dt <- input_dt
  plot_dt[, Pos := Pos+as.numeric(region)-1]
  plot_dt_orig <- copy(plot_dt)
  if(use_sigma){
    plot_dt[, toxicity_low := toxicity - 1.96*sigma]
    plot_dt[, toxicity_high := toxicity + 1.96*sigma]
  }
  else{
    plot_dt[, Pos := factor(Pos, levels = unique(Pos[order(Pos, decreasing = F)]))]
    plot_dt_summarise <- data.table::as.data.table(plyr::ddply(plot_dt, plyr::.(Pos), plyr::summarize, toxicity = mean(toxicity)))
    plot_dt_summarise[, toxicity_low := plyr::ddply(plot_dt, plyr::.(Pos), plyr::summarize, toxicity_low = mean(toxicity)-1.96*sd(toxicity)/sqrt(length(toxicity)))[,"toxicity_low"]]
    plot_dt_summarise[, toxicity_high := plyr::ddply(plot_dt, plyr::.(Pos), plyr::summarize, toxicity_high = mean(toxicity)+1.96*sd(toxicity)/sqrt(length(toxicity)))[,"toxicity_high"]]
    plot_dt_summarise[, Pos := as.numeric(as.character(Pos))]
    plot_dt <- plot_dt_summarise
  }
  d <- ggplot2::ggplot(plot_dt, ggplot2::aes(Pos, toxicity))
  #Plot Eisenberg segments?
  if(plot_rect){
    d <- d +
      ggplot2::geom_rect(ggplot2::aes(xmin=312-0.5, xmax=317+0.5, ymin=min(plot_dt[,toxicity_low]), ymax=plot_dt[,max(toxicity_high)]), fill = "grey", linetype = 3, colour = "darkgrey", alpha = 0) +
      ggplot2::geom_rect(ggplot2::aes(xmin=300-0.5, xmax=306+0.5, ymin=min(plot_dt[,toxicity_low]), ymax=plot_dt[,max(toxicity_high)]), fill = "#eeeeee", linetype = 3, colour = "darkgrey", alpha = 0) +
      ggplot2::geom_rect(ggplot2::aes(xmin=321-0.5, xmax=326+0.5, ymin=min(plot_dt[,toxicity_low]), ymax=plot_dt[,max(toxicity_high)]), fill = "#eeeeee", linetype = 3, colour = "darkgrey", alpha = 0) +
      ggplot2::geom_rect(ggplot2::aes(xmin=328-0.5, xmax=333+0.5, ymin=min(plot_dt[,toxicity_low]), ymax=plot_dt[,max(toxicity_high)]), fill = "#eeeeee", linetype = 3, colour = "darkgrey", alpha = 0) +
      ggplot2::geom_rect(ggplot2::aes(xmin=333-0.5, xmax=343+0.5, ymin=min(plot_dt[,toxicity_low]), ymax=plot_dt[,max(toxicity_high)]), fill = "#eeeeee", linetype = 3, colour = "darkgrey", alpha = 0) +
      ggplot2::geom_rect(ggplot2::aes(xmin=370-0.5, xmax=375+0.5, ymin=min(plot_dt[,toxicity_low]), ymax=plot_dt[,max(toxicity_high)]), fill = "#eeeeee", linetype = 3, colour = "darkgrey", alpha = 0) +
      # ggplot2::geom_rect(ggplot2::aes(xmin=396-0.5, xmax=402+0.5, ymin=min(plot_dt$toxicity_low), ymax=max(plot_dt$toxicity_high)), fill = "#eeeeee", linetype = 1, colour = "darkgrey", alpha = 0) +
      ggplot2::annotate("text", label = "LARKS", x = (312+317)/2, y = plot_dt[,max(toxicity_high)]+label_offset, size = 4) + 
      ggplot2::annotate("text", label = "SZ", x = (300+306)/2, y = plot_dt[,max(toxicity_high)]+label_offset, size = 4) + 
      ggplot2::annotate("text", label = "SZ", x = (321+326)/2, y = plot_dt[,max(toxicity_high)]+label_offset, size = 4) + 
      ggplot2::annotate("text", label = "SZ", x = (328+333)/2, y = plot_dt[,max(toxicity_high)]+label_offset, size = 4) + 
      ggplot2::annotate("text", label = "SZ", x = (333+343)/2, y = plot_dt[,max(toxicity_high)]+label_offset, size = 4) + 
      ggplot2::annotate("text", label = "SZ", x = (370+375)/2, y = plot_dt[,max(toxicity_high)]+label_offset, size = 4)
  }
  #Plot toxicity hotspot?
  if(plot_hotspot){
    mean_toxicity <- plot_dt_orig[,mean(toxicity)]
    pos_hotspot <- plot_dt_orig[,.(hotspot = mean(toxicity)>mean_toxicity),by=Pos][hotspot==T,Pos]
    d <- d + 
      ggplot2::geom_rect(
        ggplot2::aes(xmin=min(pos_hotspot)-0.5, xmax=max(pos_hotspot)+0.5, ymin=min(plot_dt[,toxicity_low]), ymax=plot_dt[,max(toxicity_high)]), 
        fill = "grey", 
        linetype = 0) +
      ggplot2::geom_hline(yintercept = mean_toxicity, linetype = 2)
  }
  d <- d + ggplot2::geom_linerange(ggplot2::aes(ymin=toxicity_low, ymax=toxicity_high)) +
    ggplot2::geom_smooth(data = plot_dt_orig, method = "loess", span = span, colour = "black", se = F) +
    ggplot2::geom_point() +
    # ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    # ggplot2::geom_smooth(method = lm) +
    ggplot2::geom_vline(xintercept = 332-0.5, linetype = 2) +
    ggplot2::ylab(ylab) +
    ggplot2::xlab("TDP-43 amino acid position") +
    ggplot2::theme_classic() + 
    ggplot2::coord_cartesian(xlim = c(min(plot_dt[,Pos]), max(plot_dt[,Pos])))
  d
  ggplot2::ggsave(file=output_file, width=12, height=3)
}
