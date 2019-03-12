
tau_specificity_score <- function(x, 
  min_length=3){

  ### variables 
  # x: vector of scores
  # min_length: minimum number of scores (otherwise returns NA)

  x <- x[!is.na(x)]
  if(length(x)>=min_length){
    if(max(x)==0){
      return(NA)
    }else{
      return(sum(1-x/max(x))/(length(x)-1))
    }
  }else{
    return(NA)
  }
}

#' tardbpdms__guenther_structure_propensity_plot
#'
#' Plot kernel structure propensity results.
#'
#' @param kernel_scores (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param pdb_names character vector of PDB ids
#' @param score_type name of score
#' @param pdb_list list of lists 
#' @param Nsamples number of randomisations for score significance specificity test
#' @param position_offset offset position
#' @param plot_ylim y-axis limits
#' @param width plot width (default:6)
#' @param height plot height (default:6)
#'
#' @return Nothing
#' @export
#' @import data.table
tardbpdms__guenther_structure_propensity_plot <- function(
	kernel_scores,
  outpath,
  colour_scheme,
  pdb_names,
  score_type,
  pdb_list,
  Nsamples=10000,
  position_offset = 1,
  plot_ylim = c(0, 4.5),
  width = 10,
  height = 5
  ){

  plot_df <- do.call("rbind", lapply(lapply(kernel_scores[pdb_names], as.data.frame), '[', c(paste0(score_type, "_kernel_p"), "Pos")))
  colnames(plot_df) <- c("kernel_p", "Pos")
  #Remove NAs
  plot_df <- plot_df[!is.na(plot_df$kernel_p),]
  plot_df$neglogp <- (-log10(plot_df$kernel_p))
  plot_df$pdb_structure <- sapply(strsplit(rownames(plot_df), "\\."), '[', 1)
  #Start position
  plot_df$start_position <- unlist(sapply(pdb_list[plot_df$pdb_structure], '[', "idx_start"))
  #Reformat pdb structure name
  plot_df$pdb_structure <- paste0(plot_df$pdb_structure, ": ", 
    unlist(sapply(pdb_list[plot_df$pdb_structure], '[', "idx_start")), "-",
    unlist(sapply(pdb_list[plot_df$pdb_structure], '[', "aa_seq")), "-",
    unlist(sapply(pdb_list[plot_df$pdb_structure], '[', "idx_start")) + nchar(unlist(sapply(pdb_list[plot_df$pdb_structure], '[', "aa_seq")))-1)
  plot_df$position <- plot_df$Pos+position_offset-1
  #Subset to positions with data for all kernels
  pos_tab <- table(plot_df$position)
  plot_df <- plot_df[plot_df$position %in% names(pos_tab)[pos_tab==length(pdb_names)],]
  #Kernel smooth score significance specificity
  tspec_obs <- tapply(1-plot_df$kernel_p, plot_df$position, tau_specificity_score)
  tspec_rand <- list()
  #Kernel smooth score significance specificity - randomisation test
  kernel_p_rand_list <-  list()
  plot_dt <- data.table(plot_df)
  set.seed(1)
  for(i in unique(plot_dt[,pdb_structure])){
    kernel_p_rand_list[[i]] = matrix(sample(plot_dt[pdb_structure==i,c(kernel_p)],(plot_dt[pdb_structure==i,.N])*Nsamples,replace = T),nrow = plot_dt[pdb_structure==i,.N],ncol=Nsamples)
  }
  kernel_p_rand <- do.call("rbind", kernel_p_rand_list)
  tspec_rand <- apply(kernel_p_rand, 2, function(x){tapply(1-x, plot_df$position, tau_specificity_score)})
  tspec_sig <- (apply((tspec_rand - matrix(tspec_obs, nrow = length(tspec_obs), ncol = Nsamples, byrow = F))>=0, 1, sum)+1)/(Nsamples+1)
  #Plot
  plot_df_spec <- data.frame(
  	kernel_p = tspec_sig, 
  	Pos = NA, 
  	neglogp = -log10(tspec_sig), 
  	pdb_structure = "Kernel specificity", 
  	position = as.numeric(names(tspec_sig)), 
  	start_position = NA, 
  	test = "Kernel specificity", stringsAsFactors = F)
  plot_df$test <- "Kernel score"
  plot_df <- rbind(plot_df, plot_df_spec)
  #Convert to factor
  plot_df$pdb_structure <- factor(plot_df$pdb_structure, levels=unique(plot_df$pdb_structure[order(plot_df$start_position, decreasing = F)]))
  #Data frame for start positions
  plot_df_start <- plot_df[!duplicated(plot_df$start_position),]
  plot_df_start <- plot_df_start[rep(1:dim(plot_df_start)[1], each = 2),]
  plot_df_start$test <- rep(unique(plot_df_start$test), dim(plot_df_start)[1]/2)
  #Plot with significance of specificity
  plot_cols <- unlist(colour_scheme, use.names = FALSE)[1:length(levels(plot_df[,"pdb_structure"]))]
  names(plot_cols)[1:length(levels(plot_df[,"pdb_structure"]))] <- levels(plot_df[,"pdb_structure"])
  plot_cols["Kernel specificity"] <- "black"
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(position, neglogp, color = pdb_structure)) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) + 
    ggplot2::geom_vline(data = plot_df_start, ggplot2::aes(xintercept = start_position, color = pdb_structure), linetype = 2) + 
    ggplot2::ylab("Significance, -log10(P-value)") +
    ggplot2::xlab("Start position") +
    ggplot2::coord_cartesian(ylim = plot_ylim) + 
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::scale_colour_manual(values = plot_cols)
  d <- d + ggplot2::facet_grid(test~.)
  suppressWarnings(suppressMessages(ggplot2::ggsave(file=outpath, width=width, height=height)))
}
