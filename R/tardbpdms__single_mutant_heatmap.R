
#' tardbpdms__single_mutant_heatmap
#'
#' ggplot tile heatmap wrapper.
#'
#' @param input_df matrix of heatmap values (required)
#' @param variable_name name of variable to use for heatmap cells (defaut:"fitness")
#' @param mut_dict list of mutant codes where name is used to label the corresponding cells
#' @param text_size text size (default:2.5)
#' @param output_file path to output file (required)
#' @param colour_low colour scale lower limit colour -- passed to scale_colour_gradient2 (default:"blue")
#' @param colour_high colour scale upper limit colour -- passed to scale_colour_gradient2 (default:"red")
#' @param width plot width in "units" (default:12)
#' @param height plot height in "units" (default:4)
#' @param x_annotation x-axis title (default:"position")
#' @param xaxis_angle rotation angle for x tick labels -- passed to element_text (default:90)
#' @param na_colour colour to use for NA values (default:"grey50")
#' @param na_text cell text annotation to use for NA values (default:empty string)
#'
#' @return Nothing
#' @export
#' @import data.table
tardbpdms__single_mutant_heatmap<-function(
  input_df, 
  variable_name="fitness", 
  mut_dict=NULL,#=list("F"=fALS_muts , "S"=sALS_muts), 
  text_size=2.5, 
  output_file, 
  colour_low='blue', 
  colour_high='red', 
  width=12, 
  height=4, 
  x_annotation="position", 
  xaxis_angle=90, 
  na_colour="grey50", 
  na_text=""){

  aa_obj <- Biostrings::AAString("GAVLMIFYWKRHDESTCNQP")
  aa_list <- Biostrings::AMINO_ACID_CODE[strsplit(as.character(aa_obj), NULL)[[1]]]
  aa_list["*"] <- "X"

  #Only single AA mutants
  dms_1aa <- input_df[input_df$Nmut_aa==1,]
  #Absolute position
  dms_1aa$Pos_abs <- as.numeric(substr(dms_1aa$mut_code, 2, 4))
  #WT sequence
  wt_seq <- unique(dms_1aa[order(dms_1aa[,"Pos_abs"]),c("WT_AA", "Pos_abs")])[,"WT_AA"]

  #Construct heatmap matrix
  heat_mat <- matrix(nrow = length(aa_list), ncol = max(dms_1aa$Pos_abs)-min(dms_1aa$Pos_abs)+1)
  rownames(heat_mat) <- names(aa_list)
  colnames(heat_mat) <- min(dms_1aa$Pos_abs):max(dms_1aa$Pos_abs)
  for(aa_pos in min(dms_1aa$Pos_abs):max(dms_1aa$Pos_abs)){
    for(aa_id in names(aa_list)){
      temp_index <- which(dms_1aa$Pos_abs==aa_pos & dms_1aa$Mut==aa_id)
      if(length(temp_index)==1){
        heat_mat[aa_id,as.character(aa_pos)] <- dms_1aa[temp_index,variable_name]
      }
    }
  }
  #Add reported fALS mutations
  heat_mat_text <- matrix("", nrow = length(aa_list), ncol = max(dms_1aa$Pos_abs)-min(dms_1aa$Pos_abs)+1)
  rownames(heat_mat_text) <- names(aa_list)
  colnames(heat_mat_text) <- min(dms_1aa$Pos_abs):max(dms_1aa$Pos_abs)
  #List of mutants to label supplied?
  if(!is.null(mut_dict)){
    for(mut_letter in names(mut_dict)){
      for(mut in mut_dict[[mut_letter]]){
        if(substr(mut,2,4) %in% colnames(heat_mat_text)){
          heat_mat_text[substr(mut,5,5),substr(mut,2,4)] <- paste0(heat_mat_text[substr(mut,5,5),substr(mut,2,4)], mut_letter)
        }
      }
    }
  }
  #Text alternative for NA values supplied?
  if(na_text!=""){
    heat_mat_text[is.na(heat_mat)] <- na_text
  }
  xtick_labels <- NULL
  if(x_annotation == "both"){
    xtick_labels <- paste0(colnames(heat_mat), ":", wt_seq)
  }else if(x_annotation == "sequence"){
    xtick_labels <- wt_seq
  }else if(x_annotation == "position"){
    xtick_labels <- colnames(heat_mat)
  }
  colnames(heat_mat) <- paste0(":",colnames(heat_mat))
  #Plot
  tardbpdms__tile_heatmap_wrapper(
    heat_mat, 
    input_matrix_text = heat_mat_text, 
    text_size = text_size, 
    output_file = output_file, 
    colour_low = colour_low, 
    colour_high = colour_high, 
    xlab = "TDP-43 amino acid position", 
    ylab = "Mutant AA", 
    cluster = "none", 
    width = width, 
    height = height, 
    xaxis_angle = xaxis_angle, 
    xaxis_hjust = 0.5, 
    xaxis_size = NULL, 
    xtick_labels = xtick_labels, 
    na_colour = na_colour)
}
