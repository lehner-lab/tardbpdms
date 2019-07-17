
#' tardbpdms_human_disease_mutations
#'
#' Check fitness bias of human disease mutations.
#'
#' @param toxicity_dt data.table with single mutant toxicity values (required)
#' @param outpath output path for plots and saved objects (required)
#' @param missense_AF_file table of missense mutation allele frequencies (required)
#' @param disease_mut_file table of human disease mutations and classifications (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
tardbpdms_human_disease_mutations <- function(
  toxicity_dt,
  missense_AF_file,
  disease_mut_file,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

	#Return previous results if analysis not executed
	if(!execute){
		return()
	}

  #Display status
  message(paste("\n\n*******", "running stage: tardbpdms_human_disease_mutations", "*******\n\n"))

	#Create output directory
	tardbpdms__create_dir(tardbpdms_dir = outpath)

	#Single AA mutants only (no STOPs)
	singles_dt <- copy(toxicity_dt[Nmut_aa==1])
	#Double AA mutants only (no STOPs)
	doubles_dt <- copy(toxicity_dt[Nmut_aa==2])
	doubles_dt[, toxicity := toxicity_cond]
	tox_dt <- rbind(singles_dt, doubles_dt, fill = T)

	#Toxicity hotspot positions
	mean_toxicity <- singles_dt[STOP==F,mean(abs(toxicity))]
	Pos_abs_hotspot <- singles_dt[STOP==F,.(hotspot = mean(abs(toxicity))>mean_toxicity),by=Pos_abs][hotspot==T,Pos_abs]

	#Disease mutations
	dis_mut <- read.table(disease_mut_file, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
	fALS_muts <- rownames(dis_mut)[!is.na(dis_mut$fALS)]
	sALS_muts <- rownames(dis_mut)[!is.na(dis_mut$sALS)]
	num_patients <- apply(dis_mut[,c("fALS", "sALS")], 1, sum, na.rm = T)

	#AA code translation dict
	aa_obj <- Biostrings::AAString("GAVLMIFYWKRHDESTCNQP")
	aa_list <- Biostrings::AMINO_ACID_CODE[strsplit(as.character(aa_obj), NULL)[[1]]]
	aa_list["*"] <- "X"
	aa_list_rev <- names(aa_list)
	names(aa_list_rev) <- aa_list

	#Detected human missense mutations
	miss_mut <- as.data.frame(fread(missense_AF_file))
	miss_mut <- miss_mut[nchar(miss_mut[,1])==11,]
	miss_mut[,1] <- gsub("^p.", "", miss_mut[,1])
	rownames(miss_mut) <- paste0(aa_list_rev[substr(miss_mut[,1], 1, 3)], substr(miss_mut[,1], 4, 6), aa_list_rev[substr(miss_mut[,1], 7, 9)])

	#Add mutation information
	tox_dt[, hmut_cat := "Never observed (gnomAD)"]
	tox_dt[!is.na(miss_mut[mut_code,2]), hmut_cat := "Observed (gnomAD)"]
	tox_dt[mut_code %in% rownames(dis_mut), hmut_cat := "other"]
	tox_dt[mut_code %in% fALS_muts, hmut_cat := "fALS"]
	tox_dt[mut_code %in% sALS_muts, hmut_cat := "sALS"]
	tox_dt[mut_code %in% sALS_muts & mut_code %in% fALS_muts, hmut_cat := "fALS and sALS"]
	tox_dt[, hmut_recurrent := FALSE]
	tox_dt[mut_code %in% names(num_patients)[num_patients>1], hmut_recurrent := TRUE]

	fwrite(tox_dt[mut_code %in% rownames(dis_mut),], file = "dis_mut_tab.tsv", sep = "\t")

	### Toxicity bias of TDP-43 mutations
	###########################

	#Toxicity histogram of all single and double observed versus unobserved human missense mutations
	tox_bias_all <- tox_dt[hmut_cat %in% c("fALS", "sALS", "fALS and sALS", "other"),.(pvalue = wilcox.test(toxicity)$p.value, n = .N)]
	tox_bias_fALSr <- tox_dt[hmut_cat %in% c("fALS") & hmut_recurrent==T,.(pvalue = wilcox.test(toxicity)$p.value, n = .N)]
	set.seed(1)
	plot_df <- as.data.frame(tox_dt[,c("toxicity", "Nmut_aa", "STOP", "hmut_cat", "hmut_recurrent")])
	plot_df[,"hmut_cat"] <- factor(plot_df[,"hmut_cat"], levels = c("fALS", "sALS", "fALS and sALS", "other", "Never observed (gnomAD)", "Observed (gnomAD)"))
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(toxicity, ..density..)) +
    ggplot2::geom_density() +
    ggplot2::geom_jitter(data = plot_df[!grepl("gnomAD", plot_df[,"hmut_cat"]),], ggplot2::aes(x = toxicity, y = toxicity*0-1, color = hmut_cat, shape = hmut_recurrent)) +
    ggplot2::xlab("Toxicity") +
    ggplot2::geom_vline(xintercept = 0, linetype=2) +
    ggplot2::geom_vline(xintercept = tox_dt[Nmut_aa==1 & STOP==T,median(toxicity)], linetype=2, colour = "darkgrey") +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][1:4])) +
    ggplot2::scale_shape_manual(values = c(1, 19)) + 
    ggplot2::annotate("text", label = paste0("P-value (all) = ", round(tox_bias_all[,"pvalue"], 3), " (", tox_bias_all[,"n"], ")") , x = 0.2, y = -0.5) +
    ggplot2::annotate("text", label = paste0("P-value (fALSr) = ", round(tox_bias_fALSr[,"pvalue"], 3), " (", tox_bias_fALSr[,"n"], ")") , x = 0.2, y = -1)
  ggplot2::ggsave(file=file.path(outpath, 'human_disease_mut_toxicity.pdf'), width=5, height=3)

	#Test human mutation toxicity bias
  # wilcox.test(tox_dt[hmut_cat %in% c("fALS", "sALS", "fALS and sALS"),toxicity])
  # wilcox.test(tox_dt[hmut_cat=="fALS" & hmut_recurrent==T,toxicity])
  # wilcox.test(tox_dt[hmut_cat=="fALS" & hmut_recurrent==F,toxicity])
  # wilcox.test(tox_dt[hmut_cat=="sALS" & hmut_recurrent==T,toxicity])
  # wilcox.test(tox_dt[hmut_cat=="sALS" & hmut_recurrent==F,toxicity])
  # wilcox.test(tox_dt[hmut_cat=="fALS and sALS" & hmut_recurrent==T,toxicity])
  # wilcox.test(tox_dt[hmut_cat=="fALS and sALS" & hmut_recurrent==F,toxicity])
  # wilcox.test(tox_dt[hmut_cat %in% c("fALS", "sALS", "fALS and sALS") & !(hmut_cat=="fALS" & hmut_recurrent==T),toxicity])
}

