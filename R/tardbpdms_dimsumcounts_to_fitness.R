
#' tardbpdms_dimsumcounts_to_fitness
#'
#' Estimate fitness of singles, doubles and silent AA substitution variants from DiMSum formatted output counts.
#'
#' @param counts_path path to DiMSum output (required)
#' @param outpath output path for plots and saved objects (required)
#' @param generations_path path to generations data (required)
#' @param all_reps list of replicates to retain (required)
#' @param min_mean_input_read_count minimum mean input read count for high confidence variants (default:10)
#' @param min_input_read_count_doubles minimum input read count for doubles used to derive prior for Bayesian doubles correction (default:50)
#' @param lam_d Poisson distribution for score likelihood (default:0.025)
#' @param numCores Number of available CPU cores (default:1)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
tardbpdms_dimsumcounts_to_fitness <- function(
  counts_path,
  outpath,
  generations_path,
  all_reps,
  min_mean_input_read_count = 10,
  min_input_read_count_doubles = 50,
  lam_d = 0.025,
  numCores = 1,
  execute = TRUE
  ){

  #Do nothing if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: tardbpdms_dimsumcounts_to_fitness", "*******\n\n"))

  #Create output directory
  tardbpdms__create_dir(tardbpdms_dir = outpath)

  ### Minimum input read count in any replicate (to remove very low confidence variants)
  ###########################

  #Minimum input read count
  min_input_read_count <- 5

  ### Load generations data
  ###########################

  generations <- read.table(generations_path, sep = "\t", row.names = 1, header = T)

  ### Load variant count data
  ###########################

  #load variant data from RData file
  e1 <- new.env() 
  load(counts_path)
  all_data <- get('variant_data_merge', e1)
  rm(e1)

  #Retain nucleotide variants with same length as WT and max 2 amino acids only
  all_data <- all_data[nchar(nt_seq)==nchar(all_data[WT==T,nt_seq]) & Nmut_aa<=2,]

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")

  #WT nucleotide sequences
  wt_ntseq <- all_data[WT==T,nt_seq]
  wt_ntseq_split <- strsplit(wt_ntseq,"")[[1]]
  #WT AA sequences
  wt_AAseq <- all_data[WT==T,aa_seq]
  wt_AAseq_split <- strsplit(wt_AAseq,"")[[1]]

  #Sample names
  input_samples <- colnames(all_data)[grep('^input', colnames(all_data))]
  output_samples <- colnames(all_data)[grep('^output', colnames(all_data))]

  #Add number of codons affected by mutations
  all_data[,Nmut_codons := length(unique(ceiling(which(strsplit(nt_seq,"")[[1]] != wt_ntseq_split)/3))),nt_seq]
  #Table of read counts stratified by AA and codon mutations
  all_data[,.N,.(Nmut_aa,Nmut_codons)][order(Nmut_aa,Nmut_codons)]
  #Table of read counts stratified by nucleotide mutations
  all_data[,.N,.(Nmut_nt)][order(Nmut_nt)]
  #Example histogram of input1 counts split by number of nucleotide mutations 
  d <- ggplot2::ggplot(all_data[between(Nmut_nt,1,4)  & input1_e1_s0_bNA_count > 0],
         ggplot2::aes(input1_e1_s0_bNA_count,color=factor(Nmut_nt),linetype = Nmut_codons > 0)) +
    ggplot2::geom_density() +
    ggplot2::scale_x_log10()
  ggplot2::ggsave(file.path(outpath, "1_input1_e1_s0_bNA_count_hist.pdf"), d, width = 7, height = 5)

  ### Retain only variants with >=min_input_read_count input readcounts in any biological replicate 
  ### and either purely nonsynonymous mutations
  ### or purely silent mutations
  ###########################

  #### only retain variants with at least 5 input readcounts in any of the biological replicates
  dataset <- copy(all_data[((Nmut_codons-Nmut_aa) == 0 | Nmut_aa == 0) & rowSums(all_data[,input_samples,with=F]>min_input_read_count) != 0])
  dataset[,.N,.(Nmut_nt,Nmut_aa,Nmut_codons)][order(Nmut_nt,Nmut_aa,Nmut_codons)]

  #Plot pairwise input sample count correlations for all single mutants (non-synonymous only)
  set.seed(1)
  d <- GGally::ggpairs(log10(dataset[Nmut_nt==1][sample(x = .N,min(c(.N,10000))),grep(names(dataset),pattern="input"),with=F]+1),
          upper=list(continuous = "cor"))
  ggplot2::ggsave(file.path(outpath, "2_scatterplotmatrix_singlemutants.pdf"), d, width = 10, height = 10)
  
  #Plot pairwise input sample count correlations for random sample of 10k variants (non-synonymous only)
  temp <- dataset[sample(x = .N,10000)]
  d <- GGally::ggpairs(cbind(log10(temp[,grep(names(dataset),pattern="input"),with=F]+1), Nmut_nt=as.factor(temp[,Nmut_nt])),
    columns = 1:length(grep(names(dataset),pattern="input")),
    mapping = ggplot2::aes(color = Nmut_nt),
    upper=list(continuous = "cor"))
  ggplot2::ggsave(file.path(outpath, "2_scatterplotmatrix_random10kmutants.pdf"), d, width = 10, height = 10)

  ### Aggregate counts from variants that are identical at the AA level and without synonymous mutations
  ###########################

  #Add merge_seq (aa_seq if nonsynonymous mutations, otherwise nt_seq)
  dataset[,merge_seq := aa_seq,nt_seq]
  dataset[Nmut_aa==0,merge_seq := nt_seq,nt_seq]

  #For all count columns
  idx = names(dataset)[grep(names(dataset),pattern="_count$")]
  for (i in seq_along(idx)) {
    #Aggregate counts accross identical AA variants
    dataset[,paste0(idx[i],"_agg") := sum(.SD),merge_seq,.SDcols = idx[i]]
  }
  #Retain only one row per AA variant
  dataset_syn = dataset[!duplicated(merge_seq),.SD,merge_seq,.SDcols = c("aa_seq","Nmut_nt","Nmut_aa","Nmut_codons","WT","STOP",names(dataset)[grep(names(dataset),pattern="_agg$")])]

  ### Aggregate counts for biological output replicates
  ###########################

  #Add up counts for biological output reps
  for (E in all_reps) {
    idx = names(dataset_syn)[grep(names(dataset_syn),pattern = paste0("e",E,"_s1_b"))]
    dataset_syn[,paste0("count_e",E,"_s1") := rowSums(.SD),,.SDcols = idx]
    names(dataset_syn)[grep(names(dataset_syn),pattern = paste0("input",E))] = paste0("count_e",E,"_s0")
  }
  dataset_syn = unique(dataset_syn[,.SD,merge_seq,.SDcols = c("aa_seq","Nmut_nt","Nmut_aa","Nmut_codons","WT","STOP",names(dataset_syn)[grep(names(dataset_syn),pattern="^count")])])

  ### Calculate fitness and count-based error
  ###########################

  #For all input replicates
  for (E in all_reps) {
    f_wt_corr = dataset_syn[WT == T,
                            log(.SD[,1]/.SD[,2]),,
                            .SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
    dataset_syn[,paste0("fitness",E,"_uncorr") := log(unlist(.SD[,1]) / unlist(.SD[,2])) - rep(t(f_wt_corr),nrow(dataset_syn)),
            ,.SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
    s_wt_counts = dataset_syn[WT == T,
                              1/.SD[,1] + 1/.SD[,2],,
                              .SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
    dataset_syn[,paste0("sigma",E,"_uncorr") := sqrt(1/unlist(.SD[,1]) + 1/unlist(.SD[,2]) + rep(t(s_wt_counts),nrow(dataset_syn))),
            ,.SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
  }
  #Remove unnecessary columns
  dataset_syn = dataset_syn[,.SD,merge_seq,.SDcols = c("aa_seq","Nmut_nt","Nmut_aa","Nmut_codons","WT","STOP",names(dataset_syn)[grep(names(dataset_syn),pattern="^count|^fitness|^sigma")])]

  #Result table with fitness and count-based error for all input replicates
  DT = NULL
  for (E in all_reps) {
    DT = rbind(DT, dataset_syn[,.(count_in = get(paste0("count_e", E, "_s0")),sigma = get(paste0("sigma", E, "_uncorr")),fitness = get(paste0("fitness", E, "_uncorr")),rep = paste0("rep", E),Nmut_aa)])
  }

  #Fitness density plots
  d <- ggplot2::ggplot(DT[Nmut_aa != 0],ggplot2::aes(fitness,color=factor(Nmut_aa),linetype = count_in > 10)) +
    ggplot2::geom_density()
  ggplot2::ggsave(file.path(outpath, "3_fitness_densities.pdf"), d, width = 7, height = 5)

  #Fitness vs. input count hexagonal heatmap
  d <- ggplot2::ggplot(DT,ggplot2::aes(count_in,fitness)) +
    ggplot2::geom_hex() +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap( ~ rep)
  ggplot2::ggsave(file.path(outpath, "3_fitness_count_hexbin.pdf"), d, width = 7, height = 5)

  #Sigma vs. fitness  hexagonal heatmap
  d <- ggplot2::ggplot(DT,ggplot2::aes(fitness,sigma)) +
    ggplot2::geom_hex() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap( ~ rep)
  ggplot2::ggsave(file.path(outpath, "3_sigma_fitness_hexbin.pdf"), d, width = 7, height = 5)

  ### Identify position and identity of single AA mutations (and all silent mutants)
  ###########################

  #WT
  wildtype = dataset_syn[WT==TRUE,]

  #Single AA mutants and all silent mutants
  # singles_silent = dataset_syn[Nmut_aa==1 | (is.na(WT) & Nmut_aa==0),]
  singles_silent = dataset_syn[!is.na(apply(dataset_syn[,.SD,,.SDcols=paste0("fitness", all_reps, "_uncorr")], 1, sum)) & (Nmut_aa==1 | (is.na(WT) & Nmut_aa==0)),]
  #Add position, mutant AA, WT AA and mean input count
  singles_silent[,Pos := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split),aa_seq]
  singles_silent[,Mut := strsplit(aa_seq,"")[[1]][Pos],aa_seq]
  singles_silent[,WT_AA := wt_AAseq_split[Pos],aa_seq]
  #Mean counts
  singles_silent[,mean_count := rowMeans(.SD),,.SDcols = paste0("count_e", all_reps, "_s0")]
  singles_silent[mean_count >= min_mean_input_read_count,is.reads0 := TRUE]

  #Mean input count 
  d <- ggplot2::ggplot(singles_silent[Nmut_aa==1],ggplot2::aes(mean_count)) + 
    ggplot2::geom_density() + 
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept = 400)
  ggplot2::ggsave(file.path(outpath, "4_singles_meaninputcount.pdf"), d, width = 7, height = 5)
  #lower peak is from singles with AA changes 2 nt away

  #Remove unnecessary columns and rename fitness and sigma columns
  singles_silent = singles_silent[,cbind(Pos,WT_AA,Mut,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,.SD),,.SDcols = grep("_uncorr", names(singles_silent))]
  names(singles_silent) <- gsub("_uncorr", "", names(singles_silent))
  singles_silent[,.N,.(Nmut_nt,Nmut_aa,Nmut_codons)][order(Nmut_nt,Nmut_aa,Nmut_codons)]

  ### Identify position and identity of double AA mutations
  ###########################

  #Double AA mutants
  doubles = dataset_syn[Nmut_aa==2]
  #Add position, mutant AA, WT AA and mean input count
  doubles[,Pos1 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[1],aa_seq]
  doubles[,Pos2 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[2],aa_seq]
  doubles[,Mut1 := strsplit(aa_seq,"")[[1]][Pos1],aa_seq]
  doubles[,Mut2 := strsplit(aa_seq,"")[[1]][Pos2],aa_seq]
  doubles[,WT_AA1 := wt_AAseq_split[Pos1],aa_seq]
  doubles[,WT_AA2 := wt_AAseq_split[Pos2],aa_seq]
  #Mean counts
  doubles = merge(doubles,singles_silent[,.(Pos,Mut,s1_mean_count = mean_count)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
  doubles = merge(doubles,singles_silent[,.(Pos,Mut,s2_mean_count = mean_count)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))
  doubles[,mean_count := rowMeans(.SD),,.SDcols = paste0("count_e", all_reps, "_s0")]

  doubles[,.N,s1_mean_count > 400 & s2_mean_count > 400]
  #only few doubles that have individual amino acid changes 2nt away from wild-type
  # for epistasis analysis only AAchanges 1nt away from wildtype will be useful

  #Mean input count of doubles split accoring to mean input count of singles
  d <- ggplot2::ggplot(doubles,ggplot2::aes(mean_count,..count..,color = s1_mean_count > 400 & s2_mean_count > 400)) +
    ggplot2::geom_density() + 
    ggplot2::scale_x_log10() +
    ggplot2::labs(color="high single input")
  ggplot2::ggsave(file.path(outpath, "4_doubles_meaninputcount.pdf"), d, width = 7, height = 5)

  ### Bayesian framework for fitness estimation for double mutants
  ###########################

  #Bin mean counts for replicate 1
  doubles[,counts_for_bins := .SD[[1]],,.SDcols = paste0("count_e",all_reps[1],"_s0")]
  doubles[,bin_count := findInterval(log10(counts_for_bins),seq(0.5,4,0.25))]
  doubles[,.(.N,mean(counts_for_bins)),bin_count][order(bin_count)]

  #Plot fitness densities for different mean count bins (replicate 1)
  doubles[,fitness_for_bins := .SD[[1]],,.SDcols = paste0("fitness",all_reps[1],"_uncorr")]
  d <- ggplot2::ggplot(doubles[between(bin_count,2,8)],ggplot2::aes(fitness_for_bins,..density..,color=factor(bin_count))) +
    ggplot2::geom_density(adjust=1)
  ggplot2::ggsave(file.path(outpath, "4_doubles_bayesian_framework1.pdf"), d, width = 7, height = 5)

  #Plot fitness densities for mean counts greater/less than 50 (replicate 1)
  d <- ggplot2::ggplot(doubles,ggplot2::aes(fitness_for_bins,..density..,color=counts_for_bins >= min_input_read_count_doubles)) +
    ggplot2::geom_density(adjust=1)
  ggplot2::ggsave(file.path(outpath, "4_doubles_bayesian_framework2.pdf"), d, width = 7, height = 5)
  #>> try to estimate what fitness scores are for variants with low sequence coverage
  # use double mutants with variants >= min_input_read_count_doubles counts 

  save.image(file = file.path(outpath, "Rsession1.RData"))

  ## calculate posterior double mutant fitness based on prior from single mutants
  postpois_conditioned_singleF <- function(i){  
    require(data.table)
    count_in = double_data[i,count_in]
    count_out = double_data[i,count_out]
    lam_in = exp(seq(floor(log(count_in+0.1)-max(c(0.5,1/log10(count_in+1.75)))),(log(count_in+0.1)+max(c(0.5,1/log10(count_in+1.75)))),lam_d))
    lam_out = exp(seq(floor(log(count_out+0.1)-max(c(0.5,1/log10(count_out+1.75)))),(log(count_out+0.1)+max(c(0.5,1/log10(count_out+1.75)))),lam_d))
    lam_low = range(log(lam_out))[1] - range(log(lam_in))[2]
    lam_high = range(log(lam_out))[2] - range(log(lam_in))[1]
    idx = row(matrix(NA,nrow=length(lam_out),ncol=length(lam_in))) - col(matrix(NA,nrow=length(lam_out),ncol=length(lam_in)))
    likelihood = sapply(split(outer(dpois(count_out,lambda = lam_out),dpois(count_in,lambda = lam_in)),idx),sum)
    score_prior = density(score_prior_cond[,.(fdist = sqrt((double_data[i,F1]-F1)^2+(double_data[i,F2]-F2)^2),F)][
      order(fdist)][1:Nneighbours,F],
      from = (lam_low-wt_corr),
      to = (lam_high-wt_corr),
      n = as.integer(as.character(round((lam_high-lam_low)/lam_d + 1)))) #super weird bug
    posterior = score_prior$y*likelihood
    
    moments = list()
    moments[1] = weighted.mean(x = score_prior$x,w = posterior)
    moments[2] = sqrt(sum(( moments[[1]]-score_prior$x)^2 * posterior)/
                        sum(posterior))
    return(moments)
  }

  # Setup cluster
  clust <- parallel::makeCluster(numCores) #This line will take time

  #Calculate conditional fitness and sigma
  for (E in all_reps) {
    #wildtype "correction" to calculate scores
    wt_corr <- dataset_syn[WT == T,log(unlist(.SD[,1]) / unlist(.SD[,2])),,.SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
    #data for prior calculation
    double_data <- doubles[,.(Pos1,Mut1,Pos2,Mut2,count_in = unlist(.SD[,1]),count_out = unlist(.SD[,2]),
                             F = unlist(.SD[,3])),,
                          .SDcols = c(paste0("count_e",E,"_s0"),paste0("count_e",E,"_s1"),paste0("fitness",E,"_uncorr"))]
    # double_data = merge(double_data,singles_silent[,.(Pos,Mut,F1 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
    double_data <- merge(double_data,singles_silent[!is.na(singles_silent[,paste0("fitness",E)]) & is.reads0==T,.(Pos,Mut,F1 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
    # double_data = merge(double_data,singles_silent[,.(Pos,Mut,F2 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))
    double_data <- merge(double_data,singles_silent[!is.na(singles_silent[,paste0("fitness",E)]) & is.reads0==T,.(Pos,Mut,F2 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))
    
    Nneighbours <- 500
    score_prior_cond <- double_data[count_in >= min_input_read_count_doubles & F > -Inf & F1 > -Inf & F2 > -Inf]

    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("double_data","lam_d","wt_corr","score_prior_cond","Nneighbours"), envir = environment())

    #posterior fitness conditioned on single fitness
    t=proc.time()
    helper <- parallel::parSapply(clust,X = 1:nrow(double_data), postpois_conditioned_singleF)
    print(proc.time()-t)
    helper1 <- matrix(unlist(helper),nrow=2)
    double_data[,paste0("fitness",E,"_cond") := helper1[1,]]
    double_data[,paste0("sigma",E,"_cond") := helper1[2,]]
    doubles <- merge(doubles, double_data[,.SD,,.SDcols = c("Pos1", "Pos2", "Mut1", "Mut2", paste0("fitness",E,"_cond"), paste0("sigma",E,"_cond"))], by = c("Pos1", "Pos2", "Mut1", "Mut2"), all.x = T)
  }
  parallel::stopCluster(clust)

  #Scatterplot matrix - singles
  d <- GGally::ggpairs(singles_silent[Nmut_aa==1,grep(names(singles_silent),pattern="fitness"),with=F],
          upper=list(continuous = "cor"))
  ggplot2::ggsave(file.path(outpath, "4_doubles_bayesian_framework_scattermatrix_singles.pdf"), d, width = 10, height = 10)

  #Scatterplot matrix - doubles, uncorrected
  set.seed(1)
  d <- GGally::ggpairs(doubles[apply(doubles[,.SD,,.SDcols = paste0("fitness",all_reps,"_uncorr")]==(-Inf), 1, sum)==0
                  ][sample(x = .N,1000),grep(names(doubles),pattern=paste0("fitness[", all_reps_str, "]_uncorr")),with=F],
          upper=list(continuous = "cor"))
  ggplot2::ggsave(file.path(outpath, "4_doubles_bayesian_framework_scattermatrix_doubles_uncorr.pdf"), d, width = 10, height = 10)

  #Scatterplot matrix - doubles, conditional
  set.seed(1)
  d <- GGally::ggpairs(doubles[apply(doubles[,.SD,,.SDcols = paste0("fitness",all_reps,"_uncorr")]==(-Inf), 1, sum)==0
                  ][sample(x = .N,1000),grep(names(doubles),pattern=paste0("fitness[", all_reps_str, "]_cond")),with=F],
          upper=list(continuous = "cor"))
  ggplot2::ggsave(file.path(outpath, "4_doubles_bayesian_framework_scattermatrix_doubles_cond.pdf"), d, width = 10, height = 10)

  #replicate 2 looks off compared to the others

  ### Normalise fitness and sigma for differences in number of generations (between biological replicates)
  ###########################

  #Mean generations per biological replicate
  sample_replicate <- unlist(sapply(strsplit(names(dataset)[grep("count", names(dataset))], "_"), '[', 2))
  names(sample_replicate) <- unlist(sapply(strsplit(names(dataset)[grep("count", names(dataset))], "_"), '[', 1))
  generations_mean <- tapply(generations$n_generations, sample_replicate[rownames(generations)], mean)

  #Normalise fitness and sigma by mean generations
  for (E in all_reps) {
    for(i in c("fitness", "sigma")){
      dataset_syn[,paste0(i, E, "_uncorr") := .SD/generations_mean[paste0("e", E)],,.SDcols = paste0(i, E, "_uncorr")]
      singles_silent[,paste0(i, E) := .SD/generations_mean[paste0("e", E)],,.SDcols = paste0(i, E)]
      doubles[,paste0(i, E, "_uncorr") := .SD/generations_mean[paste0("e", E)],,.SDcols = paste0(i, E, "_uncorr")]
      doubles[,paste0(i, E, "_cond") := .SD/generations_mean[paste0("e", E)],,.SDcols = paste0(i, E, "_cond")]
    }
  }

  save.image(file = file.path(outpath, "Rsession2.RData"))

  ### Merge fitness values
  ###########################

  #first check potential for replicate error
  dataset_syn[,var_fitness := rowSums((rowMeans(.SD[,1:nchar(all_reps_str)],na.rm=T) - .SD[,1:nchar(all_reps_str)])^2/(.SD[,(nchar(all_reps_str)+1):(2*nchar(all_reps_str))]^2),na.rm=T) / 
                rowSums(1/(.SD[,(nchar(all_reps_str)+1):(2*nchar(all_reps_str))]^2 ),na.rm=T),
              ,.SDcols = c(grep(names(dataset_syn),pattern=paste0("fitness[", all_reps_str, "]")),grep(names(dataset_syn),pattern=paste0("sigma[", all_reps_str, "]")))]
  dataset_syn[,avg_sigma := rowMeans(.SD[,(nchar(all_reps_str)+1):(2*nchar(all_reps_str))],na.rm=T),
              ,.SDcols = c(grep(names(dataset_syn),pattern=paste0("fitness[", all_reps_str, "]")),grep(names(dataset_syn),pattern=paste0("sigma[", all_reps_str, "]")))]
  dataset_syn[,isNA := rowSums(is.na(.SD)),,.SDcols = grep(names(dataset_syn),pattern=paste0("fitness[", all_reps_str, "]"))]
  dataset_syn[isNA==0 & var_fitness != Inf & avg_sigma != Inf,avg_sigma_fit := loess(var_fitness ~ avg_sigma,span=0.75)$fitted]

  replicate_error = dataset_syn[,min(sqrt(avg_sigma_fit),na.rm=T)]
  print(replicate_error) #roughly 1% replicate error

  #Average sigma versus fitness replicate error
  d <- ggplot2::ggplot(dataset_syn[isNA==0],ggplot2::aes(avg_sigma)) +
    ggplot2::geom_point(ggplot2::aes(y=sqrt(var_fitness))) +
    ggplot2::geom_line(ggplot2::aes(y=sqrt(avg_sigma_fit)),color="red") +
    ggplot2::geom_line(ggplot2::aes(y=sqrt(avg_sigma^2 + replicate_error^2)),color="darkgreen") +
    ggplot2::geom_abline(color="yellow",linetype=2) + 
    ggplot2::scale_x_log10() + 
    ggplot2::scale_y_log10()
  ggplot2::ggsave(file.path(outpath, "5_fitness_replicateerror_vs_avgsigma.pdf"), d, width = 7, height = 5)

  #### singles
  fitness_rx = singles_silent[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]"),colnames(singles_silent))]
  sigma_rx = sqrt(singles_silent[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]"),colnames(singles_silent))]^2 + 
                    matrix(replicate_error^2,nrow = dim(fitness_rx)[1],ncol = dim(fitness_rx)[2]))
  # sigma_s2 = random_effect_model(fitness_rx,sigma_rx)
  # singles_silent[,fitness := rowSums(fitness_rx/(sigma_rx^2 + sigma_s2),na.rm=T)/rowSums(1/(sigma_rx^2 + sigma_s2),na.rm=T)]
  singles_silent[,fitness := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
  # singles_silent[,sigma := sqrt(1/rowSums(1/(sigma_rx^2+sigma_s2),na.rm=T))]
  singles_silent[,sigma := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  d <- ggplot2::ggplot(singles_silent[Nmut_aa==1],ggplot2::aes(fitness,sigma)) + 
    ggplot2::geom_hex() + 
    ggplot2::scale_y_log10() + 
    ggplot2::coord_cartesian(ylim=c(0.01,1))
  ggplot2::ggsave(file.path(outpath, "5_sigma_vs_fitness_singles.pdf"), d, width = 5, height = 5)

  #### doubles
  #uncorrected fitness
  fitness_rx = doubles[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]_uncorr"),colnames(doubles))]
  sigma_rx = sqrt(doubles[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]_uncorr"),colnames(doubles))]^2 + 
                    matrix(replicate_error^2,nrow = dim(fitness_rx)[1],ncol = dim(fitness_rx)[2]))
  # sigma_s2 = random_effect_model(fitness_rx,sigma_rx)
  # doubles[,fitness_uncorr := rowSums(fitness_rx/(sigma_rx^2 + sigma_s2),na.rm=T)/rowSums(1/(sigma_rx^2 + sigma_s2),na.rm=T)]
  doubles[,fitness_uncorr := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
  # doubles[,sigma_uncorr := sqrt(1/rowSums(1/(sigma_rx^2+sigma_s2),na.rm=T))]
  doubles[,sigma_uncorr := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  d <- ggplot2::ggplot(doubles,ggplot2::aes(fitness_uncorr,sigma_uncorr)) + 
    ggplot2::geom_hex() + 
    ggplot2::scale_y_log10() + 
    ggplot2::coord_cartesian(ylim=c(0.01,1))
  ggplot2::ggsave(file.path(outpath, "5_sigma_vs_fitness_doubles_uncorr.pdf"), d, width = 5, height = 5)

  #conditioned fitness
  fitness_rx = doubles[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]_cond"),colnames(doubles))]
  sigma_rx = sqrt(doubles[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]_cond"),colnames(doubles))]^2 + 
                    matrix(replicate_error^2,nrow = dim(fitness_rx)[1],ncol = dim(fitness_rx)[2]))
  # sigma_s2 = random_effect_model(fitness_rx,sigma_rx)
  # doubles[,fitness_cond := rowSums(fitness_rx/(sigma_rx^2 + sigma_s2),na.rm=T)/rowSums(1/(sigma_rx^2 + sigma_s2),na.rm=T)]
  doubles[,fitness_cond := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
  # doubles[,sigma_cond := sqrt(1/rowSums(1/(sigma_rx^2+sigma_s2),na.rm=T))]
  doubles[,sigma_cond := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  d <- ggplot2::ggplot(doubles,ggplot2::aes(fitness_cond,sigma_cond)) + 
    ggplot2::geom_hex() + 
    ggplot2::scale_y_log10() + 
    ggplot2::coord_cartesian(ylim=c(0.01,1))
  ggplot2::ggsave(file.path(outpath, "5_sigma_vs_fitness_doubles_cond.pdf"), d, width = 5, height = 5)

  #Plot to compare double fitness estimates
  p1=ggplot2::ggplot(doubles,ggplot2::aes(mean_count,fitness_uncorr)) + 
    ggplot2::geom_hex() + 
    ggplot2::scale_x_log10() +
    ggplot2::scale_fill_continuous(trans="log10")
  p2=ggplot2::ggplot(doubles,ggplot2::aes(mean_count,fitness_cond)) + 
    ggplot2::geom_hex()+ 
    ggplot2::scale_x_log10() +
    ggplot2::scale_fill_continuous(trans="log10")
  p3=ggplot2::ggplot(doubles[between(bin_count,2,8)],ggplot2::aes(fitness_uncorr,..scaled..,color=factor(bin_count))) +
    ggplot2::geom_density(adjust=1)
  p4=ggplot2::ggplot(doubles[between(bin_count,2,8)],ggplot2::aes(fitness_cond,..scaled..,color=factor(bin_count))) +
    ggplot2::geom_density(adjust=1)
  p5=ggplot2::ggplot(doubles,ggplot2::aes(fitness_uncorr,sigma_uncorr)) + 
    ggplot2::geom_hex() + 
    ggplot2::scale_y_log10() +
    ggplot2::coord_cartesian(ylim = c(0.05,1.5))
  p6=ggplot2::ggplot(doubles,ggplot2::aes(fitness_cond,sigma_cond)) + 
    ggplot2::geom_hex()+ 
    ggplot2::scale_y_log10() +
    ggplot2::coord_cartesian(ylim = c(0.05,1.5))
  ggplot2::theme_set(ggplot2::theme_minimal())
  #Plot
  d <- cowplot::plot_grid(plotlist = list(p1,p2,p3,p4,p5,p6),nrow=3)
  rm(p1,p2,p3,p4,p5,p6)
  ggplot2::ggsave(file.path(outpath, "5_doubles_fitness_estimates.pdf"), d, width = 10, height = 10)

  #Plot fitness values against each other
  set.seed(1)
  d <- GGally::ggpairs(doubles[sample(.N,1000),.(fitness_uncorr,fitness_cond)])
  ggplot2::ggsave(file.path(outpath, "5_doubles_fitness_estimates_scattermatrix.pdf"), d, width = 10, height = 10)


  #Plot sigma values against each other
  d <- GGally::ggpairs(doubles[,.(sigma_uncorr,sigma_cond)])
  ggplot2::ggsave(file.path(outpath, "5_doubles_sigma_estimates_scattermatrix.pdf"), d, width = 10, height = 10)


  ### Output replicate data files
  ###########################

  # define which variants have enough reads
  wildtype[,is.reads0 := TRUE]
  singles_silent[mean_count >= min_mean_input_read_count,is.reads0 := TRUE]
  doubles[mean_count >= min_mean_input_read_count,is.reads0 := TRUE]

  #Save objects
  save(dataset_syn, singles_silent, doubles, file = file.path(outpath, "DMS_processed_data.RData"))

  ### Output plain text files
  ###########################

  ##### finalize data.tables
  silent = singles_silent[Nmut_aa==0,.(Pos,WT_AA,Mut,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,fitness,sigma)]
  singles = singles_silent[Nmut_aa==1,.(Pos,WT_AA,Mut,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,fitness,sigma)]

  #for doubles #add single mutant fitness/sigma values to double mutant table
  doubles[,fitness1 := singles[Pos == Pos1 & Mut == Mut1,fitness],.(Pos1,Mut1)]
  doubles[,sigma1 := singles[Pos == Pos1 & Mut == Mut1,sigma],.(Pos1,Mut1)]
  doubles[,fitness2 := singles[Pos == Pos2 & Mut == Mut2,fitness],.(Pos2,Mut2)]
  doubles[,sigma2 := singles[Pos == Pos2 & Mut == Mut2,sigma],.(Pos2,Mut2)]

  doubles = doubles[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,
                       fitness1,sigma1,fitness2,sigma2,
                       fitness_uncorr,sigma_uncorr,
                       fitness_cond,sigma_cond)]

  #Exclude variants with STOP codons from downstream fitness analyses
  wildtype[,is.fitness := TRUE]
  silent[,is.fitness := TRUE]
  singles[,is.fitness := !STOP]
  doubles[,is.fitness := !STOP]

  #write data to files
  write.table(x = wildtype, file = file.path(outpath, "DMS_wildtype.txt"),
              quote = F,row.names = F, col.names = T)
  write.table(x = silent, file = file.path(outpath, "DMS_silent.txt"),
              quote = F,row.names = F, col.names = T)
  write.table(x = singles, file = file.path(outpath, "DMS_singles.txt"),
              quote = F,row.names = F, col.names = T)
  write.table(x = doubles, file = file.path(outpath, "DMS_doubles.txt"),
              quote = F,row.names = F, col.names = T)
}
