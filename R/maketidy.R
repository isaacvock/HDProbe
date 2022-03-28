#' Put mutation probing data into tidy format
#'
#' @param wt_muts Dataframe from csv containing mutation profiling data for reference sample
#' @param ko_muts Dataframe from csv containing mutation profiling data for experimental sample
#' @param nreps Numeric; number of replicates
#' @param target Logical; if TRUE, then data is treated as targeted probing data
#' @param SHAPE Logical; if TRUE, then data is treated as SHAPE data rather than DMS-seq data
#' @param totcut Numeric; read count cutoff to be considered for analysis
#' @param rate_cut Numeric; 1/(max mutation rate) allowed; used for SNP filtering
#' @importFrom magrittr %>%

#' @return a dataframe that can be passed to HDProbe
#' @export
maketidy <- function(wt_muts, ko_muts, nreps, target = FALSE, SHAPE = FALSE,
                     totcut = 500, rate_cut = 4){
  ## Filter out anything with less than totcut reads and anomalously high mutation rates
  # Also filtering for any mutations for investigation purposes

  if(nreps == 3){
    wt_muts <- wt_muts %>% dplyr::filter(R1_cov > totcut  & R2_cov > totcut & R3_cov > totcut) %>%
      dplyr::filter(R1_mut < R1_cov/rate_cut & R2_mut < R2_cov/rate_cut & R3_mut < R3_cov/rate_cut)

    ko_muts <- ko_muts %>% dplyr::filter(R1_cov > totcut  & R2_cov > totcut & R3_cov > totcut) %>%
      dplyr::filter(R1_mut < R1_cov/rate_cut & R2_mut < R2_cov/rate_cut & R3_mut < R3_cov/rate_cut)

  }else if(nreps == 2){
    wt_muts <- wt_muts %>% dplyr::filter(R1_cov > totcut  & R2_cov > totcut) %>%
      dplyr::filter(R1_mut < R1_cov/rate_cut & R2_mut < R2_cov/rate_cut)

    ko_muts <- ko_muts %>% dplyr::filter(R1_cov > totcut  & R2_cov > totcut) %>%
      dplyr::filter(R1_mut < R1_cov/rate_cut & R2_mut < R2_cov/rate_cut)

  }else{
    stop("Isaac needs to write better code")
  }



  if(target){
    wt_nts <- wt_muts[,c("pos")]

  }else{
    wt_nts <- wt_muts[,c("chr", "pos", "strand")]

  }

  ### Now for KO data




  if(target){
    ko_nts <- ko_muts[,c("pos")]
    # Find nts in both ko and wt

    both_nts <- dplyr::semi_join(ko_nts, wt_nts)


    # Filter out for only nts in both datasets
    wt_muts <- dplyr::semi_join(wt_muts, both_nts, by = c("pos") )
    ko_muts <- dplyr::semi_join(ko_muts, both_nts, by = c("pos"))


    if(!SHAPE){
      # Filter out for mutations that can occur in DMS-seq
      wt_muts <- wt_muts %>% dplyr::filter(nt %in% c("A", "C"))
      ko_muts <- ko_muts %>% dplyr::filter(nt %in% c("A", "C"))
    }



    ### Convert to tidy format and combine data

    muts <- rep(0, times = nrow(wt_muts)*2*nreps)
    covs <- muts
    R_ID <- muts
    E_ID <- muts
    strand <- muts

    nsites <- nrow(wt_muts)


    for(i in 1:nreps){

      mut_id <- grep(paste0("R", i, "_mut"), colnames(wt_muts))
      cov_id <- grep(paste0("R", i, "_cov"), colnames(wt_muts))

      muts[(2*(i-1)*nsites + 1):(2*i*nsites)]  <- c(unclass(wt_muts[,mut_id])[[1]], unclass(ko_muts[,mut_id])[[1]])
      covs[(2*(i-1)*nsites + 1):(2*i*nsites)]  <- c(unclass(wt_muts[,cov_id])[[1]], unclass(ko_muts[,cov_id])[[1]])

      R_ID[(2*(i-1)*nsites + 1):(2*i*nsites)] <- rep(i, times = nrow(wt_muts)*2)
      E_ID[(2*(i-1)*nsites + 1):(2*i*nsites)] <- rep(1:2, each = nrow(wt_muts))


    }

    Muts_df <- tibble::tibble(P_ID = rep(1:nrow(wt_muts), times = nreps*2),
                              nmuts = muts,
                              ntrials = covs,
                              R_ID = R_ID,
                              E_ID = E_ID,
                              pos = rep(wt_muts$pos, times = nreps*2))

    return(Muts_df)
  }else{
    ko_nts <- ko_muts[,c("chr", "pos", "strand")]

    # Find nts in both ko and wt

    both_nts <- dplyr::semi_join(ko_nts, wt_nts)


    # Filter out for only nts in both datasets
    wt_muts <- dplyr::semi_join(wt_muts, both_nts, by = c("chr", "pos", "strand"))
    ko_muts <- dplyr::semi_join(ko_muts, both_nts, by = c("chr", "pos", "strand"))

    if(!SHAPE){
      # Filter out for mutations that can occur in DMS-seq
      wt_muts <- wt_muts %>% dplyr::filter(nt %in% c("A", "C"))
      ko_muts <- ko_muts %>% dplyr::filter(nt %in% c("A", "C"))
    }



    ### Convert to tidy format and combine data

    muts <- rep(0, times = nrow(wt_muts)*2*nreps)
    covs <- muts
    R_ID <- muts
    E_ID <- muts
    strand <- muts

    nsites <- nrow(wt_muts)


    for(i in 1:nreps){

      mut_id <- grep(paste0("R", i, "_mut"), colnames(wt_muts))
      cov_id <- grep(paste0("R", i, "_cov"), colnames(wt_muts))

      muts[(2*(i-1)*nsites + 1):(2*i*nsites)]  <- c(unclass(wt_muts[,mut_id])[[1]], unclass(ko_muts[,mut_id])[[1]])
      covs[(2*(i-1)*nsites + 1):(2*i*nsites)]  <- c(unclass(wt_muts[,cov_id])[[1]], unclass(ko_muts[,cov_id])[[1]])

      R_ID[(2*(i-1)*nsites + 1):(2*i*nsites)] <- rep(i, times = nrow(wt_muts)*2)
      E_ID[(2*(i-1)*nsites + 1):(2*i*nsites)] <- rep(1:2, each = nrow(wt_muts))


    }

    Muts_df <- tibble::tibble(P_ID = rep(1:nrow(wt_muts), times = nreps*2),
                              nmuts = muts,
                              ntrials = covs,
                              R_ID = R_ID,
                              E_ID = E_ID,
                              chr = rep(wt_muts$chr, times = nreps*2),
                              pos = rep(wt_muts$pos, times = nreps*2),
                              strand = rep(wt_muts$strand, times = nreps*2))

    return(Muts_df)

  }




}
