#' Simulate heteroskedastic mutational probing data
#'
#' @param npos Numeric; number of positions to simulate mutational data for
#' @param nreps Numeric; number of replicates to simulate for each of 2 conditions
#' @param het_m Numeric; slope of sequencing depth-variance relationship on log-log scale
#' @param het_b Numeric; intercept of sequencing depth-variance relationship on log-log scale
#' @param nsig Numeric; number of significant changes in mutation rate
#' @param mu_e Numeric; differences in logit(mutation rate) for significant nts drawn from normal(mu_e, sig_e)
#' @param sig_e Numeric; effect size parameter (see mu_e for details)
#' @param lp_mean Numeric; average logit(mutation rate)
#' @param lp_sd Numeric; standard deviation of logit(mutation rates)
#' @param r_mean Numeric; Avg. read counts drawn from lognormal(r_mean, r_sd)
#' @param r_sd Numeric; read counts drawn from lognormal(r_mean, r_sd)
#' @param sig_T Numeric; standard deviation of logit(p) variances on log scale
#' @param phi Numeric; read counts for each sample drawn from negative binomial distribution with size = phi
#' @importFrom magrittr %>%
#'
#' @return List with 2 componenets, one of which is a dataframe that can be passed to HDProbe
#' @export
sim_muts <- function(npos, nreps = 3,
                     het_m = -0.7, het_b = 1.3,
                     nsig = round(npos/10), mu_e = 0, sig_e = 0.75,
                     lp_mean = -4.35, lp_sd = 0.8,
                     r_mean = 6.75, r_sd = 1.2, sig_T = 0.05,
                     phi = 20){

  ### Simulate read coverage mean
  cov_means <- stats::rlnorm(npos*2, meanlog = r_mean, sdlog = r_sd)

  ### Simulate reads in each replicate

  reads <- stats::rnbinom(npos*nreps*2, mu = rep(cov_means, each = nreps), size = phi)


  ### Simulate heteroskedastic mutation rates

  ## Simulate mutation rate effect sizes
  effs <- c(rep(0, times = npos), rep(0, times = npos - nsig), stats::rnorm(nsig, mean = mu_e, sd = sig_e))

  ## Avg. ctl rate for each position
  avg_lps <- stats::rnorm(npos, mean = lp_mean, sd = lp_sd)

  avg_lps <- rep(avg_lps, times = 2) + effs

  var_sds <- sqrt(10^(stats::rnorm(npos*nreps*2, mean = het_m*log10(rep(cov_means, each = nreps)) + het_b, sd = sig_T)) )

  ## Individual rate for each position and each replicate
  lps <- stats::rnorm(npos*nreps*2, mean = rep(avg_lps, each = nreps), sd = var_sds)



  ### Simulate mutational data

  nmuts <- stats::rbinom(npos*nreps*2, size = reads, prob = pracma::sigmoid(lps))


  ### Create tibble

  mut_df <- tibble::tibble(cov = reads, muts = nmuts,
                           rep_id = rep(1:nreps, times = npos*2),
                           pos_id = rep(rep(1:npos, each = nreps), times = 2),
                           exp_id = rep(1:2, each = npos*nreps))



  ##### Analyze with HDProbe
  colnames(mut_df) <- c("ntrials", "nmuts", "R_ID",
                        "P_ID", "E_ID")

  avg_truth <- tibble::tibble(logit_p = avg_lps,
                              logit_effs = effs,
                              P_ID = rep(1:npos, times = 2),
                              E_ID = rep(1:2, each = npos))

  rep_truth <- tibble::tibble(logit_p = lps,
                              P_ID = mut_df$P_ID,
                              E_ID = mut_df$E_ID,
                              R_ID = mut_df$R_ID)


  sim_list <- list(avg_truth = avg_truth,
                   rep_truth = rep_truth,
                   mut_df = mut_df,
                   var_sds = var_sds)

  return(sim_list)

}
