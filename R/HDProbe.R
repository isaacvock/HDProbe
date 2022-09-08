
#' Calculate variance with delta approximation
#'
#' @param mu Numeric; expectation value of random variable
#' @param var Numeric; variance of random variable
inv_var <- function(mu, var){

  var1 <- (exp(2*mu)/((1 + exp(mu))^4))*var
  var2_num <- exp(2*mu)*((exp(mu) - 1)^2)
  var2_den <- (exp(mu) + 1)^6

  var2 <- (var2_num/var2_den)*(var^2)

  totvar <- var1 - var2


  return(totvar)
}

#' Calculate variance with delta approximation using beta distribution
#'
#' @param alpha Numeric; shape1 parameter of beta distribution
#' @param beta Numeric; shape2 parameter of beta distribution
var_calc <- function(alpha, beta){

  EX <- alpha/(alpha + beta)
  VX <- (alpha*beta)/(((alpha + beta)^2)*(alpha + beta + 1))


  var1 <- (((1/EX) + 1/(1 - EX))^2)*VX
  var2 <- -(-1/(4*(EX^2)) + 1/(4*((1-EX)^2)))^2
  var3 <- VX^2

  totvar <- var1 - var2*var3

  return(totvar)
}

#' Calculate variance with delta approximation using frequentist model
#'
#' @param nmut Numeric; number of mutations
#' @param n Numeric; number of reads
freqvar <- function(nmut, n){

  phat <- nmut/n

  var1 <- (((1/phat) + 1/(1 - phat))^2)*((phat*(1-phat))/n)
  var2 <- -(-1/(4*(phat^2)) + 1/(4*((1-phat)^2)))^2
  var3 <- ((phat*(1-phat))/n)^2

  totvar <- var1 - var2*var3

  if(is.na(totvar)){
    totvar <- 0
  }

  return(totvar)
}

#' Estimate mutation rates and variances
#'
#' @param Muts_df a dataframe in form required by HDProbe
#' @param alpha_p Numeric; shape1 of the beta prior used for mutation rate estimation
#' @param beta_p Numeric; shape2 of the beta prior used for mutation rate estimation
#' @importFrom magrittr %>%
#' @return a dataframe of same form as Muts_df but with additional info
#' @export
EstimateRates <- function(Muts_df, alpha_p, beta_p){

  # Estimate rates
  Muts_df <- Muts_df %>% dplyr::mutate(MLE = (nmuts+alpha_p)/(ntrials+beta_p),
                                       MLE_u = (nmuts + 1)/(ntrials + 1))

  # Estimate uncertainties
  Muts_df <- Muts_df %>% dplyr::mutate(phat_var = HDProbe::var_calc(MLE*ntrials + alpha_p, ntrials - MLE*ntrials + beta_p), logit_phat = pracma::logit(MLE),
                                       phat_varu = HDProbe::var_calc(MLE_u*ntrials + 1, ntrials - MLE_u*ntrials + 1), logit_phat_u = pracma::logit(MLE_u))

}

#' Filter out low read count sites for variance trend estimation
#'
#' @param Muts_df a dataframe in form required by HDProbe
#' @importFrom magrittr %>%
#' @return a dataframe of same form as Muts_df but with additional info
#' @export
FilterSites <- function(Muts_df, nreps, cutoff = 1000){

  Filter_sites <- Muts_df %>% dplyr::mutate(check = ifelse(ntrials >= cutoff, 1, 0)) %>%
    dplyr::group_by(P_ID) %>%
    dplyr::summarise(check = sum(check)) %>% dplyr::filter(check == nreps*2) %>% dplyr::select(P_ID)

  P_filter <- Filter_sites$P_ID

  Filter_df <- Muts_df[Muts_df$P_ID %in% P_filter,]

  ## Reorder by P_ID so that I can create new P_ID
  Filter_df <- Filter_df[order(Filter_df$P_ID),]

}

#' Filter out low read count sites for variance trend estimation
#'
#' @param Filter_df a dataframe
#' @param nbin Number of read count bins
#' @importFrom magrittr %>%
#' @return a dataframe
#' @export
RepVarCalc <- function(Filter_df, nbin, alpha_p = 2, beta_p = 102){

  # Commenting out as I fix variance estimation
  # Avg_df <- Filter_df %>% dplyr::group_by(E_ID, P_ID) %>%
  #   dplyr::summarise(tot_var = stats::var(logit_phat),
  #                    phat_var = mean(phat_var),
  #                    ntrials = mean(ntrials),
  #                    lphat_mu = mean(logit_phat)) %>% dplyr::ungroup() %>%
  #   dplyr::mutate(bin_ID = as.numeric(Hmisc::cut2(ntrials, g = nbin))) %>%
  #   dplyr::mutate(var_rep = tot_var - phat_var)

  Avg_df <- Filter_df %>% dplyr::group_by(E_ID, P_ID) %>%
    dplyr::summarise(tot_var = stats::var(logit_phat),
                     tot_trials = sum(ntrials),
                     ntrials = mean(ntrials),
                     lphat_mu = mean(logit_phat)) %>% dplyr::ungroup() %>%
    dplyr::mutate(bin_ID = as.numeric(Hmisc::cut2(ntrials, g = nbin)))

  Avg_df <- Avg_df %>%
    dplyr::mutate(phat_var = HDProbe::var_calc(lphat_mu*tot_trials + alpha_p, tot_trials - lphat_mu*tot_trials + beta_p)) %>%
    dplyr::mutate(var_rep = tot_var - phat_var)

  return(Avg_df)
}

#' Estimate replicate variability trend
#'
#' @param Filter_df a dataframe
#' @param Homosked Logical; if TRUE, homoskedasticity is assumed
#' @importFrom magrittr %>%
#' @return a list of linear model fits
#' @export
VarianceTrend <- function(Filter_df, Homosked = FALSE){

  if(Homosked){
    Filter_df <- Filter_df %>% dplyr::group_by(E_ID) %>%
      dplyr::summarise(avg_var_rep = mean(var_rep),
                RV_var = stats::var(var_rep)) %>%
      dplyr::mutate(avg_var_rep = log10(avg_var_rep))


    ## Regress avg_reads vs. avg_sd
    int_vect <- rep(0, times = 2)
    slope_vect <- int_vect

    for(i in 1:2){

      int_vect[i] <- Filter_df$avg_var_rep[Filter_df$E_ID == i]
      slope_vect[i] <- 0

    }

  }else{
    min_vr <- min(Filter_df$var_rep[Filter_df$var_rep > 0])

    Filter_df <- Filter_df %>%
      dplyr::mutate(var_rep = ifelse(var_rep < 0, min_vr, var_rep))


    Filter_df <- Filter_df %>% dplyr::group_by(E_ID, bin_ID) %>%
      dplyr::summarise(avg_var_rep = mean(var_rep), avg_reads = log10(mean(ntrials))) %>%
      dplyr::filter(avg_var_rep > 0)

    Filter_df <- Filter_df %>% dplyr::filter(avg_var_rep > 0)

    Filter_df <- Filter_df %>% dplyr::mutate(avg_var_rep = log10(avg_var_rep))


    ## Regress log10(avg. reads) vs. log10(avg. variance)
    int_vect <- rep(0, times = 2)
    slope_vect <- int_vect

    for(i in 1:2){
      heterosked_lm <- stats::lm(avg_var_rep ~ avg_reads, data = Filter_df[Filter_df$E_ID == i,] )
      h_int <- summary(heterosked_lm)$coefficients[1,1]
      h_slope <- summary(heterosked_lm)$coefficients[2,1]
      #lm_list[[i]] <- c(h_int, h_slope)

      int_vect[i] <- h_int
      slope_vect[i] <- h_slope
    }
  }

  lm_list <- list(slope_vect = slope_vect,
                  int_vect = int_vect,
                  Avg_RV = Filter_df$avg_var_rep,
                  Avg_reads = Filter_df$avg_reads,
                  E_ID = Filter_df$E_ID)

  return(lm_list)



}

#' Estimate variability about replicate variability trend
#'
#' @param Avg_df a dataframe
#' @param estiamte_var Variance of log10(sample variance)
#' @param slope_vect Vector of variance trend slopes
#' @param int_vect Vector of variance trend intercepts
#' @param nreps Number of replicates
#' @importFrom magrittr %>%
#' @return a list containing a dataframe and a vector
#' @export
TrendVariance <- function(Avg_df, estimate_var,
                          slope_vect, int_vect, nreps){

  # Avg_df <- Avg_df %>% dplyr::ungroup() %>% dplyr::rowwise() %>%
  #   dplyr::mutate(var_exp =slope_vect[E_ID]*log10(ntrials) + int_vect[E_ID],
  #                 ltot_var = log10(tot_var) + (2*(tot_var^4))/((tot_var^2)*log(10)*(nreps - 1))) %>%
  #   #dplyr::mutate(var_sdiff = ((ltot_var - ( (2/((ntrials - 1)*log(10)*log(10))) + var_exp ) )^2) - estimate_var )
  #   dplyr::mutate(var_sdiff = ((ltot_var - var_exp)^2))

  min_vr <- min(Avg_df$var_rep[Avg_df$var_rep > 0])

  Avg_df <- Avg_df %>%
    dplyr::mutate(var_rep = ifelse(var_rep < 0, min_vr, var_rep))


  Avg_df <- Avg_df %>% dplyr::ungroup() %>% dplyr::rowwise() %>%
    dplyr::mutate(var_exp = slope_vect[E_ID]*log10(ntrials) + int_vect[E_ID],
                  var_sdiff = log10(var_rep) - var_exp )


  Avg_df2 <- Avg_df %>% dplyr::group_by(E_ID) %>%
    dplyr::summarise(sig_T2 = var(var_sdiff) - estimate_var) # Variance about the trend


  sig_T2 <- Avg_df2$sig_T2[order(Avg_df2$E_ID)]

  Var_list <- list(Avg_df, sig_T2)
}

#' Average mutation rate estimates and regularize variance estiamte
#'
#' @param Muts_df a dataframe
#' @param bg_pval Numeric; p-value cutoff for calling a nt's mutation rate as being no greater than the background error rate
#' @param bg_rate Numeric; background mutation rate assumed
#' @param alpha_p Numeric; shape1 of the beta prior used for mutation rate estimation
#' @param beta_p Numeric; shape2 of the beta prior used for mutation rate estimation
#' @param slope_vect Vector of variance trend slopes
#' @param int_vect Vector of variance trend intercepts
#' @param nreps Number of replicates
#' @param sig_T2 Variance about variance trend
#' @importFrom magrittr %>%
#' @return dataframe
#' @export
AvgAndReg <- function(Muts_df, bg_pval, bg_rate, alpha_p, beta_p,
                     slope_vect, int_vect, nreps, sig_T2){

  Muts_df <- Muts_df %>%
    dplyr::mutate(MLE = ifelse(stats::pbinom(nmuts, ntrials, prob = bg_rate, lower.tail = FALSE) < bg_pval, MLE, bg_rate))


  ## Going to just use the trend as an exact replicate variability estimate for now

  Muts_df <- Muts_df %>% dplyr::group_by(P_ID, E_ID, GF) %>%
    dplyr::summarise(Avg_lp = stats::weighted.mean(pracma::logit(MLE), w = ntrials), # Average logit(estimate)
                     Avg_lp_u = stats::weighted.mean(pracma::logit(MLE_u), w = ntrials),
                     Var_lp = stats::var(pracma::logit(MLE))*(nreps - 1)/nreps, # Total variance of logit(estimate)
                     avg_reads = mean(ntrials),
                     ntrials = sum(ntrials),
                     nmuts = sum(nmuts))

  # Estimate variance using Bayesian posterior
  lvar_est <- HDProbe::var_calc(Muts_df$nmuts + alpha_p, Muts_df$ntrials - Muts_df$nmuts + beta_p)

  # Total variance trend
  sig_o2 <- (10^(slope_vect[Muts_df$E_ID]*log10(Muts_df$avg_reads) + int_vect[Muts_df$E_ID] )) + lvar_est


  # Estimate variance of variance on natural scale with delta approximation
  if(any(is.infinite(sig_T2))){
    sig_Ts <- rep(Inf, times = nrow(Muts_df))

    # Prior degrees of freedom of scaled inverse chi-squared distribution (see BDA3)
    nu_o <- 0

  }else{
    sig_Ts <- (log(10)^2)*(10^(2*sig_o2))*sig_T2[Muts_df$E_ID] - (((log(10)^4)*(10^(2*sig_o2))*(sig_T2[Muts_df$E_ID]^2))/4)

    # Prior degrees of freedom of scaled inverse chi-squared distribution (see BDA3)
    nu_o <- mean((2*(sig_o2^2)/sig_Ts) + 4)

  }




  if(is.infinite(nu_o)){
    # Posterior total variance estimate
    Muts_df$rep_var <- sig_o2

  }else{
    # Posterior total variance estimate
    Muts_df$rep_var <- (nu_o*sig_o2 + nreps*Muts_df$Var_lp)/(nu_o + nreps - 1)

  }


  Reg_list <- list(Muts_df, nu_o, sig_Ts, sig_o2)

  return(Reg_list)

}

#' Test for differential mutation rates
#'
#' @param Muts_df a dataframe
#' @param lden Denominator of test statistic
#' @param nreps Number of replicates
#' @param nu_o Prior degrees of freedom
#' @importFrom magrittr %>%
#' @return dataframe
#' @export
DiffMutTest <- function(Muts_df, lden, nreps, nu_o){
  Testmut <- Muts_df %>% dplyr::ungroup() %>% dplyr::group_by(P_ID) %>%
    dplyr::summarise(l_num = Avg_lp[E_ID == 2] - Avg_lp[E_ID == 1],
                     #l_den = sqrt(sum(lvar_est + rep_var) ), # nreps + 1 to convert to min MSE variance estimate
                     reads = sum(ntrials))

  Testmut$lden <- lden

  Testmut <- Testmut %>% dplyr::mutate(lz = l_num/(lden*sqrt(2/nreps)),
                                       #lpval = 2*stats::pt(-abs(lz), df = 2*nreps - 2 + 2*nu_o),
                                       lpval = 2*stats::pnorm(-abs(lz)),
                                       lpadj = stats::p.adjust(lpval, method = "BH"))

  return(Testmut)

}


#' Perform differential mutation analysis
#'
#' @param  Muts_df A dataframe in the form required by HDProbe
#' @param nreps Integer; Number of replicates
#' @param homosked Logical; if TRUE, then the slope of the mean-variance relationship is set to 0
#' @param bg_pval Numeric; p-value cutoff for calling a nt's mutation rate as being no greater than the background error rate
#' @param bg_rate Numeric; background mutation rate assumed
#' @param alpha_p Numeric; shape1 of the beta prior used for mutation rate estimation. If alpha_p or beta_p are NULL, then this is
#' estimated using beta distribution method of moments estimator
#' @param beta_p Numeric; shape2 of the beta prior used for mutation rate estimation. If alpha_p or beta_p are NULL, then this
#' is estimated using beta distribution method of moments estimator
#' @param One_ctl Logical; if TRUE, then mutation rate comparisons are made to the global average control sample mutation
#' rate
#' @param Gene_ctl Logical; if TRUE, then mutation rate comparisons are made to the
#' gene-wide average control sample mutation rate
#' @param var_of_var Variance of variance to be used to tune regularization. If NULL, this
#' is estimated empirically (i.e., from the data)
#' @param ... Parameters that can be supplied to internal functions
#' @importFrom magrittr %>%
#' @return A list with three elements
#' @export
#'
HDProbe <- function(Muts_df, nreps, homosked = FALSE,
                    bg_pval = 1, bg_rate = 0.002,
                    alpha_p = NULL, beta_p = NULL,
                    filter_het = 1000, filter_hom = 100,
                    One_ctl = FALSE, Gene_ctl = FALSE, var_of_var = NULL, ...){

  if(is.null(alpha_p) | is.null(beta_p)){
    # Estimate alpha_p and beta_p
    E_df <- Muts_df %>% ungroup() %>%
      dplyr::mutate(mutrate = (nmuts + 1)/(ntrials + 1)) %>%
      summarise(E = mean(mutrate), V = var(mutrate))

    alpha_p <- E_df$E*(((E_df$E*(1-E_df$E))/E_df$V) - 1)
    beta_p <- (alpha_p*(1-E_df$E))/E_df$E

    rm(E_df)
  }




  # Estimate mutation rates and estimate uncertainties
  Muts_df <- EstimateRates(Muts_df, alpha_p, beta_p)


  # Filter for variance trend estimation
  Filter_df <- FilterSites(Muts_df, nreps = nreps, ...)

  ## Create new P_ID
  nsf <- nrow(Filter_df)/(nreps*2)

  if(nsf < filter_het & !homosked){
    message("Not enough nucleotides make it past read count filter; using previously identified conservative RV trend.")

    lm_list <- vector("list", length = 2)

    slope_vect <- rep(-0.6, times = 2)

    int_vect <- rep(1.4, times = 2)

  }else if(nsf < filter_hom & homosked){
    message("Not enough nucleotides make it past read count filter; using previously identified conservative RV trend")


    int_vect <- rep(log10(0.2), times = 2)
    slope_vect <- rep(0, times = 2)


  }else if(homosked){



    Filter_df$P_ID <- rep(1:nsf, each = nreps*2)


    # Number of bins for replicate variability trend fitting
    nbin <- 1

    # Estimate replicate variability for each position
    Avg_df <- RepVarCalc(Filter_df, nbin)

    # Duplicate df so that one can be edited
    Filter_df <- Avg_df

    # Estimate replicate variability trend
    lm_list <- VarianceTrend(Filter_df, Homosked = TRUE)

    rm(Filter_df)

    # Variance of the log10(sample variance)
    int_vect <- lm_list$int_vect
    slope_vect <- lm_list$slope_vect

    # Varaince of log10(sample variance)
    estimate_var <- pracma::psi(1,(nreps-1)/2)*(log10(exp(1))^2)

    Var_list <- TrendVariance(Avg_df, estimate_var, int_vect, slope_vect, nreps)


    Avg_df <- Var_list[[1]]
    sig_T2 <- Var_list[[2]]

    rm(Var_list)

  }else{

    Filter_df$P_ID <- rep(1:nsf, each = nreps*2)


    # Number of bins for replicate variability trend fitting
    nbin <- max(c(round((nsf)/1000), 10))

    # Estimate replicate variability for each position
    Avg_df <- RepVarCalc(Filter_df, nbin)

    # Duplicate df so that one can be edited
    Filter_df <- Avg_df

    # Estimate replicate variability trend
    lm_list <- VarianceTrend(Filter_df)


    rm(Filter_df)

    # Intercept and slope of replicate variability trend
    int_vect <- lm_list$int_vect
    slope_vect <- lm_list$slope_vect


    # Variance of log10(sample variance)
    estimate_var <- pracma::psi(1,(nreps-1)/2)*(log10(exp(1))^2)

    # Calc. variance about trend
    Var_list <- TrendVariance(Avg_df, estimate_var, int_vect, slope_vect, nreps)

    Avg_df <- Var_list[[1]]
    sig_T2 <- Var_list[[2]]

    rm(Var_list)

  }


  if(!is.null(var_of_var)){
    sig_T2 <- rep(var_of_var, times = length(sig_T2))
  }

  # Average over replicates and regularize variance
  Reg_list <- AvgAndReg(Muts_df, bg_pval, bg_rate, alpha_p, beta_p,
                       slope_vect, int_vect, nreps, sig_T2)

  Muts_df <- Reg_list[[1]]
  nu_o <- Reg_list[[2]]
  sig_Ts <- Reg_list[[3]]
  sig_o2 <- Reg_list[[4]]

  if(One_ctl){

    Avg_ctl <- stats::weighted.mean(Muts_df$Avg_lp[Muts_df$E_ID == 1], weights = Muts_df$ntrials[Muts_df$E_ID == 1])

    Muts_df <- Muts_df %>%
      dplyr::mutate(Avg_lp = ifelse(E_ID == 1, Avg_ctl, Avg_lp))
  }else if(Gene_ctl){

    Gene_df <- Muts_df[Muts_df$E_ID == 1,] %>% dplyr::group_by(GF) %>%
      dplyr::summarise(Avg_lp = stats::weighted.mean(Avg_lp, weights = ntrials))

    Muts_df <- Muts_df %>% dplyr::rowwise() %>%
      dplyr::mutate(Avg_lp = ifelse(E_ID == 1, Gene_df$Avg_lp[Gene_df$GF == GF], Avg_lp))

  }



  rm(Reg_list)


  # Denominator of test statistic
  lden <- sqrt( (((nreps - 1)*(Muts_df$rep_var[Muts_df$E_ID == 1]))  + ((nreps-1)*(Muts_df$rep_var[Muts_df$E_ID == 2])))/(2*nreps - 2) )



  #lden <- sqrt( (((nreps - 1)*(Muts_df$rep_var[Muts_df$E_ID == 1] + Muts_df$lvar_est[Muts_df$E_ID == 1]))  + ((nreps-1)*(Muts_df$rep_var[Muts_df$E_ID == 2]  + Muts_df$lvar_est[Muts_df$E_ID == 2]) ))/(2*nreps - 2) )

  # Perform differential mutation rate testing
  Testmut <- DiffMutTest(Muts_df, lden, nreps, nu_o)

  Results <- list(Diffmut = Testmut,
                  Mut_est = Muts_df,
                  lm_fit = list(slopes = slope_vect,
                                ints = int_vect),
                  Full_lm = lm_list,
                  nu_o = nu_o,
                  sig_Ts = sig_Ts,
                  sig_T2 = sig_T2,
                  sig_o2 = sig_o2,
                  Avg_df = Avg_df,
                  Priors = c(alpha_p = alpha_p, beta_p = beta_p))

  return(Results)

}

