
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

#' Perform differential mutation analysis
#'
#' @param  Muts_df A dataframe in the form required by HDProbe
#' @param nreps Integer; Number of replicates
#' @param homosked Logical; if TRUE, then the slope of the mean-variance relationship is set to 0
#' @param bg_pval Numeric; p-value cutoff for calling a nt's mutation rate as being no greater than the background error rate
#' @param bg_rate Numeric; background mutation rate assumed
#' @param alpha_p Numeric; shape1 of the beta prior used for mutation rate estimation
#' @param beta_p Numeric; shape2 of the beta prior used for mutation rate estiamtion
#' @importFrom magrittr %>%
#' @return A list with three elements
#' @export
#'
HDProbe <- function(Muts_df, nreps, homosked = FALSE,
                    bg_pval = 1, bg_rate = 0.002,
                    alpha_p = 2, beta_p = 102){

  Muts_df <- Muts_df %>% dplyr::mutate(MLE = (nmuts+alpha_p)/(ntrials+beta_p),
                                MLE_u = (nmuts + 1)/(ntrials + 1))



  ### Calculate logit variance

  Muts_df <- Muts_df %>% dplyr::mutate(phat_var = HDProbe::var_calc(MLE*ntrials + alpha_p, ntrials - MLE*ntrials + beta_p), logit_phat = pracma::logit(MLE),
                                phat_varu = HDProbe::var_calc(MLE_u*ntrials + 1, ntrials - MLE_u*ntrials + 1), logit_phat_u = pracma::logit(MLE_u))

  ## Look for sites to keep
  Filter_sites <- Muts_df %>% dplyr::mutate(check = ifelse(ntrials >= 1000, 1, 0)) %>%
    dplyr::group_by(P_ID) %>%
    dplyr::summarise(check = sum(check)) %>% dplyr::filter(check == nreps*2) %>% dplyr::select(P_ID)

  P_filter <- Filter_sites$P_ID

  Filter_df <- Muts_df[Muts_df$P_ID %in% P_filter,]

  ## Reorder by P_ID so that I can create new P_ID
  Filter_df <- Filter_df[order(Filter_df$P_ID),]


  ## Create new P_ID
  nsf <- nrow(Filter_df)/(nreps*2)

  if(nsf < 1000 & !homosked){
    message("Not enough nucleotides make it past read count filter; using previously identified conservative RV trend.")

    lm_list <- vector("list", length = 2)

    lm_list[[1]] <- c(1.4, -0.6)
    lm_list[[2]] <- c(1.4, -0.6)

  }else if(nsf < 100 & homosked){
    message("Not enough nucleotides make it past read count filter; using previously identified conservative RV trend")

    lm_list <- vector("list", length = 2)


    lm_list[[1]] <- c(log10(0.2), 0)
    lm_list[[2]] <- c(log10(0.2), 0)

  }else if(homosked){
    Filter_df$P_ID <- rep(1:nsf, each = nreps*2)


    ### Function to assess replicate variability trend

    Filter_df <- Filter_df %>% dplyr::group_by(E_ID, P_ID) %>%
      dplyr::summarise(tot_var = stats::var(logit_phat),
                phat_var = mean(phat_var),
                ntrials = sum(ntrials)) %>% dplyr::ungroup() %>%
      dplyr::mutate(var_rep = tot_var - phat_var)

    Filter_df <- Filter_df %>% dplyr::group_by(E_ID) %>%
      dplyr::summarise(avg_var_rep = mean(var_rep),
                RV_var = stats::var(var_rep)) %>%
      dplyr::mutate(avg_var_rep = log10(avg_var_rep))


    ## Regress avg_reads vs. avg_sd
    lm_list <- vector("list", length = 2)
    lm_var <- lm_list

    for(i in 1:2){

      lm_list[[i]] <- c(Filter_df$avg_var_rep[Filter_df$E_ID == i], 0)

      lm_var[[i]] <- Filter_df$RV_var[Filter_df$E_ID == i]
    }

    rm(Filter_df)

  }else{

    Filter_df$P_ID <- rep(1:nsf, each = nreps*2)


    ### Function to assess replicate variability trend
    nbin <- max(c(round((nsf)/1000), 10))

    Avg_df <- Filter_df %>% dplyr::group_by(E_ID, P_ID) %>%
      dplyr::summarise(tot_var = stats::var(logit_phat),
                phat_var = mean(phat_var),
                ntrials = mean(ntrials),
                lphat_mu = mean(logit_phat)) %>% dplyr::ungroup() %>%
      dplyr::mutate(bin_ID = as.numeric(Hmisc::cut2(ntrials, g = nbin))) %>%
      dplyr::mutate(var_rep = tot_var - phat_var)

    Filter_df <- Filter_df %>% dplyr::group_by(E_ID, P_ID) %>%
      dplyr::summarise(tot_var = stats::var(logit_phat),
                phat_var = mean(phat_var),
                ntrials = mean(ntrials)) %>% dplyr::ungroup() %>%
      dplyr::mutate(bin_ID = as.numeric(Hmisc::cut2(ntrials, g = nbin))) %>%
      dplyr::mutate(var_rep = tot_var - phat_var)

    Filter_df <- Filter_df %>% dplyr::group_by(E_ID, bin_ID) %>%
      dplyr::summarise(avg_var_rep = mean(var_rep), avg_reads = log10(mean(ntrials))) %>%
      dplyr::filter(avg_var_rep > 0)

    Filter_df <- Filter_df %>% dplyr::filter(avg_var_rep > 0)

    Filter_df <- Filter_df %>% dplyr::mutate(avg_var_rep = log10(avg_var_rep))


    ## Regress avg_reads vs. avg_sd
    lm_list <- vector("list", length = 2)
    int_vect <- rep(0, times = 2)
    slope_vect <- int_vect

    lm_var <- lm_list

    for(i in 1:2){
      heterosked_lm <- stats::lm(avg_var_rep ~ avg_reads, data = Filter_df[Filter_df$E_ID == i,] )
      h_int <- summary(heterosked_lm)$coefficients[1,1]
      h_slope <- summary(heterosked_lm)$coefficients[2,1]
      #lm_list[[i]] <- c(h_int, h_slope)

      int_vect[i] <- h_int
      slope_vect[i] <- h_slope

      lm_var[[i]] <- stats::var(stats::residuals(heterosked_lm))
    }

    rm(Filter_df)

    delta_est <- function(n){
      trials <- 10000

      vars <- rep(0, times = trials)
      for(i in 1:trials){
        vars[i] <- var(rnorm(n, mean = 3, sd = 1))
      }

      return(var(log10(vars)))
    }

    estimate_var <- delta_est(nreps)

    # Determine variance of variance
      # Maybe just need to delta estimate log10(tot_var) and log10(phat_var)
    Avg_df <- Avg_df %>% dplyr::ungroup() %>% dplyr::rowwise() %>%
      dplyr::mutate(var_exp =slope_vect[E_ID]*log10(ntrials) + int_vect[E_ID],
                    ltot_var = log10(tot_var) + (2*(tot_var^4))/((tot_var^2)*log(10)*(nreps - 1))) %>%
      dplyr::mutate(var_sdiff = ((ltot_var - ( (2/((ntrials - 1)*log(10)*log(10))) + var_exp ) )^2) - estimate_var )


    Avg_df2 <- Avg_df %>% dplyr::group_by(E_ID) %>%
      dplyr::summarise(sig_T2 = mean(var_sdiff)) # Variance about the trend


    sig_T2 <- Avg_df2$sig_T2[order(Avg_df2$E_ID)]

    rm(Avg_df2)

  }



  ### Average over replicates

  ## I just realized that if I want to use hierarchical modeling I can
  ## use rep_var trend + est_var as mean to shrink towards
  ## Just have to be smart about how strongly I shrink estimate

  Muts_df <- Muts_df %>%
    dplyr::mutate(MLE = ifelse(stats::pbinom(nmuts, ntrials, prob = bg_rate, lower.tail = FALSE) < bg_pval, MLE, bg_rate))


  ## Going to just use the trend as an exact replicate variability estimate for now

  Muts_df <- Muts_df %>% dplyr::group_by(P_ID, E_ID) %>%
    dplyr::summarise(Avg_lp = stats::weighted.mean(pracma::logit(MLE), w = ntrials), # Average logit(estimate)
              Avg_lp_u = stats::weighted.mean(pracma::logit(MLE_u), w = ntrials),
              Var_lp = stats::var(pracma::logit(MLE))*(nreps - 1)/nreps,
              avg_reads = mean(ntrials),
              ntrials = sum(ntrials),
              nmuts = sum(nmuts)) #%>%
    # dplyr::group_by(P_ID, E_ID) %>%
    # dplyr::mutate(rep_var = (10^(lm_list[[E_ID]][2]*log10(avg_reads) + lm_list[[E_ID]][1])),#/(nreps + 1),
    #        lvar_est = HDProbe::var_calc(nmuts + alpha_p, ntrials - nmuts + beta_p))


  lvar_est <- HDProbe::var_calc(Muts_df$nmuts + alpha_p, Muts_df$ntrials - Muts_df$nmuts + beta_p)
  sig_o2 <- (10^(slope_vect[Muts_df$E_ID]*log10(Muts_df$avg_reads) + int_vect[Muts_df$E_ID] )) + lvar_est


  # Estimate variance of variance on natural scale with delta approximation
  sig_Ts <- (log(10)^2)*(10^(2*sig_o2))*sig_T2[Muts_df$E_ID] - (((log(10)^4)*(10^(2*sig_o2))*(sig_T2[Muts_df$E_ID]^2))/4)


  par2 <- -(8 + (2*(sig_o2^2)/sig_Ts))
  cubic_soln <- purrr::map(.x = par2, .f = function(x) RConics::cubic(c(1, x, 20, -16)))
  nu_o <- purrr::map_dbl(.x = cubic_soln, .f = function(x) max(Re(x[Im(x) == 0]) ))



  Muts_df$rep_var <- (nu_o*sig_o2 + nreps*Muts_df$Var_lp)/(nu_o + nreps - 2)



  lden <- sqrt( (((nreps - 1)*(Muts_df$rep_var[Muts_df$E_ID == 1]))  + ((nreps-1)*(Muts_df$rep_var[Muts_df$E_ID == 2])))/(2*nreps - 2) )



  #lden <- sqrt( (((nreps - 1)*(Muts_df$rep_var[Muts_df$E_ID == 1] + Muts_df$lvar_est[Muts_df$E_ID == 1]))  + ((nreps-1)*(Muts_df$rep_var[Muts_df$E_ID == 2]  + Muts_df$lvar_est[Muts_df$E_ID == 2]) ))/(2*nreps - 2) )


  Testmut <- Muts_df %>% dplyr::ungroup() %>% dplyr::group_by(P_ID) %>%
    dplyr::summarise(l_num = Avg_lp[E_ID == 2] - Avg_lp[E_ID == 1],
              #l_den = sqrt(sum(lvar_est + rep_var) ), # nreps + 1 to convert to min MSE variance estimate
              reads = sum(ntrials))

  Testmut$lden <- lden

  Testmut <- Testmut %>% dplyr::mutate(lz = l_num/(lden*sqrt(2/nreps)),
                                lpval = 2*stats::pnorm(-abs(lz)),
                                lpadj = stats::p.adjust(lpval, method = "BH"))


  Results <- list(Diffmut = Testmut,
                  Mut_est = Muts_df,
                  lm_fit = list(slopes = slope_vect,
                                ints = int_vect),
                  hyperps = list(sig_o2 = sig_o2,
                                 nu_o = nu_o,
                                 lvar_est = lvar_est,
                                 sig_T2 = sig_T2),
                  Avg_df = Avg_df)

  return(Results)

}

