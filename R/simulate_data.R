#' Simulate species observations
#'
#' @param quadratic add quadratic effect to the predictor
#' @param sim_list parameters and latent variables from "simulate_system"
#'
#' @returns list of 1) dataframe of simulated species observations, 2) study design for the observatioons
#' @export
#'
#' @examples
sim_sp_obs <- function(quadratic=F,
                       sim_list) {
  species_resp <- sim_list[[1]]
  sites <- sim_list[[2]]
  sp_order <- unique(species_resp$species) # sort species order

  # set up DFs
  sim_sp_obs <- as.data.frame(matrix(nrow=length(sites$val[sites$lv == 1]),
                                     ncol=max(species_resp$species)))
  sim_sdesign <- as.data.frame(matrix(nrow=length(sites$val[sites$lv == 1]),
                                      ncol=3))
  names(sim_sdesign) <- c("group", "site", "time")
  names(sim_sp_obs) <- paste0("sp_", sp_order)

  # simulate sp obs
  # need id column for unique combo of time X site
  sites$timesite <- rep(1:(length(sites$t[sites$lv == 1])), max(sites$lv))
  for (j in sp_order) {
    for (ts in unique(sites$timesite)) {
      pos <- sites[sites$timesite==ts,]
      linpred <- 0
      for (lv in pos$lv){
        sp_vals <- species_resp[species_resp$species==j & species_resp$lv==lv,]
        linpred <- linpred + sp_vals$intercept + pos$val[pos$lv==lv]*sp_vals$slope
        if (quadratic) {
          linpred <- linpred + (pos$val[pos$lv==lv]^2)*sp_vals$quad
        }
      }
      sim_obs <- rpois(n=1,lambda = exp(linpred))
      sim_sp_obs[ts,j] <- sim_obs
      sim_sdesign[ts,] <- c(pos$group[1], pos$site[1], pos$t[1])
    }
  }
  list(sim_sp_obs, sim_sdesign)
}





#' Simulate a system for testing the OU restoration gllvm model
#'
#' @param num_lv number of latent variables
#' @param sim_resp_par list of arguments for the species loadings / intercepts
#' @param sim_drift_par list of arguments for the OU-process on the latent variables
#'
#' @returns
#' @export
#'
#' @examples
simulate_system  <- function(num_lv = 1,
                             sim_resp_par = list(FALSE), # list of arguments to the sp.resp-functions
                             sim_drift_par= list(FALSE) # list of arguments to the site-drift-functions
) {

  # add dummy parameters if number of LVs larger than no. of params specified
  while (min(length(sim_resp_par), length(sim_drift_par)) < num_lv) {
    if (length(sim_resp_par) < num_lv) {
      sim_resp_par[[length(sim_resp_par) + 1]] <- FALSE
    }
    if (length(sim_drift_par) < num_lv) {
      sim_drift_par[[length(sim_drift_par) + 1]] <- FALSE
    }
  }

  # simulate species responses (why didn't ifelse() work here??)
  if (length(sim_resp_par[[1]]) != 1){
    resp_df <- do.call(sp_params, sim_resp_par[[1]])
  } else {
    resp_df <- sp_params()
  }
  resp_df$lv <- 1 # for first latent variable


  if (num_lv > 1) { # for more than 1 lv
    for (l in 2:num_lv) {
      if (length(sim_resp_par[[l]]) != 1){
        df_lv <- do.call(sp_params, sim_resp_par[[l]])
        df_lv$lv <- l
        resp_df <- rbind(resp_df, df_lv)
      } else {
        df_lv <- sp_params()
        df_lv$lv <- l
        resp_df <- rbind(resp_df, df_lv)
      }
    }
  }

  ## simulate drifts
  if (length(sim_drift_par[[1]]) != 1) {
    site_df <- do.call(sim_site_drift, sim_drift_par[[1]])
  } else {
    site_df <- sim_site_drift()
  }
  site_df$lv <- 1
  if (num_lv > 1) { # for more than 1 lv
    for (l in 2:num_lv) {
      if (length(sim_drift_par[[l]]) != 1) {
        df_lv <- do.call(sim_site_drift, sim_drift_par[[l]])
        df_lv$lv <- l
        site_df <- rbind(site_df, df_lv)
      } else {
        df_lv <- sim_site_drift()
        df_lv$lv <- l
        site_df <- rbind(site_df, df_lv)
      }
    }
  }
  list("sp_responses"=resp_df, "site_positions"=site_df)
}



# define the species parameters for 1 latent variable
# can have both linar and unodal response to latent variable
sp_params <- function(n_sp = 10,
                      coef_ic = 5:-4,
                      coef_lin = -4:5,
                      coef_q = -0.8) {

  data.frame(
    "intercept" = coef_ic,
    "slope" = coef_lin,
    "quad" = coef_q,
    "species" = 1:n_sp)
}



# simulate OU process (with optional shock) for 1 LV
sim_site_drift <- function (n_groups = 3, # number of treatment groups (incl. ref)
                            n_s_in_g = 10, # number of sites in each group (incl. ref)
                            # can be vector of numbers, with length = n_groups
                            n_t = 20, # number of time steps
                            t_end = 20, # end time
                            pos = c(0,0,0), # start position of sites
                            mean = 0, # process mean (reverting to this)
                            #sd_0 = c(.3,.3,.3), # standard error of groups initial distribution
                            sd = c(.5,.5,.5), # standard error of groups (per step)
                            theta = c(4,4,4), # phi params for the group
                            pre_mean = c(0,0,0), # mean before shock
                            shock = c(0,0,0),
                            shock_t = c(5,5,5)
) {

  # total number of sites
  tot_sites <- sum(n_s_in_g)*n_groups

  # time vector
  t <- seq(0,t_end,length.out=n_t)

  # total observations
  tot_obs <- tot_sites*length(t)

  # make out-data-frame
  site_df <- data.frame("site" = 1:tot_sites,
                        "group"=rep(1:n_groups, each=n_s_in_g)
  )
  sim_out <- data.frame("val"=numeric(),
                        "mean"=numeric(),
                        "t" = numeric(),
                        "site" = character(),
                        "group"=character()
  )

  # set mean for time series
  if (length(mean) == 1) {
    mean <- rep(mean,n_groups)
  }

  for (g in 1:n_groups) {
    # means
    if (shock[g]==0) {
      mn <- pos[g]*exp(-1/theta[g]*t) + mean[g]
    } else {
      mn <- c(rep(pre_mean[g], length(t[t<shock_t[g]])),
              (pre_mean[g] + shock[g])*exp(-1/theta[g]*t[1:(length(t)-length(t[t<shock_t[g]]))] + mean[g])
      )
    }
    # varcov for ou
    vc <- as.matrix(exp(-1/theta[g]*dist(t)))
    diag(vc) <- 1
    vc <- (sd[g]^2)*vc
    # simulate ou trajectories for group
    if (length(n_s_in_g) == 1) {
      g_sims <- mnormt::rmnorm(n_s_in_g,mn,vc)
    } else {
      g_sims <- mnormt::rmnorm(n_s_in_g[g],mn,vc)
    }
    for (s in 1:nrow(g_sims)) {
      # add trajectories to out-data frame
      sim_out <- rbind(sim_out,
                       data.frame("val"=g_sims[s,],
                                  "mean"=mn,
                                  "t" = t,
                                  "site" = site_df$site[site_df$group==g][s],
                                  "group"=g
                       ))
    }
  }
  sim_out
}


# simulate "true" hitting time distribution for OU process
sim_site_HT <- function (n_groups = 3, # number of treatment groups (incl. ref)
                            n_s_in_g = 10, # number of sites in each group (incl. ref)
                            # can be vector of numbers, with length = n_groups
                         t = NULL,
                            n_t = 20, # number of time steps
                            t_end = 20, # end time
                            pos = c(0,0,0), # start position of sites
                            mean = 0, # process mean (reverting to this)
                            #sd_0 = c(.3,.3,.3), # standard error of groups initial distribution
                            sd = c(.5,.5,.5), # standard error of groups (per step)
                            theta = c(4,4,4), # phi params for the group
                            pre_mean = c(0,0,0), # mean before shock
                            shock = c(0,0,0),
                            shock_t = c(5,5,5),
                         N=2000
) {

  # total number of sites
  tot_sites <- sum(n_s_in_g)*n_groups

  # time vector
  t <- seq(0,t_end,length.out=n_t)


  # total observations
  tot_obs <- tot_sites*length(t)

  # make out-data-frame
  sim_out <- data.frame("group"=numeric(),
                        "value"=numeric())

  # set mean for time series
  if (length(mean) == 1) {
    mean <- rep(mean,n_groups)
  }

  for (g in 1:n_groups) {
    # means
    if (shock[g]==0) {
      mn <- pos[g]*exp(-1/theta[g]*t) + mean[g]
    } else {
      mn <- c(rep(pre_mean[g], length(t[t<shock_t[g]])),
              (pre_mean[g] + shock[g])*exp(-1/theta[g]*t[1:(length(t)-length(t[t<shock_t[g]]))] + mean[g])
      )
    }
    # varcov for ou
    vc <- as.matrix(exp(-1/theta[g]*dist(t)))
    diag(vc) <- 1
    vc <- (sd[g]^2)*vc
    # simulate ou trajectories for group
    g_sims <- mnormt::rmnorm(N,mn,vc)
    for (r in 1:nrow(g_sims)) {
      hit_index <- as.numeric(which(g_sims[r,] >= -sd[g] & g_sims[r,] <= sd[g]))[1]
      sim_out <- rbind(sim_out,
                       data.frame("group"=g,
                                  "value"=t[hit_index]))

    }
  }
  sim_out
}


