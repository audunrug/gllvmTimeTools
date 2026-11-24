#' Simulate hitting times from fitted gllvm object
#'
#' @param object fitted gllvm model object
#' @param rangeX vector of timesteps to predict for
#' @param group_names names of groups defined by the next arguments
#' @param cc_init canonical coefficient start values to predict for (indices)
#' @param cc_LTI canonical coefficient long-term intercepts (indices)
#' @param cc_ref indexes of LTIs to compare trajectory with. If FALSE, treats 0 as the reference
#' @param rho rho (indices)
#' @param N number of iterations
#' @param stdevs standard deviation for the reference range, defined by sigma_lv
#'
#' @returns dataframe with 2 columns: group = group name from group_names, value = first hitting times for each iteration X group. If no hitting time recorded, value = NA.
#' @export
#'
#' @examples
sim_HT <- function(object, #gllvm
                   rangeX, #prediction range
                   group_names, # names of treatment groups
                   cc_init, # canonical coefficient start values to predict for (indices)
                   cc_LTI=NULL, # canonical coefficient long-term intercept (indices)
                   cc_ref=0,
                   constant_rho=F,
                   rho, # rho (indices)
                   stdevs = 1,
                   N=2000 # number of runs
) {

  # variance-covariance of the data
  incl <- object$Hess$incl
  Vf <- vcov(object)
  sims <- data.frame(x=rangeX)
  #hits = c()
  df <- data.frame(group=character(), value=numeric())
  for (i in 1:N) {
    sim_full <- try(MASS::mvrnorm(1,object$TMBfn$par[incl],Vf))

    # loop through group parameters
    for (i in 1:length(group_names)) {
      cc_init_t <- sim_full[names(sim_full)=="b_lv"][cc_init[i]]
      if (is.null(cc_LTI)) {
        cc_LTI_t <- 0
      } else {
        cc_LTI_t <- sim_full[names(sim_full)=="b_lv"][cc_LTI[i]]
      }
      if (constant_rho) {
        rho_t <- object$TMBfn$par[names(sim_full)=="rho_lvc"][rho[i]]
        rho_t <- exp(-1/(exp(rho_t)))
      } else {
        rho_t <- exp(-1/(exp(sim_full[names(sim_full)=="rho_lvc"][rho[i]])))
      }

      sd_lv <- sim_full[names(sim_full)=="sigmaLV"]

      if (is.null(cc_LTI_t)) {
        temp <- cc_init_t*(rho_t^(rangeX))
      } else {
        temp <- cc_LTI_t + cc_init_t*(rho_t^(rangeX))
      }

      mvsd <- as.matrix(rho_t^dist(rangeX,diag=T,upper=T))
      diag(mvsd) <- diag(mvsd) + 1
      temp <- temp + MASS::mvrnorm(1,rep(0,length(rangeX)),sd_lv*mvsd)

      if (cc_ref !=0) {
        cc_ref <- sim_full[names(sim_full)=="b_lv"][cc_ref[i]]
      }

      hit_index <- as.numeric(which(temp >= (cc_ref-sd_lv) & temp <= (cc_ref+sd_lv)))[1]
      df <- df |> rbind(data.frame("group"=group_names[i],"value"=rangeX[hit_index]))
    }
  }
  df
}
