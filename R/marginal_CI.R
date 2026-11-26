#' Confidence intervals for marginal latent variables
#'
#' @description
#' Simple function for extracting a bootstrap confidence interval for only the marginal latent variable scores
#'
#'
#' @param object gllvm object
#' @param N numer of iterations
#' @param level significance level
#'
#' @returns
#' @export
#'
#' @examples
marginal_CI <- function(object, #gllvm
                        N=2000, # number of runs
                        level=0.05 # CI cutoff level
) {

  # NB: could be made more sophisticated in the future?
  params <- gllvm:::simulate_params_gllvm(object, N, 1, level=0)
  df <- data.frame(u=numeric())
  # make columns
  for (i in 1:mod_coral$num.lv.c) {df[[paste0("LV_",i)]] <- numeric()}

  rws <- nrow(object$y)
  for(n in 1:N){
    newobject <- gllvm:::perturb.gllvm(object, params, n)
    margs <- gllvm::getLV(newobject,type="marginal")
    df_temp <- data.frame(u=1:rws)
    for (i in 1:object$num.lv.c){
      df_temp[[paste0("LV_",i)]] <- margs[,i]
    }
    df <- df |> rbind(df_temp)
  }

  df_lw <- df |>
    group_by(u) |>
    summarise_if(is.numeric,quantile,probs=(0+(level/2)))
  df_up <- df |>
    group_by(u) |>
    summarise_if(is.numeric,quantile,probs=(1-(level/2)))

  df_out <- list("upper"=df_up, "lower"=df_lw)
  df_out
}
#
#
#
#
#
#
#
#
#
#
#
#   # variance-covariance of the data
#   incl <- object$Hess$incl
#   Vf <- vcov(object)
#   sims <- data.frame(x=rangeX)
#   #hits = c()
#   df <- data.frame(group=character(), time=numeric(), value=numeric())
#   for (i in 1:n) {
#     sim_full <- try(MASS::mvrnorm(1,object$TMBfn$par[incl],Vf))
#
#     # loop through group parameters
#     for (i in 1:length(group_names)) {
#       cc_init_t <- sim_full[names(sim_full)=="b_lv"][cc_init[i]]
#       if (is.null(cc_LTI)) {
#         cc_LTI_t <- 0
#       } else {
#         cc_LTI_t <- sim_full[names(sim_full)=="b_lv"][cc_LTI[i]]
#       }
#       rho_t <- exp(-1/(exp(sim_full[names(sim_full)=="rho_lvc"][rho[i]])))
#       #sd_lv <- sim_full[names(sim_full)=="sigmaLV"]
#       if (is.null(cc_LTI_t)) {
#         temp <- cc_init_t*(rho_t^(rangeX))
#       } else {
#         temp <- cc_LTI_t + cc_init_t*(rho_t^(rangeX))
#       }
#       df <- df |> rbind(data.frame(group=group_names[i],
#                                    time=rangeX,
#                                    value=temp))
#     }
#   # summarize by quantiles
#   df_out <- df |>
#     group_by(group,time) |>
#     summarise(lwr = quantile(value,probs=(0+(level/2))),
#               upr = quantile(value,probs=(1-(level/2))))
#   df_out
# }
