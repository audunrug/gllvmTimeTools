#' Setup model arguments for a gllvm model with stationary and/or non-stationary ou-processes on the site level
#'
#' @param X_time Numeric; time variable
#' @param X_group Factor / character; group levl
#' @param X_sampling_units factor; sampling units that are correlatied over time (e.g. permanent plots)
#' @param X_other vector or dataframe; auxilliary environmental variables to include in latent variable formula
#' @param reference numeric or character; which group, if any, is the reference
#' @param impact_time numeric or vector; is there an impact at some point in the time series (NULL=decay starts at first time point)
#' @param sd_other vector or dataframe; auxilliary study design variable(s) (for row.eff, typically)
#' @param LTIs logical; estimate group-specific long-term intercepts separate from the reference
#' @param num_lv number of latent variables in the model
#' @param X_other_lv other predictors to include as canonical coefficients
#' @param decaying which groups should exhibit exponential decays (can include reference)
#'
#' @returns list of components that correspond to arguments for the "gllvm" function.
#' @export
#'
#' @examples
#'
setup_model_arguments <- function(X_time,
                                  X_group,
                                  X_sampling_units = NULL,
                                  num_lv = 1,
                                  X_other_lv = NULL,
                                  X_other = NULL,
                                  sd_other = NULL,
                                  LTIs = FALSE, # group-specific long-term intercepts
                                  reference = NULL, # which group is the reference
                                  decaying = NULL, # which groups exhibit exponential decays (can include reference)
                                  impact_time = 0 # is there an impact at some point in the time series
                                  # could also be a vector with values for each treatment group
                                  # NB must be in alphabetical order of non-reference factors
) {

  # mini helper function
  group_by_ref <- function (groups,reference) {
    if (!is.null(reference)) {
      nonref_groups <- levels(as.factor(groups))[!(levels(as.factor(groups)) %in% reference)]
      groups <- factor(groups, levels=c(reference, nonref_groups))
    } else { # if no reference level
      groups <- factor(groups)
    }
  }

  # setup time and pre-treatment groups (non-reference groups in alphabetical order)
  gs_dec <- sort(decaying)
  gs <- sort(unique(X_group[!(factor(X_group) %in% reference)]))
  X_group_sd <- group_by_ref(X_group, reference) # save studydesign sd group to preserve covariance stationarity
  X_group <- as.character(X_group) # temporary conversion to char for pre-group addition

  if (!is.null(impact_time)) {
    if (length(impact_time)==1) {
      X_time_0 <- X_time - impact_time
      for (g in gs_dec) {
        X_group[X_group==g & X_time_0<0] <- paste0("pre_",g)}
    } else {
      # allow impacts at different times for different groups
      X_time_0s <- cbind(sapply(impact_time,y=X_time,function(x,y){y-x}))
      for (i in 1:ncol(X_time_0s)) {
        g <- gs_dec[i]
        X_group[X_group==g & X_time_0s[,i]<0] <- paste0("pre_",g)
      }
    }
  } else {X_time_0 <- X_time}

  # setup group factors
  X_group <- group_by_ref(X_group, reference)

  #if there are pre-disturbance time points for the treatment groups, add that
  # data frame of design matrix
  X_df <- data.frame(
    "time" = X_time,
    "group" = X_group,
    "group_sd" = X_group_sd)

  if (!(is.null(X_sampling_units))) {
    X_df$unit <- X_sampling_units # add sampling units to DF

    # make studyDesign matrix
    sd_df <- X_df |>
      mutate(time = factor(time),
             unit = factor(unit)) |>
      arrange(unit,group_sd,time) |>
      select(-group)


    # model matrix for time X group (minus reference)
    timeVar <- c()
    if (length(impact_time)==1) {
      for (g in gs_dec) {
        n <- paste0("g",g, "_X_time")
        X_df[,n] <- X_time_0*(X_df$group==g)
        timeVar <- c(timeVar, n)
      }
      #X_ts_mod <- model.matrix(~X_df$group_sd:as.numeric(as.character(X_time_0))) |>
      #  as.data.frame()
    }
    if (length(impact_time)>1) { # if impact time is different for different groups
      #X_ts_mod <- matrix(nrow=length(X_group), ncol=length(gs))
      for (i in 1:length(gs_dec)) {
        n <- paste0("g",gs_dec[i], "_X_time")
        X_df[,n] <- X_time_0s[,i]*(X_group==gs_dec[i])
        timeVar <- c(timeVar, n)
      }
    }

    X_df[,-c(1:4)][X_df[,-c(1:4)]<0] <- 0 # disallow negative values
    #X_ts_mod <- lapply(X_ts_mod, factor(levels=levels(sd_df$time)))

    #colnames(X_ts_mod) <- paste0("g",gs, "_X_time")
    X_df <- X_df |> select(-group_sd)

    # make distance design matrix
    distLV <- matrix(as.numeric(as.character(sd_df$time)),ncol=1)

    # make rho map
    # bit tricky if you have more than one LV
    sd_rho <- sd_df |> select(unit, group_sd) |>
      mutate(group=as.numeric(group_sd)) |>
      unique()
    if (num_lv == 1) {
      setMap <- list(rho_lvc = sd_rho$group_sd)}
    else if (num_lv > 1) {
      sd_rho_list <- c(sapply(1:num_lv, y=sd_rho$group, FUN = function(x,y){y + max(y)*(x-1)}))
      sd_rho_mat <- t(matrix(ncol=num_lv, sd_rho_list))
      setMap <- list(rho_lvc = as.numeric(c(sd_rho_mat)))
    }

    # make final model matrix
    if (LTIs) {
      X_df_out <- X_df |> select(-time,-unit)
    } else {
      X_df_out <- X_df |> select(-group,-time,-unit)
    }

    # add additional environmental predictors if any are specified
    if (!(is.null(X_other))) {X_df_out <- cbind(X_df_out,X_other)}

    # add additional SD predictors if any are specified
    if (!(is.null(sd_other))) {sd_df <- cbind(sd_df,sd_other)}

    # add formula objects
    if (LTIs) {X_terms <- colnames(X_df_out)} else {X_terms <- colnames(X_df_out)}
    form_X <- as.formula(paste("~ ", paste(X_terms, collapse = " + ")))
    form_sd <- ~ corExp(0 + time | unit)

    # return
    return (list(X=X_df_out,
                 studyDesign=sd_df,
                 lv.formula=form_X,
                 lvCor=form_sd,
                 timeVar=timeVar,
                 distLV=distLV,
                 setMap=setMap))}
}
