#*******************************************************************************
#
# Project: Medical University Innsbruck - Intrauterine insemination (IUI)
# Date:    2020-04-15
# Author:  Patrick Rockenschaub
# Purpose: Functionality to model and evaluate discrete survival models
#
#******************************************************************************* 

library(discSurv)

#' Wrapper around discsurv::dataLong, transforming a dataset with one row per
#' patient into a dataset with one row per timestep (in our case cycle).
#' 
#' @param data data.frame with discrete time date; must contain a variable 
#'   `cycle` (integer time index > 0) and `hcg` (binary outcome 0/1 indicating)
#'   the outcome at the last observed time step, were 0 is negative (=censored).
#'  
#' @return data.frame 
reshape_to_long <- function(data) {
  dataLong(
    dataSet = as.data.frame(data), 
    timeColumn = "cycle", 
    censColumn = "hcg", 
    timeAsFactor = FALSE
  )
}


#' Given a multiple imputed dataset, fit a discrete survival model to each 
#' imputed dataset (uses a simple logistic regression as the model of 
#' choice). 
#' 
#' @param imp imputations as mice::mids object
#' @param formula analysis formula
#' 
#' @return mice::mira object
fit_model_on_imputation <- function(imp, formula) {
  call <- match.call()
  if (!is.mids(imp)) {
    stop("The data must have class mids")
  }
  analyses <- as.list(seq_len(imp$m))
  for (i in seq_along(analyses)) {
    data.i <- imp %>% 
      complete(i) %>% 
      reshape_to_long()
    analyses[[i]] <- glm(formula, data = data.i, family = binomial())
  }
  object <- list(call = call, call1 = imp$call, nmis = imp$nmis, 
                 analyses = analyses)
  oldClass(object) <- c("mira", "matrix")
  object
}


#' Extract pooled model coefficients from a model fit on imputed data
#' 
#' @param fits mice::mira object containing the fits for each imputed dataset
#' 
#' @return model summary
get_coefs <- function(fits) {
  summary(pool(fits))
}



#' Exponentiate model coefficients and calcualte 95% confidence intervals
#' 
#' @param df mice::mipo.summary
#' 
#' @return `df` with added modified columns `estimate`, `lower`, and `upper`
calculate_or <- function(df) {
  df %>% 
    mutate(
      lower = estimate + qnorm(0.025) * std.error, 
      upper = estimate + qnorm(0.975) * std.error
    ) %>% 
    mutate(across(all_of(c("estimate", "lower", "upper")), ~ round(exp(.), 3)))
}



#' Calculate standard error for the C-statistic based on Hanley's formula
#' 
#' @param x estimated c-index
#' @param y observed binary outcomes
#' 
#' @return standard error as scalar
#' 
#' @note copied from Snell et al. (2021)
calculate_hanley_se <- function(x, y){
  n1 <- sum(y == 0)
  n2 <- sum(y == 1)
  Q1 <- x / (2 - x)
  Q2 <- 2 * x^2 / (1 + x)
  sqrt((x*(1-x) + (n1-1)*(Q1-x^2)+(n2-1)*(Q2-x^2)) / (n1*n2))
}


#' Wrapper around discSurv::evalCindex to calculate the c-index for a 
#' discrete survival model.
#' 
#' @param single_fit single glm model object
#' @param train data.frame used for training
#' @param test data.frame on which to evaluate
#' 
#' @return estimated c-index as scalar
calculate_cindex <- function(single_fit, train, test) {
  preds <- predict(single_fit, newdata = test)
  evalCindex(
    marker = preds, 
    newTime = test$timeInt, 
    newEvent = test$y, 
    trainTime = train$timeInt, 
    trainEvent = train$y
  )
}


#' Get the c-index for a model train on imputed data
#' 
#' @param fits mice::mira object containing the fits for each imputed dataset
#' @param imp_train mice::mids object with the imputed data used for training
#' @param imp_test mice::mids object with the imputed data to evaluate on
#' 
#' @return pooled results as tibble with columns `auc` and `auc_se`
get_cindex <- function(fits, imp_train, imp_test) {
  # Format the data for discSurv::evalCindex used in calculate_cindex()
  format_data <- function(data, which) {
    data %>% mutate(
      timeInt = cycle, 
      y = hcg
    )
  }
  
  # Calculate the c-indices
  cs <- pmap_dbl(
    .l = list(
      fits$analyses, 
      complete(imp_train, "all"),
      complete(imp_test, "all")
    ),
    .f = ~ calculate_cindex(..1, format_data(..2), format_data(..3))
  )
  
  # Pool the results using Rubin's rules
  qbar = mean(cs)
  ubar = mean(map_dbl(cs, calculate_hanley_se, complete(imp_test)$hcg) ^ 2)
  b = var(cs)
  
  tibble(c = qbar, c_se = sqrt(ubar + b))  
}


#' Calculate 95% confidence intervals for the AUC
#' 
#' @param df tbl_df with columns `c` and `c_se`
#' 
#' @return `df` with added modified columns `c`, `lower`, and `upper`
calculate_auc_ci <- function(df) {
  df %>% 
    mutate(
      lower = c + qnorm(0.025) * c_se, 
      upper = c + qnorm(0.975) * c_se
    )
}


#' Get the average Akaike Information Criterion for a set of models fit on 
#' multiple imputated data. 
#' 
#' @param fits mice::mira object containing the fits for each imputed dataset 
get_aic <- function(fits) {
  mean(map_dbl(fits$analyses, AIC))
}

#' Get the average Bayesian Information Criterion for a set of models fit on 
#' multiple imputated data. 
#' 
#' @param fits mice::mira object containing the fits for each imputed dataset 
get_bic <- function(fits) {
  mean(map_dbl(fits$analyses, BIC))
}
