#*******************************************************************************
#
# Project: Medical University Innsbruck - Intrauterine insemination (IUI)
# Date:    2020-04-15
# Author:  Patrick Rockenschaub
# Purpose: Functionality to perform multiple imputatio
#
#******************************************************************************* 

library(mice)


#' Add variables needed for the imputation process to a data.frame
#' 
#' @param data data.frame to which to add the variables; needs to contain the 
#'   variables `bmi`, `amh`, `sperm`, and `female_male`
#' @return data.frame
add_imputation_variables <- function(data) {
  data %>% 
    mutate(
      # Bring closer to normality
      log_bmi = log(bmi),
      log_amh = log(amh),
      log_sperm = log(sperm),
      
      # flag for patients were male infertility is a factor
      known_male_factor = (
        !is.na(female_male) & 
          female_male %in% c("male", "female and male")
      )
    )
  
}


#' Run one set of multiple imputations using MICE
#' 
#' @param data data.frame to be imputed
#' @param m number of resulting imputed datasets
#' @param maxit number of iterations for MICE
#' @param seed random number seed
#' @param print flag to indicate whether progress should be output to console
#' 
#' @return imputations as object of class mice::mids
run_imputation <- function(data, m = 10, maxit = 10, seed = 42, print = FALSE) {
  init <- mice(data, maxit = 0)
  
  # Define which predictors to use for active/passive imputation
  predmat <- init$predictorMatrix
  impute_vars <- c("hcg", "age", "log_bmi", "diagnosis", "log_amh", "log_sperm", "follicles")
  passive_vars <- c("bmi", "amh", "sperm")
  predmat[, !colnames(predmat) %in% impute_vars] <- 0
  predmat[!colnames(predmat) %in% c(impute_vars, passive_vars), ] <- 0
  predmat["log_sperm", "known_male_factor"] <- 1
  
  # Do not impute variables that weren't marked for active imputation using 
  # mice default methods
  methods <- init$method
  methods[!names(methods) %in% impute_vars] <- ""
  
  # Define passive imputation for those marked as passive
  methods["bmi"] <- "~ I(exp(log_bmi))"
  methods["amh"] <- "~ I(exp(log_amh))"
  methods["sperm"] <- "~ I(exp(log_sperm))"
  
  # Perform imputation
  imp <- mice(
    data, 
    m = m, 
    maxit = maxit, 
    meth = methods, 
    pred = predmat, 
    visitSequence = c(impute_vars, passive_vars),
    seed = seed,
    printFlag = print
  )
}
