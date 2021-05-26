#*******************************************************************************
#
# Project: Medical University Innsbruck - Intrauterine insemination (IUI)
# Date:    2020-04-15
# Author:  Patrick Rockenschaub
# Purpose: Fit the model
#
#******************************************************************************* 

#' ---
#' title: "Model the data"
#' subtitle: "Medical University Innsbruck - Intrauterine insemination (IUI)"
#' output: pdf_document
#' ---

#+ setup, message = FALSE 
source("00_init_workspace.R")
source(glue("{.dir_src}/imputation_helpers.R"))
source(glue("{.dir_src}/modelling_helpers.R"))
source(glue("{.dir_src}/misc_helpers.R"))

library(lspline)
library(rsample)
library(furrr)

plan(multicore)

mode <- "simple"

data <- read_rds(glue("{.dir_der}/data.rds"))
imp <- read_rds(glue("{.dir_der}/imp.rds"))


#' # Create an increasingly complex model

#+ choose-formulas
base <- y ~ timeInt

if (mode == "simple") {
  formula1 <- base %>% 
    update.formula(~ . + I(amh < 1) + I(sperm >= 5 & sperm < 15) + I(sperm < 5))
  formula2 <- formula1 %>% 
    update.formula(
      ~ . + I(!diagnosis %in% c("anovulatory", "no_known_female_inf"))
    )
  formula3 <- formula2 %>% 
    update.formula(~ . + I(age > 35) + I(bmi > 30))
} else {
  formula1 <- base %>% 
    update.formula(~ . + lspline(amh - 1, 0) + lspline(sperm - 15, 0))
  formula2 <- formula1 %>% 
    update.formula(~ . + diagnosis)
  formula3 <- formula2 %>% 
    update.formula(~ . + I(age - 30) + I(bmi - 20))
}

formulas <- list(formula1, formula2, formula3)


#+ run-models
fits <- future_map(
  .x = formulas,
  .f = fit_model_on_imputation,
  imp = imp,
  .options = furrr_options(seed = TRUE)
)

coefs <- map(fits, get_coefs)
cs <- map(fits, get_cindex, imp, imp)
aics <- map(fits, get_aic)
bics <- map(fits, get_bic)

for(i in seq_along(fits)) {
  print(coefs[[i]])
  print(
    coefs[[i]] %>%
      mutate(
        lower = estimate + qnorm(0.025) * std.error, 
        upper = estimate + qnorm(0.975) * std.error
      ) %>% 
      mutate(
        across(all_of(c("estimate", "lower", "upper")), ~ round(exp(.), 3)),
        ci = glue("{estimate} ({lower}-{upper})")
      )
  )
  print(cs[[i]]) # only apparent performance
  print(aics[[i]])
  print(bics[[i]])
}

write_rds(fits, glue("{.dir_res}/fits.rds"))



#' # Perform internal validation via optimism-adjusted bootstrap

#+ bs-sample
set.seed(999)
boot <- bootstraps(
    data %>% add_imputation_variables(), 
    times = 100, 
    strata = "hcg"
  )

#+ bs-impute
boot %<>% mutate(
  imp_train = future_map2(
    .x = splits, 
    .y = 1:length(splits),
    .f = ~ run_imputation(analysis(.x), m = 10, maxit = 20, seed = 357 * .y),
    .options = furrr_options(seed = TRUE)
  ),
  imp_test = future_map2(
    .x = imp_train,
    .y = splits, 
    .f = ~ mice.mids(.x, newdata = assessment(.y), maxit = 20),
    .options = furrr_options(seed = TRUE)
  )
)

#+ bs-fit
for (i in seq_along(formulas)) {
  boot %<>% mutate(
    "fits{i}" := future_map(
      .x = imp_train, 
      .f = fit_model_on_imputation,
      formula = formulas[[i]],
      .options = furrr_options(seed = TRUE)
    )
  )
}

#+ bs-eval
for (i in seq_along(formulas)) {
  boot %<>% mutate(
    # "Apparent" bootstrap performance
    "c_bs{i}" := future_pmap(
      .l = list(get(glue("fits{i}")), imp_train, imp_train),
      .f = get_cindex,
      .options = furrr_options(seed = TRUE)
    ), 
    # "Test" bootstrap performance
    "c_orig{i}" := future_pmap(
      .l = list(get(glue("fits{i}")), imp_train),
      .f = get_cindex,
      imp_test = imp,   # original non-bootstrapped data
      .options = furrr_options(seed = TRUE)
    )
  )
}

for (i in seq_along(formulas)) {
  print(glue("Formula {i}:"))
  opt <- map_dbl(boot[[glue("c_bs{i}")]], "c") - map_dbl(boot[[glue("c_orig{i}")]], "c")
  print(cs[[i]]["c"] - mean(opt))
  print(sd(opt))
  
  print("")
}


