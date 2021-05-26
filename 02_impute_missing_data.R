#*******************************************************************************
#
# Project: Medical University Innsbruck - Intrauterine insemination (IUI)
# Date:    2020-04-15
# Author:  Patrick Rockenschaub
# Purpose: Inspect and impute missing data using MICE
#
#******************************************************************************* 

#' ---
#' title: "Inspect and impute missing data"
#' subtitle: "Medical University Innsbruck - Intrauterine insemination (IUI)"
#' output: pdf_document
#' ---


#+ setup, message = FALSE 
source("00_init_workspace.R")
source(glue("{.dir_src}/imputation_helpers.R"))

library(nnet)

#+ load-data
data <- read_rds(glue("{.dir_der}/data.rds"))


#' # Describe data distribution and missingness

#+ overall-pattern-of-missingness
invisible(
  md.pattern(
    data[, names(data)[names(data) != "drug"]], 
    rotate.names = TRUE,
  )
)

#+ define-models
lin_reg <- function(f) glm(f, data = data, family = gaussian) %>% summary()
log_reg <- function(f) glm(f, data = data, family = binomial) %>% summary()
mnl_reg <- function(f) multinom(substitute(f), data = data)


#' ## Diagnoses
#' 

#+ describe-diagnosis
summary(data$diagnosis)

(is.na(diagnosis) ~ hcg + age + log(bmi) + log(sperm) + log(amh) + cycle) %>% 
  log_reg() # Only outcome slightly predictive of missingness

multinom(
  diagnosis ~ hcg + age + log(bmi) + log(sperm) + log(amh) + cycle,
  data = data
) %>% 
  summary()


#' ## Anti-mueller hormone (AMH)
#' 

summary(data$amh)
hist(data$amh)
hist(log(data$amh))

#' For AMH, no variable predicts whether or not the value is missing. This 
#' strengthens the case for it being missing completely at random.

(is.na(amh) ~ hcg +  age + log(bmi) + diagnosis + log(sperm) + cycle) %>% 
  log_reg() # No variable predicts missingness

#' There were however a couple of variables that predicted the value of 
#' AMH were it was missing. These included the outcome (+), higher age (-), 
#' higher BMI (-), and a diagnosis of anovulation (+).

(log(amh) ~ hcg + age + log(bmi) + diagnosis + log(sperm) + cycle) %>% 
  lin_reg() 


#' ## Sperm count
#' 

summary(data$sperm)
hist(data$sperm)
hist(log(data$sperm))

#' For sperm counts, no variable predicts whether or not the value is missing. This 
#' strengthens the case for it being missing completely at random.

(is.na(sperm) ~ hcg + age + log(bmi) + diagnosis + log(amh) + cycle) %>% 
  log_reg() 

#' There were however a couple of variables that predicted the value of 
#' spermm count were it was missing. These included the outcome (+), higher 
#' age of the female patient (+), higher BMI (-), and a female diagnosis of 
#' endometriosis or tubal obstruction (+). All other diagnosis also showed a 
#' positive but non-significant relationship compared to those patients where
#' no female infertility is known. This reflects the large proportion of 
#' patients were male infertility is the underlying problem, as confirmed by 
#' additional analyses below.

(log(sperm) ~ hcg + age + log(bmi) + diagnosis + log(amh) + cycle) %>% 
  lin_reg()  # Positive: Age, hcgpositive and all diagnoses compared to no known reason

(log(sperm) ~ hcg + age + log(bmi) + log(amh) + cycle + 
    I(diagnosis != "no_known_female_inf")) %>% 
  lin_reg()

(log(sperm) ~ hcg + age + log(bmi)  + log(amh) + cycle + 
  I(diagnosis != "no_known_female_inf") + 
  I(!is.na(female_male) & female_male %in% c("male", "female and male"))) %>% 
  lin_reg()


#' ## BMI
#' 

summary(data$bmi)
hist(data$bmi)
hist(log(data$bmi))

#' For sperm counts, no variable predicts whether or not the value is missing. This 
#' strengthens the case for it being missing completely at random.

(is.na(bmi) ~ hcg + age + log(sperm) + diagnosis + log(amh) + cycle) %>% 
  log_reg()

#' For BMI, a diagnosis of anovulation (+) and AMH (-) were associated with 
#' the observed value.

(log(bmi) ~ hcg + age + log(sperm) + diagnosis + log(amh) + cycle) %>% 
  lin_reg() 


# Follicles
summary(data$follicles)
hist(data$follicles)

#' For follicles, there was only a weak association with cycle number.
(is.na(follicles) ~ hcg + age + log(sperm) + diagnosis + log(amh) + cycle) %>% 
  log_reg()

#' and no association with any of the other variables
(follicles > 1 ~ hcg + age + log(sperm) + diagnosis + log(amh) + log(bmi) +  cycle) %>% 
  log_reg() 



#' # Impute missing values
#' 
#' Before imputation, we create separate variables for AMH, sperm count, and 
#' BMI that are log-transformed versions of the original variables in order 
#' to make more normally distributed. These logged versions will be imputed 
#' via the MICE algorithm (using PMM) and the original variables will then 
#' be passively imputed via an inverse transformation (= exp()). In addition, 
#' a new variable known_male_factor is created with is only used in the 
#' imputation of sperm count.
#' 

#+ transform-for-imputation
data <- data %>% 
  add_imputation_variables()

#' Using the newly added variables, we perform multiple imputation of 10 
#' imputed datasets with 20 iterations of the algorithm for convergence.

#+ impute
imp <- run_imputation(data, m = 10, maxit = 20, print = TRUE, seed = 42)

#' After the algorithm finishes, we inspect the means and standard deviations
#' of each actively imputed variable for convergence and ensure ourselves that 
#' the imputed values appear sensible.

#+ inspect-convergence
plot(imp)

#+ inspect-imputed-value
bwplot(imp, log_bmi + log_sperm + log_amh ~ .imp)
densityplot(imp, .imp ~ log_bmi + log_sperm + log_amh)


#' Finally, we save the imputed datasets.

#+ save-data
write_rds(imp, glue("{.dir_der}/imp.rds"))

