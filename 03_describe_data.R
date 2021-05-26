#*******************************************************************************
#
# Project: Medical University Innsbruck - Intrauterine insemination (IUI)
# Date:    2020-04-15
# Author:  Patrick Rockenschaub
# Purpose: Perform descriptive and exploratory analysis of the study data
#
#******************************************************************************* 

#' ---
#' title: "Describe the data"
#' subtitle: "Medical University Innsbruck - Intrauterine insemination (IUI)"
#' output: pdf_document
#' ---
#' 
#+ setup, message = FALSE 
source("00_init_workspace.R")
source(glue("{.dir_src}/modelling_helpers.R"))
source(glue("{.dir_src}/misc_helpers.R"))

library(tableone)
format <- "markdown"

#+ load-data
data <- read_rds(glue("{.dir_der}/data.rds"))

#' # Table 1
#' 
#' We begin by creating the patient characteristics table shown in the published
#' article. That table is using the cleaned (but not imputed dataset) and 
#' summarises medians and counts of the data.
#' 

#+ create-table-one
covars <- c("age", "bmi", "diagnosis", "amh", "amh_cat", "female_male", "sperm", "sperm_cat", "drug", "cycle")

tbl1 <- CreateTableOne(
  covars, "hcg", 
  data = data %>% 
    mutate(
      cycle = pmin(cycle, 4),
      amh_cat = factor(amh < 1, c(TRUE, FALSE), c("<1", ">1")),
      sperm_cat = cut(sperm, c(-Inf, 5, 15, Inf), c("0-5", "5-15", ">15"))
    ), 
  factorVars = "cycle",
  addOverall = TRUE
)

#+ print-table-one
tbl1 %>% 
  print(
    missing = TRUE, 
    nonnormal = c("bmi", "amh", "sperm", "cycle"),
    showAllLevels = TRUE
  ) %>% 
  kableone(
    format = format,
    booktabs = TRUE
  )



#' # Investigate relationship with the outcome
#' 
#' So far, the data contained one row per patient. This row contained the 
#' insemination cycle for which the data was collected, which could either 
#' result in pregnancy (= outcome) or not. For each patient, either the first
#' cycle during which they got pregnant was chosen or their last observed 
#' cycle. This has two implications. First, patients who have a negative 
#' recorded outcome were censored after their cycle evaluated here (e.g., 
#' because they moved on to other treatments like IVF). Second, patients who
#' were evaluated after their first cycle must have had previous (unsuccessful)
#' tries. These also need to be taken into account or otherwise the analysis
#' will be biased. We therefore make those cycles explicit and bring the data 
#' into a long format (i.e., where there is one row per cycle).

#+ reshape-to-long
data_long <- data %>%
  reshape_to_long() %>% 
  add_cycle_date()

#' We now use this long format to investigate the relationship of individual 
#' predictors and the outcome.
#' 
#' ## Time
#' 
#' The success rate might have changed over the course of the study period. We
#' therefore fit a linear as well as a LOESS smoother to the data to look 
#' for any observable patterns.

#+ plot-date
plot_preg <- plot_preg %+% data_long
plot_preg %+% 
  aes(x = cycle_date) + 
  geom_smooth(method = "loess")

#' While there appears to be a slight decrease in success rate over the course
#' of the study period, the uncertainty around that decline is too large to 
#' confidently conclude a decline.

#+ fit-date
fit_time <- fit_preg(~ time_length(ymd("2013-01-01") %--% cycle_date, u = "y"))
summary(fit_time)
predict_preg(fit_time, tibble(cycle_date = ymd("2008-01-01", "2018-01-01")))


#' ## Cycle
#' 
#' The success rate may further differ by how often an IUI was tried already. 
#' Since this is a discrete problem with very low numbers of tries > 4, we use 
#' only linear regression (and LOESS, which would be too flexible and overfit)
#' to investigate any trends.

#+ tab-cycle
tabulate_preg(timeInt)

#+ plot-cycle
plot_preg %+% aes(x = timeInt)

#' Looking only at the outcomes, there hardly seem to be any pregnancies after
#' the fourth cycle. Other authors have used this argument to claim a reduced
#' effectiveness of IUI after that point. However, this would be an erroneus 
#' conclusion, since the number of patients that proceed to have more than 
#' 4 IUIs is extremely small and we would therefore expect only a few success. 
#' Let's investigate this issue further by first fitting a linear regression to 
#' the data. 

#+ fit-cycle
fit_cycle <- fit_preg( ~ timeInt)
summary(fit_cycle)
predict_preg(fit_cycle, tibble(timeInt = 1:9))

#' As already indicated in the plot, the estimated effect decrease over time, 
#' although in principle negative, is highly non-significant. While there might
#' be some negative effect over time, we are highly uncertain about its 
#' magnitude (and it might even be strongly positive!). We can see this further
#' by creating table 2 shown in the manuscript, which adds expected upper and 
#' lower limits for the number pregnancies that we would expect to see if the 
#' study-wide success rate of ~10% were constant across cycles.
#' 

#+ table-2-format
format_table_two <- . %>% 
  select(timeInt, total, total_percent, everything()) %>% 
  set_names(c("Cycle", "N patients", "% patients", 
              "N pregn", "% pregn", "cum-% preg",
              "5% expected", "95% expected")) %>% 
  t() %>% 
  kable(
    format = format,
    booktabs = TRUE
  )

#+ table-2
data_long %>% 
  group_by(timeInt) %>% 
  summarise(
    total = n(), 
    pregn = sum(y)
  ) %>% 
  calc_p_and_ci() %>% 
  format_table_two()

fisher.test(with(data_long, table(timeInt, y)))

#' We can see that all the values lie within the expected limits. The same 
#' is true for Khalil et al. (2001), which we argue incorrectly claimed 
#' otherwise. We do not that their observed success rate for the first cycle
#' is slightly higher than that of all other cycles, but this might be due 
#' to a selection bias or at least does not support a plateau after 3 or 4 
#' cycles.

#+ table-2-comparion
tibble(
  timeInt = c(1:7),
  total = c(893, 661, 498, 212, 105, 64, 40), 
  pregn = c(130, 72, 50, 19, 13, 5, 5)
) %>%
  calc_p_and_ci() %>% 
  format_table_two()


#' ## Predictors
#' 
#' The following analyses show similar investigations for all predictors
#' included in the later models.
#' 
#' ### Anti-mueller hormone

#+ tab-miss-amh
tabulate_preg(is.na(amh))

#+ plot-amh
plot_preg %+% 
  aes(x = log(amh)) + 
  geom_smooth(method = "loess")

#' #### Sperm count

#+ tab-miss-sperm
tabulate_preg(is.na(sperm))

#+ plot-sperm
plot_preg %+% 
  aes(x = log(sperm)) + 
  geom_smooth(method = "loess") + 
  geom_vline(xintercept = log(c(7.5, 15)))


#' ### BMI

#+ tab-miss-bmi
tabulate_preg(is.na(bmi))

#+ plot-bmi
plot_preg %+% 
  aes(x = bmi) + 
  geom_smooth(method = "loess")


#' ### Age

#+ tab-age
plot_preg %+% 
  aes(x = age) + 
  geom_smooth(method = "loess")





