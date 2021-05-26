#*******************************************************************************
#
# Project: Medical University Innsbruck - Intrauterine insemination (IUI)
# Date:    2020-04-15
# Author:  Patrick Rockenschaub
# Purpose: Convert the data collected by clinicians and stored as SPSS file
#          into a analysis-ready R file
#
#******************************************************************************* 

#' ---
#' title: "Converting data from SPSS to R"
#' subtitle: "Medical University Innsbruck - Intrauterine insemination (IUI)"
#' output: pdf_document
#' ---

#+ setup, message = FALSE 
source("00_init_workspace.R")

library(haven)

#' # Dataset collection
#' 
#' The original dataset was collected manually from the electronic patient 
#' record by one of the co-authors (Alfons Wachter). All information was 
#' stored in an SPSS dataset. The original dataset contained a large number
#' of variables, many of which are not used in this project. Column names are 
#' for the most part in German.

#' ## Load and inspect the data

#+ load-data
original_data <- glue("{.dir_raw}/original_data.sav") %>% 
  read_spss() %>% 
  as_factor()

nrow(original_data)  # 761
names(original_data)

#+ glimpse-scrambled-data
set.seed(1)
original_data %>% 
  mutate(across(everything(), sample, replace = FALSE)) %>% 
  glimpse()


#' ## Select relevant variables and rows, recode and rename
#' 
#' For the analysis presented in this study, we require the following 
#' variables:
#' * patient_id: a unique, anonymised patient identifier
#' * cycle: the insemination cycle for which data was collected
#' * cycle_date: the date when insemination was performed
#' * diagnosis: the infertility diagnosis of the female patient (endometriosis, 
#'              uterine_malformation, tubal_obstruction, anovulatory, or 
#'              no_known_female_inf)
#' * female_male: additional information on whether an infertility diagnosis
#'                was made for the women only, for her partner only, or for 
#'                both. 
#' * drug: the drug used to stimulate ovulation (gonadotropin, clomifen, both, 
#'         or none)
#' * age: patient age at insemination in years
#' * bmi: patient BMI (= weight / height in metres squared)
#' * amh: levels of anti-mueller hormone
#' * sperm: number of motile sperm inseminated
#' * hcg: was a positive hcg (=pregnancy) measured; outcome

#+ select-variables
data <- original_data %>% 
  select(
    patient_id = Patnr,
    cycle = ausgewerteteIUINr,
    cycle_date = date,
    diagnosis = Diagnose, 
    female_male = Diagnosenverteilung,
    drug = Stimulation,
    age = age,
    bmi = BMI,
    amh = AMH, 
    sperm = bewegliche_inseminierte_Sperm,
    follicles = Leitfollikel,
    hcg = resultHCG
  ) 

#' From this dataset, all patients with a missing or "other" diagnosis were 
#' removed (N = 2), as were patients with no information on age (N = 1). 
#' The data was then inspected for outlying values and cleaned
#' accordingly after discussing those cases with the clinicians on the team.

#+ exclude-based-on-diagnosis
data <- data %>% 
  filter(
    is.na(diagnosis) | diagnosis != "andere",
    !is.na(age)
  )
nrow(data) # 758

#+ crude-summary
summary(data)

#+ recode-variables
data <- data %>% 
  mutate(
    # turn outcome into binary 0/1
    hcg = as.integer(hcg == "positiv"),
    
    # rename and group infertility diagnoses as described above
    diagnosis = case_when(
        diagnosis == "Endometriose" ~ "endometriosis",
        diagnosis == "Uterusfehlbildung" ~ "uterine_malformation",
        diagnosis == "Tubenfx" ~ "tubal_obstruction",
        diagnosis == "idiopathische Sterilität der Frau" ~ "no_known_female_inf",
        diagnosis %in% c("PCOS", "hypothalamisch-hypophysäre Ovarialinsuffizien")
        ~ "anovulatory"
      ) %>% 
      factor(c("no_known_female_inf", "endometriosis", "anovulatory", 
               "tubal_obstruction", "uterine_malformation")),
    
    # recode stimulation drug and turn into a factor
    drug = drug %>% 
      as.character() %>% 
      str_to_lower() %>% 
      factor(
        c("gonadotropin", "clomifen", "clomifen + gonadotropin", "nostim")
      ) %>% 
      fct_recode( 
        both = "clomifen + gonadotropin",
        none = "nostim"
      ),
    
    # clean AMH for three patients were erroneously coded as 16000. These were
    # checked and describe values of >16.
    amh = if_else(amh == 16000, 16, amh),
    
    # fix a single patient BMI that was recorded 2.69 (typo from 1.69)
    bmi = if_else(patient_id == 581, 55 / (1.69) ^ 2, bmi),
    
    # set sperm count to missing where sperm = 0 (N = 3)
    sperm = if_else(sperm == 0, NA_real_, sperm),
    
    # set follicle count to missing where follicle = 0 (N = 4)
    follicles = if_else(follicles == 0, NA_real_, follicles),
    
    # assume idiopathic infertility for the rows were no male infertility was
    # reported
    female_male = case_when(
      !is.na(female_male) ~ female_male,
      diagnosis == "no_known_female_inf" ~ factor("none", levels(female_male)),
      TRUE ~ factor("female", levels(female_male))
    ),
    
    # fix a single anovulatory patient mistakenly classified as male inf only
    female_male = case_when(
      diagnosis != "no_known_female_inf" & female_male == "male" ~ 
        factor("female and male", levels(female_male)),
      TRUE ~ female_male
    )
  )

#+ crude-summary-after-recode
summary(data)



#' Finally, we save the cleaned dataset.

#+ save-data
write_rds(data, glue("{.dir_der}/data.rds"))






