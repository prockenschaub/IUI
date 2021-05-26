#*******************************************************************************
#
# Project: Medical University Innsbruck - Intrauterine insemination (IUI)
# Date:    2020-04-15
# Author:  Patrick Rockenschaub
# Purpose: Derive a clinical score
#
#******************************************************************************* 

#' ---
#' title: "Derive a clinical score"
#' subtitle: "Medical University Innsbruck - Intrauterine insemination (IUI)"
#' output: pdf_document
#' ---

#+ setup, message = FALSE 
source("00_init_workspace.R")
source(glue("{.dir_src}/modelling_helpers.R"))

library(mice)

imp <- read_rds(glue("{.dir_der}/imp.rds"))
fits <- read_rds(glue("{.dir_res}/fits.rds"))


#' Calculate the curve of cumulative probability of pregnancy
#' 
#' @param score integer score between 0 and 5
#' @param intercept model intercept on the scale of the linear predictor
#' @param slope slope coefficient on the scale of the linear predictor
#' @return vector
calculate_curve <- function(score, intercept, slope) {
  logistic <- function(x) 1 / (1 + exp(-x))
  
  p <- logistic(intercept + (score) * slope)
  
  cum_p <- vector("numeric", 6)
  cum_p[1] <- p
  for (i in 2: length(cum_p)) { 
    cum_p[i] <- cum_p[i - 1] + (1 - cum_p[i - 1]) * p
  }
  cum_p
}


#+ calculate-score
slope <- 0.5  # Manually chosen common denominator
intercept <- get_coefs(fits[[2]])[1, "estimate"] - slope * 5 # make 5 highest

score <- tibble(
    score = 0:5
  ) %>% 
  mutate(
    p = map(score, calculate_curve, intercept, slope),
    score = factor(score, 5:0)
  )

#+ plot-score
g_score <- score %>% 
    unnest(p) %>% 
    group_by(score) %>% 
    mutate(cycle = 1:6) %>% 
    ungroup() %>% 
    ggplot(aes(cycle, p, group = factor(score), colour = factor(score))) + 
    geom_line(size = 1.2) +
    scale_x_continuous(breaks = 1:6) + 
    scale_y_continuous(labels = scales::percent, breaks = 0:5 / 5) + 
    scale_colour_viridis_d(end = 0.8, direction = -1) + 
    coord_cartesian(ylim = c(0, 1), expand = FALSE) + 
    labs(
      x = "\nStimulation cycles",
      y = "Cumulative probability of pregnancy\n", 
      colour = "Score"
    ) +
    theme_bw() + 
    theme(
      panel.grid.minor = element_blank()
    )

g_score
ggsave("score_curves.png", width = 14, height = 14, dpi = 600, units = "cm")


#+ exact-probs
score_in_cohort <- complete(imp, "all") %>% 
  reduce(bind_rows) %>% 
  mutate(
    score = factor(
      0 + 2 * (amh >= 1) + 
        1 * (sperm >= 5) + 1 * (sperm >= 15) + 
        1 * map_lgl(diagnosis, ~ . %in% c("anovulatory", "no_known_female_inf")),
      5:0
    )
  )
  
counts <- score_in_cohort %>% 
  group_by(score) %>% 
  summarise(
    n = n(),
    cycles = mean(cycle)
  ) %>% 
  mutate(
    n = n / 10,
    p = scales::percent(n / sum(n))
  )

probs <- score %>% 
  mutate(
    map_dfr(
      p, 
      ~ as_tibble(matrix(., nrow = 1, dimnames = list(0, str_c("c", 1:6))))
  )) %>% 
  select(-p)

inner_join(counts, probs, by = "score")


#+ calculate-missed-preg
score_in_cohort %>% 
  select(patient_id, cycle, hcg, score) %>% 
  filter(
    cycle < 4,
    hcg == 0
  ) %>% 
  inner_join(
    probs %>% select(score, c1:c3), 
    by = "score"
  ) %>% 
  mutate(
    missed_p = case_when(
      cycle == 1 ~ c3, 
      cycle == 2 ~ c2,
      cycle == 3 ~ c1
    )
  ) %>% 
  summarise(exp_missed = sum(missed_p) / 10)
