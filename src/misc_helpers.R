

add_cycle_date <- function(data_long) {
  data_long %>%
    group_by(patient_id) %>% 
    mutate(
      cycle_date = cycle_date %m-% days(28 * (cycle - timeInt - 1))
    ) %>%
    ungroup()
}


tabulate_preg <- function(var) {
  var <- rlang::enexpr(var)
  data_long %>% 
    group_by(!!var, y) %>% 
    summarise(n = n(), .groups = "drop_last") %>% 
    mutate(p = n / sum(n)) %>% 
    ungroup() %>% 
    pivot_wider(
      names_from = "y",
      values_from = c("n", "p"),
      values_fill = 0
    )
} 


fit_preg <- function(form) {
  glm(update.formula(y ~ 1, form), data = data_long, family = binomial)
}

predict_preg <- function(fit, newdata) {
  predict(fit, newdata, type = "response", se = TRUE) %>% 
    as_tibble()
}

plot_preg <- ggplot(data = NULL, aes(y = y)) + 
  geom_point(
    position = position_jitter(width = 0.3, height = 0.01), 
    alpha = 0.5
  ) + 
  geom_smooth(method = "lm", colour = "green") + 
  coord_cartesian(ylim = c(-0.05, 1.05)) + 
  theme_bw()

binom_exp <- function(n, p, level = 0.95) {
  map_dbl(
    .x = c(low = 1 - level, high = 1 + level) / 2, 
    .f = qbinom, 
    size = n, 
    prob = p
  )
}

calc_p_and_ci <- . %>% 
  mutate(
    total_percent = total / total[1],
    percent = pregn / total,
    cum_percent = cumsum(pregn / total[1]),
    map_dfr(total, binom_exp, p = sum(pregn) / sum(total))
  ) %>% 
  mutate(across(contains("percent"), scales::percent, accuracy = 0.1))