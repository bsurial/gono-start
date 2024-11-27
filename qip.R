my_power <- function(incidence, effect, fup, power) {
  p <- power_Poisson(lambda1 = incidence, 
                     lambda2 = incidence * effect, 
                     t1 = fup, t2 = fup, equal.sample = TRUE, 
                     power = 0.8,
                     sig.level = 0.025, alternative = "one.sided")
  tibble(power = p$power, 
         N_control = ceiling(p$N),
         N_intervention = ceiling(p$N), 
         effect = p$lambda2/p$lambda1)
}

my_power(0.083, 0.7, 2.5, 0.8)

# Power for a range of effect sizes
map_dfr(seq(0.5, 0.9, 0.05), ~my_power(0.083, .x, 2.5, 0.8))


power_Poisson(lambda1 = 0.083, 
              lambda2 = 0.083 * 0.9, 
              n1 = 587,
              t1 = 2.5, t2 = 2.5, equal.sample = FALSE, 
              power = 0.8,
              sig.level = 0.025, alternative = "one.sided")

rbinom(100, 1, 2/3) |> 
  sum()



p_cont <- 0.083 # Estimated cessation proportion in control group (see details in study protocol)
p_int_uptake <- 0.5 * p_cont # Estimated cessation proportion in intervention group, among uptakers (see details in study protocol)
p_int_non_uptake <- 0.083 # Estimated cessation proportion in intervention group, among non-uptakers (see details in study protocol)
alpha <- 0.05 # Significance level
power <- 0.80 # Desired power, keep it fixed at min. 80%
attrition <- 0.03 # attrition/LTFU rate across both arms (see updated SHCS cohort profile publication Scherrer et al.)


uptake <- 0.5 # 70% uptake
p_int <- (p_int_uptake*uptake) + (p_int_non_uptake*(1-uptake)) # Estimated proportion in intervention group, with non-uptake integrated
effect_delta <- p_int-p_cont

# print effect size
cat("effect size delta with non-uptake incorporated in intervention:", effect_delta)




my_power <- function(incidence, effect, uptake, fup, power) {
  # Adjust incidence 2 for uptake
  
  p_cont <- incidence
  p_int_uptake <- effect * p_cont
  p_int_non_uptake <- incidence
  p_int <- (p_int_uptake*uptake) + (p_int_non_uptake*(1-uptake))
  
  p <- power_Poisson(lambda1 = p_cont, 
                     lambda2 = p_int, 
                     t1 = fup, t2 = fup, equal.sample = TRUE, 
                     power = 0.8,
                     sig.level = 0.025, alternative = "one.sided")
  tibble(power = p$power, 
         N_control = ceiling(p$N),
         N_intervention = ceiling(p$N), 
         uptake = uptake,
         effect = p_int_uptake / p_cont, 
         effect_adj = p_int / p_cont)
}

my_power(incidence = 0.083, 
         effect = 0.5, 
         uptake = 0.5, 
         fup = 2.5, 
         power = 0.8)


# Power for a range of effect sizes and uptakes, 2.5 years
map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 0.5, 
  fup = 2.5, 
  power = 0.8))

map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 0.55, 
  fup = 2.5, 
  power = 0.8))

map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 0.60, 
  fup = 2.5, 
  power = 0.8))

map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 0.65, 
  fup = 2.5, 
  power = 0.8))

map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 0.70, 
  fup = 2.5, 
  power = 0.8))


# Power for a range of effect sizes and uptakes, 3 years
map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 1, # no adjustment for uptake
  fup = 3, 
  power = 0.8))


map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 0.5, 
  fup = 3, 
  power = 0.8))

map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 0.55, 
  fup = 3, 
  power = 0.8))

map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 0.60, 
  fup = 3, 
  power = 0.8))

map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 0.65, 
  fup = 3, 
  power = 0.8))

map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 0.70, 
  fup = 3, 
  power = 0.8))



map_dfr(seq(0.5, 0.7, 0.05), ~my_power(
  incidence = 0.083, 
  effect = .x, 
  uptake = 0.65, 
  fup = 2.5, 
  power = 0.8))


# Power for a range of uptakes
map_dfr(seq(0.5, 0.9, 0.05), ~my_power(
  incidence = 0.083, 
  effect = 0.5, 
  uptake = .x, 
  fup = 2.5, 
  power = 0.8)) |> 
  ggplot(aes(x = uptake, y = N_control)) + 
  geom_line()
