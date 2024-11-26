library(future.apply)
library(PASSED)
library(tidyverse)
library(broom)
library(patchwork)
plan(multisession, workers = 11)


n_sim = 1000          # Number of simulations for each power estimate
fup_time = 2.5        # Mean study duration per participant in years
incidence = 5.41/100   # Incidence rate per 100 person-years
effect = 0.5          # Effect size, IRR
beta = 0.1            # Beta to stop the simulation


# To save some time, let's start with a parametric estimation to know
# where to start

temp <- power_Poisson(lambda1 = incidence, lambda2 = effect * incidence, 
                      t1 = fup_time, t2 = fup_time,
                      power = 0.80, equal.sample = TRUE, 
                      sig.level = 0.025, alternative = "one.sided")


sample_pois <- function(n_per_arm, incidence, fup, effect) {
  df <- tibble(id = 1:(2*n_per_arm), 
               incidence = incidence, 
               grp = rbinom(2*n_per_arm, 1, 0.5),
               fup = fup) |> 
    mutate(incidence = if_else(grp == 1, effect * incidence, incidence)) |> 
    mutate(count = rpois(2*n_per_arm, incidence * fup))
  
  m <- glm(count ~ grp + offset(log(fup)), family = poisson, data = df)
  
  coefs <- summary(m)$coefficients

  exp(coefs[2,1] + 1.96 * coefs[2,2]) < 1.0

}


# Calculate some power estimates, stop at 90% power (1 - beta)
n_per_arm <- ceiling(temp$N)-100    # Number of participants per arm
i <- 1                              # Index for dataframe
pwr <- c()                          # Vector to store power estimates
ns <- c()                           # Vector to store number of participants per arm


power <- mean(future_replicate(n_sim, sample_pois(n_per_arm, incidence, fup_time, effect)))
while (power < (1 - beta)) {
  n_per_arm <- n_per_arm + 10
  power <- mean(future_replicate(n_sim, sample_pois(n_per_arm, incidence, fup_time, effect)))
  list(n_per_arm, power)
  cat("Power with", n_per_arm * 2, "participants:", power, "\n")
  pwr[i] <- c(power)
  ns[i] <- c(n_per_arm)
  i <- i + 1
}

samplesize_calculated <- tibble(pwr, ns)

# Use a loess model to get estimate of the smoother that is closest to 80%
mod <- loess(pwr ~ ns, data = samplesize_calculated, span = 0.75)
new_ns <- tibble(ns = (first(ns)):last(ns))
new_ns$p <- predict(mod, newdata = new_ns, type = "response")

p0.8 <- new_ns |> 
  filter(p>=0.8) |> 
  slice(1)

# Generate graph
p_sim <- tibble(pwr, ns) |> 
  ggplot(aes(x = ns, y = pwr)) + 
  geom_line(linewidth = 0.25) + 
  geom_smooth(color = "darkred", span = 0.75) + 
  geom_point(aes(x = ns, y = p), color = "darkred", data = p0.8, size = 3) +
  ggrepel::geom_text_repel(aes(label = glue::glue("N per arm: {ns}"), x = ns, y = p), data = p0.8, 
                           nudge_x = -25, nudge_y = 0.05, segment.curvature = 0.1, 
                           point.padding = 1,
                           box.padding = 5,
                           arrow = arrow(length = unit(0.015, "npc")), 
                           family = "Roboto", fontface = "bold") +
  geom_hline(yintercept = 0.80, linetype = "dashed") + 
  # scale_x_continuous(breaks = seq(200, 500, 50)) +
  labs(title = "Simulation based", 
  x = "
       Number of participants per arm
       ", y = "Power") + 
  theme_minimal(base_family = "Roboto") + 
  theme(plot.caption = element_text(hjust = 0.5, 
                                    face = "italic"), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank())




my_power <- function(n, incidence, effect, fup) {
  p <- power_Poisson(n1 = n, lambda1 = incidence, 
                     lambda2 = incidence * effect, 
                     t1 = fup, t2 = fup, equal.sample = TRUE, 
                     sig.level = 0.025, alternative = "one.sided")
  tibble(power = p$power, 
         n_per_arm = p$N)
}


p_par <- map_dfr(seq(temp$N-100,temp$N+100, 1), ~my_power(n = .x, incidence = incidence, 
                                            effect = effect, fup = fup_time)) |> 
  ggplot(aes(n_per_arm, power)) +
  geom_line(color = "darkred") +
  geom_point(aes(x = temp$N, y = 0.80), color = "darkred") + 
  geom_label(aes(x = temp$N, y = 0.80, label = glue::glue("N per arm: {ceiling(temp$N)}")), nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, linetype = 2) +
  labs(title = "Parametric", 
       x = "Number of patients per arm", y = "Power") +
  theme_minimal(base_family = "Roboto") + 
  theme(plot.caption = element_text(hjust = 0.5, 
                                    face = "italic"), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank())


p_sim + p_par &
  patchwork::plot_annotation(caption =glue::glue("
                           Incidence: {incidence*100} per 100 PY, Effect size: IRR {effect}, Trial duration: {fup_time} years, One-sided alpha: 0.025
                           Each estimate is based on {n_sim} simulations. Model used: glm(count ~ grp + offset(log(fup)), data = dd, family = ´poisson´)"))

ggsave("02-samplesize_superiority_5.41.png", width = 10, height = 5, bg = "white", dpi = 300)
