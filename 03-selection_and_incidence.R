library(tidyverse)
library(haven)
library(here)        
library(lubridate)
library(broom)
library(patchwork)

# Specify folder of data
data_folder <- here("data", "2410stata")

# Load the data
pat <- read_dta(here(data_folder, "pat.dta"))
tailor <- read_dta(here(data_folder, "tailor.dta"))
std <- read_dta(here(data_folder, "std.dta"))
fup <- read_dta(here(data_folder, "fup.dta"))


# Look for MSM with active follow-up and occasional partners
last_fup <- fup |> 
  arrange(id, desc(fupdate)) |> 
  group_by(id) |> 
  slice(1) |> 
  select(id, last_fup = fupdate, last_center = center,
         last_source = source, last_physician = physician)

msm <- tailor |> 
  filter(riskgroup == "MSM")


msm_fup <- msm |> 
  left_join(last_fup, by = "id") |> 
  filter(last_fup >= dmy("01.01.2023"))

msm_fup_occas <- fup |> 
  filter(id %in% msm_fup$id) |> 
  arrange(id, desc(fupdate)) |> 
  group_by(id) |> 
  slice_head(n = 8) |> 
  summarise(any_occas = any(p_occas %in% c(1, 9)))

count_data <- msm_fup_occas |>
  left_join(last_fup |> select(id, last_source, last_center, last_physician), by = "id") 


count_data |>
  filter(any_occas == TRUE) |> 
  count(last_source) |> 
  mutate(p = n / sum(n)*100) |> 
  knitr::kable(digit = 1)

count_data |>
  filter(any_occas == TRUE, 
         last_source == 2) |> 
  count(last_physician, last_center, name = "n") |> 
  mutate(total = sum(n), 
         p = n / total * 100) |> 
  mutate(CPZH = last_physician %in% c("Hampel", "Wissel", "Ruff", "Yerguz", 
                                      "Komaromi", "Kalubi")) |> 
  arrange(desc(CPZH), last_center, desc(p)) |> 
  mutate(p_cum = cumsum(p), 
         n_cum = cumsum(n)) |> 
  knitr::kable(digit = 1)



count_data |>
  filter(any_occas == TRUE, 
         last_source == 3) |> 
  count(last_source, last_center, name = "n") |> 
  mutate(total = sum(n), 
         p = n / total * 100, 
         last_center = case_when(last_center == 10 ~ "ZH", 
                                 last_center == 20 ~ "BS", 
                                 last_center == 30 ~ "BE", 
                                 last_center == 40 ~ "GE", 
                                 last_center == 50 ~ "LAU", 
                                 last_center == 60 ~ "SG", 
                                 last_center == 70 ~ "LUG")) |> 
  knitr::kable(digit = 1)


count_data |>
  filter(any_occas == TRUE, 
         last_source == 3, last_center == 10) |> 
  count(last_source, last_center, last_physician, name = "n") |> 
  mutate(total = sum(n), 
         p = n / total * 100) |> 
  arrange(desc(p)) |> 
  mutate(p_cum = cumsum(p), 
         n_cum = cumsum(n)) |>
  knitr::kable(digit = 1)


CPZH <- c("Hampel", "Wissel", "Ruff", "Yerguz", "Komaromi", "Kalubi")
GPs <- c("Aceto", "Depmeier", "Kovari", "Klingler")

# Select MSM with active follow-up and occasional partners
# and who are followed in main cohort center, or in CPZH, or in large GPs in ZH

elig_msm <- count_data |> 
  filter(any_occas == TRUE) |> 
  filter(last_source == 1 |
           (last_source == 2 & last_physician %in% CPZH) |
           (last_source == 3 & last_center == 10 & last_physician %in% GPs)) 


nrow(elig_msm) / nrow(msm_fup_occas |> filter(any_occas == TRUE)) * 100



# Calculate incidence of STI in those individuals September 2017

temp <- fup |> 
  filter(id %in% elig_msm$id) |> 
  filter(fupdate >= dmy("01.09.2017")) |> 
  select(id, fupdate) |> 
  arrange(id, fupdate)

std <- std |> 
  arrange(id, std_date)

# Join temp2 on temp. Note that we have an std_date in temp2 which indicates
# the time point at which an std occurred. We also have fupdate in temp. We 
# want to join temp2 on temp where std_date is between 2 fupdates.


# Calculate fup since September 2017
fup_time <- temp |> 
  group_by(id) |> 
  summarise(fup_duration = as.numeric(last(fupdate) - first(fupdate))/365.25) |> 
  filter(fup_duration > 0)


# Calculate number of CT and NGO infections since September 2017
ctgo <- std |> 
  filter(std_date >= dmy("01.09.2017")) |> 
  group_by(id) |> 
  summarise(n_any_ctgo = sum(type %in% c(1, 3)), 
            n_symptomatic_ctgo = sum(type %in% c(1, 3) & std_symptoms == 1))

# Join datasets and replace NA with 0
incidence_df <- fup_time |> 
  left_join(ctgo, by = "id") |> 
  mutate(across(n_any_ctgo:n_symptomatic_ctgo, ~replace_na(., 0)))


# Calculate incidence using a poisson regression with fup_duration as offset


# --- Model for any CT or NGO infection
   m_any <- glm(n_any_ctgo ~ 1 + offset(log(fup_duration)), 
                family = poisson, data = incidence_df) 
   
   tidy(m_any, exp = TRUE, conf.int = TRUE) |> 
     knitr::kable(digits = 3)

# |term        | estimate| std.error| statistic| p.value| conf.low| conf.high|
# |:-----------|--------:|---------:|---------:|-------:|--------:|---------:|
# |(Intercept) |    0.197|      0.02|   -79.754|       0|    0.189|     0.205|

    
    
# --- Model for symptomatic CT or NGO infection
   m_sym <- glm(n_symptomatic_ctgo ~ 1 + offset(log(fup_duration)), 
                family = poisson, data = incidence_df)
   
   tidy(m_sym, exp = TRUE, conf.int = TRUE) |>
     knitr::kable(digits = 3)

# |term        | estimate| std.error| statistic| p.value| conf.low| conf.high|
# |:-----------|--------:|---------:|---------:|-------:|--------:|---------:|
# |(Intercept) |    0.083|     0.031|   -79.324|       0|    0.078|     0.089|
#    
    
    
    
    
# Now look at the screening in those individuals (variable introduced in April 2020)

screening_uptake <- fup |> 
  arrange(id, fupdate) |> 
  filter(id %in% elig_msm$id) |> 
  filter(fupdate >= dmy("01.04.2020")) |> 
  select(id, fupdate, nd_std_screening) |> 
  # cacluate the number of visits per ID
  group_by(id) |> 
  summarise(n_visits = n(),
            n_screens = sum(nd_std_screening == 1, na.rm = TRUE),
            fup_duration = as.numeric(last(fupdate) - first(fupdate))/365.25) |> 
  mutate(p_screened = n_screens / n_visits * 100)


# Summary of Screening uptake
median(screening_uptake$p_screened) # Median 22.2%
mean(screening_uptake$p_screened) # Mean 27.9%


screening_uptake |> 
  left_join(last_fup) |> 
  ggplot() + 
  gghalves::geom_half_violin(aes(y = p_screened, x = fct_rev(factor(last_center))), 
                             side = "r", 
                             fill = "grey90", 
                             width = 1.3) + 
  gghalves::geom_half_boxplot(aes(y = p_screened, x = fct_rev(factor(last_center))), 
                              side = "l", width = 0.2, 
                              fill = "grey50", ) +
  labs(x = "Center
  ", y = "
       Percentage of visits with screening") + 
  theme_minimal(base_family = "Signika Negative") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  coord_flip()

ggsave(here("figures", "03-uptake_by_center_violin.png"), 
       dpi = 300, width = 10, height = 7, bg = "white")

screening_uptake |> 
  left_join(last_fup) |> 
  ggplot(aes(x = p_screened, after_stat(count))) +
  geom_density(fill = "grey", color = "black") + 
  scale_x_continuous(breaks = seq(0, 100, 10)) + 
  labs(x = "
       Percentage of visits with screening", y = "Number of individuals
       ") + 
  theme_minimal(base_family = "Signika Negative") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggsave(here("figures", "03-uptake_overall.png"), 
       dpi = 300, width = 10, height = 7, bg = "white")


screening_uptake |> 
  left_join(last_fup) |> 
  mutate(last_source = paste("Last source =", last_source), 
         last_center = case_when(last_center == 10 ~ "ZH", 
                                 last_center == 20 ~ "BS", 
                                 last_center == 30 ~ "BE", 
                                 last_center == 40 ~ "GE", 
                                 last_center == 50 ~ "LAU", 
                                 last_center == 60 ~ "SG", 
                                 last_center == 70 ~ "LUG")) |> 
  ggplot(aes(x = p_screened)) +
  geom_density(fill = "grey", color = "black", adjust = 0.66) + 
  facet_grid(last_center ~ last_source) + 
  scale_x_continuous(breaks = seq(0, 100, 25)) + 
  labs(x = "
  Percentage of visits with screening", y = "Kernel density
       ") + 
  theme_minimal(base_family = "Signika Negative") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        legend.position = "none", 
        # rotate facet text
        strip.text.y = element_text(angle = 0, size = 10), 
        strip.text.x.top = element_text(face = "bold", size = 10))

ggsave(here("figures", "03-uptake_by_center.png"),
       dpi = 300, width = 10, height = 7, bg = "white")

