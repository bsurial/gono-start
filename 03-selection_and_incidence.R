library(tidyverse)
library(haven)
library(here)        
library(lubridate)

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


(1425 + 201 + 575) / nrow(msm_fup_occas |> filter(any_occas == TRUE)) * 100








