###############################################################################
# DATA ANALYSIS AND VISUALIZATION for Gačnik et al. (2025)
# CITATION: Gačnik J., Žagar K., Hatvani I. G., Kern Z., and Vreča P. (2025): Climate change reflected in 40-year isotopic composition trends of precipitation in Slovenia, Environmental Research
# CODE AUTHOR: JAN GAČNIK, MAY 2025
# R version: 4.4.2
###############################################################################
# DO THIS STEP ONLY THE FIRST TIME! Installation of needed packages
install.packages(c("mblm", "ggspatial", "sf", "rnaturalearth", "gsignal", "Rlibeemd", "ggpmisc", "lubridate", "slider", "scales", "trend", "segmented", "tidyverse", "readxl"))

# CHANGE THE FOLLOWING PARAMETERS:
path_scripts <- "ENTER/PATH/TO/DIRECTORY/Precipitation_isotopes_data_analysis/scripts" # Change to your directory path where the R scripts are located, use "/" and not "\"
path_data <- "ENTER/PATH/TO/DIRECTORY/Precipitation_isotopes_data_analysis/data" # Change to your directory path where the data is located, use "/" and not "\"
path_figures <- "ENTER/PATH/TO/DIRECTORY/Precipitation_isotopes_data_analysis/figures" # Change to your directory path where you want to save figures, use "/" and not "\"

###############################################################################
# CODE BELOW DOES NOT NEED ALTERING
###############################################################################
Ljubljana_filename <- "Ljubljana_1981-2024.xlsx"
Ljubljana_filename_metdata <- "Ljubljana-daily-meteodata.xlsx"
Vienna_filename <- "Vienna_1973-2024.xlsx"
Vienna_filename_metdata <- "Vienna-daily-meteodata.xlsx"
Vienna_filename_precipitation <- "Vienna_new_precip.xlsx"
Portoroz_filename <- "Portoroz_2000-2024.xlsx"
Portoroz_filename_metdata <- "Portoroz-daily-meteodata.xlsx"
Climate_filename <- "Climate_normals.xlsx"
HYSPLIT_filename <- "HYSPLIT_timeseries.xlsx"

MOI_filename <- "MOI.txt"
NAO_filename <- "NAO.csv"
figure_type <- "png"
dpi_set <- 500

# Load required packages for the current session
setwd(path_scripts)
lapply(c("mblm", "ggspatial", "sf", "rnaturalearth", "gsignal", "Rlibeemd", "ggpmisc", "lubridate", "slider", "scales", "trend", "segmented", "tidyverse", "readxl"), require, character.only = TRUE)

# Calling needed functions
source("Functions_MA&RMA.R")
source("Functions_Hilbert_Huang_transform.R")

###############################################################################
# Data reading, adjusting, and filtering
###############################################################################
Ljubljana_data <- read_xlsx(sheet = "Data Ljubljana", path = file.path(path_data, Ljubljana_filename), trim_ws = TRUE) %>%
  dplyr::select(Sample_ID, Station_ID, Name, Year, Month, P, T, RH, δ18O, δ2H, d, δ18On, δ2Hn, dn) %>%
  mutate(Site = "Ljubljana") %>%
  mutate(Period = "") %>%
  mutate(Period = case_when(
    Year >= 1981 & Year < 1992  ~ "1981\u20131991",
    Year >= 1992 & Year < 2002  ~ "1992\u20132002",
    Year >= 2002 & Year < 2013  ~ "2003\u20132013",
    Year >= 2014 ~ "2014\u20132024",
    TRUE ~ Period)) %>%
  mutate(Season = "") %>%
  mutate(Season = case_when(
    Month %in% c(12, 1, 2)  ~ "winter",
    Month %in% c(3, 4, 5)  ~ "spring",
    Month %in% c(6, 7, 8)  ~ "summer",
    Month %in% c(9, 10, 11)  ~ "autumn",
    TRUE ~ Season)) %>%
  mutate(Season = factor(Season, levels = c("winter", "spring", "summer", "autumn"))) %>%
  mutate(P_δ18O_δ2H = P) %>%
  mutate(P_δ18O_δ2H = case_when(
    is.na(δ18O) | is.na(δ2H) ~ NA,
    TRUE ~ P_δ18O_δ2H)) %>%
  mutate(rolling_δ18O = slide_dbl(δ18O, mean, na.rm = TRUE, .before = 120),
         rolling_δ2H = slide_dbl(δ2H, mean, na.rm = TRUE, .before = 120),
         rolling_d = slide_dbl(d, mean, na.rm = TRUE, .before = 120),
         rolling_T = slide_dbl(T, mean, na.rm = TRUE, .before = 120),
         rolling_P = slide_dbl(P, mean, na.rm = TRUE, .before = 120)) %>%
  mutate(weighted_rolling_δ18O = slide2_dbl(.x = δ18O, .y = P_δ18O_δ2H, .before = 120, .f = ~ sum(.x * ((.y/sum(.y, na.rm = TRUE))), na.rm = TRUE)),
         weighted_rolling_δ2H = slide2_dbl(.x = δ2H, .y = P_δ18O_δ2H, .before = 120, .f = ~ sum(.x * ((.y/sum(.y, na.rm = TRUE))), na.rm = TRUE)),
         weighted_rolling_d = slide2_dbl(.x = d, .y = P_δ18O_δ2H, .before = 120, .f = ~ sum(.x * ((.y/sum(.y, na.rm = TRUE))), na.rm = TRUE))) %>%
  mutate("count δ18O" = slide_dbl(δ18O, .before = 120, .f = ~ sum(!is.na(.))),
         "count δ2H" = slide_dbl(δ2H, .before = 120, .f = ~ sum(!is.na(.))),
         "count d" = slide_dbl(d, .before = 120, .f = ~ sum(!is.na(.)))) %>%
  mutate(Date = as.POSIXct(make_date(Year, Month)),
         time_since_start = interval(min(Date), Date) %/% months(1),
         Date = as.Date(make_date(Year, Month)))

Ljubljana_data_filtered <- dplyr::filter(Ljubljana_data, is.na(δ2H) == FALSE & is.na(δ18O) == FALSE & is.na(P) == FALSE) %>%
  mutate(time_since_start = interval(min(Date), Date) %/% months(1))

Ljubljana_data_longer <-  Ljubljana_data %>% 
  pivot_longer(cols = c(δ18O, δ2H, d, T, P), names_to = "Variable", values_to = "Value") %>%
  pivot_longer(cols = c(weighted_rolling_δ18O, weighted_rolling_δ2H, weighted_rolling_d), names_to = "Weighted_rolling_variable", values_to = "Weighted_rolling_value") %>%
  pivot_longer(cols = c(rolling_δ18O, rolling_δ2H, rolling_d, rolling_T, rolling_P), names_to = "Rolling_variable", values_to = "Rolling_value") %>%
  mutate(Variable_color = Variable) %>%
  mutate(Variable = case_when(
    Variable == "δ18O" ~ "italic(delta)^18*O",
    Variable == "δ2H" ~ "italic(delta)^2*H",
    Variable == "d" ~ "d-excess",
    TRUE ~ Variable
  )) %>%
  mutate(Variable = factor(Variable, levels = c("italic(delta)^18*O", "italic(delta)^2*H", "d-excess", "T", "P")),
         Rolling_variable = factor(Rolling_variable, levels = c("rolling_δ18O", "rolling_δ2H", "rolling_d", "rolling_T", "rolling_P")),
         Weighted_rolling_variable = factor(Weighted_rolling_variable, levels = c("weighted_rolling_δ18O", "weighted_rolling_δ2H", "weighted_rolling_d")))

Ljubljana_data_longer_month <- Ljubljana_data %>%
  group_by(Month) %>%
  summarise(mean_T = mean(T, na.rm = TRUE),
            mean_P = mean(P, na.rm = TRUE),
            weighted_mean_δ18O = sum(δ18O * ((P_δ18O_δ2H/sum(P_δ18O_δ2H, na.rm = TRUE))), na.rm = TRUE),
            weighted_mean_δ2H = sum(δ2H * ((P_δ18O_δ2H/sum(P_δ18O_δ2H, na.rm = TRUE))), na.rm = TRUE),
            weighted_mean_d = sum(d * ((P_δ18O_δ2H/sum(P_δ18O_δ2H, na.rm = TRUE))), na.rm = TRUE)) %>%
  pivot_longer(cols = c(weighted_mean_δ18O, weighted_mean_δ2H, weighted_mean_d, mean_T, mean_P), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = case_when(
    grepl("_δ18O", Variable) ~ "italic(delta)^18*O",
    grepl("_δ2H", Variable) ~ "italic(delta)^2*H",
    grepl("_d", Variable) ~ "d-excess",
    grepl("_T", Variable) ~ "T",
    grepl("_P", Variable) ~ "P",
    TRUE ~ Variable)) %>%
  mutate(Variable = factor(Variable, levels = c("italic(delta)^18*O", "italic(delta)^2*H", "d-excess", "T", "P")))

Ljubljana_data_longest <- Ljubljana_data %>%
  mutate(P_δ18O_δ2H = P) %>%
  mutate(P_δ18O_δ2H = case_when(
    is.na(δ18O) | is.na(δ2H) ~ NA,
    TRUE ~ P_δ18O_δ2H)) %>%
  mutate(weighted_rolling_δ18O = slide2_dbl(.x = δ18O, .y = P_δ18O_δ2H, .before = 120, .f = ~ sum(.x * ((.y/sum(.y, na.rm = TRUE))), na.rm = TRUE)),
         weighted_rolling_δ2H = slide2_dbl(.x = δ2H, .y = P_δ18O_δ2H, .before = 120, .f = ~ sum(.x * ((.y/sum(.y, na.rm = TRUE))), na.rm = TRUE)),
         weighted_rolling_d = slide2_dbl(.x = d, .y = P_δ18O_δ2H, .before = 120, .f = ~ sum(.x * ((.y/sum(.y, na.rm = TRUE))), na.rm = TRUE))) %>%
  mutate("count δ18O" = slide_dbl(δ18O, .before = 120, .f = ~ sum(!is.na(.))),
         "count δ2H" = slide_dbl(δ2H, .before = 120, .f = ~ sum(!is.na(.))),
         "count d" = slide_dbl(d, .before = 120, .f = ~ sum(!is.na(.)))) %>%
  pivot_longer(cols = c(δ18O, δ2H, d, `T`, `P`, rolling_δ18O, rolling_δ2H, rolling_d, rolling_T, rolling_P, weighted_rolling_δ18O, weighted_rolling_δ2H, weighted_rolling_d), names_to = "Variable", values_to = "Value") %>%
  mutate(Value = case_when(
    Year < 1991 & Site == "Ljubljana" & grepl("rolling_", Variable) ~ NA,
    TRUE ~ Value
  )) %>%
  mutate(Rolling = "") %>%
  mutate(Rolling = ifelse(grepl("weighted", Variable), "Precipitation-weighted", "Unweighted")) %>%
  mutate(Rolling = factor(Rolling, c("Unweighted", "Precipitation-weighted"))) %>%
  mutate(Variable = gsub("weighted_", "", Variable)) %>%
  mutate(Variable_color = Variable) %>%
  mutate(Variable_color = gsub("rolling_", "", Variable_color)) %>%
  mutate(Variable = gsub("rolling_", "10-year rolling ", Variable)) %>%
  mutate(Variable = case_when(
    Variable == "δ18O" ~ "italic(delta)^18*O",
    Variable == "10-year rolling δ18O" ~ "\"10-year rolling\"~italic(delta)^18*O",
    Variable == "δ2H" ~ "italic(delta)^2*H",
    Variable == "10-year rolling δ2H" ~ "\"10-year rolling\"~italic(delta)^2*H",
    Variable == "d" ~ '"d-excess"',
    Variable == "10-year rolling d" ~ '"10-year rolling d-excess"',
    Variable == "T" ~ '"T"',
    Variable == "10-year rolling T" ~ '"10-year rolling T"',
    Variable == "P" ~ '"P"',
    Variable == "10-year rolling P" ~ '"10-year rolling P"',
    TRUE ~ Variable
  )) %>%
  mutate(Variable = factor(Variable, levels = c("italic(delta)^18*O", "\"10-year rolling\"~italic(delta)^18*O", "italic(delta)^2*H", "\"10-year rolling\"~italic(delta)^2*H", '"d-excess"', '"10-year rolling d-excess"', '"T"', '"10-year rolling T"', '"P"', '"10-year rolling P"')))

Ljubljana_data_all_supplement <- Ljubljana_data %>% 
  mutate(P_weighting = P) %>%
  pivot_longer(cols = c(δ18O, δ2H, d, T, P), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  reframe(Mean = mean(Value, na.rm = TRUE),
          "Weighted mean" = sum(Value * ((P_weighting /sum(P_weighting , na.rm = TRUE))), na.rm = TRUE),
          StDev = sd(Value, na.rm = TRUE),
          Median = median(Value, na.rm = TRUE),
          Q1 = quantile(Value, 0.25, na.rm = TRUE),
          Q3 = quantile(Value, 0.75, na.rm = TRUE),
          Min = min(Value, na.rm = TRUE),
          Max = max (Value, na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 1)))

#write.csv(Ljubljana_data_all_supplement, "Ljubljana_statistics_all.csv", row.names = FALSE)

Ljubljana_data_month_supplement <- Ljubljana_data %>% 
  mutate(P_weighting = P) %>%
  pivot_longer(cols = c(δ18O, δ2H, d, T, P), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable, Month) %>%
  reframe(Mean = mean(Value, na.rm = TRUE),
          "Weighted mean" = sum(Value * ((P_weighting /sum(P_weighting , na.rm = TRUE))), na.rm = TRUE),
          StDev = sd(Value, na.rm = TRUE),
          Median = median(Value, na.rm = TRUE),
          Q1 = quantile(Value, 0.25, na.rm = TRUE),
          Q3 = quantile(Value, 0.75, na.rm = TRUE),
          Min = min(Value, na.rm = TRUE),
          Max = max (Value, na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 1)))

#write.csv(Ljubljana_data_month_supplement, "Ljubljana_statistics_monthly.csv", row.names = FALSE)

Ljubljana_data_season_supplement <- Ljubljana_data %>% 
  mutate(P_weighting = P) %>%
  pivot_longer(cols = c(δ18O, δ2H, d, T, P), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable, Season) %>%
  reframe(Mean = mean(Value, na.rm = TRUE),
          "Weighted mean" = sum(Value * ((P_weighting /sum(P_weighting , na.rm = TRUE))), na.rm = TRUE),
          StDev = sd(Value, na.rm = TRUE),
          Median = median(Value, na.rm = TRUE),
          Q1 = quantile(Value, 0.25, na.rm = TRUE),
          Q3 = quantile(Value, 0.75, na.rm = TRUE),
          Min = min(Value, na.rm = TRUE),
          Max = max (Value, na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 1)))

#write.csv(Ljubljana_data_season_supplement, "Ljubljana_statistics_seasonal.csv", row.names = FALSE)

Ljubljana_data_year_slopes <- Ljubljana_data %>%
  dplyr::filter(!is.na(δ18O) & !is.na(δ2H)) %>%
  group_by(Year) %>%
  summarise(slope_MA = RMA_weighted(δ18O, δ2H, P)$slopes$MA,
            slope_MA_weighted = RMA_weighted(δ18O, δ2H, P)$slopes$MA_weighted,
            slope_RMA = RMA_weighted(δ18O, δ2H, P)$slopes$RMA,
            slope_RMA_weighted = RMA_weighted(δ18O, δ2H, P)$slopes$RMA_weighted,
            intercept_MA = RMA_weighted(δ18O, δ2H, P)$intercepts$MA,
            intercept_MA_weighted = RMA_weighted(δ18O, δ2H, P)$intercepts$MA_weighted,
            intercept_RMA = RMA_weighted(δ18O, δ2H, P)$intercepts$RMA,
            intercept_RMA_weighted = RMA_weighted(δ18O, δ2H, P)$intercepts$RMA_weighted) %>%
  mutate(across(everything(), ~round(., digits = 2)))

#write.csv(Ljubljana_data_year_slopes, "Ljubljana_statistics_yearly_regressions.csv", row.names = FALSE)

Ljubljana_data_slopes_longer <- Ljubljana_data_year_slopes %>% 
  pivot_longer(
    cols = starts_with("slope") | starts_with("intercept"),
    names_to = c(".value", "type"),
    names_pattern = "(slope|intercept)_(.*)"
  ) %>% 
  rename(Slope = slope, Intercept = intercept) %>% 
  pivot_longer(cols = c("Slope", "Intercept"), names_to = "statistic", values_to = "value") %>% 
  mutate(facet = type) %>% 
  mutate(facet = case_when(
    grepl("RMA", facet) ~ "RMA",
    grepl("RMA", facet) == FALSE  ~ "MA",
    TRUE ~ facet
  )) %>% 
  mutate(weighting = "") %>% 
  mutate(weighting = case_when(
    grepl("_weighted", type) ~ "Precipitation-weighted",
    grepl("_weighted", type) == FALSE ~ "Unweighted",
  ))
  
Ljubljana_slopes_Sen <- Ljubljana_data_slopes_longer %>%
  group_by(statistic, type) %>%
  group_modify(~ {
    model <- mblm(value ~ Year, data = .x, repeated = FALSE)
    tibble(
      Slope = model$coefficients[[2]],
      Intercept = model$coefficients[[1]],
    )
  })

Ljubljana_yearly_Sen_and_MK <- Ljubljana_data_longer %>%
  filter(!is.na(Value)) %>%
  group_by(Year, Variable) %>%
  summarise(Value = mean(Value), .groups = "drop") %>%
  group_by(Variable) %>%
  group_modify(~{
    df <- .x 
    model <- mblm(Value ~ Year, data = df, repeated = FALSE)
    mk_inc <- mk.test(df$Value, alternative = "greater")
    mk_dec <- mk.test(df$Value, alternative = "less")
    mk_two_sided <- mk.test(df$Value, alternative = "two.sided")
    tibble(Slope = unname(coef(model)[2]),
           Intercept = unname(coef(model)[1]),
           p_increasing = mk_inc$p.value,
           p_decreasing = mk_dec$p.value,
           p_two_sided = mk_two_sided$p.value)
  }) %>%
  ungroup()

#write.csv(Ljubljana_yearly_Sen_and_MK, file = "Ljubljana_yearly_Sen_and_MK.csv", row.names = FALSE)

Ljubljana_data_slopes_longer_2  <-  Ljubljana_data_slopes_longer %>% 
  left_join(Ljubljana_slopes_Sen, by = c("statistic", "type"))

# Portorož data
Portoroz_data <- read_xlsx(sheet = "Data Portorož", path = file.path(path_data, Portoroz_filename), trim_ws = TRUE) %>%
  dplyr::select(Sample_ID, Station_ID, Name, Year, Month, P, T, RH, δ18O, δ2H, d, δ18On, δ2Hn, dn) %>%
  mutate(Site = "Portorož") %>%
  mutate(Season = "") %>%
  mutate(Season = case_when(
    Month %in% c(12, 1, 2)  ~ "winter",
    Month %in% c(3, 4, 5)  ~ "spring",
    Month %in% c(6, 7, 8)  ~ "summer",
    Month %in% c(9, 10, 11)  ~ "autumn",
    TRUE ~ Season)) %>%
  mutate(Season = factor(Season, levels = c("winter", "spring", "summer", "autumn"))) %>%
  mutate(Date = as.POSIXct(make_date(Year, Month)),
         time_since_start = interval(min(Date), Date) %/% months(1),
         Date = as.Date(make_date(Year, Month)))

# Vienna data
Vienna_data <- read_xlsx(sheet = "Data Vienna", path = file.path(path_data, Vienna_filename), trim_ws = TRUE) %>%
  dplyr::select(Sample_ID, Station_ID, Name, Year, Month, P, T, RH, δ18O, δ2H, d, δ18On, δ2Hn, dn) %>% 
  mutate(Site = "Vienna") %>%
  mutate(Season = "") %>%
  mutate(Season = case_when(
    Month %in% c(12, 1, 2)  ~ "winter",
    Month %in% c(3, 4, 5)  ~ "spring",
    Month %in% c(6, 7, 8)  ~ "summer",
    Month %in% c(9, 10, 11)  ~ "autumn",
    TRUE ~ Season)) %>%
  mutate(Season = factor(Season, levels = c("winter", "spring", "summer", "autumn"))) %>%
  mutate(Date = as.POSIXct(make_date(Year, Month)),
         time_since_start = interval(min(Date), Date) %/% months(1),
         Date = as.Date(make_date(Year, Month)))

########################################################################################################
# Hilbert - Huang transform (HHT) of different time series, ensemble empirical mode decomposition (EEMD)
########################################################################################################
# Prepare MOI data
lines <- readLines(file.path(path_data, MOI_filename))
lines_fixed <- gsub("^(\\d{4})\\.\\s*(\\d{1,2})\\.(\\d{1,2})\\.", "\\1 \\2 \\3 ", lines)
parsed <- lapply(lines_fixed, function(line) {
  nums <- regmatches(line, gregexpr("-?\\d+\\.?\\d*", line))[[1]]
  as.numeric(nums)
})
parsed <- parsed[sapply(parsed, length) == 4]

MOI_timeseries <- as.data.frame(do.call(rbind, parsed)) %>%
  rename(Year = V1, Month = V2, Day = V3, MOI = V4) %>%
  mutate(Date = make_date(year = Year, month = Month)) %>%
  group_by(Date) %>%
  summarise(MOI = mean(MOI, na.rm = TRUE))
MOI_plot <- MOI_timeseries %>%
  mutate(Variable = "MOI")

# Prepare NAO data
NAO_timeseries <- read.csv(file.path(path_data, NAO_filename), header = TRUE, check.names = FALSE, sep = ";") %>% 
  rename(Year = 1) %>% 
  pivot_longer(cols = 2:13, names_to = "Month", values_to = "NAO") %>% 
  mutate(Date = make_date(year = Year, month = Month))

# Prepare HYSPLIT data
Metdata_filtered_Ljubljana <- read_xlsx(file.path(path_data, Ljubljana_filename_metdata), skip = 4) %>%
  rename(Day = "...1") %>%
  mutate(Day = as.Date(Day)) %>%
  mutate(across(where(is.character), ~ ifelse(. == "ne", FALSE, .))) %>%
  mutate(across(where(is.character), ~ ifelse(. == "da", TRUE, .))) %>%
  mutate(across(c(`dež`, `sneg`, `dež s snegom`), ~ as.logical(.))) %>%
  dplyr::filter(., `dež` | `sneg` | `dež s snegom`)
Metdata_filtered_Portoroz <- read_xlsx(file.path(path_data, Portoroz_filename_metdata), skip = 4) %>%
  rename(Day = "...1") %>%
  mutate(Day = as.Date(Day)) %>%
  mutate(across(where(is.character), ~ ifelse(. == "ne", FALSE, .))) %>%
  mutate(across(where(is.character), ~ ifelse(. == "da", TRUE, .))) %>%
  mutate(across(c(`dež`, `sneg`, `dež s snegom`), ~ as.logical(.))) %>%
  dplyr::filter(., `dež` | `sneg` | `dež s snegom`)
Metdata_filtered_Vienna <- read_xlsx(file.path(path_data, Vienna_filename_metdata)) %>% 
  rename(Day = day) %>% 
  mutate(precipitation_TRUE_FALSE = ifelse(precipitation == -1 | precipitation == 0, FALSE, TRUE)) %>% 
  dplyr::filter(precipitation_TRUE_FALSE) %>% 
  mutate(Day = as.Date(Day))

HYSPLIT_timeseries_Ljubljana <- read_xlsx(file.path(path_data, HYSPLIT_filename), sheet = "Ljubljana") %>%
  mutate(Variable = "Mediterranean air masses") %>%
  mutate(Day = as.Date(Date)) %>%
  dplyr::filter(Day %in% Metdata_filtered_Ljubljana$Day) %>% # FILTERS OUT TRAJECTORIES ON DAYS WITHOUT RAIN
  mutate(Date = floor_date(Day, unit = "month")) %>%
  group_by(Date) %>%
  summarize(percent_in_mediterranean = mean(percent_in_mediterranean, na.rm = TRUE)) %>%
  mutate(Rolling_mediterranean = slide_dbl(percent_in_mediterranean, mean, na.rm = TRUE, .before = 120),
         Variable = "Mediterranean air masses",
         Year = year(Date),
         Month = month(Date)) %>%
  mutate(Season = case_when(
    Month %in% c(12, 1, 2) ~ "winter",
    Month %in% c(3, 4, 5) ~ "spring",
    Month %in% c(6, 7, 8) ~ "summer",
    Month %in% c(9, 10, 11) ~ "autumn",
    TRUE ~ NA
  )) %>%
  mutate(Season = factor(Season, levels = c("winter", "spring", "summer", "autumn")))

HYSPLIT_timeseries_Portoroz <- read_xlsx(file.path(path_data, HYSPLIT_filename), sheet = "Portorož") %>%
  mutate(Variable = "Mediterranean air masses") %>%
  mutate(Day = as.Date(Date)) %>%
  dplyr::filter(Day %in% Metdata_filtered_Portoroz$Day) %>% # FILTERS OUT TRAJECTORIES ON DAYS WITHOUT RAIN
  mutate(Date = floor_date(Day, unit = "month")) %>%
  group_by(Date) %>%
  summarize(percent_in_mediterranean = mean(percent_in_mediterranean, na.rm = TRUE)) %>%
  mutate(Rolling_mediterranean = slide_dbl(percent_in_mediterranean, mean, na.rm = TRUE, .before = 120),
         Variable = "Mediterranean air masses",
         Year = year(Date),
         Month = month(Date)) %>%
  mutate(Season = case_when(
    Month %in% c(12, 1, 2) ~ "winter",
    Month %in% c(3, 4, 5) ~ "spring",
    Month %in% c(6, 7, 8) ~ "summer",
    Month %in% c(9, 10, 11) ~ "autumn",
    TRUE ~ NA
  )) %>%
  mutate(Season = factor(Season, levels = c("winter", "spring", "summer", "autumn")))

HYSPLIT_timeseries_Vienna <- read_xlsx(file.path(path_data, HYSPLIT_filename), sheet = "Vienna") %>%
  mutate(Variable = "Mediterranean air masses") %>%
  mutate(Day = as.Date(Date)) %>%
  dplyr::filter(Day %in% Metdata_filtered_Vienna$Day) %>% # FILTERS OUT TRAJECTORIES ON DAYS WITHOUT RAIN
  mutate(Date = floor_date(Day, unit = "month")) %>%
  group_by(Date) %>%
  summarize(percent_in_mediterranean = mean(percent_in_mediterranean, na.rm = TRUE)) %>%
  mutate(Rolling_mediterranean = slide_dbl(percent_in_mediterranean, mean, na.rm = TRUE, .before = 120),
         Variable = "Mediterranean air masses",
         Year = year(Date),
         Month = month(Date)) %>%
  mutate(Season = case_when(
    Month %in% c(12, 1, 2) ~ "winter",
    Month %in% c(3, 4, 5) ~ "spring",
    Month %in% c(6, 7, 8) ~ "summer",
    Month %in% c(9, 10, 11) ~ "autumn",
    TRUE ~ NA
  )) %>%
  mutate(Season = factor(Season, levels = c("winter", "spring", "summer", "autumn")))

HYSPLIT_summary_Ljubljana <- HYSPLIT_timeseries_Ljubljana %>% 
  group_by(Season) %>% 
  summarise(mean_percent = mean(percent_in_mediterranean, na.rm = TRUE))
HYSPLIT_summary_Portoroz <- HYSPLIT_timeseries_Portoroz %>% 
  group_by(Season) %>% 
  summarise(mean_percent = mean(percent_in_mediterranean, na.rm = TRUE))
HYSPLIT_summary_Vienna <- HYSPLIT_timeseries_Vienna %>% 
  group_by(Season) %>% 
  summarise(mean_percent = mean(percent_in_mediterranean, na.rm = TRUE))

# Some basic HYSPLIT trend tests (increasing or decreasing, and Sens's slope)
mk.test(HYSPLIT_timeseries_Ljubljana$Rolling_mediterranean, alternative = "less")
sens.slope(HYSPLIT_timeseries_Ljubljana$Rolling_mediterranean)
mk.test(HYSPLIT_timeseries_Portoroz$Rolling_mediterranean, alternative = "less")
sens.slope(HYSPLIT_timeseries_Portoroz$Rolling_mediterranean)
mk.test(HYSPLIT_timeseries_Vienna$Rolling_mediterranean, alternative = "less")
sens.slope(HYSPLIT_timeseries_Vienna$Rolling_mediterranean)

# Prepare time series and clean NAs for each variable
Ljubljana_data_decomposition <- Ljubljana_data %>%
  left_join(dplyr::select(HYSPLIT_timeseries_Ljubljana, Date, percent_in_mediterranean), by = "Date") %>%
  mutate(percent_in_mediterranean = case_when(
    is.na(percent_in_mediterranean) ~ 0,
    TRUE ~ percent_in_mediterranean
  ))
Portoroz_data_decomposition <- Portoroz_data %>%
  left_join(dplyr::select(HYSPLIT_timeseries_Portoroz, Date, percent_in_mediterranean), by = "Date") %>%
  mutate(percent_in_mediterranean = case_when(
    is.na(percent_in_mediterranean) ~ 0,
    TRUE ~ percent_in_mediterranean
  ))
Vienna_data_decomposition <- Vienna_data %>%
  left_join(dplyr::select(HYSPLIT_timeseries_Vienna, Date, percent_in_mediterranean), by = "Date") %>%
  mutate(percent_in_mediterranean = case_when(
    is.na(percent_in_mediterranean) ~ 0,
    TRUE ~ percent_in_mediterranean
  ))

# Hilbert Huang transform (HHT) and ensemble empirical mode decomposition (EEMD)
# List of variables to be analyzed with HHT & EEMD
variables_Ljubljana <- c("δ18On", "δ2Hn", "dn", "T", "P", "percent_in_mediterranean")
variables_Portoroz <- c("δ18On", "δ2Hn", "dn", "T", "P", "percent_in_mediterranean")
variables_Vienna <- c("δ18On", "δ2Hn", "dn", "T", "P", "percent_in_mediterranean")

# Prepare data once (time vector)
time_vector_Ljubljana <- Ljubljana_data_decomposition$Date[5:length(Ljubljana_data_decomposition$Date)]
time_vector_MOI <- MOI_timeseries$Date
time_vector_NAO <- NAO_timeseries$Date
time_vector_Portoroz <- Portoroz_data_decomposition$Date[10:length(Portoroz_data_decomposition$Date)]
time_vector_Vienna <- Vienna_data_decomposition$Date[1:(length(Vienna_data_decomposition$Date)-13)]

# !!!!THESE STEPS CAN TAKE UP TO 20 MIN!!!!
# 1000 Monte Carlo simulations for each variable (significance test) and analysis for each variable
results_Ljubljana <- lapply(variables_Ljubljana, function(var) {
  message("Analyzing variable: ", var)
  analyze_variable(Ljubljana_data_decomposition[[var]][-c(1:4)], var, time_vector)
})

results_Portoroz <- lapply(variables_Portoroz, function(var) {
  message("Analyzing variable: ", var)
  analyze_variable(Portoroz_data_decomposition[[var]][-c(1:9)], var, time_vector_Portoroz)
})

results_Vienna <- lapply(variables_Vienna, function(var) {
  message("Analyzing variable: ", var)
  analyze_variable(Vienna_data_decomposition[[var]][1:(length(Vienna_data_decomposition[[var]]) - 13)], var, time_vector_Vienna)
})

results_MOI <- lapply("MOI", function(var) {
  analyze_variable(MOI_timeseries[["MOI"]], var, time_vector_MOI)
})

results_NAO <- lapply("NAO", function(var) {
  analyze_variable(NAO_timeseries[["NAO"]], var, time_vector_NAO)
})

eemd_df_list <- Map(function(res, varname) {
  eemd_res <- as.data.frame(res$eemd_result)
  n_imfs <- ncol(eemd_res) - 1
  eemd_res %>%
    mutate(time = time_vector_Ljubljana) %>%
    pivot_longer(cols = -time, names_to = "IMF", values_to = "Value") %>%
    mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs), "IMF_2_3", "Residual")),
           Variable = varname)
}, results_Ljubljana, variables_Ljubljana)

eemd_df_list_MOI <- Map(function(res, varname) {
  eemd_res <- as.data.frame(res$eemd_result)
  n_imfs <- ncol(eemd_res) - 1
  eemd_res %>%
    mutate(time = time_vector_MOI) %>%
    pivot_longer(cols = -time, names_to = "IMF", values_to = "Value") %>%
    mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs), "IMF_2_3", "Residual")),
           Variable = varname)
}, results_MOI, "MOI")

eemd_df_list_NAO <- Map(function(res, varname) {
  eemd_res <- as.data.frame(res$eemd_result)
  n_imfs <- ncol(eemd_res) - 1
  eemd_res %>%
    mutate(time = time_vector_NAO) %>%
    pivot_longer(cols = -time, names_to = "IMF", values_to = "Value") %>%
    mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs), "IMF_2_3", "Residual")),
           Variable = varname)
}, results_NAO, "NAO")

eemd_df_list_Portoroz <- Map(function(res, varname) {
  eemd_res <- as.data.frame(res$eemd_result)
  n_imfs <- ncol(eemd_res) - 1
  eemd_res %>%
    mutate(time = time_vector_Portoroz) %>%
    pivot_longer(cols = -time, names_to = "IMF", values_to = "Value") %>%
    mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs), "IMF_2_3", "Residual")),
           Variable = varname)
}, results_Portoroz, variables_Portoroz)

eemd_df_list_Vienna <- Map(function(res, varname) {
  eemd_res <- as.data.frame(res$eemd_result)
  n_imfs <- ncol(eemd_res) - 1
  eemd_res %>%
    mutate(time = time_vector_Vienna) %>%
    pivot_longer(cols = -time, names_to = "IMF", values_to = "Value") %>%
    mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs), "IMF_2_3", "Residual")),
           Variable = varname)
}, results_Vienna, variables_Vienna)


# Combine all summaries into one big table
full_summary <- do.call(rbind, c(
  lapply(results_Ljubljana, `[[`, "summary"))) %>% 
  #lapply(results_MOI, `[[`, "summary"))) %>%
  mutate(Variable = case_when(
    Variable == "d" ~ "d-excess",
    TRUE ~ Variable)
  ) %>%
  dplyr::select(Variable, IMF, Dominant_Period, p_value, everything()) %>%
  dplyr::select(-Energy) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

full_summary_test <- full_summary %>%
  dplyr::filter(Variable != "MOI")

full_summary_mediterranean <- full_summary %>%
  dplyr::filter(Variable == "percent_in_mediterranean")

full_summary_MOI <- full_summary %>%
  dplyr::filter(Variable == "MOI")

full_summary_NAO <- do.call(rbind, c(
  lapply(results_NAO, `[[`, "summary"))) %>%
  mutate(Variable = case_when(
    Variable == "d" ~ "d-excess",
    TRUE ~ Variable)
  ) %>%
  mutate(IMF = str_replace(IMF, "IMF(\\d+)", "IMF \\1")) %>%
  dplyr::select(Variable, IMF, Dominant_Period, p_value, everything()) %>%
  dplyr::select(-Energy) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

full_summary_Portoroz <- do.call(rbind, c(
  lapply(results_Portoroz, `[[`, "summary"))) %>%
  mutate(Variable = case_when(
    Variable == "d" ~ "d-excess",
    TRUE ~ Variable)
  ) %>%
  mutate(IMF = str_replace(IMF, "IMF(\\d+)", "IMF \\1")) %>%
  dplyr::select(Variable, IMF, Dominant_Period, p_value, everything()) %>%
  dplyr::select(-Energy) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

full_summary_Vienna <- do.call(rbind, c(
  lapply(results_Vienna, `[[`, "summary"))) %>%
  mutate(Variable = case_when(
    Variable == "d" ~ "d-excess",
    TRUE ~ Variable)
  ) %>%
  mutate(IMF = str_replace(IMF, "IMF(\\d+)", "IMF \\1")) %>%
  dplyr::select(Variable, IMF, Dominant_Period, p_value, everything()) %>%
  dplyr::select(-Energy) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

#write.csv(full_summary, "IMF_full_summary_Ljubljana.csv", row.names = FALSE)
#write.csv(full_summary_test, "IMF_full_summary_test.csv", row.names = FALSE)
#write.csv(full_summary_mediterranean, "IMF_full_summary_mediterranean.csv", row.names = FALSE)
#write.csv(full_summary_MOI, "IMF_full_summary_MOI.csv", row.names = FALSE)
#write.csv(full_summary_NAO, "IMF_full_summary_NAO.csv", row.names = FALSE)
#write.csv(full_summary_Portoroz, "IMF_full_summary_Portoroz.csv", row.names = FALSE)
#write.csv(full_summary_Vienna, "IMF_full_summary_Vienna.csv", row.names = FALSE)

# Name the list elements by variable and generate data frames for further plotting
names(eemd_df_list) <- variables_Ljubljana

eemd_df_T <- eemd_df_list[["T"]]
eemd_df_P <- eemd_df_list[["P"]]
eemd_df_mediterranean <- eemd_df_list[["percent_in_mediterranean"]]
eemd_df_MOI <- eemd_df_list_MOI[[1]]
eemd_df_NAO <- eemd_df_list_NAO[[1]]
n_imfs_P <- length(unique(eemd_df_P$IMF)) - 1
eemd_df_d18O <- eemd_df_list[["δ18On"]] %>%
  mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs_P), "Residual", "IMF_2_3"))) %>%
  rbind(eemd_df_T, eemd_df_P, eemd_df_mediterranean)
eemd_df_d2H <- eemd_df_list[["δ2Hn"]] %>%
  mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs_P), "Residual", "IMF_2_3"))) %>%
  rbind(eemd_df_T, eemd_df_P, eemd_df_mediterranean)
eemd_df_d <- eemd_df_list[["dn"]] %>%
  mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs_P), "Residual", "IMF_2_3"))) %>%
  rbind(eemd_df_T, eemd_df_P, eemd_df_mediterranean, eemd_df_MOI) %>%
  mutate(Variable = case_when(
    Variable == "d" ~ "d-excess",
    TRUE ~ Variable)
  )

eemd_df_T_Portoroz <- eemd_df_list_Portoroz[[4]]
eemd_df_P_Portoroz <- eemd_df_list_Portoroz[[5]]
eemd_df_d18O_Portoroz <- eemd_df_list_Portoroz[[1]]
eemd_df_d2H_Portoroz <- eemd_df_list_Portoroz[[2]]
eemd_df_d_Portoroz <- eemd_df_list_Portoroz[[3]]

#### ORTHOGONALITY CHECK ####
# Select only IMF columns from the list for the variable that you want to check (in this case [[1]] chooses the d18O variable)
orthogonality_df <- eemd_df_list_Ljubljana[[1]] %>%
  dplyr::select(time, IMF, Value) %>%
  pivot_wider(names_from = IMF, values_from = Value)
imfs <- orthogonality_df %>% dplyr::select(-time)

# Function to compute orthogonality index between two vectors
orth_index <- function(a, b) {
  sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
}

# Initialize empty matrix and fill it
n <- ncol(imfs)
oi_matrix <- matrix(0, n, n)

for (i in 1:n) {
  for (j in 1:n) {
    oi_matrix[i, j] <- orth_index(imfs[[i]], imfs[[j]])
  }
}

# Name the matrix and view it (END OF ORTHOGONALITY CHECK)
colnames(oi_matrix) <- rownames(oi_matrix) <- colnames(imfs)
print(round(oi_matrix, 3))

# Data frames with only significant IMFs
eemd_df_T_significant <- eemd_df_list[["T"]] %>%
  semi_join(full_summary %>% dplyr::filter(Variable == "T", p_value < 0.05) %>% dplyr::select(IMF),
            by = "IMF")
eemd_df_P_significant <- eemd_df_list[["P"]] %>%
  semi_join(full_summary %>% dplyr::filter(Variable == "P", p_value < 0.05) %>% dplyr::select(IMF),
            by = "IMF")
eemd_df_mediterranean_significant <- eemd_df_list[["percent_in_mediterranean"]] %>%
  semi_join(full_summary %>% dplyr::filter(Variable == "percent_in_mediterranean", p_value < 0.05) %>% dplyr::select(IMF),
            by = "IMF")
eemd_df_d18O_significant <- eemd_df_list[["δ18On"]] %>%
  mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs_P), "IMF_2_3", "Residual"))) %>%
  semi_join(full_summary %>% dplyr::filter(Variable == "δ18On", p_value < 0.05) %>% dplyr::select(IMF),
            by = "IMF")
eemd_df_d2H_significant <- eemd_df_list[["δ2Hn"]] %>%
  mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs_P), "IMF_2_3", "Residual"))) %>%
  semi_join(full_summary %>% dplyr::filter(Variable == "δ2Hn", p_value < 0.05) %>% dplyr::select(IMF),
            by = "IMF")
eemd_df_d_significant <- eemd_df_list[["dn"]] %>%
  mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs_P), "IMF_2_3", "Residual"))) %>%
  semi_join(full_summary %>% dplyr::filter(Variable == "dn", p_value < 0.05) %>% dplyr::select(IMF),
            by = "IMF") %>%
  rbind(eemd_df_T_significant, dplyr::filter(eemd_df_P_significant, IMF %in% c("IMF3", "IMF4", "IMF5"))) %>%
  mutate(Variable = case_when(
    Variable == "dn" ~ "d-excess",
    Variable == "δ2Hn" ~ "δ2H",
    Variable == "δ18On" ~ "δ18O",
    TRUE ~ Variable)
  ) %>%
  mutate(IMF = case_when(
    IMF == "IMF2" & Variable == "T" ~ "IMF3",
    TRUE ~ IMF
  )) %>%
  mutate(IMF = case_when(
    IMF == "IMF2" ~ "1-year oscillation (IMF 2)",
    IMF == "IMF3" ~ "1-year oscillation (IMF 3)",
    IMF == "IMF_2_3" ~ "1-year oscillation (IMF 2 + IMF 3)",
    IMF == "IMF4" ~ "2.6-year oscillation (IMF 4)",
    IMF == "IMF5" ~ "7.5-year oscillation (IMF 5)",
    TRUE ~ IMF
  )) %>%
  mutate(Value = case_when(
    Variable == "P" ~ Value/50,
    Variable == "T" ~ Value/5,
    TRUE ~ Value
  ))

eemd_df_significant_plot <- eemd_df_d18O_significant %>%
  rbind(eemd_df_d2H_significant, dplyr::filter(eemd_df_d, Variable == "dn" & IMF == "IMF3")) %>%
  mutate(Variable = case_when(
    Variable == "dn" ~ "d-excess",
    Variable == "δ2Hn" ~ "δ2H",
    Variable == "δ18On" ~ "δ18O",
    TRUE ~ Variable)
  ) %>%
  mutate(Facet = Variable) %>%
  rbind(mutate(eemd_df_T_significant, Facet = "δ18O"),
        mutate(eemd_df_T_significant, Facet = "δ2H"),
        mutate(eemd_df_T_significant, Facet = "d-excess")) %>%
  mutate(IMF = case_when(
    IMF == "IMF2" ~ "1-year oscillation (IMF 2)",
    IMF == "IMF3" ~ "1-year oscillation (IMF 3)",
    IMF == "IMF_2_3" ~ "1-year oscillation (IMF 2 + IMF 3)",
    IMF == "IMF4" ~ "2.6-year oscillation (IMF 4)",
    IMF == "IMF5" ~ "7.5-year oscillation (IMF 5)",
    TRUE ~ IMF
  ))
eemd_df_significant_plot_2 <- eemd_df_mediterranean_significant %>%
  rbind(eemd_df_list[["dn"]] %>%
          mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:n_imfs_P), "Residual", "IMF_2_3"))) %>%
          semi_join(full_summary %>% dplyr::filter(Variable == "d-excess", p_value < 0.05) %>% dplyr::select(IMF),
                    by = "IMF")) %>%
  mutate(Variable = case_when(
    Variable == "d" ~ "d-excess",
    Variable == "percent_in_mediterranean" ~ "% in Mediterranean",
    TRUE ~ Variable)
  ) %>%
  dplyr::filter(IMF == c("IMF 5", "IMF 6")) %>%
  mutate(IMF = case_when(
    IMF == "IMF 5" & Variable == "d-excess" ~ "IMF 6 (% Mediterranean, 15-year) & IMF 5 (d-excess, 7.5-year)",
    IMF == "IMF 6" & Variable == "% in Mediterranean" ~ "IMF 6 (% Mediterranean, 15-year) & IMF 5 (d-excess, 7.5-year)",
    TRUE ~ IMF
  ))

# Generate amplitude_data frames for each variable and assign to objects
amplitude_df_list <- Map(function(res, var) {
  create_amplitude_df(res$eemd_result, time_vector_Ljubljana, var)
}, results_Ljubljana, variables_Ljubljana)
names(amplitude_df_list) <- variables_Ljubljana

amplitude_df_list_Portoroz <- Map(function(res, var) {
  create_amplitude_df(res$eemd_result, time_vector_Portoroz, var)
}, results_Portoroz, variables_Portoroz)
names(amplitude_df_list_Portoroz) <- variables_Portoroz

amplitude_df_list_Vienna <- Map(function(res, var) {
  create_amplitude_df(res$eemd_result, time_vector_Vienna, var)
}, results_Vienna, variables_Vienna)
names(amplitude_df_list_Vienna) <- variables_Vienna

amplitude_df_d18O <- amplitude_df_list[["δ18On"]]
amplitude_df_d2H <- amplitude_df_list[["δ2Hn"]]
amplitude_df_d <- amplitude_df_list[["dn"]]
amplitude_df_T <- amplitude_df_list[["T"]]
amplitude_df_P <- amplitude_df_list[["P"]]
amplitude_df_mediterranean <- amplitude_df_list[["percent_in_mediterranean"]]

amplitude_df_d18O_Portoroz <- amplitude_df_list_Portoroz[["δ18On"]]
amplitude_df_d2H_Portoroz <- amplitude_df_list_Portoroz[["δ2Hn"]]
amplitude_df_d_Portoroz <- amplitude_df_list_Portoroz[["dn"]]
amplitude_df_T_Portoroz <- amplitude_df_list_Portoroz[["T"]]
amplitude_df_P_Portoroz <- amplitude_df_list_Portoroz[["P"]]
amplitude_df_mediterranean_Portoroz <- amplitude_df_list_Portoroz[["percent_in_mediterranean"]]

amplitude_df_d18O_Vienna <- amplitude_df_list_Vienna[["δ18On"]]
amplitude_df_d2H_Vienna <- amplitude_df_list_Vienna[["δ2Hn"]]
amplitude_df_d_Vienna <- amplitude_df_list_Vienna[["dn"]]
amplitude_df_T_Vienna <- amplitude_df_list_Vienna[["T"]]
amplitude_df_P_Vienna <- amplitude_df_list_Vienna[["P"]]
amplitude_df_mediterranean_Vienna <- amplitude_df_list_Vienna[["percent_in_mediterranean"]]

# Data frames LJUBLJANA with only significant IMFs
amplitude_df_d18O_significant <- amplitude_df_list[["δ18On"]] %>%
  dplyr::filter(IMF == "Residual" | IMF %in% (full_summary %>% dplyr::filter(Variable == "δ18On", p_value < 0.05) %>% pull(IMF)))
amplitude_df_d2H_significant <- amplitude_df_list[["δ2Hn"]] %>%
  dplyr::filter(IMF == "Residual" | IMF %in% (full_summary %>% dplyr::filter(Variable == "δ2Hn", p_value < 0.05) %>% pull(IMF)))
amplitude_df_d_significant <- amplitude_df_list[["dn"]] %>%
  mutate(Variable = case_when(
    Variable == "dn" ~ "d-excess",
    TRUE ~ Variable)
  ) %>%
  dplyr::filter(IMF == "Residual" | IMF %in% (full_summary %>% dplyr::filter(Variable == "dn", p_value <= 0.05) %>% pull(IMF)))
amplitude_df_T_significant <- amplitude_df_list[["T"]] %>%
  dplyr::filter(IMF == "Residual" | IMF %in% (full_summary %>% dplyr::filter(Variable == "T", p_value < 0.05) %>% pull(IMF)))
amplitude_df_P_significant <- amplitude_df_list[["P"]] %>%
  dplyr::filter(IMF == "Residual" | IMF %in% (full_summary %>% dplyr::filter(Variable == "P", p_value < 0.05) %>% pull(IMF)))
amplitude_df_mediterranean_significant <- amplitude_df_list[["percent_in_mediterranean"]] %>%
  dplyr::filter(IMF == "Residual" | IMF %in% (full_summary %>% dplyr::filter(Variable == "percent_in_mediterranean", p_value <= 0.05) %>% pull(IMF))) %>%
  mutate(IMF = case_when(
    IMF == "IMF3" ~ "1-year oscillation (IMF 3)",
    IMF == "IMF6" ~ "15-year oscillation (IMF 6)",
    TRUE ~ IMF
  ))

amplitude_df_plot <- amplitude_df_d18O_significant %>%
  rbind(amplitude_df_d2H_significant, amplitude_df_d_significant) %>%
  mutate(Variable = case_when(
    Variable == "dn" ~ "d-excess",
    Variable == "δ2Hn" ~ "δ2H",
    Variable == "δ18On" ~ "δ18O",
    TRUE ~ Variable)
  ) %>%
  mutate(Variable = factor(Variable, levels = c("δ18O", "δ2H", "d-excess"))) %>%
  mutate(IMF = case_when(
    IMF == "IMF2" ~ "1-year oscillation (IMF 2)",
    IMF == "IMF3" ~ "1-year oscillation (IMF 3)",
    IMF == "IMF_2_3" ~ "1-year oscillation (IMF 2 + IMF 3)",
    IMF == "IMF4" ~ "2.6-year oscillation (IMF 4)",
    IMF == "IMF5" ~ "7.5-year oscillation (IMF 5)",
    TRUE ~ IMF
  ))

# FOR POWER SPECTRUM
imf2_values <- amplitude_df_list[["dn"]] %>%
  dplyr::filter(IMF == "IMF3", Type == "IMF oscillation") %>%
  pull(Value)
spec <- spectrum(imf2_values, plot = FALSE)
spec_df <- data.frame(
  Frequency = spec$freq,
  Power = spec$spec
)
max_freq <- spec_df$Frequency[which.max(spec_df$Power)]

###############################################################################
# PLOTTING
###############################################################################
# GRAPHICAL ABSTRACT MAP & FIGURE 1
###############################################################################
# Creating and obtaining shapefiles for map plots
borders_sf <- ne_countries(scale = 10, continent = "Europe")
names_sf <- data.frame(
  Station_name = c("SVN", "AUT", "HRV", "ITA", "HUN"),
  lon = c(14.9, 15, 15.35, 12.5, 17.2),
  lat = c(46.3, 47.5, 45, 46.1, 46.8)
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
stations_sf <- data.frame(
  Station_name = c("Ljubljana", "Vienna", "Portorož"),
  lon = c(14.512352, 16.37250, 13.5799694),
  lat = c(46.065507, 48.20833, 45.5166333)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# GRAPHICAL ABSTRACT MAP
ggplot() +
  geom_sf(data = borders_sf, fill = NA) +
  geom_sf(data = dplyr::filter(borders_sf, sov_a3 == "SVN"), fill = "grey85") +
  geom_text(data = names_sf, aes(x = st_coordinates(geometry)[,1], 
                                 y = st_coordinates(geometry)[,2], 
                                 label = Station_name),
            size = 1.8) +
  geom_sf(data = stations_sf, fill = "black", size = 0.6, shape = 21, stroke = 0.4) + 
  geom_text(data = stations_sf, aes(x = st_coordinates(geometry)[,1], 
                                    y = st_coordinates(geometry)[,2], 
                                    label = Station_name),
            color = "black", size = 1.6, fontface = "bold", nudge_y = -0.15) +
  coord_sf(xlim = c(12, 17.5), ylim = c(44.5, 48.5)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank())

ggsave(filename = paste("Map_slovenia_graphical_abstract.png"), path = path_figures, device = "png", dpi = 1500, width = 60, height = 50, units = "mm")

# FIGURE 1 MAP
ggplot() +
  geom_sf(data = borders_sf, fill = NA) +
  geom_sf(data = dplyr::filter(borders_sf, sov_a3 == "SVN"), fill = "grey85") +
  geom_text(data = names_sf, aes(x = st_coordinates(geometry)[,1],
                                 y = st_coordinates(geometry)[,2],
                                 label = Station_name),
            size = 2.5) +
  ggspatial::annotation_scale(style = "ticks", location = "bl",
                              height = unit(0.12, "cm"),
                              pad_x = unit(3, "mm"),
                              pad_y = unit(3, "mm"),
                              text_cex = 0.55,
                              line_width = 0.3) +
  ggspatial::annotation_north_arrow(which_north = "true", location = "tr",
                                    pad_x = unit(-0.8, "mm"),
                                    pad_y = unit(2.5, "mm"),
                                    height = unit(1, "cm"),
                                    width = unit(1, "cm"),
                                    style = ggspatial::north_arrow_fancy_orienteering(
                                      fill = c("grey40", "white"),
                                      line_col = "grey20",
                                      text_size = 6,
                                      line_width = 0.3)) +
  geom_sf(data = stations_sf, fill = "black", size = 1, shape = 21, stroke = 0.4) + 
  geom_text(data = stations_sf, aes(x = st_coordinates(geometry)[,1], 
                                    y = st_coordinates(geometry)[,2], 
                                    label = Station_name),
            color = "black", size = 2.5, fontface = "bold", nudge_y = -0.15) +
  coord_sf(xlim = c(12, 17.5), ylim = c(44.5, 48.5)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())

ggsave(filename = paste("Map_slovenia",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 80, units = "mm")

# Reading climate data
climate_data <- read_xlsx(path = file.path(path_data, Climate_filename), sheet = 1, skip = 1) %>% 
  mutate(Site = "Ljubljana") %>% 
  rbind(read_xlsx(path = file.path(path_data, Climate_filename), sheet = 2, skip = 1) %>% 
          mutate(Site = "Portorož"),
        read_xlsx(path = file.path(path_data, Climate_filename), sheet = 3, skip = 1) %>% 
          mutate(Site = "Vienna")) %>% 
  rename("T [°C]" = `T / °C`,
         "P [mm]" = `P / mm`) %>% 
  mutate(Month = factor(Month, levels = month.name, labels = month.abb)) %>% 
  pivot_longer(cols = c ("P [mm]", "T [°C]"), names_to = "Variable", values_to = "Value") %>% 
  mutate(Variable =  factor(Variable, levels = c("T [°C]", "P [mm]")))

climate_summary <- climate_data %>% 
  group_by(Site, Variable) %>% 
  summarise(Mean = mean(Value),
            StDev = sd(Value)) %>% 
  mutate(label = sprintf("Mean = %.1f\n SD = %.1f", Mean, StDev), Month = "Jan", y = Inf)

# FIGURE 1 MONTHLY LINE PLOTS
for (station in unique(climate_data$Site)) {
  ggplot(data = dplyr::filter(climate_data, Site == station), aes(x = Month, y = Value, color = Variable, group = Variable)) +
    geom_line(linewidth = 0.3) + 
    geom_text(data = dplyr::filter(climate_summary, Site == station), aes(x = Month, y = y, label = label),
              hjust = -0.1, vjust = 1.2, inherit.aes = FALSE,
              size = 1.7) +
    facet_wrap(~ Variable, ncol = 1, scales = "free_y") + 
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_line(linewidth = 0.2),
          panel.grid.minor.y = element_blank(),
          panel.grid = element_line(linewidth = 0.2),
          legend.position = "",
          text = element_text(size = 6, color = "black"),
          strip.text = element_text(margin = margin(0.03, 0, 0.03, 0, unit = "cm")),
    )
  ggsave(filename = paste0("Climate_", station, ".", figure_type), path = path_figures, device = figure_type, dpi = 1500, width = 44, height = 33, units = "mm")
}

###############################################################################
# FIGURE 2a
###############################################################################
ggplot(data = dplyr::filter(Ljubljana_data_longest, Variable != "3H (TU)" & Variable != "10-year rolling 3H"), aes(x = Date, y = Value, color = Variable_color, linetype = Rolling)) +
  geom_line(linewidth = 0.3) +
  labs(y = expression("\u03B4"^"18"*"O, \u03B4"^"2"*"H, and d-excess: [\u2030]    T: [°C]     P: [mm]")) +
  guides(color = "none") +
  scale_x_date(
    limits = as.Date(c("1981-01-01", "2025-01-01")),
    breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) + # %b for month
  scale_y_continuous(labels = ~sub("-", "\u2212", .x)) +
  scale_color_manual(values = c("#00BA38", "#F564E3","#619CFF", "#F8766D", "#B79F00")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2, labeller = label_parsed) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 7, colour = "black"),
        strip.text = element_text(size = 7, colour = "black", margin = margin(0.02, 0, 0.02, 0, unit = "cm")),
        legend.position = "top",
        legend.box.spacing = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.key.spacing.x = unit(0.7, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        #panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggsave(filename = paste("Ljubljana_basic",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 180, height = 110, units = "mm")

###############################################################################
# FIGURE 2b
###############################################################################
ggplot(data = dplyr::filter(Ljubljana_data_longer, Variable != "3H (TU)" & is.na(Variable) == FALSE), aes(x = Month, y = Value, colour = Variable, group = Month)) +
  geom_violin(trim = FALSE, scale = "area", alpha = 0.5, linewidth = 0.3) +  # Trim = FALSE to avoid cutting off the tails
  geom_point(data = Ljubljana_data_longer_month, aes(x = Month, y = Value, group = Month), color = "black", size = 0.6) +
  labs(y = expression("\u03B4"^"18"*"O, \u03B4"^"2"*"H, and d-excess: [\u2030]    T: [°C]     P: [mm]")) +
  scale_color_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#619CFF", "#F564E3")) +
  scale_x_continuous(breaks = seq(1, 12, 1),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
                     expand = expansion(mult = 0.01, add = 0)) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1, labeller = label_parsed) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 7, colour = "black"),
        strip.text = element_text(size = 6, colour = "black", margin = margin(0.02, 0, 0.02, 0, unit = "cm")),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(filename = paste("Ljubljana_violin_month",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 110, units = "mm")

###############################################################################
# FIGURE 3
###############################################################################
ggplot(data = Ljubljana_data_slopes_longer_2, aes(x = Year, linetype = weighting)) +
  geom_line(aes(y = value), linewidth = 0.25) +
  geom_segment(data = dplyr::filter(Ljubljana_data_slopes_longer_2, weighting == "Unweighted"),
              aes(x = 1981,  xend = 2024, y = Intercept + Slope * 1981, yend = Intercept + Slope * 2024), linewidth = 0.3, color = "dodgerblue3", inherit.aes = FALSE) +
  labs(y = expression("Slope: [1/\u2030],    intercept: [\u2030]")) +
  scale_x_continuous(breaks = seq(1981, 2025, 5)) +
  scale_y_continuous(labels = ~sub("-", "\u2212", .x)) +
  scale_linetype_manual(values = c(2,1)) +
  facet_grid(statistic ~ facet, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 7, colour = "black"),
        strip.text = element_text(size = 6, colour = "black", , margin = margin(0.03, 0.03, 0.03, 0.03, unit = "cm")),
        legend.position = "top",
        legend.box.spacing = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.key.spacing.y = unit(0.7, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.margin = margin(t = 0, r = 23, b = 0, l = 0),
        legend.text = element_text(size = 6, colour = "black"),
        legend.title = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = paste("Ljubljana_slopes",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 90, units = "mm")

###############################################################################
# DECOMPOSITION PLOTS (GRAPHICAL ABSTRACT, FIGURE 4)
###############################################################################
# DECOMPOSITION GRAPHICAL ABSTRACT #1
abstract_temperature <- amplitude_df_T_significant %>% 
  dplyr::filter(IMF == "IMF2" | IMF == "Residual") %>% 
  dplyr::filter(Type == "IMF oscillation") %>% 
  mutate(Value = case_when(
    IMF == "Residual" ~ Value-20,
    TRUE ~ Value
  )) %>% 
  mutate(Variable_facet = "italic(delta)^18*O") %>% 
  mutate(IMF_facet = case_when(
    IMF == "IMF2" ~  '"1-year oscillation"',
    IMF == "Residual" ~ '"Residual trend"',
    TRUE ~ NA
  ))
amplitude_df_plot_filtered_2 <- amplitude_df_plot_filtered %>% 
  dplyr::filter(Variable == "δ18O" & Type == "IMF oscillation") %>% 
  mutate(IMF_facet = case_when(
    IMF_facet == '"1-year oscillation (IMF 2 + IMF 3)"' ~ '"1-year oscillation"',
    IMF_facet == "Residual" ~ '"Residual trend"',
    TRUE ~ IMF_facet
  )) %>% 
  rbind(abstract_temperature)

ggplot() +
  geom_line(data = amplitude_df_plot_filtered_2, aes(x = time, y = Value, group = Variable, color = Variable), linewidth = 0.2) +
  scale_x_date(
    limits = as.Date(c("1981-01-01", "2025-01-01")),
    breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) + # %b for month
  scale_color_manual(values = c("grey30", "dodgerblue3")) +
  scale_y_continuous(labels = ~sub("-", "\u2212", .x)) +
  facet_wrap(~ IMF_facet, ncol = 1, scales = "free_y", labeller = label_parsed) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, colour = "black", angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 14, colour = "black", margin = margin(0.3,0.3,0.3,0.3, "mm")),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 9.5), "pt"),
        legend.position = "",
        legend.title = element_blank()) +
  labs(y = expression("[\u2030]"))

ggsave(filename = paste("Ljubljana_decomposition_ABSTRACT",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 80, height = 60, units = "mm")

# DECOMPOSITION GRAPHICAL ABSTRACT #2
ggplot() +
  geom_line(data = mutate(Ljubljana_data, facet = "Raw trend"), aes(x = Date, y = δ18O), color = "grey30", linewidth = 0.2) +
  scale_x_date(
    limits = as.Date(c("1981-01-01", "2025-01-01")),
    breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) + # %b for month
  scale_y_continuous(labels = ~sub("-", "\u2212", .x)) +
  facet_wrap(~ facet) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 14, colour = "black", margin = margin(0.3,0.3,0.3,0.3, "mm")),
        panel.grid.minor.y = element_blank(),
        legend.position = "",
        legend.title = element_blank()) +
  labs(y = expression("[\u2030]"))

ggsave(filename = paste("Ljubljana_decomposition_ABSTRACT_2",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 80, height = 30, units = "mm")

# FIGURE 4a
amplitude_df_plot_filtered <- amplitude_df_plot %>% 
  dplyr::filter(!(IMF %in% c("1-year oscillation (IMF 2)", "1-year oscillation (IMF 3)") & Variable %in% c("δ18O", "δ2H"))) %>% 
  dplyr::filter(!(IMF == "1-year oscillation (IMF 2 + IMF 3)" & Variable == "d-excess")) %>%
  mutate(Variable_facet = Variable,
         IMF_facet = IMF) %>% 
  mutate(Variable_facet = case_when(
    Variable_facet == "δ18O" ~ "italic(delta)^18*O",
    Variable_facet == "δ2H" ~ "italic(delta)^2*H",
    Variable_facet == "d-excess" ~ '"d-excess"',
    TRUE ~ Variable_facet
  )) %>% 
  mutate(Variable_facet = factor(Variable_facet, levels = c("italic(delta)^18*O", "italic(delta)^2*H", '"d-excess"'))) %>%
  mutate(IMF_facet = case_when(
    IMF_facet == "1-year oscillation (IMF 2 + IMF 3)" ~ '"1-year oscillation (IMF 2 + IMF 3)"',
    IMF_facet == "1-year oscillation (IMF 2)" ~ '"1-year oscillation (IMF 2)"',
    IMF_facet == "1-year oscillation (IMF 3)" ~ '"1-year oscillation (IMF 3)"',
    IMF_facet == "2.6-year oscillation (IMF 4)" ~ '"2.6-year oscillation (IMF 4)"',
    IMF_facet == "7.5-year oscillation (IMF 5)" ~ '"7.5-year oscillation (IMF 5)"',
    IMF_facet == "Resiudal" ~ '"Residual"',
    TRUE ~ IMF_facet
  ))
ggplot() +
  geom_line(data = amplitude_df_plot_filtered, aes(x = time, y = Value, group = interaction(Variable, Type), linetype = Type, color = Variable_facet), linewidth = 0.2) +
  scale_color_manual(values = c("italic(delta)^18*O" = "#F8766D",
                                "italic(delta)^2*H" = "#B79F00",
                                '"d-excess"' = "#00BA38"), labels = scales::parse_format()) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_x_date(
    limits = as.Date(c("1981-01-01", "2025-01-01")),
    breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) + # %b for month
  scale_y_continuous(labels = ~sub("-", "\u2212", .x)) +
  guides(color = guide_legend(override.aes = list(linewidth = 0.4))) +
  facet_wrap(Variable_facet ~ IMF_facet, ncol = 2, scales = "free_y", labeller = label_parsed) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, colour = "black", angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 7, colour = "black"),
        strip.text = element_text(size = 6, colour = "black", margin = margin(0.3,0.3,0.3,0.3, "mm")),
        panel.grid.minor.y = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.spacing = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(size = 6, colour = "black", margin = margin(t = 0, r = 0, b = 0, l = 0.5)),
        legend.key.spacing.x = unit(0.3, units = "cm"), 
        legend.key.width = unit(0.3, units = "cm"),
        legend.spacing.x =  unit(1.5, units = "cm"),
        legend.title = element_blank()) +
  labs(y = expression("[\u2030]"))

ggsave(filename = paste("Ljubljana_decomposition_significant",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 110, height = 110, units = "mm")

# FIGURE 4b
ggplot(spec_df, aes(x = Frequency*12, y = Power)) +
  geom_line(linewidth = 0.2) +
  geom_vline(xintercept = max_freq*12, color = "#00BA38", linewidth = 0.25) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6)) + 
  annotate("text", x = max_freq*12, y = max(spec_df$Power), label = sprintf("Maximum frequency: %.1f", max_freq*12), vjust = 0.7, hjust = -0.1, color = "#00BA38", size = 2) +
  labs(title = "Power Spectrum of IMF 3 (d-excess)", x = "Frequency (cycles per year)", y = "Spectral Density") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, colour = "black"),
        axis.title.x = element_text(size = 7, colour = "black"),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 7, colour = "black"),
        strip.text = element_text(size = 6, colour = "black", margin = margin(0.5,0.5,0.5,0.5, "mm")),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 7.5),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.spacing = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(size = 6, colour = "black", margin = margin(t = 0, r = 0, b = 0, l = 0.5)),
        legend.key.spacing.x = unit(0.3, units = "cm"), 
        legend.key.width = unit(0.2, units = "cm"),
        legend.spacing.x =  unit(1.5, units = "cm"),
        legend.title = element_blank()) 

ggsave(filename = paste("Ljubljana_power_spectrum_d_IMF3",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 70, height = 40, units = "mm")

# FIGURE 4c
eemd_df_significant_plot_filtered <- eemd_df_significant_plot %>%
  dplyr::filter(!(IMF %in% c("1-year oscillation (IMF 2)", "1-year oscillation (IMF 3)") & Variable %in% c("δ18O", "δ2H"))) %>% 
  dplyr::filter(!(IMF == "1-year oscillation (IMF 2 + IMF 3)" & Variable == "T")) %>%
  dplyr::filter(IMF != "Residual") %>%
  mutate(Facet = case_when(
    Facet == "δ18O" ~ "italic(delta)^18*O",
    Facet == "δ2H" ~ "italic(delta)^2*H",
    Facet == "d-excess" ~ '"d-excess"',
    TRUE ~ Facet
  )) %>%
  mutate(Variable_color = Variable) %>% 
  mutate(Variable_color = case_when(
    Variable_color == "δ18O" ~ "italic(delta)^18*O",
    Variable_color == "δ2H" ~ "italic(delta)^2*H",
    Variable_color == "d-excess" ~ '"d-excess"',
    Variable_color == "T" ~ '"T"',
    TRUE ~ Variable_color
  )) %>% 
  mutate(Variable_color = factor(Variable_color, levels = c("italic(delta)^18*O", "italic(delta)^2*H", '"d-excess"', '"T"'))) %>% 
  mutate(Facet = factor(Facet, levels = c("italic(delta)^18*O", "italic(delta)^2*H", '"d-excess"')))
ggplot() +
  geom_line(data = eemd_df_significant_plot_filtered, aes(x = time, y = Value, color = Variable_color), linewidth = 0.2) +
  scale_color_manual(values = c("italic(delta)^18*O" = "#F8766D",
                                "italic(delta)^2*H" = "#B79F00",
                                '"d-excess"' = "#00BA38",
                                '"T"' = "grey30"), labels = scales::parse_format()) + 
  scale_x_date(
    limits = as.Date(c("1981-01-01", "2025-01-01")),
    breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) + # %b for month
  scale_y_continuous(labels = ~sub("-", "\u2212", .x)) +
  guides(color = guide_legend(override.aes = list(linewidth = 0.4))) +
  facet_wrap (~ Facet, ncol = 1, scales = "free_y", labeller = label_parsed) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, colour = "black", angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 7, colour = "black"),
        strip.text = element_text(size = 6, colour = "black", margin = margin(0.2,0.2,0.2,0.2, "mm")),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"),
        legend.spacing.y = unit(0, units = "cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.spacing = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(size = 6, colour = "black", margin = margin(t = 0, r = 0, b = 0, l = 0.5)),
        legend.key.spacing.x = unit(0.3, units = "cm"), 
        legend.key.width = unit(0.3, units = "cm"),
        legend.key.height = unit(0.4, units = "cm"),
        legend.spacing.x =  unit(1.5, units = "cm"),
        legend.title = element_blank()) +
  labs(y = expression(italic(delta)^18 * O * ", " * italic(delta)^2 * H * " and d-excess [\u2030], T [°C]"))

ggsave(filename = paste("Ljubljana_decomposition_significant_TEMPERATURE",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 70, height = 65, units = "mm")

###############################################################################
# HYSPLIT PLOTS (FIGURE 5 & supplemental figures)
###############################################################################
# Figure 5a, upper plot
ggplot() +
  geom_line(data = HYSPLIT_timeseries_Ljubljana, aes(x = Date, y = percent_in_mediterranean), color = "salmon4", linewidth = 0.3) +
  scale_x_date(
    limits = as.Date(c("1981-05-01", "2025-01-01")),
    breaks = seq(as.Date("1981-05-01"), as.Date("2025-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) +
  facet_wrap(~ Variable) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, colour = "black", angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 7, colour = "black", margin = margin(0.5,0.5,0.5,0.5, "mm")),
        legend.position = "none") +
  labs(y = expression("%"))

ggsave(filename = paste("Ljubljana_HYSPLIT_timetrend",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 40, units = "mm")

# Figure 5a, middle plot
ggplot() +
  geom_line(data = dplyr::mutate(filter(HYSPLIT_timeseries_Ljubljana, Year > 1990), Variable = "10-year rolling Mediterranean air masses"), aes(x = Date, y = Rolling_mediterranean), color = "salmon4", linewidth = 0.3) +
  scale_x_date(
    limits = as.Date(c("1981-05-01", "2025-01-01")),
    breaks = seq(as.Date("1981-05-01"), as.Date("2025-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, colour = "black", angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7, colour = "black", margin = margin(0.5,0.5,0.5,0.5, "mm")),
        legend.position = "none") +
  labs(y = expression("%"))

ggsave(filename = paste("Ljubljana_HYSPLIT_rolling_timetrend",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 40, units = "mm")

# Figure 5a, bottom plot
ggplot(data = HYSPLIT_timeseries_Ljubljana, aes(x = Season, y = percent_in_mediterranean, colour = Variable, group = Season), linewidth = 0.3) +
  geom_violin(trim = FALSE, scale = "area", alpha = 0.5, width = 0.3, linewidth = 0.3, position = position_nudge(0.165)) +  # Trim = FALSE to avoid cutting off the tails
  geom_boxplot(outlier.size = 0.5, linewidth = 0.3, width = 0.3, position = position_nudge(-0.165)) +
  geom_point(data =HYSPLIT_summary_Ljubljana, aes(x = Season, y = mean_percent, group = Season), color = "black", size = 0.8, position = position_nudge(0.155)) +
  scale_color_manual(values = c("salmon4")) +
  labs(y = expression("%")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, colour = "black", hjust = 0.5, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 7, colour = "black", margin = margin(0.5,0.5,0.5,0.5, "mm")),
        legend.position = "none")

ggsave(filename = paste("Ljubljana_HYSPLIT_seasonality",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 40, units = "mm")

# Figure 5b, bottom plot
ggplot() +
  geom_line(data = dplyr::filter(amplitude_df_mediterranean_significant, IMF != "IMF_2_3"), aes(x = time, y = Value, group = interaction(Variable, Type), linetype = Type), linewidth = 0.2, color = "black") +
  scale_linetype_manual(values = c(2, 1)) +
  scale_x_date(
    limits = as.Date(c("1981-01-01", "2025-01-01")),
    breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) + # %b for month
  facet_wrap(~ IMF, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, colour = "black", angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 7, colour = "black"),
        strip.text = element_text(size = 6, colour = "black", margin = margin(0.5,0.5,0.5,0.5, "mm")),
        panel.grid.minor.y = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.spacing = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(size = 6, colour = "black", margin = margin(t = 0, r = 0, b = 0, l = 0.5)),
        legend.key.spacing.x = unit(0.3, units = "cm"), 
        legend.key.width = unit(0.2, units = "cm"),
        legend.spacing.x =  unit(1.5, units = "cm"),
        legend.title = element_blank()) +
  labs(y = expression("[%]"))

ggsave(filename = paste("Ljubljana_decomposition_significant_HYSPLIT",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 120, units = "mm")

# Portorož HYSPLIT supplemental Figure a, upper plot
ggplot() +
  geom_line(data = HYSPLIT_timeseries_Portoroz, aes(x = Date, y = percent_in_mediterranean), color = "salmon4", linewidth = 0.3) +
  scale_x_date(
    limits = as.Date(c("2000-01-01", "2025-01-01")),
    breaks = seq(as.Date("2000-01-01"), as.Date("2025-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("2000-01-01"), as.Date("2025-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) +
  facet_wrap(~ Variable) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, colour = "black", angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 7, colour = "black", margin = margin(0.5,0.5,0.5,0.5, "mm")),
        legend.position = "none") +
  labs(y = expression("%"))

ggsave(filename = paste("Portoroz_HYSPLIT_timetrend",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 40, units = "mm")

# Portorož HYSPLIT supplemental Figure a, bottom plot
ggplot(data = HYSPLIT_timeseries_Portoroz, aes(x = Season, y = percent_in_mediterranean, colour = Variable, group = Season), linewidth = 0.3) +
  geom_violin(trim = FALSE, scale = "area", alpha = 0.5, width = 0.3, linewidth = 0.3, position = position_nudge(0.165)) +  # Trim = FALSE to avoid cutting off the tails
  geom_boxplot(outlier.size = 0.5, linewidth = 0.3, width = 0.3, position = position_nudge(-0.165)) +
  geom_point(data = HYSPLIT_summary_Portoroz, aes(x = Season, y = mean_percent, group = Season), color = "black", size = 0.8, position = position_nudge(0.155)) +
  scale_color_manual(values = c("salmon4")) +
  labs(y = expression("%")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, colour = "black", hjust = 0.5, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 7, colour = "black", margin = margin(0.5,0.5,0.5,0.5, "mm")),
        legend.position = "none")

ggsave(filename = paste("Portoroz_HYSPLIT_seasonality",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 40, units = "mm")

# Portorož HYSPLIT supplemental Figure b
ggplot() +
  geom_line(data = dplyr::filter(amplitude_df_mediterranean_Portoroz, IMF == "IMF6" | IMF ==  "Residual"), aes(x = time, y = Value, group = interaction(Variable, Type), linetype = Type), linewidth = 0.2, color = "black") +
  scale_linetype_manual(values = c(2, 1)) +
  scale_x_date(
    limits = as.Date(c("2000-01-01", "2025-01-01")),
    breaks = seq(as.Date("2000-01-01"), as.Date("2025-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("2000-01-01"), as.Date("2025-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) + # %b for month
  facet_wrap(~ IMF, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, colour = "black", angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 7, colour = "black"),
        strip.text = element_text(size = 6, colour = "black", margin = margin(0.5,0.5,0.5,0.5, "mm")),
        panel.grid.minor.y = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.spacing = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(size = 6, colour = "black", margin = margin(t = 0, r = 0, b = 0, l = 0.5)),
        legend.key.spacing.x = unit(0.3, units = "cm"), 
        legend.key.width = unit(0.2, units = "cm"),
        legend.spacing.x =  unit(1.5, units = "cm"),
        legend.title = element_blank()) +
  labs(y = expression("[%]"))

ggsave(filename = paste("Portoroz_decomposition_significant_HYSPLIT",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 85, units = "mm")

# Vienna HYSPLIT supplemental Figure a, upper plot
ggplot() +
  geom_line(data = HYSPLIT_timeseries_Vienna, aes(x = Date, y = percent_in_mediterranean), color = "salmon4", linewidth = 0.3) +
  scale_x_date(
    limits = as.Date(c("1973-01-01", "2024-01-01")),
    breaks = seq(as.Date("1973-01-01"), as.Date("2024-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1973-01-01"), as.Date("2024-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) +
  facet_wrap(~ Variable) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, colour = "black", angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 7, colour = "black", margin = margin(0.5,0.5,0.5,0.5, "mm")),
        legend.position = "none") +
  labs(y = expression("%"))

ggsave(filename = paste("Vienna_HYSPLIT_timetrend",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 40, units = "mm")

# Vienna HYSPLIT supplemental Figure a, bottom plot
ggplot(data = HYSPLIT_timeseries_Vienna, aes(x = Season, y = percent_in_mediterranean, colour = Variable, group = Season), linewidth = 0.3) +
  geom_violin(trim = FALSE, scale = "area", alpha = 0.5, width = 0.3, linewidth = 0.3, position = position_nudge(0.165)) +  # Trim = FALSE to avoid cutting off the tails
  geom_boxplot(outlier.size = 0.5, linewidth = 0.3, width = 0.3, position = position_nudge(-0.165)) +
  geom_point(data = HYSPLIT_summary_Vienna, aes(x = Season, y = mean_percent, group = Season), color = "black", size = 0.8, position = position_nudge(0.155)) +
  scale_color_manual(values = c("salmon4")) +
  labs(y = expression("%")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, colour = "black", hjust = 0.5, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 7, colour = "black", margin = margin(0.5,0.5,0.5,0.5, "mm")),
        legend.position = "none")

ggsave(filename = paste("Vienna_HYSPLIT_seasonality",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 40, units = "mm")

# Vienna HYSPLIT supplemental Figure b
ggplot() +
  geom_line(data = dplyr::filter(amplitude_df_mediterranean_Vienna, IMF == "IMF3" | IMF ==  "Residual"), aes(x = time, y = Value, group = interaction(Variable, Type), linetype = Type), linewidth = 0.2, color = "black") +
  scale_linetype_manual(values = c(2, 1)) +
  scale_x_date(
    limits = as.Date(c("1973-01-01", "2024-01-01")),
    breaks = seq(as.Date("1973-01-01"), as.Date("2024-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1973-01-01"), as.Date("2024-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) + # %b for month
  facet_wrap(~ IMF, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, colour = "black", angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 7, colour = "black"),
        strip.text = element_text(size = 6, colour = "black", margin = margin(0.5,0.5,0.5,0.5, "mm")),
        panel.grid.minor.y = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.spacing = margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(size = 6, colour = "black", margin = margin(t = 0, r = 0, b = 0, l = 0.5)),
        legend.key.spacing.x = unit(0.3, units = "cm"), 
        legend.key.width = unit(0.2, units = "cm"),
        legend.spacing.x =  unit(1.5, units = "cm"),
        legend.title = element_blank()) +
  labs(y = expression("[%]"))

ggsave(filename = paste("Vienna_decomposition_significant_HYSPLIT",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1500, width = 90, height = 85, units = "mm")

###############################################################################
# Plots for comparison of decomposition residual model fits (supplemental figures)
###############################################################################
# Select residuals and shape the data into appropriate format
d18O_res <- eemd_df_list[[1]] %>% 
  dplyr::filter(IMF == "Residual") %>% 
  mutate(time_num = as.numeric(time) - as.numeric(min(time)))
d2H_res <- eemd_df_list[[2]] %>% 
  dplyr::filter(IMF == "Residual") %>% 
  mutate(time_num = as.numeric(time) - as.numeric(min(time)))
d_exc_res <- eemd_df_list[[3]] %>% 
  dplyr::filter(IMF == "Residual") %>% 
  mutate(time_num = as.numeric(time) - as.numeric(min(time)))

# Fit linear and exponential models to the data
lm_lin_d18O <- lm(Value ~ time_num, data = d18O_res)
fit_exp_d18O <- nls(
  Value ~ a * exp(b * time_num) + c,
  data = d18O_res,
  start = list(a = diff(range(d18O_res$Value)),
               b = 1e-4,
               c = min(d18O_res$Value))
)
lm_lin_d2H <- lm(Value ~ time_num, data = d2H_res)
fit_exp_d2H <- nls(
  Value ~ a * exp(b * time_num) + c,
  data = d2H_res,
  start = list(a = diff(range(d2H_res$Value)),
               b = 1e-4,
               c = min(d2H_res$Value))
)
lm_lin_d_exc <- lm(Value ~ time_num, data = d_exc_res)
fit_exp_d_exc <- nls(
  Value ~ a * exp(b * time_num) + c,
  data = d_exc_res,
  start = list(a = diff(range(d_exc_res$Value)),
               b = 1e-4,
               c = min(d_exc_res$Value))
)

d18O_res$pred_linear <- predict(lm_lin_d18O)
d18O_res$pred_exp <- predict(fit_exp_d18O)
d2H_res$pred_linear <- predict(lm_lin_d2H)
d2H_res$pred_exp <- predict(fit_exp_d2H)
d_exc_res$pred_linear <- predict(lm_lin_d_exc)
d_exc_res$pred_exp <- predict(fit_exp_d_exc)

# Compare AIC (model fitness criterion)
AIC(lm_lin_d18O, fit_exp_d18O)
AIC(lm_lin_d2H, fit_exp_d2H)
AIC(lm_lin_d_exc, fit_exp_d_exc)

# Plot, δ18O
ggplot(d18O_res, aes(x = time, y = Value)) +
  geom_line(aes(color = "Decomposition residual"), size = 0.7) +
  geom_line(aes(y = pred_linear, color = "Linear fit"), size = 0.7) +
  geom_line(aes(y = pred_exp, color = "Exponential fit"), size = 0.7) +
  scale_color_manual(values = c("Decomposition residual" = "black", "Linear fit" = "blue", "Exponential fit" = "red")) +
  labs(x = "",
       color = "",
       y = expression("δ"^18*"O [‰]")
  ) +
  scale_y_continuous(labels = ~sub("-", "\u2212", .x)) +
  scale_x_date(
    limits = as.Date(c("1981-01-01", "2024-01-01")),
    breaks = seq(as.Date("1981-01-01"), as.Date("2024-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1981-01-01"), as.Date("2024-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.3, 0.8))

ggsave(filename = paste("Ljubljana_d18O_residual_models",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 700, width = 90, height = 85, units = "mm")

# Plot, δ2H
ggplot(d2H_res, aes(x = time, y = Value)) +
  geom_line(aes(color = "Decomposition residual"), size = 0.7) +
  geom_line(aes(y = pred_linear, color = "Linear fit"), size = 0.7) +
  geom_line(aes(y = pred_exp, color = "Exponential fit"), size = 0.7) +
  scale_color_manual(values = c("Decomposition residual" = "black", "Linear fit" = "blue", "Exponential fit" = "red")) +
  labs(x = "",
       color = "",
       y = expression("δ"^2*"H [‰]")
  ) +
  scale_y_continuous(labels = ~sub("-", "\u2212", .x)) +
  scale_x_date(
    limits = as.Date(c("1981-01-01", "2024-01-01")),
    breaks = seq(as.Date("1981-01-01"), as.Date("2024-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1981-01-01"), as.Date("2024-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.3, 0.8))

ggsave(filename = paste("Ljubljana_d2H_residual_models",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 700, width = 90, height = 85, units = "mm")

# Plot, d-excess
ggplot(d_exc_res, aes(x = time, y = Value)) +
  geom_line(aes(color = "Decomposition residual"), size = 0.7) +
  geom_line(aes(y = pred_linear, color = "Linear fit"), size = 0.7) +
  geom_line(aes(y = pred_exp, color = "Exponential fit"), size = 0.7) +
  scale_color_manual(values = c("Decomposition residual" = "black", "Linear fit" = "blue", "Exponential fit" = "red")) +
  labs(x = "",
       color = "",
       y = expression("d-excess [‰]")
  ) +
  scale_y_continuous(labels = ~sub("-", "\u2212", .x)) +
  scale_x_date(
    limits = as.Date(c("1981-01-01", "2024-01-01")),
    breaks = seq(as.Date("1981-01-01"), as.Date("2024-01-01"), by = "5 year"),
    minor_breaks = seq(as.Date("1981-01-01"), as.Date("2024-01-01"), by = "1 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.3, 0.8))

ggsave(filename = paste("Ljubljana_d_excess_residual_models",".",figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 700, width = 90, height = 85, units = "mm")

###############################################################################
# DECOMPOSITION CORRELATION BETWEEN IMFs, LAG CHECK, etc.
###############################################################################
# LAG CHECK and correlation for IMF 2 MOI
MOI_IMF2 <- dplyr::filter(eemd_df_MOI, IMF == "IMF2" & time > as.Date("1981-04-01"))$Value
d_IMF3 <- dplyr::filter(eemd_df_d, IMF == "IMF3" & Variable == "dn" & year(time) < 2023)$Value # TAKE CARE, MANY FILTERS!
d_IMF4 <- dplyr::filter(eemd_df_d, IMF == "IMF4" & Variable == "dn" & year(time) < 2023)$Value # TAKE CARE, MANY FILTERS!
d_IMF5 <- dplyr::filter(eemd_df_d, IMF == "IMF5" & Variable == "dn" & year(time) < 2023)$Value # TAKE CARE, MANY FILTERS!

cor(MOI_IMF2, d_IMF3, use = "pairwise.complete.obs")
cor.test(MOI_IMF2, d_IMF3, use = "pairwise.complete.obs")
ccf_result <- ccf(MOI_IMF2, d_IMF3, lag.max = 20, plot = TRUE)
max_lag <- ccf_result$lag[which.max(abs(ccf_result$acf))]
print(max_lag)
max_corr <- ccf_result$acf[which.max(abs(ccf_result$acf))]
print(max_corr)

cor(MOI_IMF2, d_IMF4, use = "pairwise.complete.obs")
cor.test(MOI_IMF2, d_IMF4, use = "complete.obs")
ccf_result <- ccf(MOI_IMF2, d_IMF4, lag.max = 20, plot = TRUE)

cor(MOI_IMF2, d_IMF5, use = "pairwise.complete.obs")
cor.test(MOI_IMF2, d_IMF5, use = "complete.obs")
ccf_result <- ccf(MOI_IMF2, d_IMF5, lag.max = 20, plot = TRUE)

# LAG CHECK and correlation for IMF 4 & IMF 7 NAO for LJUBLJANA !!!
NAO_IMF4_Ljubljana <- dplyr::filter(eemd_df_NAO, IMF == "IMF4" & time >= as.Date("1981-05-01"))$Value # Date filter to match Ljubljana!
NAO_IMF7_Ljubljana <- dplyr::filter(eemd_df_NAO, IMF == "IMF7" & time >= as.Date("1981-05-01"))$Value # Date filter to match Ljubljana!
d_IMF4_Ljubljana <- dplyr::filter(eemd_df_d, IMF == "IMF4" & Variable == "dn"& year(time) < 2024)$Value # Date filter to match NAO (only till 2023)!
d_IMF5_Ljubljana <- dplyr::filter(eemd_df_d, IMF == "IMF5" & Variable == "dn"& year(time) < 2024)$Value # Date filter to match NAO (only till 2023)!
NAO_IMF4_Vienna <- dplyr::filter(eemd_df_NAO, IMF == "IMF4" & time >= as.Date("1973-01-01") & time <= as.Date("2023-11-01"))$Value # Date filter to match Ljubljana!
NAO_IMF7_Vienna <- dplyr::filter(eemd_df_NAO, IMF == "IMF7" & time >= as.Date("1973-01-01") & time <= as.Date("2023-11-01"))$Value # Date filter to match Ljubljana!
d_IMF4_Vienna <- dplyr::filter(eemd_df_list_Vienna[[3]], IMF == "IMF4")$Value 
d_IMF5_Vienna <- dplyr::filter(eemd_df_list_Vienna[[3]], IMF == "IMF5")$Value


cor(NAO_IMF4_Ljubljana, d_IMF4_Ljubljana, use = "complete.obs")
cor.test(NAO_IMF4_Ljubljana, d_IMF4_Ljubljana, use = "complete.obs")
result_ccf <- ccf(NAO_IMF4_Ljubljana, d_IMF4_Ljubljana, lag.max = 20, plot = TRUE)
max_lag <- result_ccf$lag[which.max(abs(result_ccf$acf))]
print(max_lag)
max_corr <- result_ccf$acf[which.max(abs(result_ccf$acf))]
print(max_corr)

cor(NAO_IMF4_Ljubljana, d_IMF5_Ljubljana, use = "complete.obs")
cor.test(NAO_IMF4_Ljubljana, d_IMF5_Ljubljana, use = "complete.obs")
result_ccf <-ccf(NAO_IMF4_Ljubljana, d_IMF5_Ljubljana, lag.max = 20, plot = TRUE)

cor(NAO_IMF7_Ljubljana, d_IMF4_Ljubljana, use = "complete.obs")
cor.test(NAO_IMF7_Ljubljana, d_IMF4_Ljubljana, use = "complete.obs")
result_ccf <-ccf(NAO_IMF7_Ljubljana, d_IMF4_Ljubljana, lag.max = 20, plot = TRUE)

cor(NAO_IMF7_Ljubljana, d_IMF5_Ljubljana, use = "complete.obs")
cor.test(NAO_IMF7_Ljubljana, d_IMF5_Ljubljana, use = "complete.obs")
result_ccf <-ccf(NAO_IMF7_Ljubljana, d_IMF5_Ljubljana, lag.max = 20, plot = TRUE)

cor(NAO_IMF4_Vienna, d_IMF4_Vienna, use = "complete.obs")
cor.test(NAO_IMF4_Vienna, d_IMF4_Vienna, use = "complete.obs")
result_ccf <- ccf(NAO_IMF4_Vienna, d_IMF4_Vienna, lag.max = 20, plot = TRUE)
max_lag <- result_ccf$lag[which.max(abs(result_ccf$acf))]
print(max_lag)
max_corr <- result_ccf$acf[which.max(abs(result_ccf$acf))]
print(max_corr)

cor(NAO_IMF4_Vienna, d_IMF5_Vienna, use = "complete.obs")
cor.test(NAO_IMF4_Vienna, d_IMF5_Vienna, use = "complete.obs")
result_ccf <- ccf(NAO_IMF4_Vienna, d_IMF5_Vienna, lag.max = 20, plot = TRUE)

cor(NAO_IMF7_Vienna, d_IMF4_Vienna, use = "complete.obs")
cor.test(NAO_IMF7_Vienna, d_IMF4_Vienna, use = "complete.obs")
result_ccf <- ccf(NAO_IMF7_Vienna, d_IMF4_Vienna, lag.max = 20, plot = TRUE)

cor(NAO_IMF7_Vienna, d_IMF5_Vienna, use = "complete.obs")
cor.test(NAO_IMF7_Vienna, d_IMF5_Vienna, use = "complete.obs")
result_ccf <- ccf(NAO_IMF7_Vienna, d_IMF5_Vienna, lag.max = 20, plot = TRUE)

# LAG CHECK and correlation for d-excess and HYSPLIT, Ljubljana
time_vector <- unique(eemd_df_mediterranean$time)
med_IMF3 <- dplyr::filter(eemd_df_mediterranean, IMF == "IMF3")$Value
med_IMF6 <- dplyr::filter(eemd_df_mediterranean, IMF == "IMF6")$Value
d_IMF3 <- dplyr::filter(eemd_df_d, IMF == "IMF3" & Variable == "dn")$Value
d_IMF4 <- dplyr::filter(eemd_df_d, IMF == "IMF4" & Variable == "dn")$Value
d_IMF5 <- dplyr::filter(eemd_df_d, IMF == "IMF5" & Variable == "dn")$Value

cor(med_IMF3, d_IMF3, use = "complete.obs")
cor.test(med_IMF3, d_IMF3, use = "complete.obs")
ccf(med_IMF3, d_IMF3, lag.max = 20, plot = TRUE)

cor(med_IMF3, d_IMF4, use = "complete.obs")
cor.test(med_IMF3, d_IMF4, use = "complete.obs")
ccf(med_IMF3, d_IMF4, lag.max = 20, plot = TRUE)

cor(med_IMF3, d_IMF5, use = "complete.obs")
cor.test(med_IMF3, d_IMF5, use = "complete.obs")
ccf(med_IMF3, d_IMF5, lag.max = 20, plot = TRUE)

cor(med_IMF6, d_IMF3, use = "complete.obs")
cor.test(med_IMF6, d_IMF3, use = "complete.obs")
ccf(med_IMF6, d_IMF3, lag.max = 20, plot = TRUE)

cor(med_IMF6, d_IMF4, use = "complete.obs")
cor.test(med_IMF6, d_IMF4, use = "complete.obs")
ccf(med_IMF6, d_IMF4, lag.max = 20, plot = TRUE)

cor(med_IMF6, d_IMF5, use = "complete.obs")
cor.test(med_IMF6, d_IMF5, use = "complete.obs") # THIS ONE IS SIGNIFICANT
ccf(med_IMF6, d_IMF5, lag.max = 40, plot = TRUE)

# LAG CHECK and correlation for d and HYSPLIT, Vienna
time_vector <- unique(eemd_df_list_Vienna[[6]]$time)
med_IMF3_Vienna <- dplyr::filter(eemd_df_list_Vienna[[6]], IMF == "IMF3")$Value
d_IMF3_Vienna <- dplyr::filter(eemd_df_list_Vienna[[3]], IMF == "IMF3")$Value
d_IMF4_Vienna <- dplyr::filter(eemd_df_list_Vienna[[3]], IMF == "IMF4")$Value
d_IMF5_Vienna <- dplyr::filter(eemd_df_list_Vienna[[3]], IMF == "IMF5")$Value

cor(med_IMF3_Vienna, d_IMF3_Vienna, use = "complete.obs")
cor.test(med_IMF3_Vienna, d_IMF3_Vienna, use = "complete.obs")
result_ccf <- ccf(med_IMF3_Vienna, d_IMF3_Vienna, lag.max = 20, plot = TRUE)
max_lag <- result_ccf$lag[which.max(abs(result_ccf$acf))]
print(max_lag)
max_corr <- result_ccf$acf[which.max(abs(result_ccf$acf))]
print(max_corr)

cor(med_IMF3_Vienna, d_IMF4_Vienna, use = "complete.obs")
cor.test(med_IMF3_Vienna, d_IMF4_Vienna, use = "complete.obs")
result_ccf <- ccf(med_IMF3_Vienna, d_IMF4_Vienna, lag.max = 20, plot = TRUE)

cor(med_IMF3_Vienna, d_IMF5_Vienna, use = "complete.obs")
cor.test(med_IMF3_Vienna, d_IMF5_Vienna, use = "complete.obs")
result_ccf <- ccf(med_IMF3_Vienna, d_IMF5_Vienna, lag.max = 20, plot = TRUE)

# LAG CHECK and correlation for IMF 2 Temperature
eemd_df_T <- eemd_df_list[["T"]]
time_vector <- unique(eemd_df_T$time)
T_IMF2 <- dplyr::filter(eemd_df_T, IMF == "IMF2")$Value
d_IMF3 <- dplyr::filter(eemd_df_d, IMF == "IMF3" & Variable == "dn")$Value
d_IMF4 <- dplyr::filter(eemd_df_d, IMF == "IMF4" & Variable == "dn")$Value
d_IMF5 <- dplyr::filter(eemd_df_d, IMF == "IMF5" & Variable == "dn")$Value

cor(T_IMF2, d_IMF3, use = "complete.obs")
cor.test(T_IMF2, d_IMF3, use = "complete.obs")
result_ccf <- ccf(T_IMF2, d_IMF3, lag.max = 20, plot = TRUE)
max_lag <- result_ccf$lag[which.max(abs(result_ccf$acf))]
print(max_lag)
max_corr <- result_ccf$acf[which.max(abs(result_ccf$acf))]
print(max_corr)

shifted_d_IMF3 <- c(rep(NA, 2), d_IMF3) # ADJUSTED BY LAG = 3 MONTHS
shifted_T_IMF2 <- c(T_IMF2, rep(NA, 2))
cor(shifted_T_IMF2, shifted_d_IMF3, use = "pairwise.complete.obs")

cor(T_IMF2, d_IMF4, use = "complete.obs")
cor.test(T_IMF2, d_IMF4, use = "complete.obs")
ccf(T_IMF2, d_IMF4, lag.max = 20, plot = TRUE)

cor(T_IMF2, d_IMF5, use = "complete.obs")
cor.test(T_IMF2, d_IMF5, use = "complete.obs")
ccf(T_IMF2, d_IMF5, lag.max = 20, plot = TRUE)

# PORTOROŽ!!! LAG CHECK and correlation for IMF 2 Temperature
eemd_df_d18O_Portoroz <- eemd_df_list_Portoroz[[1]]
eemd_df_d2H_Portoroz <- eemd_df_list_Portoroz[[2]]
eemd_df_dP_Portoroz <- eemd_df_list_Portoroz[[5]]
eemd_df_T_Portoroz <- eemd_df_list_Portoroz[[4]]
eemd_df_d_Portoroz <- eemd_df_list_Portoroz[[3]]
time_vector <- unique(eemd_df_T_Portoroz$time)

d18O_IMF2 <- dplyr::filter(eemd_df_d18O_Portoroz, IMF == "IMF2")$Value
d18O_res<- dplyr::filter(eemd_df_d18O_Portoroz, IMF == "Residual")$Value
d2H_IMF2 <- dplyr::filter(eemd_df_d2H_Portoroz, IMF == "IMF2")$Value
d2H_res<- dplyr::filter(eemd_df_d2H_Portoroz, IMF == "Residual")$Value
T_IMF2 <- dplyr::filter(eemd_df_T_Portoroz, IMF == "IMF2")$Value
T_res<- dplyr::filter(eemd_df_T_Portoroz, IMF == "Residual")$Value
d_IMF3 <- dplyr::filter(eemd_df_d_Portoroz, IMF == "IMF3" & Variable == "dn")$Value
d_res<- dplyr::filter(eemd_df_d_Portoroz, IMF == "Residual")$Value
P_res<- dplyr::filter(eemd_df_P_Portoroz, IMF == "Residual")$Value

cor(T_IMF2, d18O_IMF2, use = "complete.obs")
cor(T_IMF2, d2H_IMF2, use = "complete.obs")
cor(T_IMF2, d_IMF3, use = "complete.obs")
cor(T_res, d18O_res, use = "complete.obs")
cor(P_res, d18O_res, use = "complete.obs")
cor(T_res, d2H_res, use = "complete.obs")
cor(P_res, d2H_res, use = "complete.obs")
cor(T_res, d_res, use = "complete.obs")
cor(P_res, d_res, use = "complete.obs")

cor.test(T_IMF2, d_IMF3, use = "complete.obs")
result_ccf <- ccf(T_IMF2, d_IMF3, lag.max = 20, plot = TRUE)
max_lag <- result_ccf$lag[which.max(abs(result_ccf$acf))]
print(max_lag)
max_corr <- result_ccf$acf[which.max(abs(result_ccf$acf))]
print(max_corr)

shifted_d_IMF3 <- c(rep(NA, 2), d_IMF3) # ADJUSTED BY LAG = 3 MONTHS
shifted_T_IMF2 <- c(T_IMF2, rep(NA, 2))
cor(shifted_T_IMF2, shifted_d_IMF3, use = "pairwise.complete.obs")

# VIENNA!!! LAG CHECK and correlation for IMF 2 Temperature
eemd_df_d18O_Vienna <- eemd_df_list_Vienna[[1]]
eemd_df_d2H_Vienna <- eemd_df_list_Vienna[[2]]
eemd_df_P_Vienna <- eemd_df_list_Vienna[[5]]
eemd_df_T_Vienna <- eemd_df_list_Vienna[[4]]
eemd_df_d_Vienna <- eemd_df_list_Vienna[[3]]
time_vector <- unique(eemd_df_T_Vienna$time)

d18O_IMF_2_3 <- dplyr::filter(eemd_df_d18O_Vienna, IMF == "IMF_2_3")$Value
d18O_res<- dplyr::filter(eemd_df_d18O_Vienna, IMF == "Residual")$Value
d2H_IMF_2_3 <- dplyr::filter(eemd_df_d2H_Vienna, IMF == "IMF_2_3")$Value
d2H_res<- dplyr::filter(eemd_df_d2H_Vienna, IMF == "Residual")$Value
T_IMF2 <- dplyr::filter(eemd_df_T_Vienna, IMF == "IMF2")$Value
T_res<- dplyr::filter(eemd_df_T_Vienna, IMF == "Residual")$Value
d_IMF3 <- dplyr::filter(eemd_df_d_Vienna, IMF == "IMF3" & Variable == "dn")$Value
d_res<- dplyr::filter(eemd_df_d_Vienna, IMF == "Residual")$Value
P_IMF3 <- dplyr::filter(eemd_df_P_Vienna, IMF == "IMF3")$Value
P_res <- dplyr::filter(eemd_df_P_Vienna, IMF == "Residual")$Value

cor(T_IMF2, d18O_IMF_2_3, use = "complete.obs")
cor(P_IMF3, d18O_IMF_2_3, use = "complete.obs")
cor(T_IMF2, d2H_IMF_2_3, use = "complete.obs")
cor(P_IMF3, d2H_IMF_2_3, use = "complete.obs")
cor(T_IMF2, d_IMF3, use = "complete.obs")
cor(P_IMF3, d_IMF3, use = "complete.obs")
cor(T_res, d18O_res, use = "complete.obs")
cor(P_res, d18O_res, use = "complete.obs")
cor(T_res, d2H_res, use = "complete.obs")
cor(P_res, d2H_res, use = "complete.obs")
cor(T_res, d_res, use = "complete.obs")
cor(P_res, d_res, use = "complete.obs")

cor.test(T_IMF2, d_IMF3, use = "complete.obs")
result_ccf <- ccf(T_IMF2, d_IMF3, lag.max = 20, plot = TRUE)
max_lag <- result_ccf$lag[which.max(abs(result_ccf$acf))]
print(max_lag)
max_corr <- result_ccf$acf[which.max(abs(result_ccf$acf))]
print(max_corr)

shifted_d_IMF3 <- c(rep(NA, 3), d_IMF3) # ADJUSTED BY LAG = 3 MONTHS
shifted_T_IMF2 <- c(T_IMF2, rep(NA, 3))
cor(shifted_T_IMF2, shifted_d_IMF3, use = "pairwise.complete.obs")

###############################################################################
# OTHER STAT TESTS
###############################################################################
# Basic time series statistical tests
stat_data <- Ljubljana_data %>%
  dplyr::filter(Year > 1990) %>%
  dplyr::filter(!is.na(rolling_δ2H) & !is.na(rolling_δ18O) & !is.na(rolling_d))

mk.test(stat_data$rolling_δ2H, alternative = "greater")
sens.slope(stat_data$rolling_δ2H)
mk.test(stat_data$rolling_δ18O, alternative = "greater")
sens.slope(stat_data$rolling_δ18O)
mk.test(stat_data$rolling_d, alternative = "greater")
sens.slope(stat_data$rolling_d)

stat_data_weighted <- Ljubljana_data %>%
  dplyr::filter(Year > 1990) %>%
  dplyr::filter(!is.na(weighted_rolling_δ2H) & !is.na(weighted_rolling_δ18O) & !is.na(weighted_rolling_d))

mk.test(stat_data_weighted$weighted_rolling_δ2H, alternative = "greater")
sens.slope(stat_data_weighted$weighted_rolling_δ2H)
mk.test(stat_data_weighted$weighted_rolling_δ18O, alternative = "greater")
sens.slope(stat_data_weighted$weighted_rolling_δ18O)
mk.test(stat_data_weighted$weighted_rolling_d, alternative = "greater")
sens.slope(stat_data_weighted$weighted_rolling_d)

# For months
kruskal.test(data = Ljubljana_data, δ2H ~ Month)
d2H <- as.data.frame(pairwise.wilcox.test(x = Ljubljana_data$δ2H, g = Ljubljana_data$Month)$p.value) %>%
  mutate(Variable = "d2H", .before = 1)

kruskal.test(data = Ljubljana_data, δ18O ~ Month)
d18O <- as.data.frame(pairwise.wilcox.test(x = Ljubljana_data$δ18O, g = Ljubljana_data$Month)$p.value) %>%
  mutate(Variable = "d18O", .before = 1)

kruskal.test(data = Ljubljana_data, d ~ Month)
d_excess <- as.data.frame(pairwise.wilcox.test(x = Ljubljana_data$d, g = Ljubljana_data$Month)$p.value) %>%
  mutate(Variable = "d-excess", .before = 1)

stat_months <- rbind(d2H, d18O, d_excess)

stat_months <- rbind(d2H, d18O, d_excess) %>%
  mutate(across(everything(), ~ case_when(
    is.numeric(.) & .< 0.001 ~ "***",
    is.numeric(.) & . < 0.01 & . > 0.001 ~ "**",
    is.numeric(.) & . < 0.05 & . > 0.01 ~ "*",
    is.numeric(.) & . > 0.05 ~ "ns",
    is.na(.) ~ "",
    TRUE ~ "")))

# write.csv(stat_months, file = "stat_months.csv", row.names = FALSE)

# For seasons
kruskal.test(data = Ljubljana_data, δ2H ~ Season)
pairwise.wilcox.test(x = Ljubljana_data$δ2H, g = Ljubljana_data$Season)

kruskal.test(data = Ljubljana_data, δ18O ~ Season)
pairwise.wilcox.test(x = Ljubljana_data$δ18O, g = Ljubljana_data$Season)

kruskal.test(data = Ljubljana_data, d ~ Season)
pairwise.wilcox.test(x = Ljubljana_data$d, g = Ljubljana_data$Season)

kruskal.test(data = Ljubljana_data, `3H (TU)` ~ Season)
pairwise.wilcox.test(x = Ljubljana_data$`3H (TU)`, g = Ljubljana_data$Season)

# For yearly slopes & intercepts of of MA & RMA regressions
mk.test(Ljubljana_data_year_slopes$slope_MA, alternative = "less")
sens.slope(Ljubljana_data_year_slopes$slope_MA)
mk.test(Ljubljana_data_year_slopes$slope_MA_weighted, alternative = "less")
sens.slope(Ljubljana_data_year_slopes$slope_MA_weighted)
mk.test(Ljubljana_data_year_slopes$slope_RMA, alternative = "less")
sens.slope(Ljubljana_data_year_slopes$slope_RMA)
mk.test(Ljubljana_data_year_slopes$slope_RMA_weighted, alternative = "less")
sens.slope(Ljubljana_data_year_slopes$slope_RMA_weighted)

mk.test(Ljubljana_data_year_slopes$intercept_MA, alternative = "less")
sens.slope(Ljubljana_data_year_slopes$intercept_MA)
mk.test(Ljubljana_data_year_slopes$intercept_MA_weighted, alternative = "less")
sens.slope(Ljubljana_data_year_slopes$intercept_MA_weighted)
mk.test(Ljubljana_data_year_slopes$intercept_RMA, alternative = "less")
sens.slope(Ljubljana_data_year_slopes$intercept_RMA)
mk.test(Ljubljana_data_year_slopes$intercept_RMA_weighted, alternative = "less")
sens.slope(Ljubljana_data_year_slopes$intercept_RMA_weighted)