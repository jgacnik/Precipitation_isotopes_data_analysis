###############################################################################
# FIGURES FOR THE SZGG CONTRIBUTION (simplified figures, Ljubljana station)
# Based on Gačnik J., Žagar K., Hatvani I. G., Kern Z., and Vreča P. (2025): Climate change reflected in 40-year isotopic composition trends of precipitation in Slovenia, Environmental Research
# CODE AUTHOR: JAN GAČNIK, NOVEMBER 2026
# R version: 4.4.2
###############################################################################
# DO THIS STEP ONLY THE FIRST TIME! Installation of needed packages
install.packages(c("mblm", "Rlibeemd", "trend", "patchwork", "tidyverse", "readxl"))

# CHANGE THE FOLLOWING PARAMETERS:
path_scripts <- "ENTER/PATH/TO/DIRECTORY/Precipitation_isotopes_data_analysis/scripts" # Change to your directory path where the R scripts are located, use "/" and not "\"
path_data <- "ENTER/PATH/TO/DIRECTORY/Precipitation_isotopes_data_analysis/data" # Change to your directory path where the data is located, use "/" and not "\"
path_figures <- "ENTER/PATH/TO/DIRECTORY/Precipitation_isotopes_data_analysis/figures" # Change to your directory path where you want to save figures, use "/" and not "\"

###############################################################################
# CODE BELOW DOES NOT NEED ALTERING
###############################################################################
Ljubljana_filename <- "Ljubljana_1981-2024.xlsx"
Ljubljana_filename_metdata <- "Ljubljana-daily-meteodata.xlsx"
HYSPLIT_filename <- "HYSPLIT_timeseries.xlsx"
figure_type <- "png"

# Load required packages for the current session
setwd(path_scripts)
lapply(c("mblm", "Rlibeemd", "trend", "patchwork", "tidyverse", "readxl"), require, character.only = TRUE)

# Calling needed functions
source("Functions_MA&RMA.R")

# Slovenian month labels and figure width (SZGG text width is 21 - 3 - 2.5 = 15.5 cm)
months_slo <- c("jan", "feb", "mar", "apr", "maj", "jun", "jul", "avg", "sep", "okt", "nov", "dec")
width_full <- 155

###############################################################################
# Data reading, adjusting, and filtering
###############################################################################
Ljubljana_data <- read_xlsx(sheet = "Data Ljubljana", path = file.path(path_data, Ljubljana_filename), trim_ws = TRUE) %>%
  dplyr::select(Year, Month, P, T, δ18O, δ2H, d, δ18On) %>%
  mutate(Period = case_when(
    Year >= 1981 & Year <= 1991 ~ "1981–1991",
    Year >= 1992 & Year <= 2002 ~ "1992–2002",
    Year >= 2003 & Year <= 2013 ~ "2003–2013",
    Year >= 2014 ~ "2014–2024")) %>%
  mutate(Date = make_date(Year, Month))

# Precipitation used for weighting only when both isotopes are available
Ljubljana_data_filtered <- dplyr::filter(Ljubljana_data, is.na(δ2H) == FALSE & is.na(δ18O) == FALSE)

###############################################################################
# SZGG FIGURE 1 (raw monthly data)
###############################################################################
Ljubljana_data_raw <- Ljubljana_data %>%
  pivot_longer(cols = c(δ18O, δ2H, d, T, P), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable_color = Variable) %>%
  mutate(Variable = case_when(
    Variable == "δ18O" ~ "italic(delta)^18*O~'[‰]'",
    Variable == "δ2H" ~ "italic(delta)^2*H~'[‰]'",
    Variable == "d" ~ "atop('devterijev', 'presežek [‰]')",
    Variable == "T" ~ "'T [°C]'",
    Variable == "P" ~ "'P [mm]'")) %>%
  mutate(Variable = factor(Variable, levels = c("italic(delta)^18*O~'[‰]'", "italic(delta)^2*H~'[‰]'", "atop('devterijev', 'presežek [‰]')", "'T [°C]'", "'P [mm]'")),
         Variable_color = factor(Variable_color, levels = c("δ18O", "δ2H", "d", "T", "P")))

ggplot(data = Ljubljana_data_raw, aes(x = Date, y = Value, color = Variable_color, fill = Variable_color)) +
  geom_line(data = dplyr::filter(Ljubljana_data_raw, Variable_color != "P"), linewidth = 0.3) +
  geom_col(data = dplyr::filter(Ljubljana_data_raw, Variable_color == "P"), width = 25, color = NA) +
  scale_color_manual(values = c("grey30", "grey30", "grey30", "#D7301F", "#2166AC")) +
  scale_fill_manual(values = c("grey30", "grey30", "grey30", "#D7301F", "#2166AC")) +
  scale_x_date(
    limits = as.Date(c("1981-01-01", "2025-01-01")),
    breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "5 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) +
  scale_y_continuous(labels = ~sub("-", "−", sub(".", ",", .x, fixed = TRUE))) +
  facet_wrap(~ Variable, ncol = 1, scales = "free_y", labeller = label_parsed, strip.position = "left") +
  theme_bw(base_family = "serif") +
  theme(axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_blank(),
        strip.text = element_text(size = 9, colour = "black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(filename = paste("SZGG_raw_data", ".", figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1000, width = width_full, height = 140, units = "mm")

###############################################################################
# SZGG FIGURE 2 (simplified EEMD of δ18O)
###############################################################################
# Same EEMD settings as in the article (without the Monte Carlo significance test)
time_vector_Ljubljana <- Ljubljana_data$Date[5:length(Ljubljana_data$Date)]
eemd_d18O <- as.data.frame(eemd(Ljubljana_data$δ18On[-c(1:4)], noise_strength = 0.2, ensemble_size = 100, S_number = 4))
colnames(eemd_d18O) <- gsub(" ", "", colnames(eemd_d18O))
n_imfs <- ncol(eemd_d18O) - 1

eemd_df_SZGG <- data.frame(time = time_vector_Ljubljana,
                           "izmerjeno" = Ljubljana_data$δ18On[-c(1:4)],
                           "letni cikel" = eemd_d18O$IMF2 + eemd_d18O$IMF3,
                           "večletna nihanja" = rowSums(eemd_d18O[, paste0("IMF", 4:n_imfs)]),
                           "dolgoročni trend" = eemd_d18O$Residual,
                           "kratkoročni šum" = eemd_d18O$IMF1,
                           check.names = FALSE) %>%
  pivot_longer(cols = -time, names_to = "IMF", values_to = "Value") %>%
  mutate(IMF = factor(IMF, levels = c("izmerjeno", "letni cikel", "večletna nihanja", "dolgoročni trend", "kratkoročni šum")))

# Change of the residual trend over the whole period (for the figure caption)
eemd_d18O$Residual[nrow(eemd_d18O)] - eemd_d18O$Residual[1]

ggplot(data = eemd_df_SZGG, aes(x = time, y = Value, color = IMF == "izmerjeno")) +
  geom_line(linewidth = 0.25) +
  scale_color_manual(values = c("dodgerblue3", "grey50")) +
  scale_x_date(
    limits = as.Date(c("1981-01-01", "2025-01-01")),
    breaks = seq(as.Date("1981-01-01"), as.Date("2025-01-01"), by = "5 year"),
    date_labels = "%Y",
    expand = expansion(mult = 0.01, add = 0)) +
  scale_y_continuous(labels = ~sub("-", "−", sub(".", ",", .x, fixed = TRUE))) +
  facet_wrap(~ IMF, ncol = 1, scales = "free_y") +
  theme_bw(base_family = "serif") +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.y = element_text(size = 9, colour = "black"),
        strip.text = element_text(size = 9, colour = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "mm")),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  labs(y = expression(italic(delta)^18*O~"[‰]"))

ggsave(filename = paste("SZGG_EEMD_d18O", ".", figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1000, width = 85, height = 120, units = "mm")

###############################################################################
# SZGG FIGURE 3 (annual means with Sen's slope and Mann-Kendall test)
###############################################################################
# Annual means of monthly values, same as Ljubljana_yearly_Sen_and_MK in the main script
Ljubljana_yearly <- Ljubljana_data %>%
  pivot_longer(cols = c(δ18O, d, T, P), names_to = "Variable", values_to = "Value") %>%
  dplyr::filter(!is.na(Value)) %>%
  group_by(Year, Variable) %>%
  summarise(Value = mean(Value), .groups = "drop") %>%
  mutate(Variable = case_when(
    Variable == "δ18O" ~ "italic(delta)^18*O~'[‰]'",
    Variable == "d" ~ "'devterijev presežek [‰]'",
    Variable == "T" ~ "'T [°C]'",
    Variable == "P" ~ "'P [mm]'")) %>%
  mutate(Variable = factor(Variable, levels = c("italic(delta)^18*O~'[‰]'", "'devterijev presežek [‰]'", "'T [°C]'", "'P [mm]'")))

Ljubljana_yearly_Sen_and_MK <- Ljubljana_yearly %>%
  group_by(Variable) %>%
  group_modify(~{
    model <- mblm(Value ~ Year, data = .x, repeated = FALSE)
    tibble(Slope = unname(coef(model)[2]),
           Intercept = unname(coef(model)[1]),
           p_two_sided = mk.test(.x$Value, alternative = "two.sided")$p.value)
  }) %>%
  ungroup() %>%
  mutate(label = paste0(sprintf("%+.2f", 10 * Slope), " / 10 let, p ",
                        ifelse(p_two_sided < 0.001, "< 0.001", paste("=", sprintf("%.3f", p_two_sided)))),
         label = gsub(".", ",", label, fixed = TRUE))

ggplot(data = Ljubljana_yearly, aes(x = Year, y = Value)) +
  geom_line(linewidth = 0.25, color = "grey70") +
  geom_point(size = 0.6) +
  geom_abline(data = Ljubljana_yearly_Sen_and_MK, aes(intercept = Intercept, slope = Slope), linewidth = 0.4, color = "dodgerblue3") +
  geom_text(data = Ljubljana_yearly_Sen_and_MK, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.4, size = 3, family = "serif", color = "dodgerblue3") +
  scale_x_continuous(breaks = seq(1981, 2025, 10)) +
  scale_y_continuous(labels = ~sub("-", "−", sub(".", ",", .x, fixed = TRUE)), expand = expansion(mult = c(0.05, 0.2))) +
  facet_wrap(~ Variable, ncol = 2, scales = "free_y", labeller = label_parsed) +
  theme_bw(base_family = "serif") +
  theme(axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_blank(),
        strip.text = element_text(size = 9, colour = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "mm")),
        panel.grid.minor = element_blank())

ggsave(filename = paste("SZGG_yearly_Sen", ".", figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1000, width = width_full, height = 90, units = "mm")

###############################################################################
# SZGG FIGURE 4 (comparison of periods: a) seasonality, b) LMWL, c) δ18O vs. temperature)
###############################################################################
# a) Precipitation-weighted monthly means per period
Ljubljana_data_period_month <- Ljubljana_data_filtered %>%
  group_by(Period, Month) %>%
  summarise(weighted_mean_δ18O = weighted.mean(δ18O, P),
            weighted_mean_d = weighted.mean(d, P), .groups = "drop") %>%
  pivot_longer(cols = c(weighted_mean_δ18O, weighted_mean_d), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = case_when(
    Variable == "weighted_mean_δ18O" ~ "italic(delta)^18*O~'[‰]'",
    Variable == "weighted_mean_d" ~ "'devterijev presežek [‰]'")) %>%
  mutate(Variable = factor(Variable, levels = c("italic(delta)^18*O~'[‰]'", "'devterijev presežek [‰]'")))

plot_seasonality <- ggplot(data = Ljubljana_data_period_month, aes(x = Month, y = Value, color = Period)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 0.6) +
  scale_color_viridis_d(end = 0.9, direction = -1) +
  guides(color = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = 1:12, labels = months_slo) +
  scale_y_continuous(labels = ~sub("-", "−", sub(".", ",", .x, fixed = TRUE))) +
  facet_wrap(~ Variable, ncol = 1, scales = "free_y", labeller = label_parsed) +
  theme_bw(base_family = "serif") +
  theme(axis.text = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(size = 9, colour = "black", margin = margin(0.3, 0.3, 0.3, 0.3, "mm")),
        panel.grid.minor = element_blank()) +
  labs(title = "a)")

# b) LMWL per period, precipitation-weighted RMA (equations written in the plot)
LMWL_periods <- Ljubljana_data_filtered %>%
  group_by(Period) %>%
  summarise(slope_RMA_weighted = RMA_weighted(δ18O, δ2H, P)$slopes$RMA_weighted,
            intercept_RMA_weighted = RMA_weighted(δ18O, δ2H, P)$intercepts$RMA_weighted) %>%
  mutate(label = gsub(".", ",", sprintf("δ²H = %.2f δ¹⁸O + %.1f", slope_RMA_weighted, intercept_RMA_weighted), fixed = TRUE),
         y = -10 - 11 * (row_number() - 1))

plot_LMWL <- ggplot(data = Ljubljana_data_filtered, aes(x = δ18O, y = δ2H, color = Period)) +
  geom_abline(intercept = 10, slope = 8, linetype = 2, linewidth = 0.3, color = "grey40") + # GMWL
  geom_point(size = 0.5, alpha = 0.5) +
  geom_abline(data = LMWL_periods, aes(intercept = intercept_RMA_weighted, slope = slope_RMA_weighted, color = Period), linewidth = 0.5) +
  geom_text(data = LMWL_periods, aes(x = -20.5, y = y, label = label), hjust = 0, size = 2.6, family = "serif") +
  scale_color_viridis_d(end = 0.9, direction = -1, guide = "none") +
  scale_x_continuous(labels = ~sub("-", "−", .x)) +
  scale_y_continuous(labels = ~sub("-", "−", .x)) +
  theme_bw(base_family = "serif") +
  theme(axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_text(size = 9, colour = "black"),
        panel.grid.minor = element_blank()) +
  labs(title = "b)", x = expression(italic(delta)^18*O~"[‰]"), y = expression(italic(delta)^2*H~"[‰]"))

# c) δ18O vs. temperature per period
plot_d18O_T <- ggplot(data = Ljubljana_data_filtered, aes(x = T, y = δ18O, color = Period)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
  scale_color_viridis_d(end = 0.9, direction = -1, guide = "none") +
  scale_x_continuous(labels = ~sub("-", "−", .x)) +
  scale_y_continuous(labels = ~sub("-", "−", .x)) +
  theme_bw(base_family = "serif") +
  theme(axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_text(size = 9, colour = "black"),
        panel.grid.minor = element_blank()) +
  labs(title = "c)", x = "T [°C]", y = expression(italic(delta)^18*O~"[‰]"))

# Combined figure with one legend for all panels
(plot_seasonality | (free(plot_LMWL, side = "t") / plot_d18O_T)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8, family = "serif"),
        legend.key.spacing.x = unit(0.4, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        plot.title = element_text(size = 14, face = "bold", family = "serif", margin = margin(t = 0, r = 0, b = 5, l = 0)),
        plot.title.position = "panel")

ggsave(filename = paste("SZGG_periods", ".", figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1000, width = width_full, height = 140, units = "mm")

###############################################################################
# SZGG TABLE 1 (mean values per period and for the whole record)
###############################################################################
# Isotopes are precipitation-weighted (only months with both isotopes), T and P use all months
Ljubljana_table_isotopes <- Ljubljana_data_filtered %>%
  rbind(mutate(Ljubljana_data_filtered, Period = "1981–2024")) %>%
  group_by(Period) %>%
  summarise("število mesečnih vzorcev" = sprintf("%d", n()),
            "δ18O [‰]" = sprintf("%.2f", weighted.mean(δ18O, P)),
            "δ2H [‰]" = sprintf("%.1f", weighted.mean(δ2H, P)),
            "devterijev presežek [‰]" = sprintf("%.1f", weighted.mean(d, P)),
            "naklon LMWL" = sprintf("%.2f", RMA_weighted(δ18O, δ2H, P)$slopes$RMA_weighted),
            "presečišče LMWL [‰]" = sprintf("%.1f", RMA_weighted(δ18O, δ2H, P)$intercepts$RMA_weighted))

Ljubljana_table_climate <- Ljubljana_data %>%
  rbind(mutate(Ljubljana_data, Period = "1981–2024")) %>%
  group_by(Period) %>%
  summarise("T [°C]" = sprintf("%.1f", mean(T, na.rm = TRUE)),
            "P [mm/leto]" = sprintf("%.0f", sum(P, na.rm = TRUE) / n_distinct(Year)))

# Values as text with decimal comma and minus sign, whole record as the last column
Ljubljana_table_periods <- Ljubljana_table_isotopes %>%
  left_join(Ljubljana_table_climate, by = "Period") %>%
  pivot_longer(cols = -Period, names_to = "Spremenljivka", values_to = "Value") %>%
  mutate(Value = sub("-", "−", sub(".", ",", Value, fixed = TRUE))) %>%
  mutate(Period = factor(Period, levels = c("1981–1991", "1992–2002", "2003–2013", "2014–2024", "1981–2024"))) %>%
  arrange(Period) %>%
  pivot_wider(names_from = Period, values_from = Value)

print(Ljubljana_table_periods)

# Saved as comma-delimited CSV (values with decimal comma are quoted)
write_excel_csv(Ljubljana_table_periods, file.path(path_figures, "SZGG_table_periods.csv"))

###############################################################################
# SZGG FIGURE 5 (share of Mediterranean air masses per season and period)
###############################################################################
Metdata_filtered_Ljubljana <- read_xlsx(file.path(path_data, Ljubljana_filename_metdata), skip = 4) %>%
  rename(Day = "...1") %>%
  mutate(Day = as.Date(Day)) %>%
  mutate(across(where(is.character), ~ ifelse(. == "ne", FALSE, .))) %>%
  mutate(across(where(is.character), ~ ifelse(. == "da", TRUE, .))) %>%
  mutate(across(c(`dež`, `sneg`, `dež s snegom`), ~ as.logical(.))) %>%
  dplyr::filter(., `dež` | `sneg` | `dež s snegom`)

HYSPLIT_periods_Ljubljana <- read_xlsx(file.path(path_data, HYSPLIT_filename), sheet = "Ljubljana") %>%
  mutate(Day = as.Date(Date)) %>%
  dplyr::filter(Day %in% Metdata_filtered_Ljubljana$Day) %>% # FILTERS OUT TRAJECTORIES ON DAYS WITHOUT RAIN
  mutate(Date = floor_date(Day, unit = "month")) %>%
  group_by(Date) %>%
  summarize(percent_in_mediterranean = mean(percent_in_mediterranean, na.rm = TRUE)) %>%
  mutate(Year = year(Date),
         Month = month(Date)) %>%
  mutate(Season = case_when(
    Month %in% c(12, 1, 2) ~ "zima",
    Month %in% c(3, 4, 5) ~ "pomlad",
    Month %in% c(6, 7, 8) ~ "poletje",
    Month %in% c(9, 10, 11) ~ "jesen")) %>%
  mutate(Season = factor(Season, levels = c("zima", "pomlad", "poletje", "jesen"))) %>%
  mutate(Period = case_when(
    Year >= 1981 & Year <= 1991 ~ "1981–1991",
    Year >= 1992 & Year <= 2002 ~ "1992–2002",
    Year >= 2003 & Year <= 2013 ~ "2003–2013",
    Year >= 2014 ~ "2014–2024")) %>%
  group_by(Period, Season) %>%
  summarise(mean_percent = mean(percent_in_mediterranean, na.rm = TRUE), .groups = "drop")

ggplot(data = HYSPLIT_periods_Ljubljana, aes(x = Season, y = mean_percent, fill = Period)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_viridis_d(end = 0.9, direction = -1) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw(base_family = "serif") +
  theme(axis.text = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9, colour = "black"),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.spacing.x = unit(0.4, "cm"),
        legend.key.spacing.y = unit(0.1, "cm"),
        legend.text = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 1)),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  labs(y = "delež sredozemskih zračnih mas [%]")

ggsave(filename = paste("SZGG_HYSPLIT_periods", ".", figure_type, sep = ""), path = path_figures, device = figure_type, dpi = 1000, width = 85, height = 60, units = "mm")
