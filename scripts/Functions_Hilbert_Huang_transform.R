# EEMD analysis and calculation of related variables
analyze_variable <- function(time_series_clean, varname, time_vector, n_sims = 1000) {
  # EEMD
  eemd_res <- eemd(time_series_clean, noise_strength = 0.2, ensemble_size = 100, S_number = 4)
  n_imfs <- ncol(as.data.frame(eemd_res)) - 1
  imf_names <- paste0("IMF", 1:n_imfs)
  imf_names_extended <- c(imf_names, "IMF_2_3", "Residual")
  
  # Add combined IMF 2+3 
  eemd_res_df <- as.data.frame(eemd_res)  # clean names like "IMF 1"
  colnames(eemd_res_df) <- gsub(" ", "", colnames(eemd_res_df))  # Optional: remove spaces
  
  combined_imf_2_3 <- eemd_res_df[, "IMF2"] + eemd_res_df[, "IMF3"]
  eemd_res_extended <- eemd_res_df %>%
    add_column(IMF_2_3 = combined_imf_2_3, .before = "Residual")
  
  # Dominant periods
  get_dominant_period <- function(imf) {
    spec <- spectrum(imf, plot = FALSE)
    dominant_freq <- spec$freq[which.max(spec$spec)]
    return(1 / dominant_freq)
  }
  dominant_periods <- sapply(imf_names_extended, function(name) {
    get_dominant_period(eemd_res_extended[[name]])
  })
  
  period_df <- data.frame(IMF = imf_names_extended, Dominant_Period = dominant_periods)
  
  # Amplitudes from Hilbert transform
  amplitude_list <- lapply(imf_names_extended, function(name) {
    Mod(hilbert(as.numeric(eemd_res_extended[[name]])))
  })
  
  mean_amps <- sapply(amplitude_list, mean)
  max_amps <- sapply(amplitude_list, max)
  median_amps <- sapply(amplitude_list, median)
  
  amplitude_df <- data.frame(IMF = imf_names_extended,
                             Mean_Amplitude = mean_amps,
                             Max_Amplitude = max_amps,
                             Median_Amplitude = median_amps)
  
  # Significance test by Monte Carlo
  energy_original <- sapply(imf_names_extended, function(name) var(eemd_res_extended[[name]]))
  energy_surrogates <- matrix(NA, nrow = n_sims, ncol = length(imf_names_extended))
  set.seed(123)
  
  for (sim in 1:n_sims) {
    surrogate_noise <- rnorm(length(time_series_clean), mean = 0, sd = sd(time_series_clean, na.rm = TRUE))
    surrogate_eemd <- eemd(surrogate_noise, noise_strength = 0.2, ensemble_size = 100)
    
    n_imfs_surr <- ncol(as.data.frame(surrogate_eemd)) - 1
    surrogate_combined <- surrogate_eemd[, 2] + surrogate_eemd[, 3]
    residual_surr <- surrogate_eemd[, n_imfs_surr + 1]
    
    surrogate_extended <- cbind(surrogate_eemd[, 1:n_imfs], surrogate_combined, residual_surr)
    
    colnames(surrogate_extended) <- c(paste0("IMF", 1:n_imfs_surr), "IMF_2_3", "Residual")
    
    energy_surrogates[sim, ] <- apply(surrogate_extended, 2, var)
  }
  
  p_values <- sapply(1:length(imf_names_extended), function(i) {
    mean(energy_surrogates[, i] >= energy_original[i])
  })
  
  significance_df <- data.frame(IMF = imf_names_extended,
                                Energy = energy_original,
                                p_value = p_values)
  
  # Combine all summaries
  summary_df <- period_df %>%
    left_join(amplitude_df, by = "IMF") %>%
    left_join(significance_df, by = "IMF") %>%
    mutate(Variable = varname)
  
  return(list(summary = summary_df, eemd_result = eemd_res_extended))
}

# Function to create combined amplitude + oscillation long dataframe per variable
create_amplitude_df <- function(eemd_res, time_vector, varname) {
  n_imfs <- ncol(as.data.frame(eemd_res)) - 1
  
  # Calculate Hilbert amplitudes for each IMF
  amplitude_list <- lapply(1:n_imfs, function(i) {
    Mod(hilbert(as.numeric(eemd_res[, i])))
  })
  
  amplitude_df <- as.data.frame(do.call(cbind, amplitude_list))
  colnames(amplitude_df) <- c(paste0("IMF", 1:(n_imfs-1)), "IMF_2_3")
  amplitude_df <- amplitude_df %>%
    mutate(time = time_vector) %>%
    pivot_longer(cols = -time, names_to = "IMF", values_to = "Value") %>%
    mutate(Variable = varname,
           Type = "Hilbert amplitude")
  
  # Prepare IMF oscillation long format
  eemd_df <- as.data.frame(eemd_res) %>%
    mutate(time = time_vector) %>%
    pivot_longer(cols = -time, names_to = "IMF", values_to = "Value") %>%
    mutate(IMF = factor(IMF, levels = c(paste0("IMF", 1:(n_imfs-1)), "Residual", "IMF_2_3")),
           Variable = varname,
           Type = "IMF oscillation")
  
  # Combine both
  combined_df <- bind_rows(amplitude_df, eemd_df)
  return(combined_df)
}