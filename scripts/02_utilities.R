# Utilities: Functions

# Logit function
# Input p in (0,1)
logit <- function(p) {
  p <- pmax(pmin(p, 0.999), 0.001)
  return(log(p / (1 - p)))
}

# Expit function
# Input x in real line
expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

# Min-max transformation: Re-scales to [0,1]
# Input vector v
rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE, finite = TRUE)
  return( (x - rng[1]) / (rng[2] - rng[1]) )
}

# Define the bootstrap function for OLS
OLS_sub <- function(data,
                    Qform = paste("Y ~ A +", paste(
                      colnames(data)[grepl("^W", colnames(data))], collapse = " + ")))
{
  # Fit Qform to data
  Qfit <- glm2(Qform, data = data)
  
  # Predict
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  Qfit_0W <- predict(Qfit, newdata = data.frame(A = 0, data[covariates]))
  Qfit_1W <- predict(Qfit, newdata = data.frame(A = 1, data[covariates]))
  
  # Evaluate estimate
  estimate <- mean(Qfit_1W - Qfit_0W)
  return(estimate)
}

OLS_stat <- function(data, indices,
                     Qform = paste("Y ~ A +", paste(
                       colnames(data)[grepl("^W", colnames(data))], collapse = " + "))) 
{
  
  sample_data <- data[indices,]
  return(OLS_sub(sample_data, Qform = Qform))
}

# Old bootstrap function
bootstrap_variance <- function(data, statistic, alpha = 0.05, R = 1000) {
  l <- length(data)
  bootstrap_estimates <- numeric(R)
  
  for (i in 1:R) {
    resampled_data <- sample_n(as.data.frame(data), l, replace = TRUE)
    bootstrap_estimates[i] <- statistic(resampled_data)
  }
  
  return(
    tibble(
      var = var(bootstrap_estimates),
      ci_lower = quantile(bootstrap_estimates, probs = alpha/2, na.rm = TRUE),
      ci_upper = quantile(bootstrap_estimates, probs = 1-alpha/2, na.rm = TRUE)
      )
  )
}


# Function to evaluate estimators
evaluate_estimators <- function(data, estimator_list) {
  results <- list()
  
  for (name in names(estimator_list)) {
    estimator <- estimator_list[[name]]
    
    result <- estimator(data)
    
    results[[name]] <- data.frame(
      Initial.Est = ifelse("est_init" %in% names(result), result$est_init, NA),
      Estimate = ifelse("est" %in% names(result), result$est, NA),
      Standard.Error = ifelse("se" %in% names(result), result$se, NA),
      CI_Lower = ifelse("ci_lower" %in% names(result), result$ci_lower, NA),
      CI_Upper = ifelse("ci_upper" %in% names(result), result$ci_upper, NA),
      Bias = ifelse("bias" %in% names(result), result$bias, NA),
      Bias.to.se = ifelse("bias.to.se" %in% names(result), result$bias.to.se, NA)
      #MSE = if ("mse" %in% names(result)) result$mse else NA
    )
  }
  
  results_df <- do.call(rbind, results)
  results_df$Algorithm <- rownames(results_df)
  rownames(results_df) <- NULL 
  results_df <- results_df[, c("Algorithm", "Initial.Est", "Estimate", "Standard.Error", "CI_Lower", "CI_Upper", "Bias", "Bias.to.se")]#, "MSE")]
  return(results_df)
}

# Evaluate estimators
MC_est <- function(k = 100, sample_no = 1000,
                   simulation = sim1a, est_list = estimators)
{
  
  simulation_results <- replicate(
    n = k, 
    expr = evaluate_estimators(simulation(sample_no), estimator_list = est_list), 
    simplify = FALSE
  )
  
  combined_results <- do.call(rbind, simulation_results)
  combined_results[,2:8] <- signif(combined_results[,2:8], digits = 3)
  return(combined_results)
}




# Aggregate results, creating MC statistics
MC_res <- function(combined_results, true_value = 1) {
  k <- dim(combined_results)[1]
  combined_results %>% 
    group_by(Algorithm) %>%
    summarize(
      Mean_Est = mean(Estimate, na.rm = TRUE),
      Mean_Bias = mean(Bias, na.rm = TRUE),
      Mean_Var = mean(Standard.Error^2, na.rm = TRUE),
      Sample_Var = (1 / k) * sum((Estimate - mean(Estimate))^2, na.rm = TRUE),
      Coverage = mean((true_value >= CI_Lower) & (true_value <= CI_Upper), na.rm = TRUE),
      Mean_Bias.to.se = mean(Bias.to.se, na.rm = TRUE),
      Mean_MSE = Mean_Bias^2 + Sample_Var,
      Mean_CI_Lower = mean(CI_Lower, na.rm = TRUE),
      Mean_CI_Upper = mean(CI_Upper, na.rm = TRUE),
      CI_2.5th = quantile(Estimate, probs = 0.025, na.rm = TRUE),
      CI_97.5th = quantile(Estimate, probs = 0.975, na.rm = TRUE)
      #Sample_CI_Lo = Mean_Estimate - 1.96*sqrt(Sample_Var), # TYPO: This should be Estimate, not Mean_Estimate
      #Sample_CI_Up = Mean_Estimate + 1.96*sqrt(Sample_Var), # TYPO: This should be Estimate, not Mean_Estimate
      #Mean_MSE = mean(MSE, na.rm = TRUE),
      #Sample_Coverage = mean((true_value >= Estimate - 1.96*sqrt(Sample_Var)) & 
      #                         (true_value <= Estimate + 1.96*sqrt(Sample_Var)), na.rm = TRUE) # Fix Sample_CI here!
    ) %>% 
    mutate(across(-c(Algorithm), ~ signif(., digits = 3)))
}

## Plot histogram per estimator
# res = output of MC_est()
hist_per_est <- function(res, xlimits = range(res$Estimate), true_value = 1){
  res %>% ggplot(aes(x = Estimate, fill = Algorithm)) +
    geom_density(alpha = 0.4, color = "black") +
    scale_x_continuous(limits = xlimits) +
    facet_wrap(~Algorithm) + #, scales = "free") +
    labs(title = "Distribution of ATE Estimates by Algorithm",
         x = "ATE Estimate",
         y = "Density") +
    geom_vline(xintercept = true_value, color = "black", linetype = "dashed") +
    #stat_function(fun = dnorm,
    #              args = list(mean = mean(res$Estimate), sd = sd(res$Estimate)), 
    #              color = "red", size = 1) +
    #facet_wrap(~Algorithm, scales = "free") +
    theme_minimal()
}

hist_per_est_CI <- function(res, xlimits = c(min(res$CI_Lower),max(res$CI_Upper))) {
  res %>% 
    dplyr::select(Algorithm, CI_Lower, CI_Upper) %>%
    tidyr::pivot_longer(cols = c(CI_Lower, CI_Upper), 
                        names_to = "CI_Bound", 
                        values_to = "CI_Value") %>% 
    ggplot(aes(x = CI_Value, fill = CI_Bound, color = CI_Bound)) +
    geom_density(alpha = 0.4) +
    scale_x_continuous(limits = xlimits) +
    facet_wrap(~Algorithm) + #, scales = "free_y") +
    labs(title = "Distribution of Confidence Interval Bounds by Algorithm",
         x = "CI Bound Value",
         y = "Density",
         fill = "CI Bound",
         color = "CI Bound") +
    geom_vline(xintercept = true_value, color = "black", linetype = "dashed") +
    theme_minimal()
}

## Plot histograms on top of each other
# res = output of MC_est()
hist_all <- function(res, xlimits = range(res$Estimate), true_value = 1){
  res %>% ggplot(aes(x = Estimate, fill = Algorithm)) +
  geom_density(alpha = 0.6, na.rm = TRUE) +
    scale_x_continuous(limits = xlimits) +
    labs(title = "Distribution of ATE Estimates by Algorithm",
       x = "ATE Estimate",
       y = "Density") +
    geom_vline(xintercept = true_value, color = "black", linetype = "dashed") +
    theme_minimal()
}

# Plot estimate ± error bar, ordering by CI width in descending order
plot_est_CI <- function(res, custom_title = "Results",
                        xlimits = c(min(res$CI_Lower), max(res$CI_Upper)),
                        true_value = 1) {
  
  res <- res %>%
    mutate(CI_Width = CI_Upper - CI_Lower) %>%
    arrange(desc(CI_Width)) 
  
  res$Algorithm <- factor(res$Algorithm, levels = unique(res$Algorithm))
  
  ggplot(res, aes(x = Algorithm, y = Estimate)) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, show.legend = FALSE) +
    geom_point(size = 3) +
    geom_hline(yintercept = true_value, color = "black", linetype = "dashed") +
    scale_y_continuous(limits = xlimits) +
    coord_flip() +
    labs(
      title = custom_title,
      x = "",
      y = "ATE estimates"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
}