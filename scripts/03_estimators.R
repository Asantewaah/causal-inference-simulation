### List of estimator functions

source("./scripts/02_utilities.R")

## Difference-in-means estimator
# Data: Tibble
# Alpha: Allowable type I error
estimate_DM <- function(data, alpha = 0.05, true_value = 1) {
  data %>% group_by(A) %>% 
    summarise(
      mean_outcome = mean(Y, na.rm = TRUE),
      var_norm = var(Y) / n()
    ) %>% 
    summarise(
      estimator = "DM",
      Qform = NA,
      gform = NA,
      est_init = NA,
      est = mean_outcome[2] - mean_outcome[1],
      se = sqrt(var_norm[2] + var_norm[1]),
      ci_lower = est - qnorm(1-alpha) * se,
      ci_upper = est + qnorm(1-alpha) * se,
      bias = est - true_value,
      bias.to.se = bias / se
    ) %>% 
    as_tibble()
}

## Ordinary Least Squares estimator
# Data: Tibble
# Alpha: Allowable type I error
estimate_OLS <- function(data, alpha = 0.05, true_value = 1,
                         Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + ")))
                         )
  {
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  # OLS formula Y ~ A + W1 + ... + Wk
  #formula_OLS <- as.formula(paste("Y ~ A +", paste(covariates, collapse = " + "))) #, "+ A:(", paste(covariates, collapse = " + "), ")"))
  model <- lm(Qform, data = data)
  coefs_A <- coef(summary(model))["A",]
  return(
    tibble(
      estimator = "OR",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = NA,
      est_init = NA,
      est = coefs_A["Estimate"],
      se = coefs_A["Std. Error"],
      ci_lower = confint(model, "A", level = 1-alpha)[1],
      ci_upper = confint(model, "A", level = 1-alpha)[2],
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}


## Ordinary Least Squares substitution estimator
# Data: Tibble
# Alpha: Allowable type I error
estimate_OLS_sub <- function(data, alpha = 0.05, true_value = 1,
                             Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + ")))
                             )
  {
  
  estimate <- OLS_sub(data, Qform) 
  boot_res <- boot(data,
                   statistic = function(data, indices){ OLS_stat(data, indices, Qform)},
                   R = 1000)

  return(
    tibble(
      estimator = "OR_sub",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = NA,
      est_init = NA,
      est = estimate,
      se = sd(boot_res$t),
      ci_lower = quantile(boot_res$t, probs = alpha/2, na.rm = TRUE),
      ci_upper = quantile(boot_res$t, probs = 1-alpha/2, na.rm = TRUE),
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}


## Inverse Propensity Weighting estimator (survey package)
# Data: Tibble
# Alpha: Allowable type I error
# gform: Propensity score formula, e.g., A ~ W1*W2
estimate_IPW <- function(data, alpha = 0.05, true_value = 1,
                         gform = as.formula(paste("A ~", paste(
                             colnames(data)[grepl("^W", colnames(data))], collapse = " + "))
                             )
                         ) {
  # covariates <- colnames(data)[grepl("^W", colnames(data))]
  propensity_model <- glm2(gform,
                           family = binomial(link = "logit"),
                           data = data)
  
  data$weights <- ifelse(data$A == 1, 
                         1 / predict(propensity_model, type = "response", newdata = data),
                         1 / (1 - predict(propensity_model, type = "response", newdata = data))
                         )
  
  ipw_design <- svydesign(ids = ~1, weights = ~weights, data = data)
  
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  formula_IPW <- as.formula("Y ~ A")
  #formula_IPW <- as.formula(paste("Y ~ A +", paste(
  #  colnames(data)[grepl("^W", colnames(data))], collapse = " + ")))
  ipw_result <- svyglm(formula_IPW, design = ipw_design)
  
  coef_ipw <- coef(ipw_result)
  se_ipw <- sqrt(diag(vcov(ipw_result)))
  ci_ipw <- coef_ipw["A"] + c(-qnorm(1-alpha), qnorm(1-alpha)) * se_ipw["A"]

  return(
    tibble(
      estimator = "IPW",
      Qform = NA,
      gform = paste("A ~",as.character(gform)[3]),
      est_init = NA,
      est = coef_ipw["A"],
      se = se_ipw["A"],
      ci_lower = ci_ipw[1],
      ci_upper = ci_ipw[2],
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}


## Inverse Propensity Weighting estimator (by hand)
# Data: Tibble
# Alpha: Allowable type I error
estimate_IPW_hand <- function(data, alpha = 0.05, true_value = 1,
                              gform = as.formula(paste("A ~", paste(
                                colnames(data)[grepl("^W", colnames(data))], collapse = " + "))
                              )
) {
  # covariates <- colnames(data)[grepl("^W", colnames(data))]
  propensity_model <- glm2(gform,
                           family = binomial(link = "logit"),
                           data = data)
  
  data$weights <- ifelse(data$A == 1, 
                         1 / predict(propensity_model, type = "response", newdata = data),
                         1 / (1 - predict(propensity_model, type = "response", newdata = data))
  )
  
  est_hand <- data %>% summarise(mean(Y * (2*A-1) * weights)) %>% as.numeric()
  
  return(
    tibble(
      estimator = "IPW_hand",
      est = est_hand,
      se = NA,
      ci_lower = NA,
      ci_upper = NA,
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}


## TMLE estimator
# Data: Tibble
# Alpha: Allowable type I error
# Qform
# gform
# TMLE specs: Weighted, single_H, bounded_continuous_outcome
estimate_TMLE <- function(data, alpha = 0.05, true_value = 1,
                          Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + "))),
                          gform = as.formula(paste("A ~", paste(
                            colnames(data)[grepl("^W", colnames(data))], collapse = " + "))),
                          weight = TRUE, use_single_H = FALSE, transforms = TRUE) {
  
  if(transforms){
    a <- min(data$Y)
    b <- max(data$Y)
    data$Y <- (data$Y - a) / (b - a)
  }
  
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  
  # Outcome regression
  Q_init <- glm(Qform, family = gaussian, data = data)
  data$Q0W <- predict(Q_init, newdata = data.frame(A = 0, data[covariates]))
  data$Q1W <- predict(Q_init, newdata = data.frame(A = 1, data[covariates]))
  data$Q <- predict(Q_init, type = "response")
  
  # Treatment mechanism
  g_init <- glm(gform, family = binomial, data = data)
  data$g1W <- predict(g_init, type = "response")
  
  data$k_Q1 <- if (weight) 1 else 1 / data$g1W
  data$k_Q0 <- if (weight) 1 else 1 / (1 - data$g1W)
  
  data$H <-if (weight) data$A - (1 - data$A) else (data$A / data$g1W) - ((1 - data$A) / (1 - data$g1W))
  data$weights_single_H <- if (weight) ifelse(data$A == 1, 1 / data$g1W, 1 / (1 - data$g1W)) else 1
  data$H1 <- if(weight) data$A else data$A / data$g1W
  data$H0 <- if(weight) (1  - data$A) else ((1 - data$A) / (1 - data$g1W))
  data$weights_H1 <- if (weight) 1 / data$g1W else 1
  data$weights_H0 <- if (weight) 1 / (1 - data$g1W) else 1
  
  if (use_single_H) {
    if (transforms) {
      second_stage <- glm2(data$Y ~ -1 + H + offset(logit(data$Q)), 
                           family = quasibinomial, data = data, weights = weights_single_H)
      epsilon <- coef(second_stage)[1]
      Q_star1 <- expit(logit(data$Q1W) + epsilon * data$k_Q1)
      Q_star0 <- expit(logit(data$Q0W) - epsilon * data$k_Q0)
    } else {
      second_stage <- glm2(data$Y ~ -1 + H + offset(data$Q), 
                           family = gaussian, data = data, weights = weights_single_H)
      epsilon <- coef(second_stage)[1]
      Q_star1 <- data$Q1W + epsilon * data$k_Q1
      Q_star0 <- data$Q0W - epsilon * data$k_Q0
    }
  } else {
    if (transforms) {
      second_stage_1 <- glm2(data$Y ~ -1 + H1 + offset(logit(data$Q)), 
                             family = quasibinomial, data = data, weights = weights_H1)
      second_stage_0 <- glm2(data$Y ~ -1 + H0 + offset(logit(data$Q)), 
                             family = quasibinomial, data = data, weights = weights_H0)
      epsilon_1 <- coef(second_stage_1)[1]
      epsilon_0 <- coef(second_stage_0)[1]
      Q_star1 <- expit(logit(data$Q1W) + epsilon_1 * data$k_Q1)
      Q_star0 <- expit(logit(data$Q0W) + epsilon_0 * data$k_Q0)
    } else {
      second_stage_1 <- glm2(data$Y ~ -1 + H1 + offset(data$Q), 
                             family = gaussian, data = data, weights = weights_H1)
      second_stage_0 <- glm2(data$Y ~ -1 + H0 + offset(data$Q), 
                             family = gaussian, data = data, weights = weights_H0)
      epsilon_1 <- coef(second_stage_1)[1]
      epsilon_0 <- coef(second_stage_0)[1]
      Q_star1 <- data$Q1W + epsilon_1 * data$k_Q1
      Q_star0 <- data$Q0W + epsilon_0 * data$k_Q0
    }
  }
  
  if (transforms) {
    ATE_pre <- mean(data$Q1W - data$Q0W) * (b - a)
    ATE <- mean(Q_star1 - Q_star0) * (b - a)
  }else{
    ATE_pre <- mean(data$Q1W - data$Q0W)
    ATE <- mean(Q_star1 - Q_star0)
  }
  
  if  (transforms) {
    ATE_n <- ATE / (b - a)
  }else{
    ATE_n <- ATE 
  }
  
  data$H_n <- (data$A / data$g1W) - ((1 - data$A) / (1 - data$g1W))
  data$H_1 <- data$A / data$g1W
  data$H_0 <- (1 - data$A) / (1 - data$g1W)
  
  if(use_single_H){
    IF <- data$H_n * (data$Y - (data$A * Q_star1 + (1 - data$A) * Q_star0)) + (Q_star1 - Q_star0) - ATE_n
    IF_1 <- mean(data$H_n * (data$Y - (data$A * Q_star1 + (1 - data$A) * Q_star0)))
    IF_0 <- mean((Q_star1 - Q_star0) - ATE_n)
  }else{
    IF <- data$H_1* (data$Y - Q_star1) - data$H_0 * (data$Y - Q_star0) + (Q_star1 - Q_star0) - ATE_n
    IF_1 <- mean(data$H_1* (data$Y - Q_star1) - data$H_0 * (data$Y - Q_star0))
    IF_0 <- mean((Q_star1 - Q_star0) - ATE_n)
  }
  
  n <- nrow(data)
  SE <- if(transforms) sqrt(var(IF) / n) * (b - a) else sqrt(var(IF) / n)
  z <- qnorm(1 - alpha / 2)
  v <- if(transforms) (var(IF) / n) * (b - a)^2 else var(IF) / n #var(IF)/n THIS WAS A TYPO
  CI_lower <- ATE - z * SE
  CI_upper <- ATE + z * SE
  #MSE <- bias^2 + v
  
  #return(list(estimate = ATE, var = v, bias = bias, se = SE, mse = MSE, ci_lower = CI_lower, ci_upper = CI_upper))
  
  return(
    tibble(
      estimator = "TMLE",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = paste("A ~",as.character(gform)[3]),
      est_init = ATE_pre,
      est = ATE,
      se = SE,
      ci_lower = CI_lower,
      ci_upper = CI_upper,
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}


## OSE estimator (One-Step Estimator(OSE)/Augmented inverse probability weighting (AIPW))
# Data: Tibble
# Alpha: Allowable type I error
# Qform
# gform

estimate_OSE <- function(data, alpha = 0.05, true_value = 1,
                         Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + "))),
                         gform = as.formula(paste("A ~", paste(
                           colnames(data)[grepl("^W", colnames(data))], collapse = " + ")))) {
  
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  
  Q_init <- glm(Qform, family = gaussian, data = data)
  data$Q0W <- predict(Q_init, newdata = data.frame(A = 0, data[covariates]))
  data$Q1W <- predict(Q_init, newdata = data.frame(A = 1, data[covariates]))
  
  g_init <- glm(gform, family = binomial, data = data)
  data$g1W <- predict(g_init, type = "response")
  data$g0W <- 1 - data$g1W
  
  data$H1 <- data$A / data$g1W
  data$H0 <- (1 - data$A) / data$g0W
  
  data$IF <- (data$H1 - data$H0) * (data$Y - (data$A * data$Q1W + (1 - data$A) * data$Q0W)) + (data$Q1W - data$Q0W)
  
  ATE <- mean(data$IF)
  SE <- sd(data$IF) / sqrt(nrow(data))
  v <- var(data$IF)/nrow(data)
  CI_lower <- ATE - 1.96 * SE
  CI_upper <- ATE + 1.96 * SE
  #MSE <- bias^2 + v
  
  return(
    tibble(
      estimator = "OSE/A-IPW",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = paste("A ~",as.character(gform)[3]),
      est = ATE,
      se = SE,
      ci_lower = CI_lower,
      ci_upper = CI_upper,
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}

## CTMLE estimator : Greedy
# Data: Tibble
# Alpha: Allowable type I error
# Qform
# gform
# CTMLE specs: Weighted, double_H, bounded_continuous_outcome

estimate_GREEDY_CTMLE <- function(data, alpha = 0.05, true_value = 1,
                         Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + ")))) {
  
  a <- min(data$Y)
  b <- max(data$Y)
  data$Y <- (data$Y - a) / (b - a)
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  
  # Outcome regression
  
  Q_init <- glm(Qform, family = gaussian, data = data)
  data$Q0W <- predict(Q_init, newdata = data.frame(A = 0, data[covariates]))
  data$Q1W <- predict(Q_init, newdata = data.frame(A = 1, data[covariates]))
  data$Q <- predict(Q_init, type = "response")
  
  #Variable selection
  
  selected_vars <- c()
  remaining_vars <- covariates
  min_loss <- Inf
  best_g <- NULL
  
  for (k in 1:length(covariates)) {
    losses <- c()
    g_models <- list()
    
    for (cov_var in remaining_vars) {
      vars_to_use <- c(selected_vars, cov_var)
      gform  <-  as.formula(paste("A ~", paste(vars_to_use, collapse = " + ")))
      g_model <- glm2(gform, family = binomial, data = data)
      g1W <- predict(g_model, type = "response")
      
      data$H1 <- data$A 
      data$H0 <- (1 - data$A) 
      data$weights_H1 <-  1 / g1W 
      data$weights_H0 <-  1 / (1 - g1W) 
      
      second_stage_1 <- glm2(data$Y ~ -1 + H1 + offset(logit(data$Q)), family = quasibinomial, data = data, weights = weights_H1)
      second_stage_0 <- glm2(data$Y ~ -1 + H0 + offset(logit(data$Q)), family = quasibinomial, data = data, weights = weights_H0)
      epsilon_1 <- coef(second_stage_1)[1]
      epsilon_0 <- coef(second_stage_0)[1]
      Q_star1 <- expit(logit(data$Q1W) + epsilon_1 )
      Q_star0 <- expit(logit(data$Q0W) + epsilon_0 )
      
      loss <- mean((data$Y - (data$A * Q_star1 + (1 - data$A) * Q_star0))^2)
      losses <- c(losses, loss)
      
      g_models[[cov_var]] <- list(g_model = g_model, loss = loss, g1W = g1W,
                                  epsilon_1 = epsilon_1, epsilon_0 = epsilon_0,
                                  Q_star1 = Q_star1, Q_star0 = Q_star0)
      
      #cat("\nIteration:", k, "| Testing Combination:", paste(vars_to_use, collapse = " + "), "| Loss:", loss, "\n")
    }
    best_var <- remaining_vars[which.min(losses)]
    best_model <- g_models[[best_var]]
    
    if (best_model$loss < min_loss) {
      min_loss <- best_model$loss
      selected_vars <- c(selected_vars, best_var)
      remaining_vars <- setdiff(remaining_vars, best_var)
      best_g <- best_model$g_model
      best_g1W <- best_model$g1W
      epsilon_1 <- best_model$epsilon_1
      epsilon_0 <- best_model$epsilon_0
      Q_star1 <- best_model$Q_star1
      Q_star0 <- best_model$Q_star0
      best_combination <- selected_vars
      
      #cat("Selected combination so far:", paste(selected_vars, collapse = " + "), "| Updated min loss:", min_loss, "\n")
    } else {
      #cat("No improvement in loss, stopping selection.\n")
      break
    }
  }
  
  #cat("\nBest combination chosen:", paste(best_combination, collapse = " + "), "\n")
  
  ATE_pre <- mean(data$Q1W - data$Q0W) * (b - a)
  ATE <- mean(Q_star1 - Q_star0) * (b - a)
  ATE_n <- ATE/(b - a)
  
  data$H_1 <- data$A / best_g1W
  data$H_0 <- (1 - data$A) / (1 - best_g1W)
  IF <- (data$A / best_g1W) * (data$Y - Q_star1) - ((1 - data$A) / (1 - best_g1W)) * (data$Y - Q_star0) + (Q_star1 - Q_star0) - ATE_n
  
  IF_1 <- mean((data$A / best_g1W) * (data$Y - Q_star1) - ((1 - data$A) / (1 - best_g1W)) * (data$Y - Q_star0))
  IF_0 <- mean((Q_star1 - Q_star0) - ATE_n)
  SE <- sqrt(var(IF) / nrow(data))* (b - a)
  z <- qnorm(1 - alpha / 2)
  v <- var(IF)/nrow(data)
  CI_lower <- ATE - z * SE
  CI_upper <- ATE + z * SE
  #MSE <- bias^2 + v
  
  return(
    tibble(
      estimator = "C-TMLE: GREEDY",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = ifelse(length(selected_vars) > 0, 
                     paste("A ~", paste(best_combination, collapse = " + ")), "A ~ 1"),
      est_init = ATE_pre,
      est = ATE,
      se = SE,
      ci_lower = CI_lower,
      ci_upper = CI_upper,
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}

estimate_GREEDY_CTMLE1 <- function(data, alpha = 0.05, true_value = 1,
                                  Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + ")))) {
  
  a <- min(data$Y)
  b <- max(data$Y)
  data$Y <- (data$Y - a) / (b - a)
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  X <- as.matrix(data[covariates])
  
  # Outcome regression - keeping this the same
  Q_init <- glm(Qform, family = gaussian, data = data)
  data$Q0W <- predict(Q_init, newdata = data.frame(A = 0, data[covariates]))
  data$Q1W <- predict(Q_init, newdata = data.frame(A = 1, data[covariates]))
  data$Q <- predict(Q_init, type = "response")
  
  # Variable selection using greedy approach but with glmnet for each model
  selected_vars <- c()
  remaining_vars <- covariates
  min_loss <- Inf
  best_g <- NULL
  
  for (k in 1:length(covariates)) {
    losses <- c()
    g_models <- list()
    
    for (cov_var in remaining_vars) {
      vars_to_use <- c(selected_vars, cov_var)
      
      # Create design matrix for selected variables
      X_subset <- as.matrix(data[, vars_to_use, drop = FALSE])  # Explicitly select columns
      
      # Handle case where X_subset might have only one column or be empty
      if (ncol(X_subset) < 2) {
        if (ncol(X_subset) == 1){
          # If only one variable is selected, add a dummy variable
          X_subset <- cbind(X_subset, rep(1, nrow(data)))
          colnames(X_subset)[2] <- "dummy" # Ensure the dummy has a valid name
        } else {
          # No variables selected, use an intercept-only model
          g_model <- glm(data$A ~ 1, family = "binomial", data = data)
          g1W <- predict(g_model, type = "response")
          
          # Store model and results
          g_models[[cov_var]] <- list(g_model = g_model, loss = Inf, g1W = g1W,
                                      epsilon_1 = 0, epsilon_0 = 0, Q_star1 = data$Q1W, 
                                      Q_star0 = data$Q0W, vars = vars_to_use)
          losses <- c(losses, Inf)
          next
        }
      }
      
      # Using glmnet with cross-validation for lambda selection
      lambda_seq <- 10^seq(-3, 1, length.out = 50)
      cv_lasso <- tryCatch({
        cv.glmnet(X_subset, data$A, family = "binomial", alpha = 1, lambda = lambda_seq, nfolds = 5)
      }, error = function(e) {
        # If glmnet fails, return a default error object
        message("glmnet error: ", e$message)
        return(NULL)
      })
      
      if (is.null(cv_lasso)) {
        losses <- c(losses, Inf)
        g_models[[cov_var]] <- list(g_model = NULL, loss = Inf, g1W = rep(0.5, nrow(data)),
                                    epsilon_1 = 0, epsilon_0 = 0, Q_star1 = data$Q1W, 
                                    Q_star0 = data$Q0W, vars = vars_to_use)
        next
      }
      
      best_lambda <- cv_lasso$lambda.min
      
      # Fit model with best lambda
      g_model <- glmnet(X_subset, data$A, family = "binomial", alpha = 1, lambda = best_lambda)
      g1W <- predict(g_model, newx = X_subset, type = "response")
      
      # Apply bounds to avoid extreme propensity scores
      g1W <- pmax(pmin(g1W, 0.99), 0.01)
      
      data$H1 <- data$A 
      data$H0 <- (1 - data$A) 
      data$weights_H1 <- 1 / g1W 
      data$weights_H0 <- 1 / (1 - g1W) 
      
      second_stage_1 <- glm2(data$Y ~ -1 + H1 + offset(logit(data$Q)), family = quasibinomial, data = data, weights = weights_H1)
      second_stage_0 <- glm2(data$Y ~ -1 + H0 + offset(logit(data$Q)), family = quasibinomial, data = data, weights = weights_H0)
      epsilon_1 <- coef(second_stage_1)[1]
      epsilon_0 <- coef(second_stage_0)[1]
      Q_star1 <- expit(logit(data$Q1W) + epsilon_1)
      Q_star0 <- expit(logit(data$Q0W) + epsilon_0)
      
      loss <- mean((data$Y - (data$A * Q_star1 + (1 - data$A) * Q_star0))^2)
      losses <- c(losses, loss)
      
      g_models[[cov_var]] <- list(g_model = g_model, loss = loss, g1W = g1W,
                                  epsilon_1 = epsilon_1, epsilon_0 = epsilon_0,
                                  Q_star1 = Q_star1, Q_star0 = Q_star0,
                                  lambda = best_lambda, vars = vars_to_use)
    }
    
    # Select the best variable
    best_var <- remaining_vars[which.min(losses)]
    best_model <- g_models[[best_var]]
    
    if (best_model$loss < min_loss) {
      min_loss <- best_model$loss
      selected_vars <- c(selected_vars, best_var)
      remaining_vars <- setdiff(remaining_vars, best_var)
      best_g <- best_model$g_model
      best_g1W <- best_model$g1W
      epsilon_1 <- best_model$epsilon_1
      epsilon_0 <- best_model$epsilon_0
      Q_star1 <- best_model$Q_star1
      Q_star0 <- best_model$Q_star0
      best_combination <- selected_vars
      #best_lambda <- best_model$lambda
    } else {
      break
    }
  }
  
  # For the final model, double-check with complete LASSO (optional alternative approach)
  # This fits a full LASSO model with all variables and compares it to the greedy approach
  lambda_seq <- 10^seq(-3, 1, length.out = 50)
  cv_full_lasso <- cv.glmnet(X, data$A, family = "binomial", alpha = 1, lambda = lambda_seq, nfolds = 5)
  best_full_lambda <- cv_full_lasso$lambda.min
  full_g_model <- glmnet(X, data$A, family = "binomial", alpha = 1, lambda = best_full_lambda)
  full_g1W <- predict(full_g_model, newx = X, type = "response")
  full_g1W <- pmax(pmin(full_g1W, 0.99), 0.01)
  
  data$H1 <- data$A 
  data$H0 <- (1 - data$A) 
  data$weights_H1 <- 1 / full_g1W 
  data$weights_H0 <- 1 / (1 - full_g1W) 
  
  full_second_stage_1 <- glm2(data$Y ~ -1 + H1 + offset(logit(data$Q)), family = quasibinomial, data = data, weights = weights_H1)
  full_second_stage_0 <- glm2(data$Y ~ -1 + H0 + offset(logit(data$Q)), family = quasibinomial, data = data, weights = weights_H0)
  full_epsilon_1 <- coef(full_second_stage_1)[1]
  full_epsilon_0 <- coef(full_second_stage_0)[1]
  full_Q_star1 <- expit(logit(data$Q1W) + full_epsilon_1)
  full_Q_star0 <- expit(logit(data$Q0W) + full_epsilon_0)
  
  full_loss <- mean((data$Y - (data$A * full_Q_star1 + (1 - data$A) * full_Q_star0))^2)
  
  # Choose the approach with the lower loss
  if (full_loss < min_loss) {
    best_g1W <- full_g1W
    Q_star1 <- full_Q_star1
    Q_star0 <- full_Q_star0
    best_combination <- covariates # All variables used in LASSO
    best_g <- full_g_model
  }
  
  # Calculate treatment effect estimates
  ATE_pre <- mean(data$Q1W - data$Q0W) * (b - a)
  ATE <- mean(Q_star1 - Q_star0) * (b - a)
  ATE_n <- ATE/(b - a)
  
  # Influence function for statistical inference
  data$H_1 <- data$A / best_g1W
  data$H_0 <- (1 - data$A) / (1 - best_g1W)
  IF <- (data$A / best_g1W) * (data$Y - Q_star1) - ((1 - data$A) / (1 - best_g1W)) * (data$Y - Q_star0) + (Q_star1 - Q_star0) - ATE_n
  
  SE <- sqrt(var(IF) / nrow(data)) * (b - a)
  z <- qnorm(1 - alpha / 2)
  CI_lower <- ATE - z * SE
  CI_upper <- ATE + z * SE
  
  # Use the non-zero coefficients for reporting the model formula
  if (is.null(best_g) || full_loss < min_loss) {
    # If using the full model
    coef_names <- names(which(coef(cv_full_lasso, s = "lambda.min")[-1] != 0))
    gform_vars <- if (length(coef_names) > 0) paste(coef_names, collapse = " + ") else "1"
  } else {
    # If using the greedy approach
    gform_vars <- paste(best_combination, collapse = " + ")
  }
  
  return(
    tibble(
      estimator = "C-TMLE: GREEDY+LASSO",
      Qform = paste("Y ~", as.character(Qform)[3]),
      gform = paste("A ~", gform_vars),
      est_init = ATE_pre,
      est = ATE,
      se = SE,
      ci_lower = CI_lower,
      ci_upper = CI_upper,
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}


## CTMLE estimator : SuperLearner
# Data: Tibble
# Alpha: Allowable type I error
# Qform
# gform
# CTMLE specs:4 model in SL library

estimate_SL_CTMLE <- function(data, true_value = 1, alpha = 0.05, Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + "))),
                              gform = as.formula(paste("A ~", paste(
                                colnames(data)[grepl("^W", colnames(data))], collapse = " + "))),max_iter = 3) {
  
  a <- min(data$Y)
  b <- max(data$Y)
  data$Y <- (data$Y - a) / (b - a)
  
  #SuperLearner library
  SL.library <- c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.gam")
  
  
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  Q_init <- glm(Qform, family = gaussian, data = data)
  data$Q0W <- predict(Q_init, newdata = data.frame(A = 0, data[covariates]))
  data$Q1W <- predict(Q_init, newdata = data.frame(A = 1, data[covariates]))
  data$Q <- predict(Q_init, type = "response")
  
  # Initial propensity score (g)
  g_init <- glm(gform, family = binomial, data = data)
  data$g1W <- predict(g_init, type = "response")
  
  # Collaborative TMLE Iterative Selection
  for (i in 1:max_iter) {
    
    data$H1 <- data$A / data$g1W
    data$H0 <- (1 - data$A) / (1 - data$g1W)
    
    fluctuation_model1 <- glm(Y ~ -1 + H1 + offset(logit(data$Q)), family = quasibinomial, data = data)
    fluctuation_model0 <- glm(Y ~ -1 + H0 + offset(logit(data$Q)), family = quasibinomial, data = data)
    epsilon_1 <- coef(fluctuation_model1)[1]
    epsilon_0 <- coef(fluctuation_model0)[1]
    
    data$Q_star1 <- expit(logit(data$Q0W) + epsilon_1 * (1 / (1 - data$g1W)))
    data$Q_star0 <- expit(logit(data$Q1W) + epsilon_0 * (1 / data$g1W))
    
    g_superlearner <- SuperLearner(Y = data$A, X = data[, covariates], family = binomial, SL.library = SL.library)
    data$g1W <- g_superlearner$SL.predict
  }
  
  ATE_pre <- mean(data$Q1W - data$Q0W) * (b - a)
  ATE <- mean(data$Q1W - data$Q0W) * (b - a)
  ATE_n <- ATE / (b - a)
  
  IF <- (data$A / data$g1W) * (data$Y - Q_star1) - ((1 - data$A) / (1 - data$g1W)) * (data$Y - Q_star0) + (Q_star1 - Q_star0) - ATE_n
  SE <- sqrt(var(IF) / nrow(data)) * (b - a)
  z <- qnorm(1 - alpha / 2)
  CI_lower <- ATE - z * SE
  CI_upper <- ATE + z * SE
  #MSE <- bias^2 + v
  
  return(
    tibble(
      estimator = "C-TMLE: SUPERLEARNER",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = NA,
      est_init = ATE_pre,
      est = ATE,
      se = SE,
      ci_lower = CI_lower,
      ci_upper = CI_upper,
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
  
}

## C-TMLE estimator : C-TMLE 0
# Data: Tibble
# Alpha: Allowable type I error
# Qform
# gform
# Predefined CTMLE specs: Weighted, double_H, bounded_continuous_outcome

estimate_CTMLE0 <- function(data, true_value = 1, alpha = 0.05,
                            Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + "))),
                            gform = as.formula(paste("A ~", paste(
                              colnames(data)[grepl("^W", colnames(data))], collapse = " + ")))) {
  
  a <- min(data$Y)
  b <- max(data$Y)
  data$Y <- (data$Y - a) / (b - a)
  
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  
  # Outcome regression
  Q_init <- glm(Qform, family = gaussian, data = data)
  data$Q0W <- predict(Q_init, newdata = data.frame(A = 0, data[covariates]))
  data$Q1W <- predict(Q_init, newdata = data.frame(A = 1, data[covariates]))
  data$Q <- predict(Q_init, type = "response")
  
  # Treatment mechanism
  g_init <- glm(gform, family = binomial, data = data)
  data$g1W <- predict(g_init, type = "response")
  
  # derivative of the H wrt lambda (H̃_{gn,λk})
  eps <- 1e-4
  data$g1W_grad <- predict(update(g_init, data = data, start = coef(g_init) + eps), type = "response")
  data$H_grad <- ((1 - data$A) / (1 - data$g1W))^2 * (data$g1W_grad - data$g1W) + (data$A / data$g1W)^2 * (data$g1W_grad - data$g1W)
  
  data$H1 <- data$A 
  data$H0 <- (1 - data$A) 
  data$weights_H1 <-  1 / data$g1W 
  data$weights_H0 <-  1 / (1 - data$g1W) 
  
  if (length(covariates) == 1){
    second_stage_1 <- glm2(data$Y ~ -1 + H1 + offset(logit(data$Q)), family = quasibinomial, data = data, weights = weights_H1)
    second_stage_0 <- glm2(data$Y ~ -1 + H0 + offset(logit(data$Q)), family = quasibinomial, data = data, weights = weights_H0)
    epsilon_1 <- coef(second_stage_1)[1]
    epsilon_0 <- coef(second_stage_0)[1]
    Q_star1 <- expit(logit(data$Q1W) + epsilon_1 )
    Q_star0 <- expit(logit(data$Q0W) + epsilon_0 )
  } else {
    second_stage_1 <- glm2(data$Y ~ -1 + H1 + H_grad + offset(logit(data$Q)), family = quasibinomial, data = data, weights = weights_H1)
    second_stage_0 <- glm2(data$Y ~ -1 + H0 + H_grad + offset(logit(data$Q)), family = quasibinomial, data = data, weights = weights_H0)
    epsilon_1 <- coef(second_stage_1)[1]
    epsilon_grad_1 <- coef(second_stage_1)[2]
    epsilon_0 <- coef(second_stage_0)[1]
    epsilon_grad_0 <- coef(second_stage_0)[2]
    Q_star1 <- expit(logit(data$Q1W) + epsilon_1 + epsilon_grad_1 * data$H_grad)
    Q_star0 <- expit(logit(data$Q0W) + epsilon_0 + epsilon_grad_0 * data$H_grad)
  }
  
  
  ATE_pre <- mean(data$Q1W - data$Q0W) * (b - a)
  ATE <- mean(Q_star1 - Q_star0) * (b - a)
  ATE_n <- ATE/(b - a)
  
  IF <- (data$A / data$g1W ) * (data$Y - Q_star1) - ((1 - data$A) / (1 - data$g1W )) * (data$Y - Q_star0) + (Q_star1 - Q_star0) - ATE_n
  
  IF_1 <- mean((data$A / data$g1W ) * (data$Y - Q_star1) - ((1 - data$A) / (1 - data$g1W )) * (data$Y - Q_star0))
  IF_0 <- mean((Q_star1 - Q_star0) - ATE_n)
  SE <- sqrt(var(IF) / nrow(data))* (b - a)
  z <- qnorm(1 - alpha / 2)
  v <- var(IF)/nrow(data)
  CI_lower <- ATE - z * SE
  CI_upper <- ATE + z * SE
  #MSE <- bias^2 + v
  
  return(
    tibble(
      estimator = "C-TMLE0",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = paste("A ~",as.character(gform)[3]),
      est_init = ATE_pre,
      est = ATE,
      se = SE,
      ci_lower = CI_lower,
      ci_upper = CI_upper,
      bias = est - true_value,
      bias.to.se = bias / se
    )
  )
}

## C-TMLE estimator : LASSO/ C-TMLE 1
# Data: Tibble
# Alpha: Allowable type I error
# Qform
# gform
# Predefined CTMLE specs: Weighted, double_H, bounded_continuous_outcome
estimate_CTMLE1 <- function(data, true_value = 1, alpha = 0.05,
                            Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + ")))) {
  a <- min(data$Y)
  b <- max(data$Y)
  data$Y <- (data$Y - a) / (b - a)
  
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  X <- as.matrix(data[covariates])
  
  #initial outcome model
  Q_init <- glm(Qform, family = gaussian, data = data)
  data$Q0W <- predict(Q_init, newdata = data.frame(A = 0, data[covariates]))
  data$Q1W <- predict(Q_init, newdata = data.frame(A = 1, data[covariates]))
  data$Q <- predict(Q_init, type = "response")
  
  #LASSO regression estimation of propensity score
  lambda_seq <- 10^seq(-3, 1, length.out = 50)
  cv_lasso <- cv.glmnet(X, data$A, family = "binomial", alpha = 1, lambda = lambda_seq, nfolds = 5)
  
  best_lambda <- cv_lasso$lambda.min
  g_init <- glmnet(X, data$A, family = "binomial", alpha = 1, lambda = best_lambda)
  data$g1W <- predict(g_init, newx = X, type = "response")
  
  data$H1 <- data$A / data$g1W
  data$H0 <- (1 - data$A) / (1 - data$g1W)
  
  epsilon_converged <- FALSE
  max_iter <- 100
  iter <- 0
  epsilon_1 <- 0
  epsilon_0 <- 0
  
  while (!epsilon_converged && iter < max_iter) {
    second_stage_1 <- glm2(data$Y ~ -1 + H1 + offset(logit(data$Q)), family = quasibinomial, data = data)
    second_stage_0 <- glm2(data$Y ~ -1 + H0 + offset(logit(data$Q)), family = quasibinomial, data = data)
    
    new_epsilon_1 <- coef(second_stage_1)[1]
    new_epsilon_0 <- coef(second_stage_0)[1]
    
    if (abs(new_epsilon_1 - epsilon_1) < 1e-6 && abs(new_epsilon_0 - epsilon_0) < 1e-6) {
      epsilon_converged <- TRUE
    }
    
    epsilon_1 <- new_epsilon_1
    epsilon_0 <- new_epsilon_0
    iter <- iter + 1
  }
  
  Q_star1 <- expit(logit(data$Q1W) + epsilon_1)
  Q_star0 <- expit(logit(data$Q0W) + epsilon_0)
  
  ATE_pre <- mean(data$Q1W - data$Q0W) * (b - a)
  ATE <- mean(Q_star1 - Q_star0) * (b - a)
  ATE_n <- ATE / (b - a)
  
  IF <- (data$A / data$g1W) * (data$Y - Q_star1) - ((1 - data$A) / (1 - data$g1W)) * (data$Y - Q_star0) + (Q_star1 - Q_star0) - ATE_n
  SE <- sqrt(var(IF) / nrow(data)) * (b - a)
  z <- qnorm(1 - alpha / 2)
  CI_lower <- ATE - z * SE
  CI_upper <- ATE + z * SE
  
  return(
    tibble(
      estimator = "C-TMLE1",
      Qform = paste("Y ~", as.character(Qform)[3]),
      gform = NA,
      est_init = unname(ATE_pre),
      est = unname(ATE),
      se = unname(SE),
      ci_lower = unname(CI_lower),
      ci_upper = unname(CI_upper),
      bias = unname(ATE - true_value),
      bias.to.se = unname(bias / se)
    )
  )
}


estimate_CTMLE1_0 <- function(data, true_value = 1, alpha = 0.05,
                            Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + ")))) {
  a <- min(data$Y)
  b <- max(data$Y)
  data$Y <- (data$Y - a) / (b - a)
  
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  X <- as.matrix(data[covariates])
  
  # Step 1: Initial outcome model Q_n
  Q_init <- glm(Qform, family = gaussian, data = data)
  data$Q0W <- predict(Q_init, newdata = data.frame(A = 0, data[covariates]))
  data$Q1W <- predict(Q_init, newdata = data.frame(A = 1, data[covariates]))
  data$Q <- predict(Q_init, type = "response")
  
  # Step 2: Create sequence of propensity score estimators with different penalties
  K <- 50  # Number of lambda values
  lambda_seq <- 10^seq(1, -3, length.out = K)  # Monotonically decreasing
  
  # Get lambda_CV using cross-validation as recommended in text
  cv_lasso <- cv.glmnet(X, data$A, family = "binomial", alpha = 1, lambda = lambda_seq, nfolds = 5)
  lambda_CV <- cv_lasso$lambda.min
  
  # Filter lambda sequence to start from lambda_CV
  lambda_seq <- lambda_seq[lambda_seq <= lambda_CV]
  K <- length(lambda_seq)
  
  # Initialize storage for models and results
  g_models <- list()
  Q_star_models <- list()
  empirical_loss <- numeric(K)
  
  # Step 3: Build sequence of estimators Q*_n,λ for each g_n,λ
  for (k in 1:K) {
    lambda_k <- lambda_seq[k]
    
    # Fit propensity score model with current lambda
    g_models[[k]] <- glmnet(X, data$A, family = "binomial", alpha = 1, lambda = lambda_k)
    g1W_k <- predict(g_models[[k]], newx = X, type = "response")
    
    # Clever covariates based on current g
    H1_k <- data$A / g1W_k
    H0_k <- (1 - data$A) / (1 - g1W_k)
    
    # Targeted fluctuation for current lambda
    epsilon_converged <- FALSE
    max_iter <- 100
    iter <- 0
    epsilon_1 <- 0
    epsilon_0 <- 0
    
    while (!epsilon_converged && iter < max_iter) {
      second_stage_1 <- glm2(data$Y ~ -1 + H1_k + offset(logit(data$Q)), family = quasibinomial, data = data)
      second_stage_0 <- glm2(data$Y ~ -1 + H0_k + offset(logit(data$Q)), family = quasibinomial, data = data)
      
      new_epsilon_1 <- coef(second_stage_1)[1]
      new_epsilon_0 <- coef(second_stage_0)[1]
      
      if (abs(new_epsilon_1 - epsilon_1) < 1e-6 && abs(new_epsilon_0 - epsilon_0) < 1e-6) {
        epsilon_converged <- TRUE
      }
      
      epsilon_1 <- new_epsilon_1
      epsilon_0 <- new_epsilon_0
      iter <- iter + 1
    }
    
    # Store Q*_n,λ for current lambda
    Q_star1_k <- expit(logit(data$Q1W) + epsilon_1)
    Q_star0_k <- expit(logit(data$Q0W) + epsilon_0)
    Q_star_models[[k]] <- list(Q_star1 = Q_star1_k, Q_star0 = Q_star0_k, 
                               epsilon_1 = epsilon_1, epsilon_0 = epsilon_0)
    
    # Calculate empirical loss for CV selection
    predicted_Y <- data$Q + epsilon_1 * H1_k + epsilon_0 * H0_k
    empirical_loss[k] <- mean((data$Y - predicted_Y)^2)
  }
  
  # Step 4: Select the best lambda_ctmle by minimizing empirical loss
  best_index <- which.min(empirical_loss)
  lambda_ctmle <- lambda_seq[best_index]
  Q_n_ctmle <- Q_star_models[[best_index]]
  
  # Step 5: Fluctuate selected initial estimate with each g_n,λ for λK < λ < λctmle
  refinement_indices <- which(lambda_seq < lambda_ctmle)
  
  if (length(refinement_indices) > 0) {
    refined_Q_star_models <- list()
    refined_empirical_loss <- numeric(length(refinement_indices))
    
    for (j in 1:length(refinement_indices)) {
      idx <- refinement_indices[j]
      lambda_j <- lambda_seq[idx]
      
      # Get propensity score for this lambda
      g1W_j <- predict(g_models[[idx]], newx = X, type = "response")
      
      # Clever covariates
      H1_j <- data$A / g1W_j
      H0_j <- (1 - data$A) / (1 - g1W_j)
      
      # Targeted fluctuation with initial estimate Q_n_ctmle
      epsilon_converged <- FALSE
      max_iter <- 100
      iter <- 0
      epsilon_1 <- 0
      epsilon_0 <- 0
      
      while (!epsilon_converged && iter < max_iter) {
        second_stage_1 <- glm2(data$Y ~ -1 + H1_j + offset(logit(data$Q)), 
                               family = quasibinomial, data = data)
        second_stage_0 <- glm2(data$Y ~ -1 + H0_j + offset(logit(data$Q)), 
                               family = quasibinomial, data = data)
        
        new_epsilon_1 <- coef(second_stage_1)[1]
        new_epsilon_0 <- coef(second_stage_0)[1]
        
        if (abs(new_epsilon_1 - epsilon_1) < 1e-6 && abs(new_epsilon_0 - epsilon_0) < 1e-6) {
          epsilon_converged <- TRUE
        }
        
        epsilon_1 <- new_epsilon_1
        epsilon_0 <- new_epsilon_0
        iter <- iter + 1
      }
      
      # Store refined Q*_n,λ
      Q_star1_j <- expit(logit(data$Q1W) + epsilon_1)
      Q_star0_j <- expit(logit(data$Q0W) + epsilon_0)
      refined_Q_star_models[[j]] <- list(Q_star1 = Q_star1_j, Q_star0 = Q_star0_j,
                                         epsilon_1 = epsilon_1, epsilon_0 = epsilon_0)
      
      # Calculate empirical loss
      predicted_Y <- data$Q + epsilon_1 * H1_j + epsilon_0 * H0_j
      refined_empirical_loss[j] <- mean((data$Y - predicted_Y)^2)
    }
    
    # Step 6: Choose final estimate that minimizes empirical loss
    final_index <- which.min(refined_empirical_loss)
    final_model <- refined_Q_star_models[[final_index]]
  } else {
    # If there are no λ < λctmle, use the best model directly
    final_model <- Q_star_models[[best_index]]
  }
  
  # Extract final results
  Q_star1 <- final_model$Q_star1
  Q_star0 <- final_model$Q_star0
  
  # Transform ATE back to original scale
  ATE_pre <- mean(data$Q1W - data$Q0W) * (b - a)
  ATE <- mean(Q_star1 - Q_star0) * (b - a)
  ATE_n <- ATE / (b - a)
  
  # Get final g estimate for inference
  final_g_index <- ifelse(length(refinement_indices) > 0, 
                          refinement_indices[final_index], best_index)
  final_g1W <- predict(g_models[[final_g_index]], newx = X, type = "response")
  
  # Calculate influence function (solves critical equation 5 from the text)
  IF <- (data$A / final_g1W) * (data$Y - Q_star1) - 
    ((1 - data$A) / (1 - final_g1W)) * (data$Y - Q_star0) + 
    (Q_star1 - Q_star0) - ATE_n
  
  SE <- sqrt(var(IF) / nrow(data)) * (b - a)
  z <- qnorm(1 - alpha / 2)
  CI_lower <- ATE - z * SE
  CI_upper <- ATE + z * SE
  
  return(
    tibble(
      estimator = "New C-TMLE1",
      Qform = paste("Y ~", as.character(Qform)[3]),
      gform = NA,
      est_init = unname(ATE_pre),
      est = unname(ATE),
      se = unname(SE),
      ci_lower = unname(CI_lower),
      ci_upper = unname(CI_upper),
      bias = unname(ATE - true_value),
      bias.to.se = unname(bias / se)
    )
  )
}

