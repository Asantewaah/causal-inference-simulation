### List of high-dimensional estimator functions(glmnet-based TMLE and C-TMLE)

source(here::here("scripts", "01_data_generating_processes.R"))
source(here::here("scripts", "02_utilities.R"))

estimate_TMLE_glmnet <- function(data, alpha = 0.05, true_value = 1,
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
  X <- model.matrix(Qform, data = data)
  Q_init <- cv.glmnet(X, data$Y, alpha = 1)  # Using elastic net
  data$Q0W <- predict(Q_init, newx = model.matrix(as.formula(paste("~ A +", paste(covariates, collapse = " + "))), 
                                                  data = data.frame(A = 0, data[covariates])), s = "lambda.min")
  data$Q1W <- predict(Q_init, newx = model.matrix(as.formula(paste("~ A +", paste(covariates, collapse = " + "))), 
                                                  data = data.frame(A = 1, data[covariates])), s = "lambda.min")
  data$Q <- predict(Q_init, newx = X, s = "lambda.min")
  
  # Treatment mechanism
  X_g <- model.matrix(gform, data = data)
  g_model <- cv.glmnet(X_g, data$A, family = "binomial", alpha = 0.5)  # Using elastic net
  data$g1W <- predict(g_model, newx = X_g, type = "response", s = "lambda.min")
  
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
  
  #return(list(estimate = ATE, var = v, se = SE, mse = MSE, ci_lower = CI_lower, ci_upper = CI_upper))
  
  #hist(data$g1W, breaks = 30 , main = 'GLMNET TMLE Propensity Score')
  
  return(
    tibble(
      estimator = "GLMNET TMLE",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = paste("A ~",as.character(gform)[3]),
      est_init = unname(ATE_pre),
      est = unname(ATE),
      se = unname(SE),
      ci_lower = unname(CI_lower),
      ci_upper = unname(CI_upper),
      bias = unname(est - true_value),
      bias.to.se = unname(bias / se)
    )
  )
}


estimate_GREEDY_CTMLE_glmnet <- function(data, alpha = 0.05, true_value = 1,
                                  Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + ")))) {
  
  a <- min(data$Y)
  b <- max(data$Y)
  data$Y <- (data$Y - a) / (b - a)
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  
  # Outcome regression
  X <- model.matrix(Qform, data = data)
  Q_init <- cv.glmnet(X, data$Y, alpha = 1)  
  data$Q0W <- predict(Q_init, newx = model.matrix(as.formula(paste("~ A +", paste(covariates, collapse = " + "))), 
                                                  data = data.frame(A = 0, data[covariates])), s = "lambda.min")
  data$Q1W <- predict(Q_init, newx = model.matrix(as.formula(paste("~ A +", paste(covariates, collapse = " + "))), 
                                                  data = data.frame(A = 1, data[covariates])), s = "lambda.min")
  data$Q <- predict(Q_init, newx = X, s = "lambda.min")
  
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
      X_g <- model.matrix(as.formula(paste("~ ", paste(vars_to_use, collapse = " + "))), data = data)
      g_model <- cv.glmnet(X_g, data$A, family = "binomial", alpha = 1)  
      g1W <- predict(g_model, newx = X_g, type = "response", s = "lambda.min")
      
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
  
  #hist(best_g1W, breaks = 30 , main = 'GLMNET GREEDY C-TMLE Propensity Score')
  
  return(
    tibble(
      estimator = "C-TMLE: GREEDY",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = ifelse(length(selected_vars) > 0, 
                     paste("A ~", paste(best_combination, collapse = " + ")), "A ~ 1"),
      est_init = unname(ATE_pre),
      est = unname(ATE),
      se = unname(SE),
      ci_lower = unname(CI_lower),
      ci_upper = unname(CI_upper),
      bias = unname(est - true_value),
      bias.to.se = unname(bias / se)
    )
  )
}
estimate_CTMLE1_glmnet <- function(data, true_value = 1, alpha = 0.05,
                                   Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + ")))) {
  a <- min(data$Y)
  b <- max(data$Y)
  data$Y <- (data$Y - a) / (b - a)
  
  V <- 10
  
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  X <- as.matrix(data[covariates])
  n <- nrow(data)
  
  # the initial outcome model Q_n
  Q_X <- model.matrix(Qform, data = data)
  Q_init <- cv.glmnet(Q_X, data$Y, alpha = 1)
  data$Q0W <- predict(Q_init, newx = model.matrix(~ A + ., data = data.frame(A = 0, data[covariates])), s = "lambda.min")
  data$Q1W <- predict(Q_init, newx = model.matrix(~ A + ., data = data.frame(A = 1, data[covariates])), s = "lambda.min")
  data$Q <- predict(Q_init, newx = Q_X, s = "lambda.min")
  
  # Propensity score model: build lambda sequence starting at lambda.min
  ps_cv <- cv.glmnet(X, data$A, family = "binomial", alpha = 1)
  lambda_cv <- ps_cv$lambda.min
  
  #lambda sequence (from λ_CV to more regularized)
  lambda_seq <- ps_cv$glmnet.fit$lambda[ps_cv$glmnet.fit$lambda >= lambda_cv]
  lambda_seq <- sort(lambda_seq, decreasing=TRUE)
  
  #cross-validation folds
  folds <- caret::createFolds(data$A, k=V, list=TRUE)
  
  
  #collaborative model selection
  best_lambda <- NULL
  best_loss <- Inf
  current_Q <- list(Q=data$Q, Q0W=data$Q0W, Q1W=data$Q1W)
  
  for(lambda in lambda_seq) {
    
    cv_loss <- 0
    
    for(v in 1:V){ 
      
      train_idx <- setdiff(1:n, folds[[v]])
      valid_idx <- folds[[v]]
      
      g_fit <- glmnet(X[train_idx,], data$A[train_idx], 
                      family="binomial", lambda=lambda)
      g1W_valid <- as.numeric(predict(g_fit, X[valid_idx,], type="response"))
      
      
      # Collaborative targeting on training data
      Q_valid <- current_Q$Q[valid_idx]
      Q_valid0 <- current_Q$Q0W[valid_idx]
      Q_valid1 <- current_Q$Q1W[valid_idx]
      
      data$H1_k <- data$A[valid_idx] 
      data$H0_k <- (1 - data$A[valid_idx]) 
      data$weights_H1_k <-  1 / g1W_valid 
      data$weights_H0_k <-  1 / (1 - g1W_valid)
      
      # Targeted fluctuation (one-step update)
      second_stage_1 <- glm2(data$Y[valid_idx] ~ -1 + H1_k[valid_idx] + offset(logit(Q_valid1)), family = quasibinomial, data = data, weights = weights_H1_k[valid_idx])
      second_stage_0 <- glm2(data$Y[valid_idx] ~ -1 + H0_k[valid_idx] + offset(logit(Q_valid0)), family = quasibinomial, data = data, weights = weights_H0_k[valid_idx])
      epsilon_1 <- coef(second_stage_1)[1]
      epsilon_0 <- coef(second_stage_0)[1]
      
      Q_star1_k <- expit(logit(Q_valid1) + epsilon_1)
      Q_star0_k <- expit(logit(Q_valid0) + epsilon_0)
      
      
      # loss
      loss <-mean((data$Y[valid_idx] - (data$A[valid_idx] * Q_star1_k + (1 - data$A[valid_idx]) * Q_star0_k))^2)
      cv_loss <- cv_loss + loss
    }
    
    # Update best lambda
    if(cv_loss < best_loss) {
      best_loss <- cv_loss
      best_lambda <- lambda
    }
  }
  
  
  # Select best lambda by minimizing empirical loss (collaborative selection)
  final_g_model <-  glmnet(X, data$A, family="binomial", lambda=best_lambda)
  final_g1W <- as.numeric(predict(final_g_model, X, type="response"))
  
  
  #final iterative step
  
  Q_final1 <- current_Q$Q1W
  Q_final0 <- current_Q$Q0W
  
  data$H1_j <- data$A 
  data$H0_j <- (1 - data$A) 
  data$weights_H1_j <-  1 / final_g1W 
  data$weights_H0_j <-  1 / (1 - final_g1W)
  
  # Targeted fluctuation (one-step update)
  second_stage_1 <- glm2(data$Y ~ -1 + H1_j + offset(logit(Q_final1)), family = quasibinomial, data = data, weights = weights_H1_j)
  second_stage_0 <- glm2(data$Y ~ -1 + H0_j + offset(logit(Q_final0)), family = quasibinomial, data = data, weights = weights_H0_j)
  epsilon_1 <- coef(second_stage_1)[1]
  epsilon_0 <- coef(second_stage_0)[1]
  
  Q_star1 <- expit(logit(Q_final1) + epsilon_1)
  Q_star0 <- expit(logit(Q_final0) + epsilon_0)
  
  
  # Estimate ATE and inference
  ATE_pre <- mean(data$Q1W - data$Q0W) * (b - a)
  ATE <- mean(Q_star1 - Q_star0) * (b - a)
  ATE_n <- ATE / (b - a)
  
  IF <- (data$A / final_g1W) * (data$Y - Q_star1) - 
    ((1 - data$A) / (1 - final_g1W)) * (data$Y - Q_star0) + 
    (Q_star1 - Q_star0) - ATE_n
  SE <- sqrt(var(IF) / nrow(data)) * (b - a)
  z <- qnorm(1 - alpha / 2)
  CI_lower <- ATE - z * SE
  CI_upper <- ATE + z * SE
  
  #hist(final_g1W, breaks = 30, main = 'GLMNET C-TMLE1 (Lasso & Collaborative Selection) Propensity Score')
  
  tibble::tibble(
    estimator = "GLMNET C-TMLE1",
    Qform = paste("Y ~", as.character(Qform)[3]),
    gform = NA,
    est_init = unname(ATE_pre),
    est = unname(ATE),
    se = unname(SE),
    ci_lower = unname(CI_lower),
    ci_upper = unname(CI_upper),
    bias = unname(ATE - true_value),
    bias.to.se = unname((ATE - true_value) / SE)
  )
}



estimate_CTMLE0_glmnet <- function(data, true_value = 1, alpha = 0.05,
                            Qform = as.formula(paste("Y ~ A +", paste(covariates, collapse = " + "))),
                            gform = as.formula(paste("A ~", paste(
                              colnames(data)[grepl("^W", colnames(data))], collapse = " + ")))) {
  
  a <- min(data$Y)
  b <- max(data$Y)
  data$Y <- (data$Y - a) / (b - a)
  
  covariates <- colnames(data)[grepl("^W", colnames(data))]
  
  # Outcome regression
  X <- model.matrix(Qform, data = data)
  Q_init <- cv.glmnet(X, data$Y, alpha = 1)  # Using elastic net
  data$Q0W <- predict(Q_init, newx = model.matrix(as.formula(paste("~ A +", paste(covariates, collapse = " + "))), 
                                                  data = data.frame(A = 0, data[covariates])), s = "lambda.min")
  data$Q1W <- predict(Q_init, newx = model.matrix(as.formula(paste("~ A +", paste(covariates, collapse = " + "))), 
                                                  data = data.frame(A = 1, data[covariates])), s = "lambda.min")
  data$Q <- predict(Q_init, newx = X, s = "lambda.min")
  
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
  
  #hist(data$g1W, breaks = 30, main = 'GLMNET C-TMLE0 Propensity Score')
  
  return(
    tibble(
      estimator = "GLMNET C-TMLE0",
      Qform = paste("Y ~",as.character(Qform)[3]),
      gform = paste("A ~",as.character(gform)[3]),
      est_init = unname(ATE_pre),
      est = unname(ATE),
      se = unname(SE),
      ci_lower = unname(CI_lower),
      ci_upper = unname(CI_upper),
      bias = unname(est - true_value),
      bias.to.se = unname(bias / se)
    )
  )
}


