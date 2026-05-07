t<-proc.time()
# Needed libraries
library(survival)
library(survminer)
library(MASS)
library(e1071)
library(pROC)
library(caTools)
library(caret)
library(rpart)
library(neuralnet)
library(randomForest)
library(randomForestSRC)
library(xgboost)
library(splines)
library(readr)
library(timeROC)
library(rpart.plot)

# --- Function: smsurv (Weighted Breslow Estimator) ---
smsurv <- function(Time, Status, X, beta, w, model) {
  death_point <- sort(unique(subset(Time, Status == 1)))
  
  if (length(death_point) == 0) {
    return(list(survival = rep(1, length(Time)), death_times = NULL, s0_vals = NULL))
  }
  
  if(is.null(beta) || length(beta) < (ncol(X) - 1)) {
    beta <- c(beta, rep(0, (ncol(X) - 1) - length(beta)))
  }
  
  .eps <- 1e-10
  
  if (model == 'ph') {
    coxexp <- exp(X[, -1, drop = FALSE] %*% beta)
  }
  
  lambda <- numeric(length(death_point))
  event <- numeric(length(death_point))
  
  for (i in 1:length(death_point)) {
    event[i] <- sum(Status * as.numeric(Time == death_point[i]))
    if (model == 'ph') temp <- sum(as.numeric(Time >= death_point[i]) * w * drop(coxexp))
    if (model == 'aft') temp <- sum(as.numeric(Time >= death_point[i]) * w)
    
    if (!is.na(temp) && temp > .eps) {
      lambda[i] <- event[i] / temp
    } else {
      lambda[i] <- 0
    }
  }
  
  s0_at_deaths <- exp(-cumsum(lambda))
  
  survival <- sapply(Time, function(t) {
    if (t < min(death_point)) return(1)
    s0_at_deaths[max(which(death_point <= t))]
  })
  
  survival <- pmin(pmax(survival, .eps), 1 - .eps)
  
  return(list(survival = survival, death_times = death_point, s0_vals = s0_at_deaths))
}




# LOGIT 
em.Logit.Pois <- function(Time, Status, Time1, Status1,
                          X, X1, Z, Z1,
                          b, beta, s0, s01,
                          emmax, eps) {
  
  n <- length(Status)
  s  <- s0
  s1 <- s01
  
  convergence <- 1000
  i <- 1
  
  while (convergence > eps && i <= emmax) {
    UN <- matrix(exp(Z %*% b)/(1+exp(Z %*% b)),ncol=1)
    PRED <- matrix(exp(Z1 %*% b)/(1+exp(Z1 %*% b)),ncol=1)
    
    
    survival  <- drop(s^(exp(X[, -1, drop=FALSE] %*% beta)))
    survival1 <- drop(s1^(exp(X1[, -1, drop=FALSE] %*% beta)))
    
    # E-STEP: Using matched probabilities for Train (UN) and Test (PRED)
    M <- Status - (survival * log(pmax(1 - UN, 1e-6)))
    M1 <- Status1 - (survival1 * log(pmax(1 - PRED, 1e-6)))
    
    # Incidence Part (M-STEP)
    Q1 <- function(par) {
      bb <- par
      u_temp <- matrix(exp(Z %*% bb)/(1+exp(Z %*% bb)),ncol=1)
      loglik <- sum(M*log(-log(pmax(1-u_temp, 1e-10)))) + sum(log(pmax(1-u_temp, 1e-10)))
      return(-loglik)
    }
    
    update_b = optim(par=b,fn=Q1,method="Nelder-Mead")$par 
    
    # Latency Part (M-STEP)
    fit_w <- coxph(Surv(Time, Status) ~ X[, -1, drop = FALSE] + offset(log(pmax(M, 1e-4))),
                   subset = M != 0, method = "breslow")
    update_beta <- fit_w$coef
    
    update_s <- smsurv(Time, Status, X, update_beta, w = M, model = "ph")$survival
    update_s1 <- smsurv(Time1, Status1, X1, update_beta, w = M1, model = "ph")$survival
    
    convergence <- sum((update_b - b)^2) + sum((update_beta - beta)^2)
    
    b <- update_b
    beta <- update_beta 
    s <- update_s
    s1 <- update_s1
    
    i <- i+1
  }
  
  # Final calculation for returns
  UN <- matrix(exp(Z %*% b)/(1+exp(Z %*% b)),ncol=1)
  PRED <- matrix(exp(Z1 %*% b)/(1+exp(Z1 %*% b)),ncol=1)
  survival  <- drop(s^(exp(X[, -1, drop=FALSE] %*% beta)))
  survival1 <- drop(s1^(exp(X1[, -1, drop=FALSE] %*% beta)))
  Sp <- (1 - UN)^(1 - survival)
  Sp.pred <- (1 - PRED)^(1 - survival1)
  S1 = (Sp-(1-UN))/UN 
  S1.pred = (Sp.pred-(1-PRED))/PRED
  
  return(list(
    b = b,
    latencyfit = beta,
    UN = UN,
    PRED = PRED,
    Sp = Sp,
    Sp.pred = Sp.pred,
    S1 = S1,
    S1.pred = S1.pred,
    s0 = s,
    s01 = s1,
    S = survival,
    S.pred = survival1,
    tau = convergence
  ))
}

smcure.Logit.Pois <- function(train, test,  Var = F, emmax = 1000, eps = 1e-3, nboot = 100){
  
  Time <- train$t; Status <- train$d
  Time1 <- test$t; Status1 <- test$d
  
  X <- model.matrix(~ x2 + x3 + x4 + x6, data = train)
  X1 <- model.matrix(~ x2 + x3 + x4 + x6, data = test)
  Z <- model.matrix(~ x2 + x3 + x4 + x6, data = train)
  Z1 <- model.matrix(~ x2 + x3 + x4 + x6, data = test)
  
  #Scaling age and tumor thickness
  cols_scale <- colnames(X) %in% c("x2", "x6")
  scale_params <- list(
    x2 = list(mean = mean(X[, "x2"]), sd = sd(X[, "x2"])),
    x6 = list(mean = mean(X[, "x6"]), sd = sd(X[, "x6"]))
  )
  X[, cols_scale] <- scale(X[, cols_scale])
  X1[, "x2"] <- (X1[, "x2"] - scale_params$x2$mean) / scale_params$x2$sd
  X1[, "x6"] <- (X1[, "x6"] - scale_params$x6$mean) / scale_params$x6$sd
  
  Z <- X; Z1 <- X1
  
  
  
  coxfit_train <- coxph(Surv(Time, Status) ~ 1, data = train)
  coxfit_test <- coxph(Surv(Time1, Status1) ~ 1, data = test)
  
  out.data <- basehaz(coxfit_train, centered = FALSE)  # columns: hazard, time
  out.data1 <- basehaz( coxfit_test, centered = FALSE)  # columns: hazard, time
  t_grid  <- out.data$time
  t_grid1  <- out.data1$time
  idx_tr <- pmax(1, findInterval(Time,  t_grid))
  idx_te <- pmax(1, findInterval(Time1, t_grid1))
  S0_grid <- exp(-out.data$hazard)
  S01_grid <- exp(-out.data1$hazard)
  
  
  # initial full versions
  s0_full  <- exp(-out.data[,1])
  s01_full <- exp(-out.data1[,1])
  
  # fallback indexed versions
  s0_idx  <- S0_grid[idx_tr]
  s01_idx <- S01_grid[idx_te]
  
  # choose based on length match
  if (length(s0_full) == length(Time)) {
    s0_init <- s0_full
  } else {
    s0_init <- s0_idx
  }
  
  if (length(s01_full) == length(Time1)) {
    s01_init <- s01_full
  } else {
    s01_init <- s01_idx
  } 
  
  
  
  
  
  
  if(any(is.na(s0_init))) s0_init[is.na(s0_init)] <- min(s0_init, na.rm = TRUE)
  if(any(is.na(s01_init))) s01_init[is.na(s01_init)] <- min(s01_init, na.rm = TRUE)
  
  
  b_start <- rep(0, ncol(Z))
  
  
  beta_start <- coxfit_train$coefficients
  if (is.null(beta_start) || length(beta_start) != (ncol(X) - 1)) {
    beta_start <- rep(0, ncol(X) - 1)
  } else {
    beta_start <- as.numeric(beta_start)
  }
  
  
  emfit <- em.Logit.Pois(Time, Status, Time1, Status1,
                         X, X1, Z, Z1,
                         b_start, beta_start, s0_init, s01_init,
                         emmax, eps) 
  b.est <- emfit$b
  beta.est <- emfit$latencyfit
  
  if (Var) {
    n_latency <- length(beta.est)
    n_incidence <- length(b.est)
    
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = n_latency)
    incidence_boot <- matrix(NA_real_, nrow = nboot, ncol = n_incidence)
    
    cat("Starting bootstrap (", nboot, " iterations)...\n")
    
    for (i in 1:nboot) {
      boot_idx <- sample(seq_along(Status), length(Status), replace = TRUE)
      
      Time_b   <- Time[boot_idx]
      Status_b <- Status[boot_idx]
      X_b      <- X[boot_idx, , drop = FALSE]
      Z_b      <- Z[boot_idx, , drop = FALSE]
      
      boot_train_df <- data.frame(
        t  = Time_b,
        d  = Status_b,
        x2 = X_b[, "x2"],
        x3 = X_b[, "x3"],
        x4 = X_b[, "x4"],
        x6 = X_b[, "x6"]
      )
      
      coxfit_b <- try(
        coxph(Surv(t, d) ~ x2 + x3 + x4 + x6, data = boot_train_df),
        silent = TRUE
      )
      if (inherits(coxfit_b, "try-error")) next
      
      bh_b <- basehaz(coxfit_b, centered = FALSE)
      S0_b_grid <- exp(-bh_b$hazard)
      t_grid_b  <- bh_b$time
      
      idx_b <- pmax(1, findInterval(Time_b, t_grid_b))
      s0_b  <- S0_b_grid[idx_b]
      
      if (any(is.na(s0_b))) {
        s0_b[is.na(s0_b)] <- min(s0_b, na.rm = TRUE)
      }
      
      bootfit <- try(
        em.Logit.Pois(Time_b, Status_b, Time1, Status1,
                      X_b, X1, Z_b, Z1,
                      b.est, beta.est, s0_b, s01_init, emmax, eps),
        silent = TRUE
      )
      
      if (!inherits(bootfit, "try-error")) {
        latency_boot[i, ] <- bootfit$latencyfit
        incidence_boot[i, ] <- bootfit$b
      }
    }
    
    emfit$latency_se <- apply(latency_boot, 2, sd, na.rm = TRUE)
    emfit$incidence_se <- apply(incidence_boot, 2, sd, na.rm = TRUE)
    
    emfit$latency_p <- 2 * (1 - pnorm(abs(emfit$latencyfit / emfit$latency_se)))
    emfit$incidence_p <- 2 * (1 - pnorm(abs(emfit$b / emfit$incidence_se)))
  }
  
  return(emfit)
}


# SPLINE
em.Spline.Pois <- function(Time, Status, Time1, Status1,
                           X, X1, Z, Z1,
                           b, beta, s0, s01,
                           emmax, eps) {
  
  n <- length(Status)
  s  <- s0
  s1 <- s01
  
  convergence <- 1000
  i <- 1
  
  while (convergence > eps && i <= emmax) {
    UN <- matrix(exp(Z %*% b)/(1+exp(Z %*% b)),ncol=1)
    PRED <- matrix(exp(Z1 %*% b)/(1+exp(Z1 %*% b)),ncol=1)
    
    if(is.null(beta) || length(beta) != (ncol(X) - 1)) {
      beta <- rep(0, ncol(X) - 1)
    }
    
    survival  <- drop(s^(exp(X[, -1, drop=FALSE] %*% beta)))
    survival1 <- drop(s1^(exp(X1[, -1, drop=FALSE] %*% beta)))
    
    # E-STEP: Using matched probabilities for Train (UN) and Test (PRED)
    M <- Status - (survival * log(pmax(1 - UN, 1e-6)))
    M1 <- Status1 - (survival1 * log(pmax(1 - PRED, 1e-6)))
    
    # Incidence Part (M-STEP)
    Q1 <- function(par) {
      bb <- par
      u_temp <- matrix(exp(Z %*% bb)/(1+exp(Z %*% bb)),ncol=1)
      loglik <- sum(M*log(-log(pmax(1-u_temp, 1e-10)))) + sum(log(pmax(1-u_temp, 1e-10)))
      return(-loglik)
    }
    
    update_b = optim(par=b,fn=Q1,method="Nelder-Mead")$par 
    
    # Latency Part (M-STEP)
    fit_w <- coxph(Surv(Time, Status) ~ X[, -1, drop = FALSE] + offset(log(pmax(M, 1e-4))),
                   subset = M > 0, method = "breslow")
    update_beta <- fit_w$coef
    
    update_s <- smsurv(Time, Status, X, update_beta, w = M, model = "ph")$survival
    update_s1 <- smsurv(Time1, Status1, X1, update_beta, w = M1, model = "ph")$survival
    
    convergence <- sum((update_b - b)^2) + sum((update_beta - beta)^2)
    
    b <- update_b
    beta <- update_beta 
    s <- update_s
    s1 <- update_s1
    
    i <- i+1
  }
  
  # Final calculation for returns
  UN <- matrix(exp(Z %*% b)/(1+exp(Z %*% b)),ncol=1)
  PRED <- matrix(exp(Z1 %*% b)/(1+exp(Z1 %*% b)),ncol=1)
  survival  <- drop(s^(exp(X[, -1, drop=FALSE] %*% beta)))
  survival1 <- drop(s1^(exp(X1[, -1, drop=FALSE] %*% beta)))
  Sp <- (1 - UN)^(1 - survival)
  Sp.pred <- (1 - PRED)^(1 - survival1)
  S1 = (Sp-(1-UN))/UN 
  S1.pred = (Sp.pred-(1-PRED))/PRED
  
  return(list(
    b = b,
    latencyfit = beta,
    UN = UN,
    PRED = PRED,
    Sp = Sp,
    Sp.pred = Sp.pred,
    S1 = S1,
    S1.pred = S1.pred,
    s0 = s,
    s01 = s1,
    S = survival,
    S.pred = survival1,
    tau = convergence
  ))
}

smcure.Spline.Pois <- function(train, test,  Var = TRUE, emmax = 1000, eps = 1e-3, nboot = 100){
  
  Time <- train$t; Status <- train$d
  Time1 <- test$t; Status1 <- test$d
  
  X <- model.matrix(~ x2 + x3 + x4 + x6, data = train)
  X1 <- model.matrix(~ x2 + x3 + x4 + x6, data = test)
  
  #Scaling age and tumor thickness
  cols_scale <- colnames(X) %in% c("x2", "x6")
  scale_params <- list(
    x2 = list(mean = mean(X[, "x2"]), sd = sd(X[, "x2"])),
    x6 = list(mean = mean(X[, "x6"]), sd = sd(X[, "x6"]))
  )
  X[, cols_scale] <- scale(X[, cols_scale])
  X1[, "x2"] <- (X1[, "x2"] - scale_params$x2$mean) / scale_params$x2$sd
  X1[, "x6"] <- (X1[, "x6"] - scale_params$x6$mean) / scale_params$x6$sd
  
  Z <- X; Z1 <- X1
  
  
  
  coxfit_train <- coxph(Surv(Time, Status) ~ 1, data = train)
  coxfit_test <- coxph(Surv(Time1, Status1) ~ 1, data = test)
  
  out.data <- basehaz(coxfit_train, centered = FALSE)  # columns: hazard, time
  out.data1 <- basehaz( coxfit_test, centered = FALSE)  # columns: hazard, time
  t_grid  <- out.data$time
  t_grid1  <- out.data1$time
  idx_tr <- pmax(1, findInterval(Time,  t_grid))
  idx_te <- pmax(1, findInterval(Time1, t_grid1))
  S0_grid <- exp(-out.data$hazard)
  S01_grid <- exp(-out.data1$hazard)
  
  
  # initial full versions
  s0_full  <- exp(-out.data[,1])
  s01_full <- exp(-out.data1[,1])
  
  # fallback indexed versions
  s0_idx  <- S0_grid[idx_tr]
  s01_idx <- S01_grid[idx_te]
  
  # choose based on length match
  if (length(s0_full) == length(Time)) {
    s0_init <- s0_full
  } else {
    s0_init <- s0_idx
  }
  
  if (length(s01_full) == length(Time1)) {
    s01_init <- s01_full
  } else {
    s01_init <- s01_idx
  } 
  
  
  
  
  if(any(is.na(s0_init))) s0_init[is.na(s0_init)] <- min(s0_init, na.rm = TRUE)
  if(any(is.na(s01_init))) s01_init[is.na(s01_init)] <- min(s01_init, na.rm = TRUE)
  
  # --- FIXED Spline CONSTRUCTION ---
  # We only apply Splines to continuous variables x2 and x6.
  # Categorical dummies (x3, x4) remain linear.
  
  # df = 2 is the minimum for a natural Spline to be valid.
  df_cont <- 2
  
  # Basis for Training
  ns_x2 <- splines::ns(train$x2, df = df_cont)
  ns_x6 <- splines::ns(train$x6, df = df_cont)
  
   Z <- cbind(
    1,
    ns_x2,
    ns_x6,
    train$x3,
    train$x4
  )
  
  # Basis for Testing
  Z1 <- cbind(
    1,
    predict(ns_x2, test$x2),
    predict(ns_x6, test$x6),
    test$x3,
    test$x4
  )
  
  b_start <- rep(0, ncol(Z))
  beta_start <- coxfit_train$coefficients 
  
  emfit <- em.Spline.Pois(Time, Status, Time1, Status1, X, X1, Z, Z1, b_start, beta_start, s0_init, s01_init, emmax, eps)
  
  if (Var) {
    cat("Starting bootstrap (", nboot, " iterations)...\n")
    
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = length(emfit$latencyfit))
    incidence_boot <- matrix(NA_real_, nrow = nboot, ncol = length(emfit$b))
    
    for (i in 1:nboot) {
      boot_idx <- sample(seq_along(Status), length(Status), replace = TRUE)
      
      Time_b   <- Time[boot_idx]
      Status_b <- Status[boot_idx]
      X_b <- X[boot_idx, , drop = FALSE]
      Z_b <- Z[boot_idx, , drop = FALSE]
      
      boot_train_df <- data.frame(
        t  = Time_b,
        d  = Status_b,
        x2 = X_b[, "x2"],
        x3 = X_b[, "x3"],
        x4 = X_b[, "x4"],
        x6 = X_b[, "x6"]
      )
      
      coxfit_b <- try(
        coxph(Surv(t, d) ~ x2 + x3 + x4 + x6, data = boot_train_df),
        silent = TRUE
      )
      if (inherits(coxfit_b, "try-error")) next
      
      bh_b <- basehaz(coxfit_b, centered = FALSE)
      S0_b_grid <- exp(-bh_b$hazard)
      t_grid_b  <- bh_b$time
      
      idx_b <- pmax(1, findInterval(Time_b, t_grid_b))
      s0_b  <- S0_b_grid[idx_b]
      
      if (any(is.na(s0_b))) {
        s0_b[is.na(s0_b)] <- min(s0_b, na.rm = TRUE)
      }
      
      bootfit <- try(
        em.Spline.Pois(Time_b, Status_b, Time1, Status1,
                       X_b, X1, Z_b, Z1,
                       emfit$b, emfit$latencyfit,
                       s0_b, s01_init, emmax, eps),
        silent = TRUE
      )
      
      if (!inherits(bootfit, "try-error")) {
        latency_boot[i, ] <- bootfit$latencyfit
        incidence_boot[i, ] <- bootfit$b
      }
    }
    
    emfit$latency_se <- apply(latency_boot, 2, sd, na.rm = TRUE)
    emfit$incidence_se <- apply(incidence_boot, 2, sd, na.rm = TRUE)
    
    emfit$latency_p <- 2 * (1 - pnorm(abs(emfit$latencyfit / emfit$latency_se)))
    emfit$incidence_p <- 2 * (1 - pnorm(abs(emfit$b / emfit$incidence_se)))
  }
  
  
  Status_tr <- as.integer(Status)
  Status_te <- as.integer(Status1)
  
  UN <- as.numeric(emfit$UN)
  PRED <- as.numeric(emfit$PRED)
  
  S <- as.numeric(emfit$S)       # susceptible survival component at train times
  S_pred <- as.numeric(emfit$S.pred)  # susceptible survival component at test times
  
  w_post_tr <- Status_tr + ((1 - Status_tr) * (1-(1-UN)^S))
  w_post_te <- Status_te + (1 - Status_te) *  (1-(1-PRED)^S_pred)
  
  w_post_tr <- pmin(pmax(w_post_tr, 1e-6), 1 - 1e-6)
  w_post <- pmin(pmax(w_post_te, 1e-6), 1 - 1e-6)
  
  emfit$w_post_tr <- w_post_tr
  emfit$w_post <- w_post
  return(emfit)
}




# DECISION TREE
em.DT.Pois <- function(Time, Status, Time1, Status1,
                       X, X1, Z, Z1,
                       beta, s0, s01,
                       uncureprob, uncurepred,
                       emmax, eps, best_params) {
  
  n <- length(Status)
  m <- length(Status1)
  s <- s0
  s1 <- s01
  UN   <- uncureprob
  PRED <- uncurepred
  
  
  
  convergence <- 1000
  i <- 1
  
  while (convergence > eps && i <= emmax) {
    
    if(is.null(beta) || length(beta) != (ncol(X) - 1)) {
      beta <- rep(0, ncol(X) - 1)
    }
    
    # Latency survival update
    survival  <- drop(s^(exp(X[, -1, drop = FALSE] %*% beta)))
    survival1 <- drop(s1^(exp(X1[, -1, drop = FALSE] %*% beta)))
    
    # E-STEP: Calculate posterior probability of being uncured
    w_prob <- Status + (1-Status)*(1-((1-UN)^(survival)))
    w_prob <- pmin(pmax(w_prob, 1e-6), 1 - 1e-6)
    M  <- Status - (survival * log(pmax(1 - UN, 1e-6)))
    M1 <- Status1 - (survival1 * log(pmax(1 - PRED, 1e-6)))
    
    # M-STEP: Incidence Part using Data Augmentation (K=5)
    K <- 5
    V_matrix <- matrix(rbinom(n * K, size = 1, prob = rep(w_prob, each = K)), nrow = n, byrow = TRUE)
    
    cure_preds <- matrix(NA, nrow = n, ncol = K)
    pred_preds <- matrix(NA, nrow = m, ncol = K)
    for (k in 1:K) {
      yk <- as.factor(V_matrix[, k])
      yk <- factor(V_matrix[, k], levels = c(0,1), labels = c("cured","uncured"))
      
      mod_data <- data.frame(Z[, -1, drop = FALSE])
      mod_data$yk <- yk
      
      mod <- rpart(yk ~ ., data = mod_data, method = "class", control = rpart.control(
        cp = best_params$cp,
        minsplit = best_params$minsplit,
        maxdepth = best_params$maxdepth
      ))
      
      probs_train <- predict(mod, newdata = as.data.frame(Z[, -1, drop = FALSE]), type = "prob")
      probs_test <- predict(mod, newdata = as.data.frame(Z1[, -1, drop = FALSE]), type = "prob")
      
      cure_preds[, k] <- probs_train[, "uncured"]
      pred_preds[, k] <- probs_test[, "uncured"]
    }
    
    update_cureb <- rowMeans(cure_preds, na.rm = TRUE)
    update_pred  <- rowMeans(pred_preds, na.rm = TRUE)
    # M-step for beta (Latency)
    fit_w <- coxph(Surv(Time, Status) ~ X[, -1, drop = FALSE] + offset(log(pmax(M, 1e-4))), 
                   subset = M != 0, method = "breslow")
    update_beta <- fit_w$coef
    
    # Update baseline survival
    update_s  <- smsurv(Time, Status, X, beta, w = M, model = "ph")$survival
    update_s1 <- smsurv(Time1, Status1, X1, beta, w = M1, model = "ph")$survival
    
    convergence <- sum(c(update_beta-beta,mean(update_cureb)-mean(UN),mean(update_s)-mean(s))^2)
    
    UN   <- update_cureb
    PRED <- update_pred
    beta <- update_beta
    s    <- update_s
    s1   <- update_s1
    i    <- i + 1
  }
  
  Sp      = (1-UN)^(1-s)
  Sp.pred = (1-PRED)^(1-s1) 
  S1      = (Sp-(1-UN))/pmax(UN, 1e-9)
  S1.pred = (Sp.pred-(1-PRED))/pmax(PRED, 1e-9)
  
  return(list(latencyfit = beta, UN = UN, PRED = PRED,
              Sp = Sp, Sp.pred = Sp.pred, S1 = S1, S1.pred = S1.pred, S.pred = s1,
              s0 = s, S = s, tau = convergence,best_params = best_params))
}

smcure.DT.Pois <- function(train, test, Var = TRUE, emmax = 1000, eps = 1e-3, nboot = 100) {
  Time <- train$t; Status <- train$d
  Time1 <- test$t; Status1 <- test$d
  X <- model.matrix(~ x2 + x3 + x4 + x6, data = train)
  X1 <- model.matrix(~ x2 + x3 + x4 + x6, data = test)
  # Standardize
  cols_scale <- colnames(X) %in% c("x2", "x6")
  train_scaled <- scale(X[, cols_scale])
  X[, cols_scale] <- train_scaled
  X1[, cols_scale] <- scale(X1[, cols_scale], center = attr(train_scaled, "scaled:center"), scale = attr(train_scaled, "scaled:scale"))
  
  Z <- X; Z1 <- X1
  
  
  coxfit_train <- coxph(Surv(Time, Status) ~ 1, data = train)
  coxfit_test <- coxph(Surv(Time1, Status1) ~ 1, data = test)
  
  
  out.data <- basehaz(coxfit_train, centered = FALSE)  # columns: hazard, time
  S0_grid <- exp(-out.data[,1])
  
  
  out.data1 <- basehaz( coxfit_test, centered = FALSE)  # columns: hazard, time
  S01_grid <- exp(-out.data1[,1])
  
  
  
  
  s0_init  <- S0_grid # length(Time)
  s01_init <- S01_grid  # length(Time1)
  
  
  # --- Pre-tuning of Decision Tree  ---
  nw <- factor(Status, levels = c(0,1), labels = c("cured","uncured"))
  ZDT <- as.data.frame(Z[, -1, drop = FALSE])
  K <- 10; set.seed(1)
  pos <- which(nw == "uncured"); neg <- which(nw == "cured")
  fpos <- split(sample(pos), rep(1:K, length.out = length(pos)))
  fneg <- split(sample(neg), rep(1:K, length.out = length(neg)))
  folds <- lapply(1:K, function(k) sort(c(fpos[[k]], fneg[[k]])))
  auc_fast <- function(y, p){ y <- as.integer(y=="uncured"); n1<-sum(y==1); n0<-sum(y==0); if(n1==0||n0==0) return(NA_real_); r<-rank(p,ties.method="average"); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
  
  
  
  ms_grid <- c(2, 3, 4)         
  md_grid <- c(28, 30, 32, 35)    
  cp_grid <- c(0, 0.00001, 0.00005) 
  
  
  best_auc <- -Inf
  best_params <- list(cp = 0.0001, minsplit = 2, maxdepth = 30)
  
  for (cp_val in cp_grid) {
    for (ms_val in ms_grid) {
      for (md_val in md_grid) {
        
        y_all <- c(); p_all <- c()
        
        for (k in 1:K) {
          vl <- folds[[k]]
          tr <- setdiff(seq_len(nrow(ZDT)), vl)
          
          if (length(unique(nw[tr])) < 2 || length(unique(nw[vl])) < 2) next
          
          fit <- try(
            rpart::rpart(
              nw ~ .,
              data = data.frame(ZDT, nw = nw)[tr, ],
              method = "class",
              control = rpart::rpart.control(
                cp = cp_val,
                minsplit = ms_val,
                maxdepth = md_val
              )
            ),
            silent = TRUE
          )
          if (inherits(fit, "try-error")) next
          
          pv <- try(
            predict(fit, newdata = ZDT[vl, , drop = FALSE], type = "prob")[, "uncured"],
            silent = TRUE
          )
          if (inherits(pv, "try-error")) next
          
          y_all <- c(y_all, nw[vl])
          p_all <- c(p_all, pv)
        }
        
        auc_val <- if (length(y_all) > 1) auc_fast(y_all, p_all) else NA_real_
        
        if (!is.na(auc_val) && auc_val > best_auc) {
          best_auc <- auc_val
          best_params <- list(
            cp = cp_val,
            minsplit = ms_val,
            maxdepth = md_val
          )
        }
      }
    }
  }
  
  initial_mod <- rpart::rpart(nw ~ ., data = data.frame(ZDT, nw = nw), method = "class", control = rpart::rpart.control(cp = best_params$cp))
  
  uncureprob <- predict(initial_mod, newdata = as.data.frame(Z[, -1, drop = FALSE]), type = "prob")[, "uncured"]
  uncurepred <- predict(initial_mod, newdata = as.data.frame(Z1[, -1, drop = FALSE]), type = "prob")[, "uncured"]
  
  
  
  
  
  emfit <- em.DT.Pois(Time, Status, Time1, Status1,
                      X, X1, Z, Z1,
                      coxfit_train$coefficients,
                      s0_init, s01_init,
                      uncureprob, uncurepred,
                      emmax, eps, best_params)
  emfit$best_params <- best_params
  if (Var) {
    cat("Starting bootstrap (", nboot, " iterations)...\n")
    latency_boot <- matrix(NA_real_, nboot, length(emfit$latencyfit))
    for (i in 1:nboot) {
      boot_idx <- sample(1:nrow(train), replace = TRUE)
      Time_b   <- Time[boot_idx]
      Status_b <- Status[boot_idx]
      X_b      <- X[boot_idx, , drop = FALSE]
      Z_b      <- Z[boot_idx, , drop = FALSE]
      
      boot_train_df <- data.frame(
        t  = Time_b,
        d  = Status_b,
        x2 = X_b[, "x2"],
        x3 = X_b[, "x3"],
        x4 = X_b[, "x4"],
        x6 = X_b[, "x6"]
      )
      
      coxfit_b <- try(
        coxph(Surv(t, d) ~ x2 + x3 + x4 + x6, data = boot_train_df),
        silent = TRUE)
      
      if (inherits(coxfit_b, "try-error")) next
      
      bh_b <- basehaz(coxfit_b, centered = FALSE)
      S0_b_grid <- exp(-bh_b$hazard)
      t_grid_b  <- bh_b$time
      
      idx_b <- pmax(1, findInterval(Time_b, t_grid_b))
      s0_b  <- S0_b_grid[idx_b]
      
      if (any(is.na(s0_b))) {
        s0_b[is.na(s0_b)] <- min(s0_b, na.rm = TRUE)
      }
      
      b_fit <- try(
        em.DT.Pois(Time_b, Status_b, Time1, Status1, 
                   X_b, X1, Z_b, Z1, 
                   emfit$latencyfit, s0_b, s01_init, 
                   uncureprob[boot_idx], uncurepred,
                   emmax, eps, best_params),
        silent = TRUE)
      if (!inherits(b_fit, "try-error")) latency_boot[i, ] <- b_fit$latencyfit
    }
    cat("Successful DT bootstrap fits:", sum(complete.cases(latency_boot)), "out of", nboot, "\n")
    emfit$latency_se <- apply(latency_boot, 2, sd, na.rm = TRUE)
    emfit$latency_p <- 2 * (1 - pnorm(abs(emfit$latencyfit / emfit$latency_se)))
  }
  # ==========================================
  # FINAL REPRESENTATIVE TREE
  # ==========================================
  Status_tr <- as.integer(Status)
  UN_final  <- as.numeric(emfit$UN)
  S_final   <- as.numeric(emfit$S)
  
  w_post_tr <- Status_tr + (1 - Status_tr) * (1 - (1 - UN_final)^S_final)
  w_post_tr <- pmin(pmax(w_post_tr, 1e-6), 1 - 1e-6)
  
  final_class <- factor(
    ifelse(w_post_tr >= 0.5, "uncured", "cured"),
    levels = c("cured", "uncured")
  )
  
  final_tree_data <- data.frame(
    x2 = Z[, "x2"],
    x3 = Z[, "x3"],
    x4 = Z[, "x4"],
    x6 = Z[, "x6"],
    final_class = final_class
  )
  
  final_tree <- rpart::rpart(
    final_class ~ x2 + x3 + x4 + x6,
    data = final_tree_data,
    method = "class",
    weights = w_post_tr,
    model = TRUE,
    control = rpart::rpart.control(
      cp = best_params$cp,
      minsplit = best_params$minsplit,
      maxdepth = best_params$maxdepth
    )
  )
  
  emfit$final_tree <- final_tree
  emfit$var_importance <- final_tree$variable.importance
  emfit$w_post_tr <- w_post_tr
  
  return(emfit)
}

build_consensus_tree <- function(DT_fit, train_df, B = 200, seed = 123) {
  set.seed(seed)
  
  # Rebuild scaled predictors exactly like in smcure.DT.Pois
  X <- model.matrix(~ x2 + x3 + x4 + x6, data = train_df)
  cols_scale <- colnames(X) %in% c("x2", "x6")
  train_scaled <- scale(X[, cols_scale])
  X[, cols_scale] <- train_scaled
  
  Z <- X
  
  # Posterior uncure probabilities from final EM fit
  w_post_tr <- as.numeric(DT_fit$w_post_tr)
  w_post_tr <- pmin(pmax(w_post_tr, 1e-6), 1 - 1e-6)
  
  # Final hard classes for summary tree building
  y_final <- factor(
    ifelse(w_post_tr >= 0.5, "uncured", "cured"),
    levels = c("cured", "uncured")
  )
  
  tree_data <- data.frame(
    x2 = Z[, "x2"],
    x3 = Z[, "x3"],
    x4 = Z[, "x4"],
    x6 = Z[, "x6"],
    y = y_final,
    w = w_post_tr
  )
  
  vars <- c("x2", "x3", "x4", "x6")
  trees <- vector("list", B)
  root_vars <- character(B)
  var_imp_mat <- matrix(0, nrow = B, ncol = length(vars))
  colnames(var_imp_mat) <- vars
  
  for (b in 1:B) {
    idx <- sample(seq_len(nrow(tree_data)), replace = TRUE)
    boot_data <- tree_data[idx, , drop = FALSE]
    
    fit_b <- rpart::rpart(
      y ~ x2 + x3 + x4 + x6,
      data = boot_data,
      weights = boot_data$w,
      method = "class",
      model = TRUE,
      control = rpart::rpart.control(
        cp = DT_fit$best_params$cp,
        minsplit = DT_fit$best_params$minsplit,
        maxdepth = DT_fit$best_params$maxdepth
      )
    )
    
    trees[[b]] <- fit_b
    root_vars[b] <- as.character(fit_b$frame$var[1])
    
    imp <- fit_b$variable.importance
    if (!is.null(imp)) {
      tmp <- setNames(rep(0, length(vars)), vars)
      tmp[names(imp)] <- imp
      var_imp_mat[b, ] <- tmp
    }
  }
  
  root_freq <- sort(prop.table(table(root_vars[root_vars != "<leaf>"])), decreasing = TRUE)
  mean_imp <- sort(colMeans(var_imp_mat, na.rm = TRUE), decreasing = TRUE)
  sd_imp <- apply(var_imp_mat, 2, sd, na.rm = TRUE)
  
  # Consensus summary tree
  consensus_tree <- rpart::rpart(
    y ~ x2 + x3 + x4 + x6,
    data = tree_data,
    weights = tree_data$w,
    method = "class",
    model = TRUE,
    control = rpart::rpart.control(
      cp = DT_fit$best_params$cp,
      minsplit = DT_fit$best_params$minsplit,
      maxdepth = DT_fit$best_params$maxdepth
    )
  )
  
  return(list(
    consensus_tree = consensus_tree,
    bootstrap_trees = trees,
    root_freq = root_freq,
    mean_importance = mean_imp,
    sd_importance = sd_imp,
    importance_matrix = var_imp_mat
  ))
}






# --- Data Execution ---
Melanoma <- read.table("~/Desktop/Desktop/melanoma-data-1.txt", header = T)
#Melanoma <- read.table("/cloud/project/melanoma-data-1.txt", header = T)

Melanoma$x3 <- as.numeric(Melanoma$x3)
Melanoma$x4 <- as.numeric(Melanoma$x4)

#Melanoma$x3 <- as.factor( Melanoma$x3)
#Melanoma$x4 <- as.factor(Melanoma$x4)

set.seed(2026)

idx <- createDataPartition(Melanoma$d, p = 0.7, list = FALSE)
train_df <- Melanoma[idx, ]
test_df  <- Melanoma[-idx, ]



methods <- c("Logit","Spline","DT")




plot_test_rocs_ptcm <- function(methods, train_df, test_df,
                                B =500, emmax = 1000, eps = 1e-3,
                                seed_fit = 2025, seed_imp = 790,
                                legend_pos = "bottomright") {
  
  fpr_grid <- seq(0, 1, length.out = 200)
  
  kept <- character(0)
  mean_tpr_list <- list()
  fits <- list()
  auc_mean <- c()
  valid_runs_vec <- c()
  mean_fpr_list <- list()
  
  for (m in methods) {
    cat("\n====================\nMethod:", m, "\n====================\n")
    
    func <- get(paste0("smcure.", m, ".Pois"))
    
    # Fit on TRAIN, predict on TEST
    set.seed(seed_fit)
    out<- try(func(train = train_df, test =test_df, Var = T, emmax = emmax, eps = eps, nboot=200),
              silent = TRUE)
    
    if (inherits(out, "try-error")) {
      cat("FAILED during fit for", m, ":\n", as.character(out), "\n")
      next
    }
    
    
    
    fits[[m]] <- out
    
    
    if (m == "DT") {
      
      plot_tree <- out$final_tree
      
      # Replace variable names in the tree structure
      plot_tree$frame$var <- gsub("x2", "Age", plot_tree$frame$var)
      plot_tree$frame$var <- gsub("x3", "Nod-Cat", plot_tree$frame$var)
      plot_tree$frame$var <- gsub("x4", "Gender", plot_tree$frame$var)
      plot_tree$frame$var <- gsub("x6", "Tum-Thick", plot_tree$frame$var)
      
      # Also update categorical split labels if present
      if (!is.null(attr(plot_tree, "xlevels"))) {
        attr(plot_tree, "xlevels") <- setNames(
          attr(plot_tree, "xlevels"),
          gsub(
            "x2", "Age",
            gsub(
              "x3", "Nod-Cat",
              gsub(
                "x4", "Gender",
                gsub("x6", "Tum-Thick", names(attr(plot_tree, "xlevels")))
              )
            )
          )
        )
      }
      
      # Plot final representative decision tree
      rpart.plot::rpart.plot(
        plot_tree,
        type = 2,
        extra = 104,
        fallen.leaves = TRUE,
        roundint = FALSE,
        main = "",
        cex = 0.7,
        left = FALSE
      )
      
      # Variable importance
      imp <- out$var_importance
      print(imp)
      
      names(imp) <- c(
        x2 = "Age",
        x3 = "Nod-Cat",
        x4 = "Gender",
        x6 = "Tum-Thick",
        TumThick = "Tum-Thick",
        Age = "Age",
        NodCat = "Nod-Cat",
        Gender = "Gender"
      )[names(imp)]
      
      barplot(
        sort(imp, decreasing = TRUE),
        ylab = "Importance",
        las = 2
      )
    }
 
    
    
    
    # term_names <- c("x2", "x3", "x4", "x6")
    
    if (!is.null(out$latencyfit) && !is.null(out$latency_se) && !is.null(out$latency_p)) {
      
      print(data.frame(
        term    = names(out$latencyfit),
        beta    = sprintf("%.4f", as.numeric(out$latencyfit)),
        se      = sprintf("%.4f", as.numeric(out$latency_se)),
        p_value = sprintf("%.4f", as.numeric(out$latency_p))
      ))
    }
    
    # ----- Posterior Pr(uncured | observed) on TEST (PTCM) -----
    Status_te <- as.integer(test_df$d)
    PRED   <- as.numeric(out$PRED)
    S_pred <- as.numeric(out$S.pred)
    
    w_post_te <- Status_te + (1 - Status_te) *  (1-(1-PRED)^S_pred)
    w_post <- pmin(pmax(w_post_te, 1e-6), 1 - 1e-6)
    
    # We want ROC for π(x) = Pr(uncured), so:
    pi_hat  <-PRED
    pi_hat  <- pmin(pmax(pi_hat, 1e-6), 1 - 1e-6)
    pi_post <- w_post
    
    
    set.seed(seed_imp)
    thr_grid <- seq(0, 1, length.out = length(Status_te))  
    
    
    V_test_ext <- matrix(
      rbinom(length(w_post) * B, size = 1, prob = rep(w_post, each = B)),
      nrow = length(w_post), byrow = TRUE
    )
    
    
    
    tpr_mat <- matrix(NA_real_, nrow = B, ncol = length(thr_grid))
    fpr_mat <- matrix(NA_real_, nrow = B, ncol = length(thr_grid))
    aucs <- rep(NA_real_, B)
    
    
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    
    for (b in 1:B) {
      C <- V_test_ext[, b]
      if (length(unique(C)) < 2) {
        utils::setTxtProgressBar(pb, b)
        next
      }
      
      r <- pROC::roc(response = C, predictor = pi_hat, quiet = TRUE, direction = "<")
      aucs[b] <- as.numeric(pROC::auc(r))
      
      # stable curve points
      cc <- pROC::coords(
        r,
        x = thr_grid,
        input = "threshold",
        ret = c("threshold", "sensitivity", "specificity"),
        transpose = FALSE
      )
      
      tpr_mat[b, ] <- cc[, "sensitivity"]
      fpr_mat[b, ] <- 1 - cc[, "specificity"]
      
      utils::setTxtProgressBar(pb, b)
    }
    close(pb)
    
    
    # cat(sprintf("Done %s: TRAIN mean AUC(π) = %.4f\n", m, mean(aucs_tr, na.rm = TRUE)))
    mean_fpr <- colMeans(fpr_mat, na.rm = TRUE)
    mean_tpr <- colMeans(tpr_mat, na.rm = TRUE)
    
    # --- GUARD AGAINST INVALID ROC (ALL NA / NON-FINITE) ---
    if (!any(is.finite(mean_fpr)) || !any(is.finite(mean_tpr))) {
      cat("Skipping", m, ": no valid ROC points\n")
      next
    }
    
    mean_tpr_list[[m]] <- mean_tpr
    mean_fpr_list[[m]] <- mean_fpr
    
    valid_runs <- sum(is.finite(aucs))
    cat("\nValid ROC runs:", valid_runs, "out of", B, "\n")
    
    kept <- c(kept, m)
    mean_tpr_list[[m]] <- colMeans(tpr_mat, na.rm = TRUE)
    auc_mean[m] <- mean(aucs, na.rm = TRUE)
    valid_runs_vec[m] <- valid_runs
    
    
    cat(sprintf("Done %s: TEST mean AUC(π) = %.4f\n", m, auc_mean[m]))
  }
  
  if (length(kept) == 0) stop("No methods succeeded.")
  
  cols <- grDevices::rainbow(length(kept))
  ltys <- rep(1:6, length.out = length(kept))
  
  plot(mean_fpr_list[[kept[1]]], mean_tpr_list[[kept[1]]],
       type = "l", col = cols[1], lty = ltys[1], lwd = 2,
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = "ROC")
  
  if (length(kept) > 1) {
    for (k in 2:length(kept)) {
      lines(mean_fpr_list[[kept[k]]], mean_tpr_list[[kept[k]]], col = cols[k], lty = ltys[k], lwd = 2)
    }
  }
  
  legend_labels <- sprintf("%s (AUC=%.4f)", kept, auc_mean[kept])
  legend(legend_pos, legend = legend_labels,
         col = cols, lty = ltys, lwd = 2, bty = "n", cex = 0.85)
  
  invisible(list(methods = kept, auc = auc_mean, valid = valid_runs_vec,
                 fpr = fpr_grid, tpr = mean_tpr_list))
  
  
   
}


res_test <- plot_test_rocs_ptcm(methods, train_df, test_df, B = 500, emmax = 1000, eps = 1e-3)


time<-proc.time()-t
time




