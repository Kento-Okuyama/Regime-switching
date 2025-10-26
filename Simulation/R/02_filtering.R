# R/02_filtering.R (Simulation Version - 2-Factor/4-Free-Loadings Compatible)

library(progress)
library(torch)

filtering <- function(seed, N, Nt, O1, O2, U1, y1, y2, DO, init, maxIter, n_train, patience, min_delta = 1e-3) {
  set.seed(seed + init)
  lEpsilon <- 1e-3
  ceil <- 1e6
  sEpsilon <- 1e-8
  stopCrit <- 1e-4
  lr <- 1e-1
  min_lr <- 1e-3 
  wd <- 1e-2
  epsilon <- 1e-8
  betas <- c(.9, .999)
  const <- (2 * pi)^(-O1 / 2)
  
  # --- Convert data to torch tensors ---
  y1 <- torch_tensor(y1); y2 <- torch_tensor(y2)
  
  # --- Parameter Initialization ---
  B11 <- torch_tensor(rnorm(U1,0,1), requires_grad=FALSE)
  log_B12_delta <- torch_tensor(log(abs(rnorm(U1,0,1))), requires_grad=FALSE)
  B21d <- torch_tensor(rnorm(U1,0,1), requires_grad=FALSE)
  B22d <- torch_tensor(rnorm(U1,0,1), requires_grad=FALSE)
  B31 <- torch_tensor(rnorm(U1,0,1), requires_grad=FALSE)
  B32 <- torch_tensor(rnorm(U1,0,1), requires_grad=FALSE)
  gamma1 <- torch_tensor(rnorm(1,0,1), requires_grad=FALSE)
  gamma2 <- torch_tensor(rnorm(U1,0,1), requires_grad=FALSE)
  gamma3 <- torch_tensor(rnorm(1,0,1), requires_grad=FALSE)
  gamma4 <- torch_tensor(rnorm(U1,0,1), requires_grad=FALSE)
  log_Qd <- torch_tensor(log(rgamma(U1,9,4)), requires_grad=FALSE)
  log_R1d <- torch_tensor(log(rgamma(O1,9,4)), requires_grad=FALSE)
  logit_tP_SB <- torch_tensor(qlogis(runif(1,0,0.1)), requires_grad=FALSE)
  
  # Change the number of free factor loadings to 4
  Lmd1f1 <- torch_tensor(abs(rnorm(1,0,1)), requires_grad=FALSE)
  Lmd1f2 <- torch_tensor(abs(rnorm(1,0,1)), requires_grad=FALSE)
  Lmd1f3 <- torch_tensor(abs(rnorm(1,0,1)), requires_grad=FALSE)
  Lmd1f4 <- torch_tensor(abs(rnorm(1,0,1)), requires_grad=FALSE)
  
  # Initialize y2 measurement model parameters
  Lmd2 <- c(1.0, abs(rnorm(2, 0, 1))) 
  R2d <- rgamma(3,shape=9,rate=4) 
  Lmd2f <- torch_tensor(Lmd2[2:3], requires_grad = FALSE) 
  R2d[R2d < 1e-4] <- 1e-4
  log_R2d <- torch_tensor(log(R2d), requires_grad = FALSE) 
  
  # --- Initialize variables for the optimization loop ---
  sumLik_best <- -Inf
  iter <- 1
  patience_counter <- 0
  iter_best <- NA
  theta_best <- NULL
  mEta_train_end <- NULL
  mPr_train_end <- NULL
  mP_train_end <- NULL
  mPr_pred_best <- NULL
  mEta_best <- NULL
  mP_best <- NULL
  
  pb <- progress_bar$new(format = "  optim step :current/:total [:bar] - ETA: :eta | sumLik: :sumlik", total = maxIter, width = 80, clear = FALSE)
  pb$tick(0, tokens = list(sumlik = "N/A"))
  
  # Update the list of parameters subject to optimization
  theta <- list(B11=B11, log_B12_delta=log_B12_delta, B21d=B21d, B22d=B22d, B31=B31, B32=B32,
                Lmd1f1=Lmd1f1, Lmd1f2=Lmd1f2, Lmd1f3=Lmd1f3, Lmd1f4=Lmd1f4, 
                Lmd2f=Lmd2f, log_Qd=log_Qd, log_R1d=log_R1d, log_R2d=log_R2d,
                gamma1=gamma1, gamma2=gamma2, gamma3=gamma3, gamma4=gamma4, logit_tP_SB=logit_tP_SB)
  
  lapply(theta, function(p) p$requires_grad_(TRUE))
  optimizer <- torch::optim_adamw(params=theta, lr=lr, betas=betas, eps=epsilon, weight_decay=wd)
  
  while (iter <= MAX_ITER) {
    
    # --- Calculate eta2 scores from current parameters (Bartlett's factor scores) ---
    R2d <- log_R2d$exp()
    R2_inv <- torch_inverse(R2d$diag()) 
    Lmd2 <- torch_cat(c(torch_tensor(1.0), Lmd2f))$unsqueeze(1)
    Lmd2_t_R2_inv <- Lmd2$matmul(R2_inv)
    factor_weights <- linalg_solve(
      Lmd2_t_R2_inv$matmul(Lmd2$transpose(1,2)),
      Lmd2_t_R2_inv 
    )
    eta2 <- y2$matmul(factor_weights$transpose(1,2))$squeeze()
    
    jEta <- torch_full(c(N, Nt, 2, 2, U1), 0)
    jP <- torch_full(c(N, Nt, 2, 2, U1, U1), 0)
    jV <- torch_full(c(N, Nt, 2, 2, O1), NaN)
    jF <- torch_full(c(N, Nt, 2, 2, O1, O1), NaN)
    jEta2 <- torch_full(c(N, Nt, 2, 2, U1), 0)
    jP2 <- torch_full(c(N, Nt, 2, 2, U1, U1), 0)
    mEta <- torch_full(c(N, Nt+1, 2, U1), 0)
    mP <- torch_full(c(N, Nt+1, 2, U1, U1), NaN)
    W <- torch_full(c(N, Nt, 2, 2), NaN)
    jPr <- torch_full(c(N, Nt, 2, 2), 0)
    mLik <- torch_full(c(N, Nt), NaN)
    jPr2 <- torch_full(c(N, Nt, 2, 2), 0)
    mPr_filtered <- torch_full(c(N, Nt+1, 2), NaN)
    mPr_pred <- torch_full(c(N, Nt+1, 2), NaN)
    jLik <- torch_full(c(N, Nt, 2, 2), 0)
    tPr <- torch_full(c(N, Nt, 2, 2), NaN)
    KG <- torch_full(c(N, Nt, 2, 2, U1, O1), 0)
    I_KGLmd <- torch_full(c(N, Nt, 2, 2, U1, U1), NaN)
    subEta <- torch_full(c(N, Nt, 2, 2, U1), NaN)
    eta1_pred <- torch_full(c(N, Nt, U1), NaN)
    
    mEta <- torch_full(c(N, Nt+1, 2, U1), 0)
    mP <- torch_full(c(N, Nt+1, 2, U1, U1), NaN)
    mP[,1,,,] <- torch_eye(U1)
    mPr_filtered <- torch_full(c(N, Nt+1, 2), NaN)
    mPr_filtered[,1,1] <- 1 - sEpsilon
    mPr_filtered[,1,2] <- sEpsilon
    
    # Construct the new measurement model (Lmd1)
    Lmd1 <- torch_full(c(O1, U1), 0)
    # New Factor 1: Afraid to Fail (Observed variables 1-3)
    Lmd1[1,1] <- 1
    Lmd1[2,1] <- Lmd1f1
    Lmd1[3,1] <- Lmd1f2
    # New Factor 2: No Positive Affect (Observed variables 4-6)
    Lmd1[4,2] <- 1
    Lmd1[5,2] <- Lmd1f3
    Lmd1[6,2] <- Lmd1f4
    
    B12 <- B11 + log_B12_delta$exp()
    B1 <- torch_cat(c(B11, B12))$reshape(c(2, U1))
    B21 <- B21d$diag()
    B22 <- B22d$diag()
    Lmd1T <- Lmd1$transpose(1, 2) 
    Q <- log_Qd$clamp(min = -15)$exp()$diag()
    R1 <- log_R1d$clamp(min = -15)$exp()$diag() 
    
    B2 <- torch_cat(c(B21, B22))$reshape(c(2, U1, U1))
    B3 <- torch_cat(c(B31, B32))$reshape(c(2, U1))
    tP_SB <- torch_sigmoid(logit_tP_SB) * 0.1    
    tPr[,,1,2] <- tP_SB
    tPr[,,2,2] <- 1 - tP_SB   
    
    for (t in 1:Nt) {
      jEta[,t,,,] <- B1$unsqueeze(1)$unsqueeze(-2) + mEta[,t,,]$clone()$unsqueeze(2)$matmul(B2$unsqueeze(1)) + eta2$unsqueeze(-1)$unsqueeze(-1)$unsqueeze(-1) * B3$unsqueeze(1)$unsqueeze(-2)
      jP[,t,,,,] <- B2$unsqueeze(1)$unsqueeze(3)$matmul(mP[,t,,,]$clone()$unsqueeze(2))$matmul(B2$unsqueeze(1)$unsqueeze(3)) + Q$unsqueeze(1)$expand(c(1, 2, 2, -1, -1))
      jP[,t,,,,] <- jP[,t,,,,]$clone()$clip(max=ceil) + 1e-6 * torch_eye(U1)
      jV[,t,,,] <- y1[,t,]$clone()$unsqueeze(-2)$unsqueeze(-2) - jEta[,t,,,]$clone()$matmul(Lmd1T$unsqueeze(1)) # LmdT -> Lmd1T
      jF[,t,,,,] <- Lmd1$unsqueeze(1)$matmul(jP[,t,,,,]$clone()$matmul(Lmd1T$unsqueeze(1))) + R1$unsqueeze(1)$expand(c(1, 2, 2, -1, -1)) # Lmd->Lmd1, LmdT->Lmd1T, R->R1
      jF[,t,,,,] <- jF[,t,,,,]$clone()$clip(max=ceil) + 1e-6 * torch_eye(O1)
      C <- jP[,t,,,,]$clone()$matmul(Lmd1T$unsqueeze(1)); F_val <- jF[,t,,,,]$clone() # LmdT -> Lmd1T
      C_t <- C$transpose(4,5); F_t <- F_val$transpose(4,5); XT <- linalg_solve(F_t, C_t)
      KG[,t,,,,] <- XT$transpose(4,5)
      jEta2[,t,,,] <- jEta[,t,,,] + KG[,t,,,,]$clone()$matmul(jV[,t,,,]$clone()$unsqueeze(-1))$squeeze()
      I_KGLmd[,t,,,,] <- torch_eye(U1)$unsqueeze(1)$expand(c(1,2,2,-1,-1)) - KG[,t,,,,]$clone()$matmul(Lmd1$unsqueeze(1)) # Lmd -> Lmd1
      jP2[,t,,,,] <- I_KGLmd[,t,,,,]$clone()$matmul(jP[,t,,,,]$clone()$matmul(I_KGLmd[,t,,,,]$clone()$transpose(4, 5))) + KG[,t,,,,]$clone()$matmul(R1$unsqueeze(1)$expand(c(1, 2, 2, -1, -1)))$matmul(KG[,t,,,,]$clone()$transpose(4, 5)) # R -> R1
      log_det_jF <- torch_slogdet(jF[,t,,,,]$clone())[[2]]
      v_unsqueezed <- jV[,t,,,]$clone()$unsqueeze(-1)
      solved_x <- linalg_solve(jF[,t,,,,]$clone(), v_unsqueezed)
      quadratic_term <- -.5 * v_unsqueezed$transpose(4,5)$matmul(solved_x)$squeeze()$squeeze()
      jLik[,t,,] <- log(sEpsilon + const) - .5 * log_det_jF + quadratic_term
      jLik[,t,,]$clip_(min=-ceil, max=ceil)
      eta1_pred[,t,] <- mPr_filtered[,t,1]$clone()$unsqueeze(-1) * mEta[,t,1,]$clone() + mPr_filtered[,t,2]$clone()$unsqueeze(-1) * mEta[,t,2,]$clone()
      interaction_term <- (eta1_pred[,t,]$clone() * gamma4)$sum(dim=2) * eta2$squeeze()
      tPr[,t,1,1] <- (gamma1 + eta1_pred[,t,]$clone()$matmul(gamma2) + eta2$squeeze() * gamma3 + interaction_term)$sigmoid()$clip(min=sEpsilon, max=1-sEpsilon)
      tPr[,t,2,1] <- 1 - tPr[,t,1,1]
      jPr[,t,,] <- tPr[,t,,]$clone() * mPr_filtered[,t,]$clone()$unsqueeze(-2)
      log_joint_lik <- jLik[,t,,]$clone() + jPr[,t,,]$clone()$log()
      log_mLik <- torch_logsumexp(log_joint_lik, dim = c(2, 3))
      mLik[,t] <- log_mLik
      jPr2[,t,,] <- (log_joint_lik - log_mLik$unsqueeze(-1)$unsqueeze(-1))$exp()
      mPr_pred[,t+1,] <- jPr[,t,,]$sum(3)$clip(min=sEpsilon, max=1-sEpsilon)
      mPr_filtered[,t+1,] <- jPr2[,t,,]$sum(3)$clip(min=sEpsilon, max=1-sEpsilon)
      W[,t,,] <- jPr2[,t,,]$clone() / (mPr_filtered[,t+1,]$clone()$unsqueeze(-1) + sEpsilon)
      W[,t,,]$clip_(min = sEpsilon, max = 1 - sEpsilon)
      mEta[,t+1,,] <- (W[,t,,]$clone()$unsqueeze(-1) * jEta2[,t,,,]$clone())$sum(3) 
      subEta[,t,,,] <- mEta[,t+1,,]$unsqueeze(3) - jEta2[,t,,,]
      mP[,t+1,,,] <- (W[,t,,]$clone()$unsqueeze(-1)$unsqueeze(-1) * (jP2[,t,,,,] + subEta[,t,,,]$clone()$unsqueeze(-1)$matmul(subEta[,t,,,]$clone()$unsqueeze(-2))))$sum(3)
      mP[,t+1,,,] <- mP[,t+1,,,]$clone() + 1e-6 * torch_eye(U1)
    } 
    
    loss <- -mLik[, 1:n_train]$nanmean()
    
    if (!is.finite(as.numeric(loss))) {
      cat('   error in calculating the sum likelihood. Halving learning rate and continuing.\n')
      if (!is.null(theta_best)) { theta <- lapply(theta_best, function(p) p$clone()$detach()) }
      lr <- lr * 0.5; if (lr < min_lr) break;
      iter <- iter + 1; next
    }
    
    sumLik <- -as.numeric(loss)
    pb$tick(tokens = list(sumlik = sprintf("%.2f", sumLik)))
    
    if (sumLik > sumLik_best + min_delta) {
      sumLik_best <- sumLik; iter_best <- iter
      theta_best <- lapply(theta, function(p) p$clone()$detach())
      mEta_train_end <- mEta[, n_train + 1, , ]$clone()$detach()
      mPr_train_end <- mPr_filtered[, n_train + 1, ]$clone()$detach()
      mP_train_end <- mP[, n_train + 1, , , ]$clone()$detach()
      patience_counter <- 0
      mPr_pred_best <- mPr_pred$clone()$detach()
      mEta_best <- mEta$clone()$detach()
      mP_best <- mP$clone()$detach()
    } else {
      patience_counter <- patience_counter + 1
    }
    
    if (patience_counter >= patience) { 
      current_lr <- optimizer$param_groups[[1]]$lr
      if (current_lr > min_lr) {
        new_lr <- current_lr * 0.5
        optimizer$param_groups[[1]]$lr <- new_lr
        cat(paste0("\n   -> Patience reached. Reducing learning rate to: ", sprintf("%.6f", new_lr), "\n"))
        if (!is.null(theta_best)) {
          cat("   -> Restoring model to the best parameters...\n")
          for (name in names(theta)) {
            theta[[name]]$detach_()$copy_(theta_best[[name]])
            theta[[name]]$requires_grad_(TRUE)
          }
        }
        patience_counter <- 0; iter <- iter + 1; next 
      } else {
        pb$message('   Early stopping triggered.'); break
      }
    } else {
      optimizer$zero_grad()
      loss$backward()
      torch::nn_utils_clip_grad_norm_(theta, max_norm = 10.0)
      optimizer$step()
    }
    
    iter <- iter + 1; gc() 
  }
  
  if (!is.null(theta_best)) {
    theta_best <- lapply(theta_best, as.array)
    mEta_train_end <- as.array(mEta_train_end)
    mPr_train_end <- as.array(mPr_train_end)
    mP_train_end <- as.array(mP_train_end)
    mPr_pred_best <- as.array(mPr_pred_best)
    mEta_best <- as.array(mEta_best)
    mP_best <- as.array(mP_best)
  }
  
  filter <- list(iter_best = iter_best, sumLik_best = sumLik_best, theta_best = theta_best,
                 mEta_train_end = mEta_train_end, mPr_train_end = mPr_train_end, mP_train_end = mP_train_end,
                 mPr_pred_best = mPr_pred_best, mEta_best = mEta_best, mP_best = mP_best)
  gc(); return(filter)
}