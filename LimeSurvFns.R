##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Project: Functions for LIME in survival framework
## Created: Mar 26, 2022
## Author: Krithika Suresh
## Description: Modifies "lime" code found in https://github.com/thomasp85/lime/find/master
## Modifies "RMTL" code found in https://github.com/transbioZI/RMTL/blob/master/R/MTL.R to incorporate weights
## Implements smoothed estimation across the time points 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(glmnet)
library(RMTL)

# Function for prediction from a Cox model --------------------------------
predictCox <- function(model, dat, w){
  sf <- survfit(model, newdata = dat)
  pred_surv <- summary(sf, times=w)$surv
  return(pred_surv)
}

# Function for prediction from an RSF model -------------------------------
predictRSF <- function (object, newdata, times, ...) {
  N <- NROW(newdata)
  class(object) <- c("rfsrc", "grow")
  S <- randomForestSRC::predict.rfsrc(object, newdata = newdata)$survival
  if(N == 1) S <- matrix(S, nrow = 1)
  Time <- object$time.interest
  p <- cbind(1, S)[, 1 + prodlim::sindex(Time, times),drop = FALSE]
  if(NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed") #prediction fails if predictions cannot be made for every subject
  return(c(p))
}

# Function for prediction from a TDCox model ------------------------------
# model <- explainer$model
# dat <- explainer$preprocess(case_perm)
# pred_time=pred_times[1]
predictTDCox <- function(model, dat, w) {
  sfi <- model$sfi
  bet <- model$coef
  dat <- dat %>% mutate_all(~as.numeric(as.character(.x)))
  nt <- nrow(dat)
  Fw <- NULL
  
  for(i in 1:nt)
  {
    sfi$hazard<-sfi$haz0* exp(c(bet[1:(length(bet)-1)]%*%t(dat[i,])) + 
                                bet[length(bet)]*dat$x1[i]*log(sfi$time))
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time, sfi$Haz, w, subst=0)
    Fw <- rbind(Fw,exp(-tmp))
  }
  return(data.frame("pred_yr"=Fw))
}

# Prediction functions ----------------------------------------------------
predict_model.rfsrc <- function(x, newdata, pred_time, cause, ...) {
  pred <- 1-predictRisk(x, newdata = newdata, times = pred_time, cause=cause)
  pred <- data.frame("pred_yr"=as.numeric(pred))
  return(pred)
}
model_type.rfsrc <- function(x, ...) {
  return("classification")
}

predict_model.coxph <- function(x, newdata, pred_time, ...) {
  pred <- 1-predictRisk(x, newdata = newdata, times = pred_time)
  pred <- data.frame("pred_yr"=as.numeric(pred))
  return(pred)
}
model_type.coxph <- function(x, ...) {
  return("classification")
}

predict_model.tdcox <- function(x, newdata, pred_time, ...) {
  pred <- predictTDCox(x, newdata, pred_time)
  pred <- data.frame("pred_yr"=unlist(pred), row.names=NULL)
  return(pred)
}
model_type.tdcox <- function(x, ...) {
  return("classification")
}

predict_model.LearnerSurv <- function(x, newdata, pred_time, ...){
  newdata$eventtime <- rnorm(nrow(newdata))
  newdata$status <- sample(c(0,1), size=nrow(newdata), replace=TRUE)
  task.ND <- TaskSurv$new(id = "right_censored", backend = newdata, time="eventtime", event="status", type="right")
  prediction_km <- model_km$predict(task.ND)
  pred = x$predict(task.ND)
  pod = po("distrcompose", param_vals = list(form = "ph", overwrite = FALSE))
  p = pod$predict(list(base = prediction_km, pred = pred))$output
  ret <- t(p$distr$survival(pred_time))
  rownames(ret) <- NULL
  colnames(ret) <- "pred_yr"
  return(data.frame(ret))
}
model_type.LearnerSurv <- function(x, ...) {
  return("classification")
}

predict_model.FGR <- function(x, newdata, pred_time, cause, ...) {
  cif_yr <- predictRisk(x, newdata=newdata, times=pred_time, cause=cause)
  return(data.frame("pred_yr"=cif_yr))
}
model_type.FGR <- function(x, ...) {
  return("classification")
}

# Explanation data frame --------------------------------------------------
explain_dat <- function(x, explainer, labels = NULL, n_labels = NULL,
                        n_features, n_permutations = 5000,
                        feature_select = 'auto', dist_fun = 'gower',
                        kernel_width = NULL, gower_pow = 1, pred_times, scale=FALSE, seed, cause, form=NULL) {
  m_type <- "classification"
  o_type <- "prob"
  
  if (is.null(kernel_width)) {
    kernel_width <- sqrt(ncol(x)) * 0.75
  }
  kernel <- lime:::exp_kernel(kernel_width)
  
  set.seed(seed)
  #Code to create permuted data sets (want to keep this the same across the time points)
  case_perm <- lime:::permute_cases.data.frame(x, n_permutations, explainer$feature_distribution,
                                               explainer$bin_continuous, explainer$bin_cuts,
                                               explainer$use_density)
  
  case_res_full <- NULL
  case_res_full_t <- list()
  case_res <- NULL 
  
  pred_add <- predict_model(explainer$model, explainer$preprocess(case_perm), pred_time=pred_times, type=o_type, cause=cause)
  case_res_full <- pred_add
  mat_pred_add <- matrix(pred_add$pred_yr, nrow=n_permutations)
  # case_res_full_t <- as.list(pred_add)
  for(s in 1:ncol(mat_pred_add)) {
    case_res_full_t[[s]] <- data.frame("pred_yr"=mat_pred_add[,s])
  }
  case_res <- mat_pred_add[1,]
  
  case_res_full <- lime:::set_labels(case_res_full, explainer$model)
  
  case_ind <- split(seq_len(nrow(case_perm)), rep(seq_len(nrow(x)), each = n_permutations))
  
  if (dist_fun == 'gower') {
    sim <- 1 - (gower::gower_dist(case_perm[1, , drop = FALSE], case_perm[, ])) ^ gower_pow
  }
  
  perms <- lime:::numerify(case_perm, explainer$feature_type, 
                           explainer$bin_continuous, explainer$bin_cuts)
  
  if (dist_fun != 'gower') {
    sim <- kernel(c(0, dist(lime:::feature_scale(perms, explainer$feature_distribution, explainer$feature_type, explainer$bin_continuous),
                            method = dist_fun)[seq_len(n_permutations-1)]))
  }
  
  if(scale==TRUE) {
    perms <- scale(perms)
  }
  perms <- as.matrix(perms)
  
  if(!is.null(form)) {
    perms <- model.matrix(as.formula(form), data=data.frame(perms))[,-1]
    case_perm <- as.data.frame(model.matrix(as.formula(form), data=data.frame(case_perm))[,-1])
    names(case_perm) <- colnames(perms)
  }
  
  #Outputs for RMTL 
  x_list <- rep(list(as.matrix(cbind("(Intercept)"=1,perms))), length(pred_times))
  y_list <- case_res_full_t
  weights_list = rep(list(sim), length(pred_times))
  # weights_list = rep(list(rep(1,length(sim))), length(pred_times))
  
  #Outputs for Quad
  x_orig = do.call("rbind", rep(list(as.matrix(perms)), length(pred_times)))
  t_full <- rep(pred_times, each = n_permutations)
  y_full = case_res_full[, , drop = FALSE]
  weights_full = rep(sim, length(pred_times))
  feature_method = feature_select
  #append time-dependent measures of all of the variables
  mat_torig <- cbind(t_full, t_full^2)
  colnames(mat_torig) <- c("t", "t2")
  mat_t <- apply(x_orig, 2, function(x) x*t_full)
  colnames(mat_t) <- paste0("t_", colnames(mat_t))
  mat_t2 <- apply(x_orig, 2, function(x) x*t_full^2)
  colnames(mat_t2) <- paste0("t2_", colnames(mat_t2))
  x_full <- cbind(x_orig, mat_torig, mat_t, mat_t2)
  
  return(list("explainer"=explainer,
              "dat"=x,
              "perms"=perms, 
              "case_perm"=case_perm, 
              "case_res"=unlist(case_res),
              "pred_times"=pred_times,
              "x_list"=x_list, 
              "y_list"=y_list, 
              "weights_list"=weights_list,
              "x_full"=x_full,
              "y_full"=y_full,
              "weights_full"=weights_full,
              "n_permutations"=n_permutations))
}

# Explanation models -----------------------------------------------------
model_permutations_glmnet <- function(x_full, y_full, weights_full, feature_method="auto", method, n_permutations=NULL, pred_times=NULL) {
  labels <- names(y_full)
  x_full <- x_full[, colSums(is.na(x_full)) == 0 & apply(x_full, 2, var) != 0, drop = FALSE]
  label <- labels[1]
  
  features <- lime:::select_features(feature_method, x_full, y_full[[label]], weights_full, 1000)
  
  shuffle_order <- sample(length(y_full[[label]])) # glm is sensitive to the order of the examples
  
  fit <- glmnet(x_full[shuffle_order, features],
                y_full[[label]][shuffle_order],
                weights = weights_full[shuffle_order], alpha = 0, lambda = 2 / length(y_full[[label]]))
  
  # modification of Kvilseth’s coefficient of determination (Willett, Singer; American Statistician)
  # "Another Cautionary Note about R2: Its Use in Weighted Least-Squares Regression Analysis"
  r2 <- NULL
  resid <- y_full$pred_yr-predict(fit,newx=x_full)
  if(method=="LM") {
    r2 <- 1-sum(resid^2)/(sum(y_full^2)-nrow(y_full)*mean(t(y_full))^2)
  } else if(method=="Quad") {
    y_temp <- split(y_full, factor(sort(rank(row.names(y_full))%%5)))
    resid_temp <- data.frame(resid)
    resid_temp <- split(resid_temp, factor(sort(rank(row.names(resid_temp))%%5)))
    for(r in 1:length(pred_times)) {
      r2_temp <- 1-sum(resid_temp[[r]]^2)/(sum(y_temp[[r]]^2)-nrow(y_temp[[r]])*mean(t(y_temp[[r]]))^2)
      r2 <- c(r2, r2_temp)
    }
  }
  
  # #CSI: Coefficients Stability Index 
  lam <- fit$lambda
  
  coefs <- coef(fit)
  
  res <- list(
    "feature" = rownames(coefs),
    "feature_weight" = matrix(coefs),
    "model_r2" = r2,
    "model_lam" = lam, 
    stringsAsFactors = FALSE
  )
}

model_permutations_RMTLweighted <- function(x_list, y_list, weights_list, Regularization="Lasso", G) {
  if (all(weights_list[[1]][-1] == 0)) {
    stop('All permutations have no similarity to the original observation. Try setting bin_continuous to TRUE and/or increase kernel_size', call. = FALSE)
  }
  
  cv_vals <- NULL 
  Lam2 <- 0 
 
  #Standardized x covariates (necessary to improve performance)
  # https://stackoverflow.com/questions/59846325/confusion-about-standardize-option-of-glmnet-package-in-r
  x_list_temp <- x_list[[1]]
  x_list_sd <- c(as.numeric((data.frame(x_list_temp) %>% mutate_all(~sd(.x) %>% as.vector))[1,-1]))
  x_list_mean <- c(as.numeric((data.frame(x_list_temp) %>% mutate_all(~mean(.x) %>% as.vector))[1,-1]))
  x_list_temp <- data.frame(x_list_temp) %>% mutate_all(~(scale(.) %>% as.vector))
  x_list_temp[,1] <- 1 
  x_list_scale <- rep(list(x_list_temp), length(pred_times))
  
  #Weighted, scaled
  cvfit_scaleW <- cvMTLweighted(x_list_scale, y_list, weights_list,
                                type="Regression", Regularization=Regularization,
                                Lam2=Lam2, G=G)
  
  fit_scaleW <- MTLweighted(X=x_list_scale, Y=y_list, weights_list,
                            type="Regression", Regularization=Regularization,
                            Lam1=cvfit_scaleW$Lam1.min, Lam2=Lam2, G=G)
  
  coefs <- data.frame(fit_scaleW$W) %>% mutate_all(~.x/c(1,x_list_sd))
  coefs[1,] <- fit_scaleW$W[1,] + fit_scaleW$C + colSums(data.frame(fit_scaleW$W)[-1,] %>% mutate_all(~.x*-x_list_mean/x_list_sd))
  
  temp <- lapply(c(1:ncol(coefs)), function(x)
    as.matrix(x_list_scale[[x]]) %*% fit_scaleW$W[,x] + fit_scaleW$C[x])
  r2 <- NULL
  for(r in 1:length(y_list)) {
    y_temp <- y_list[[r]]
    pred_temp <- temp[[r]]
    resid_temp <- y_temp-pred_temp
    r2_temp <- 1-sum(resid_temp^2)/(sum(y_temp^2)-nrow(y_temp)*mean(t(y_temp))^2)
    r2 <- c(r2, r2_temp)
  }
  
  # #CSI: Coefficients Stability Index 
  lam <- fit_scaleW$Lam1
 
  res <- list(
    feature = colnames(x_list[[1]]),
    feature_weight = coefs,
    model_r2 = r2,
    model_lam = lam,
    stringsAsFactors = FALSE
  )
}

# Explanation functions ---------------------------------------------------
explain_LM <- function(explain_dat_model) {
  dat <- explain_dat_model$dat
  perms <- explain_dat_model$perms
  pred_times <- explain_dat_model$pred_times
  case_perm <- explain_dat_model$case_perm
  #prediction probabilities at each time point
  case_res <- explain_dat_model$case_res
  explainer <- explain_dat_model$explainer 
  
  res <- data.frame("feature" = colnames(perms))
  res_sub <- data.frame("feature" = c("(Intercept)", colnames(perms)))
  time_sub <- NULL
  r2_sub <- NULL
  lam_sub <- NULL
  for(i in 1:length(pred_times)) {
    val <- c(1, perms[1,]) #first row corresponds to true covariates
    time_temp <- model_permutations_glmnet(explain_dat_model$x_list[[i]], 
                                           explain_dat_model$y_list[[i]], 
                                           explain_dat_model$weights_list[[i]],
                                           method="LM")
    time_sub <- cbind(time_sub, val*time_temp$feature_weight)
    r2_sub <- cbind(r2_sub, time_temp$model_r2)
    lam_sub <- cbind(lam_sub, time_temp$model_lam)
  }
  colnames(time_sub) <- paste0("value_t_",pred_times)
  colnames(r2_sub) <- paste0("r2_t_",pred_times)
  colnames(lam_sub) <- paste0("lam_t_",pred_times)
  
  pred_sub <- matrix(c(rep(case_res, each=nrow(res_sub))), nrow=nrow(res_sub), ncol=length(pred_times))
  colnames(pred_sub) <- paste0("pred_t_", pred_times)
  
  #Predictions from model (glmnet) = sum of coefficient estimates
  model_pred_sub <- matrix(rep(colSums(time_sub), nrow(res_sub)), ncol=length(pred_times), nrow=nrow(res_sub),
                           byrow = TRUE)
  colnames(model_pred_sub) <- paste0("explpred_t_", pred_times)
  
  res_sub$feature_value <- c(1, unlist(lapply(case_perm[1, res$feature], as.character)))
  temp_factor <- describe_level(colnames(perms), case_perm[1,], explainer$feature_type, explainer$bin_continuous, explainer$bin_cuts, explainer$feature_distribution)
  res_sub$feature_desc <- c("(Intercept)", temp_factor$feature_desc)
  res_sub$feature_levels <- c(NA, temp_factor$feature_levels) 
  res_sub$case <- rownames(dat)
  res_new <- cbind(res_sub, pred_sub, time_sub, model_pred_sub, r2_sub, lam_sub)
  return(res_new)  
}

explain_Quad <- function(explain_dat_model) {
  dat <- explain_dat_model$dat
  perms <- explain_dat_model$perms
  nperms <- explain_dat_model$n_permutations
  pred_times <- explain_dat_model$pred_times
  case_perm <- explain_dat_model$case_perm
  case_res <- explain_dat_model$case_res
  explainer <- explain_dat_model$explainer 
  
  res <- model_permutations_glmnet(explain_dat_model$x_full, 
                                   explain_dat_model$y_full, 
                                   explain_dat_model$weights_full, 
                                   method="Quad",n_permutations=nperms, pred_times=pred_times)
  
  res_sub <- cbind(as.numeric(perms[1,]), 
                   res$feature_weight[res$feature%in%colnames(perms)],
                   res$feature_weight[grep("^t_",res$feature)],
                   res$feature_weight[grep("^t2_",res$feature)])
  res_sub <- rbind(c(1, res$feature_weight[res$feature%in%c("(Intercept)","t","t2")]),
                   res_sub)
  
  #R2 is the same for all timepoints (from a single model)
  r2_sub <- data.frame(t(res$model_r2))
  colnames(r2_sub) <- paste0("r2_t_", pred_times)
  
  lam_sub <- matrix(rep(res$model_lam[1], nrow(res_sub)*length(pred_times)), ncol=length(pred_times), nrow=nrow(res_sub))
  colnames(lam_sub) <- paste0("lam_t_", pred_times)
  
  time_sub <- t(apply(res_sub, 1, function(x) {
    val <- rep(x[1], length(pred_times))
    coef_est <- x[2:4]
    coef_dat <- model.matrix(~val+ val:pred_times + val:I(pred_times^2)-1)
    ret <- coef_dat %*% coef_est
    return(ret)
  }))
  colnames(time_sub) <- paste0("value_t_",pred_times)
  
  pred_sub <- matrix(c(rep(case_res, each=nrow(res_sub))), nrow=nrow(res_sub), ncol=length(pred_times))
  colnames(pred_sub) <- paste0("pred_t_", pred_times)
  
  #Predictions from model (glmnet) = sum of coefficient estimates
  model_pred_sub <- matrix(rep(colSums(time_sub), nrow(res_sub)), ncol=length(pred_times), nrow=nrow(res_sub),
                           byrow = TRUE)
  colnames(model_pred_sub) <- paste0("explpred_t_", pred_times)
  
  res_sub2 <- res$feature[which(res$feature%in%colnames(perms))]
  res_new <- data.frame("feature"=c("(Intercept)", res_sub2))
  res_new$feature_value <- c(1, unlist(lapply(case_perm[1, which(colnames(perms)%in%res$feature)], as.character)))
  temp_factor <- describe_level(res_new$feature, case_perm[1,], explainer$feature_type, explainer$bin_continuous, explainer$bin_cuts, explainer$feature_distribution)
  res_new$feature_desc <- temp_factor$feature_desc
  res_new$feature_levels <- temp_factor$feature_levels
  res_new$case <- rownames(dat)
  res_new <- cbind(res_new, pred_sub, time_sub, model_pred_sub, r2_sub, lam_sub)
  return(res_new)  
}

explain_RMTLweighted <- function(explain_dat_model, Regularization="Lasso", G=NULL) {
  dat <- explain_dat_model$dat
  perms <- explain_dat_model$perms
  pred_times <- explain_dat_model$pred_times
  case_perm <- explain_dat_model$case_perm
  case_res <- explain_dat_model$case_res
  explainer <- explain_dat_model$explainer 
  
  res <- model_permutations_RMTLweighted(explain_dat_model$x_list, 
                                         explain_dat_model$y_list, 
                                         explain_dat_model$weights_list,
                                         Regularization=Regularization, 
                                         G=G)

  time_sub <- matrix(unlist(apply(res$feature_weight,2,function(x) x*c(1,as.numeric(perms[1,])))), nrow=nrow(res$feature_weight), ncol=length(pred_times))
  colnames(time_sub) <- paste0("value_t_",pred_times)
  
  pred_sub <- matrix(c(rep(case_res, each=nrow(res$feature_weight))), nrow=nrow(res$feature_weight), ncol=length(pred_times))
  colnames(pred_sub) <- paste0("pred_t_", pred_times)
  
  #R2 is the same for all timepoints (from a single model)
  r2_sub <- data.frame(t(res$model_r2))
  colnames(r2_sub) <- paste0("r2_t_", pred_times)
  lam_sub <- matrix(rep(res$model_lam, nrow(res$feature_weight)*length(pred_times)), ncol=length(pred_times), nrow=nrow(res$feature_weight))
  colnames(lam_sub) <- paste0("lam_t_", pred_times)
  
  #Predictions from model (glmnet) = sum of coefficient estimates
  model_pred_sub <- matrix(rep(colSums(time_sub), nrow(res$feature_weight)), ncol=length(pred_times), nrow=nrow(res$feature_weight),
                           byrow = TRUE)
  colnames(model_pred_sub) <- paste0("explpred_t_", pred_times)
  
  res_new <- data.frame("feature"=res$feature)
  res_new$feature_value <- c(1, unlist(lapply(case_perm[1, res$feature[-1]], as.character)))
  temp_factor <- describe_level(res_new$feature, case_perm[1,], explainer$feature_type, explainer$bin_continuous, explainer$bin_cuts, explainer$feature_distribution)
  res_new$feature_desc <- temp_factor$feature_desc
  res_new$feature_levels <- temp_factor$feature_levels
  res_new$case <- rownames(dat)
  res_new <- cbind(res_new, pred_sub, time_sub, model_pred_sub, r2_sub, lam_sub)
  return(res_new)
}


# RMTL functions ----------------------------------------------------------
#Modified functions: Least square
LS_grad_eval_KS <- function( w, c, x, y, d){
  grad_w <-  t(x) %*% (d*(x %*% w + c - y)) / nrow(x)
  grad_c <-  mean(d*(x %*% w + c -y))
  funcVal <- 0.5 * mean(d*(y - x %*% w -c)^2)
  return(list(grad_w, grad_c, funcVal))
}

LS_funcVal_eval_KS <- function ( w, c, x, y, d){
  return(0.5 * mean(d*(y - x %*% w -c)^2))
}


#least-square solver for regression
LS_Lasso_KS <- function (X, Y, lam1, lam2, opts, weights){
  # private functions
  gradVal_eval <- function (W, C){
    r <- lapply(c(1:task_num), function(x)
      LS_grad_eval_KS(W[, x], C[x], X[[x]], Y[[x]], weights[[x]]))
    grad_W <- sapply(r, function(x)x[[1]]) + 2* lam2 * W
    grad_C <- sapply(r, function(x)x[[2]])
    funcVal = sum(sapply(r, function(x)x[[3]])) + lam2 * norm(W, 'f')^2
    return(list(grad_W, grad_C, funcVal))
  }
  
  funVal_eval <- function (W, C){
    return(sum(sapply(c(1:task_num), function(x)
      LS_funcVal_eval_KS(W[, x], C[x], X[[x]], Y[[x]], weights[[x]]))) +
        lam2 * norm(W, 'f')^2)
  }
  
  nonsmooth_eval <- function (W, lam1){
    return(lam1*sum(abs(W)))
  }
  
  # Main algorithm
  task_num <- length (X);
  dimension = dim(X[[1]])[2];
  Obj <- vector(); 
  
  #initialize a starting point
  if(opts$init==0){
    W0 <- matrix(0, nrow=dimension, ncol=task_num);
    C0 <- rep(0, task_num);
  }else if(opts$init==1){
    W0 <- opts$W0
    C0 <- opts$C0
  }    
  
  bFlag <- 0; 
  Wz <- W0;
  Cz <- C0;
  Wz_old <- W0;
  Cz_old <- C0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]
    
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      Wzp <- RMTL:::l1_projection(Ws - gWs/gamma, lam1 / gamma);
      Czp <- Cs - gCs/gamma;
      Fzp <- funVal_eval(Wzp, Czp);
      
      delta_Wzp <- Wzp - Ws;
      delta_Czp <- Czp - Cs;
      nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
      nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
      r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
      
      Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
        sum(delta_Czp * gCs) + gamma * r_sum
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
    }
    
    Wz_old = Wz;
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    Obj = c(Obj, Fzp + nonsmooth_eval(Wz, lam1));
    
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol)
        break;
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
  }
  
  W = Wzp;
  C = Czp;
  return(list(W=W, C=C, Obj=Obj))
}

LS_L21_KS <- function (X, Y, lam1, lam2, opts, weights){
  #------------------------------------------------------
  # private functions
  gradVal_eval <- function (W, C){
    r <- lapply(c(1:task_num), function(x)
      LS_grad_eval_KS(W[, x], C[x], X[[x]], Y[[x]], weights[[x]]))
    grad_W <- sapply(r, function(x)x[[1]]) + 2* lam2 * W
    grad_C <- sapply(r, function(x)x[[2]])
    funcVal = sum(sapply(r, function(x)x[[3]])) + lam2 * norm(W, 'f')^2
    return(list(grad_W, grad_C, funcVal))
  }
  
  funVal_eval <- function (W, C){
    return(sum(sapply(c(1:task_num), function(x)
      LS_funcVal_eval_KS(W[, x], C[x], X[[x]], Y[[x]], weights[[x]]))) +
        lam2 * norm(W, 'f')^2)
  }
  
  nonsmooth_eval <- function (W, lam1){
    return(sum(sqrt(rowSums(W^2)))*lam1)
  }
  #------------------------------------------------------
  
  task_num <- length (X);
  dimension = dim(X[[1]])[2];
  Obj <- vector(); 
  
  #initialize a starting point
  if(opts$init==0){
    W0 <- matrix(0, nrow=dimension, ncol=task_num);
    C0 <- rep(0, task_num);
  }else if(opts$init==1){
    W0 <- opts$W0
    C0 <- opts$C0
  }    
  
  bFlag <- 0; 
  Wz <- W0;
  Cz <- C0;
  Wz_old <- W0;
  Cz_old <- C0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]
    
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      Wzp <- L21_projection(Ws - gWs/gamma, lam1 / gamma);
      Czp <- Cs - gCs/gamma;
      Fzp <- funVal_eval  (Wzp, Czp);
      
      delta_Wzp <- Wzp - Ws;
      delta_Czp <- Czp - Cs;
      nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
      nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
      r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
      
      Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
        sum(delta_Czp * gCs) + gamma * r_sum
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
      
    }
    
    Wz_old = Wz;
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    
    Obj = c(Obj, Fzp + nonsmooth_eval(Wz, lam1));
    
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol)
        break;
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
  }
  
  W = Wzp;
  C = Czp;
  return(list(W=W, C=C, Obj=Obj))
}

L21_projection <- function (W, lambda ){
  thresfold <- sqrt(rowSums(W^2))
  zeros <- which(thresfold==0)              
  temp <- 1 - lambda/thresfold
  temp <- ifelse(temp<0, 0, temp)
  Wp = matrix(rep(temp, ncol(W)), nrow=length(temp))*W
  Wp[zeros,] <- 0
  return(Wp)
}

LS_Graph_KS <- function (X, Y, G, lam1, lam2, opts, weights){
  #------------------------------------------------------
  # private functions
  gradVal_eval <- function (W, C){
    r <- lapply(c(1:task_num), function(x)
      LS_grad_eval_KS(W[, x], C[x], X[[x]], Y[[x]], weights[[x]]))
    grad_W <- sapply(r, function(x)x[[1]]) + 2*lam1*W %*% GGt + 2* lam2*W
    grad_C <- sapply(r, function(x)x[[2]])
    funcVal <- sum(sapply(r, function(x)x[[3]])) +
      lam1*norm(W%*%G,'f')^2 + lam2*norm(W,'f')^2
    return(list(grad_W, grad_C, funcVal))
  }
  
  funVal_eval <- function (W, C){
    return(sum(sapply(c(1:task_num),
                      function(x) LS_funcVal_eval_KS(W[, x], C[x], X[[x]], Y[[x]], weights[[x]]))) +
             lam1*norm(W%*%G,'f')^2 + lam2*norm(W,'f')^2)
  }
  #-------------------------------------------------------    
  
  # Main algorithm
  task_num <- length (X);
  dimension = dim(X[[1]])[2];
  Obj <- vector(); 
  
  #precomputation
  GGt <- G %*% t(G)
  
  #initialize a starting point
  if(opts$init==0){
    W0 <- matrix(0, nrow=dimension, ncol=task_num);
    C0 <- rep(0, task_num);
  }else if(opts$init==1){
    W0 <- opts$W0
    C0 <- opts$C0
  }    
  
  bFlag <- 0; 
  Wz <- W0;
  Cz <- C0;
  Wz_old <- W0;
  Cz_old <- C0;
  
  t <- 1;
  t_old <- 0;
  iter <- 0;
  gamma <- 1;
  gamma_inc <- 2;
  
  while (iter < opts$maxIter){
    alpha <- (t_old - 1) /t;
    
    Ws <- (1 + alpha) * Wz - alpha * Wz_old;
    Cs <- (1 + alpha) * Cz - alpha * Cz_old;
    
    # compute function value and gradients of the search point
    r <- gradVal_eval(Ws, Cs);
    gWs <- r[[1]]
    gCs <- r[[2]]
    Fs <- r[[3]]
    
    
    # the Armijo Goldstein line search scheme
    while (TRUE){
      Wzp <- Ws - gWs/gamma;
      Czp <- Cs - gCs/gamma;
      Fzp <- funVal_eval  (Wzp, Czp);
      
      delta_Wzp <- Wzp - Ws;
      delta_Czp <- Czp - Cs;
      nrm_delta_Wzp <- norm(delta_Wzp, 'f')^2;
      nrm_delta_Czp <- sum(delta_Czp * delta_Czp);
      r_sum <- (nrm_delta_Wzp+nrm_delta_Czp)/2;
      
      
      Fzp_gamma = Fs + sum(delta_Wzp* gWs) + 
        sum(delta_Czp * gCs) + gamma * r_sum
      
      if (r_sum <=1e-20){
        bFlag=1; 
        break;
      }
      
      if (Fzp <= Fzp_gamma) break else {gamma = gamma * gamma_inc}
      
    }
    
    Wz_old = Wz; 
    Cz_old = Cz;
    Wz = Wzp;
    Cz = Czp;
    
    Obj = c(Obj, Fzp );
    
    #test stop condition.
    if (bFlag) break;
    if (iter>=2){
      if (abs( Obj[length(Obj)] - Obj[length(Obj)-1] ) <= opts$tol)
        break;
    }
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
  }
  
  W = Wzp;
  C = Czp;
  return(list(W=W, C=C, Obj=Obj))
}

MTLweighted <- function(X, Y, weights, type="Regression", Regularization="L21",
                        Lam1=0.1, Lam1_seq=NULL, Lam2=0,
                        opts=list(init=0,  tol=10^-3,
                                  maxIter=1000), G=NULL, k=2)
{
  
  #test vilidity of input data
  if (!missing(X) & !missing(Y)){
    if (all(sapply(X, class)!="matrix")){
      X <- lapply(X, function(x){as.matrix(x)})
    }
    if (all(sapply(Y, class)!="matrix")){
      Y <- lapply(Y, function(x){as.matrix(x)})
    }
  }else{
    stop("data X or Y does not exist")
  }
  
  #Assume mean-regularized multi-task learning (Evgeniou and Pontil 2004) 
  #Change to consider other specifications
  #G t x t ;Gi,j={(t−1)/t if i=j; 1/t others
  # t <- length(Y)
  # G <- matrix(1/t, nrow=t, ncol=t) 
  # diag(G) <- (t-1)/t
  
  #test the validity of problem type
  method <- "LS"
  
  #test the validity of regularization 
  allRegularizations <- c("L21", "Lasso", "Trace", "Graph")
  if (is.element(Regularization, allRegularizations)){
    method <- paste0(method, "_", Regularization,"_KS")
  }else{
    stop("Regularization is not recognizable")}
  
  #test validity of Lam1 and Lam2
  if (Lam1<0) {stop("Lam1 must be positive")}
  if (Lam2<0) {stop("Lam2 must be positive")}
  
  #collect arguments 
  args <- list(X=X, Y=Y, lam1=Lam1, lam2=Lam2, opts=opts, weights=weights)
  
  if (Regularization=="CMTL"){
    if (k>0){
      args$k <- k
    }else(stop("for CMTL, k must be positive interger"))
  }
  if (Regularization=="Graph"){
    if(!is.null(G)){
      args$G <- G
    }else{stop("graph matrix G is not provided")}
  }
  
  # #call solver
  # #sparse routine
  # if( !is.null(Lam1_seq) & length(Lam1_seq)>0){
  #   #with warm start
  #   opt <- opts
  #   for (x in Lam1_seq){
  #     args$lam1 <- x
  #     m <- do.call(method, args)
  #     opt$init <- 1
  #     opt$W0 <- m$W
  #     opt$C0 <- m$C
  #     args$opts <- opt
  #     if (x<=Lam1) break
  #   }
  # } else {
  #   #without warm start
  #   m <- do.call(method, args)
  # }
  
  if (any(Regularization==c("L21", "Lasso", "Trace"))){
    #sparse routine
    if( !is.null(Lam1_seq) & length(Lam1_seq)>0){
      #with warm start
      opt <- opts
      for (x in Lam1_seq){
        args$lam1 <- x
        m <- do.call(method, args)
        opt$init <- 1
        opt$W0 <- m$W
        opt$C0 <- m$C
        args$opts <- opt
        if (x<=Lam1) break
      }
    } else {
      #without warm start
      m <- do.call(method, args)
    }
  } else if(any(Regularization==c("Graph", "CMTL"))){
    m <- do.call(method, args)
  }
  
  m$call <- match.call()
  m$Lam1 <- args$lam1
  m$Lam2 <- args$Lam2
  m$opts <- args$opts
  m$dim <- sapply(X, function(x)dim(x))
  m$type=type
  m$Regularization=Regularization
  m$method=method
  class(m) <- "MTL"
  return(m)
}

cvMTLweighted <- function(X, Y, weights, type="Classification", Regularization="L21",
                          Lam1_seq=10^seq(1,-4, -1), Lam2=0, G=NULL, k=2,
                          opts=list(init=0, tol=10^-3, maxIter=1000),
                          stratify=FALSE, nfolds=5, ncores=2, parallel=FALSE){
  #test vilidity of input data
  if (!missing(X) & !missing(Y)){
    if (all(sapply(X, class)!="matrix")){
      X <- lapply(X, function(x){as.matrix(x)})
    }
    if (all(sapply(Y, class)!="matrix")){
      Y <- lapply(Y, function(x){as.matrix(x)})
    }
  }else{
    stop("data X or Y doesnot exists")
  }
  task_num <- length(X)
  if(stratify & type=="Regression"){
    stop("stratified CV is not applicable to regression")}
  cvPar <- RMTL:::getCVPartition(Y, nfolds, stratify)
  
  #cv
  if (!parallel){
    cvm <- rep(0, length(Lam1_seq))
    for (i in 1:nfolds){
      cv_Xtr <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[1]][[x]], ])
      cv_Ytr <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[1]][[x]]])
      cv_weightstr <- lapply(c(1:task_num),
                             function(x) weights[[x]][cvPar[[i]][[1]][[x]]])
      cv_Xte <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[2]][[x]], ])
      cv_Yte <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[2]][[x]]])
      opt <- opts
      for (p_idx in 1: length(Lam1_seq)){
        m <- MTLweighted(X=cv_Xtr, Y=cv_Ytr, weights=cv_weightstr, type=type,
                         Regularization=Regularization, Lam1=Lam1_seq[p_idx],
                         Lam2=Lam2, opts=opt, k=k, G=G)
        #non sparse model training
        if (!is.element(Regularization, c("Graph", "CMTL"))){
          opt$init <- 1
          opt$W0 <- m$W
          opt$C0 <- m$C
        }
        cv_err <- calcError(m, newX=cv_Xte, newY=cv_Yte)
        cvm[p_idx] = cvm[p_idx]+cv_err
      }
    }
    cvm = cvm/nfolds
  } else {
    requireNamespace('doParallel')
    requireNamespace('foreach')
    doParallel::registerDoParallel(ncores)
    cvm <- foreach::foreach(i = 1:nfolds, .combine="cbind") %dopar%{
      cv_Xtr <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[1]][[x]], ])
      cv_Ytr <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[1]][[x]]])
      cv_weightstr <- lapply(c(1:task_num),
                             function(x) weights[[x]][cvPar[[i]][[1]][[x]]])
      cv_Xte <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[2]][[x]], ])
      cv_Yte <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[2]][[x]]])
      opt <- opts
      cvVec=rep(0, length(Lam1_seq))
      for (p_idx in 1: length(Lam1_seq)){
        m <- MTLweighted(X=cv_Xtr, Y=cv_Ytr, weights=cv_weightstr, type=type,
                         Regularization=Regularization, Lam1=Lam1_seq[p_idx],
                         Lam2=Lam2, opts=opt, k=k, G=G)
        #non sparse model training
        if (!is.element(Regularization, c("Graph", "CMTL")) ){
          opt$init <- 1
          opt$W0 <- m$W
          opt$C0 <- m$C
        }
        cv_err <- calcError(m, newX=cv_Xte, newY=cv_Yte)
        cvVec[p_idx] <- cv_err
      }
      return(cvVec)
    }
    cvm <- rowMeans(cvm)
  }
  
  best_idx <- which(cvm==min(cvm))[1]
  cv <- list(Lam1_seq=Lam1_seq, Lam1.min=Lam1_seq[best_idx],
             Lam2=Lam2, cvm=cvm)
  class(cv) <- "cvMTL"
  return(cv)
}


# Describing factor and level ---------------------------------------------
create_factor <- function(feature, dat, type, bin_continuous, bin_cuts) 
{
  temp <- sapply(feature, function(f) {
    if (type[[f]] == "logical") {
      paste0(f, " is ", tolower(as.character(unique(dat[[f]]))))
    }
    else if (type[[f]] %in% c("character", "factor")) {
      ret <- paste0(f, " = ", as.character(levels(dat[[f]])))
      levels <- 1:length(ret)
    }
    else if (bin_continuous) {
      cuts <- bin_cuts[[f]]
      cuts[1] <- -Inf
      cuts[length(cuts)] <- Inf
      cuts <- trimws(format(cuts, digits = 3))
      ret <- c(paste0(f, " <= ", cuts[2]))
      for(i in 2:(length(cuts)-2)) {
        ret <- c(ret, paste0(cuts[i], " < ", f, " <= ", cuts[i + 1]))
      }
      ret <- c(ret, paste0(cuts[length(cuts)-1], " < ", f))
      levels <- 1:length(ret)
    }
    else {
      ret <- f
      levels <- 1
    }
    return(list(ret, levels))
  })
  return(data.frame("feature_desc"=unlist(temp[1,]), "feature_levels"=unlist(temp[2,])))
}

describe_level <- function(feature, dat, type, bin_continuous, bin_cuts, feature_distribution) 
{
  temp <- sapply(feature, function(f) {
    if(f=="(Intercept)") {
      ret <- f
      levels <- NA
    } else if (!f%in%names(type)) {
      ret <- f
      levels <- NA 
    }
    else if (type[[f]] == "logical") {
      ret <- paste0(f, " is ", tolower(as.character(unique(dat[[f]]))))
      levels <- ifelse(tolower(as.character(unique(dat[[f]])))=="false", 1, 2)
    }
    else if (type[[f]] %in% c("character", "factor")) {
      ret <- paste0(f, " = ", as.character(dat[[f]]))
      levels <- which(names(feature_distribution[[f]])==dat[[f]])
    }
    else if (type[[f]] == "null") {
      ret <- f
      levels <- NA
    }
    else if (bin_continuous) {
      cuts <- bin_cuts[[f]]
      cuts[1] <- -Inf
      cuts[length(cuts)] <- Inf
      bin <- cut(dat[[f]], unique(cuts), labels = FALSE, 
                 include.lowest = TRUE)
      cuts <- trimws(format(cuts, digits = 3))
      if (bin == 1) {
        ret <- paste0(f, " <= ", cuts[bin + 1])
      } else if (bin == length(cuts) - 1) {
        ret <- paste0(cuts[bin], " < ", f)
      } else {
        ret <- paste0(cuts[bin], " < ", f, " <= ", cuts[bin + 1])
      }
      levels <- cut(as.numeric(dat[[f]]), breaks=unique(cuts), labels=FALSE, include.lowest=TRUE)
    }
    else {
      ret <- f
      levels <- NA
    }
    return(list(ret, levels))
  })
  return(data.frame("feature_desc"=unlist(temp[1,]), "feature_levels"=unlist(temp[2,])))
}

# Regularization input functions ------------------------------------------
explain_RMTLLasso <- function(explain_dat_model) {
  explain_RMTLweighted(explain_dat_model, Regularization="Lasso", G=NULL)
}

explain_RMTLL21 <- function(explain_dat_model) {
  explain_RMTLweighted(explain_dat_model, Regularization="L21", G=NULL)
}

explain_RMTLMean <- function(explain_dat_model) {
  t <- length(explain_dat_model$pred_times)
  
  G_mean <- matrix(1/t, nrow=t, ncol=t)
  diag(G_mean) <- (t-1)/t
  
  explain_RMTLweighted(explain_dat_model, Regularization="Graph", G=G_mean)
}

explain_RMTLSmooth <- function(explain_dat_model) {
  t <- length(explain_dat_model$pred_times)
  
  G_smooth <- matrix(0, nrow=t, ncol=t+1)
  diag(G_smooth) <- 1
  for(i in 1:t) {
    G_smooth[i, i+1] <- -1
  }
  explain_RMTLweighted(explain_dat_model, Regularization="Graph", G=G_smooth)
}