#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Project: Code for implementing LIME in survival framework (Setting 1 in simulation)
## Created: Mar 26, 2022
## Author: Krithika Suresh
## Description: To demonstrate the use of functions in "LimeSurvFns.R" and computation of performance metrics
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list=ls())

library(lime) #lime
library(survival) #coxph
library(simsurv) #simsurv
library(dynpred) #evalstep
library(randomForestSRC) #rfsrc
library(riskRegression) #Score
library(philentropy) #distance
library(dplyr)
library(tidyr) 

# Load functions ----------------------------------------------------------
dir0 <- getwd() #Set wd 
dir0 <- "/Users/ksuresh/OneDrive - Michigan Medicine/sureshk/Projects/Grants/ACSIRG_KS/Code/TimeDepSimulation"
source(paste0(dir0, "/LimeFnsSim.R"))

sim <- 1 #setting simulation seed 

# Setting 1 ---------------------------------------------------------------
# No time-varying covariate effects, no correlated covariates 
n <- 1100
maxt <- 3
pred_times <- c(0.5, 1, 1.5, 2, 2.5)
true_param <- list("gamma"=1.2, "lambda"=0.5, "betas"=c(0.5, -0.5, 0.4, 0.3, -0.3, rep(0,5)))
method <- c("cox","rsf") 
vis <- c("LM", "Quad", "RMTLSmooth")
mods <- paste(rep(method, each = length(vis)), vis, sep = "_")
dat_imp <- expand.grid("method"=method, "vis"=vis)
#Set to >0 (e.g. 20) if want to compute stability metrics (increases computation time)
stability_reps <- 0

set.seed(sim)
cov <- data.frame(id = 1:n, 
                  x1 = rbinom(n, 1, 0.3), 
                  x2 = rbinom(n, 1, 0.3),
                  x3 = rbinom(n, 1, 0.2), 
                  x4 = rbinom(n, 1, 0.5), 
                  x5 = rbinom(n, 1, 0.5), 
                  y1 = rbinom(n, 1, 0.1), 
                  y2 = rbinom(n, 1, 0.2), 
                  y3 = rbinom(n, 1, 0.3), 
                  y4 = rbinom(n, 1, 0.5),
                  y5 = rbinom(n, 1, 0.8))

dat <- simsurv(dist = "weibull", lambdas = true_param$lambda, gammas = true_param$gamma,
               betas = c(x1 = true_param$betas[1], x2 = true_param$betas[2], x3 = true_param$betas[3], 
                         x4 = true_param$betas[4], x5 = true_param$betas[5],
                         y1 = 0, y2 = 0, y3 = 0, y4 = 0, y5 =0), x=cov, maxt = maxt)
cens <- rexp(n, 1/6)
dat$status[cens < dat$eventtime] <- 0
dat$eventtime <- pmin(dat$eventtime, cens)

cov <- cov %>% dplyr::mutate(across(!id, ~factor(.x)))
dat <- merge(dat, cov)

dat <- subset(dat, select = -c(id))

# Create explanation model ------------------------------------------------
expl.dat <- dat[1:100,]
model.dat <- dat[101:nrow(dat),]

train.dat <- model.dat
test.dat <- expl.dat

# Fit various prediction models -------------------------------------------
#Cox
model_cox <- coxph(Surv(eventtime, status)~ ., train.dat, x=TRUE)
#RSF
model_rsf <- rfsrc(Surv(eventtime, status) ~ ., train.dat, ntree = 500, seed=sim, importance="permute")

# Compare model performance -----------------------------------------------
#Need to provide the probabilty of having the event OR model object
metrics <- Score(list("Cox"=model_cox,
                      "RSF"=model_rsf),
                 formula=Surv(eventtime, status)~1, data=test.dat, times=pred_times, summary="ibs")

AUC <- metrics$AUC$score
BS <- metrics$Brier$score

dat_coef <- data.frame("cox"=model_cox$coefficients,
                       "cox.sq"=model_cox$coefficients^2, 
                       "rsf"=model_rsf$importance, 
                       "rsf_rel"=model_rsf$importance/max(model_rsf$importance))

# Fit explainers ----------------------------------------------------
explainer_cox <- lime(dplyr::select(expl.dat,-c("eventtime","status")), model_cox, 
                      bin_continuous = TRUE)
explainer_rsf <- lime(dplyr::select(expl.dat,-c("eventtime","status")), model_rsf, 
                      bin_continuous = TRUE)

# Compute explanations ----------------------------------------------------
for(m in method) {
  explainer_temp <- eval(parse(text=paste0("explainer_",m)))
  for(v in vis) {
    assign(paste0("explanation_",m,"_",v), NULL) 
    assign(paste0("stabilitySR_",m,"_",v), NULL)
    assign(paste0("stabilityE_",m,"_",v), NULL)
    
    for(j in 1:nrow(expl.dat)) {
      assign(paste0("L_stability_",m,"_",v), vector(mode='list', length=length(pred_times)))
      
      for(s in 0:(stability_reps-1)) {
        
        explain_temp <- explain_dat(dplyr::select(expl.dat[j,],-c("eventtime","status")),
                                    explainer_temp,
                                    n_permutations=5000,
                                    n_features=100,
                                    n_labels=1,
                                    pred_times = pred_times, 
                                    scale=FALSE,
                                    seed=sim+s
        )
        
        assign(paste0("explain_",m,"_",v), eval(parse(text=paste0("explain_",v,"(explain_temp)"))))
        
        if(s==0) {
          assign(paste0("explanation_",m,"_",v), rbind(eval(parse(text=paste0("explanation_",m,"_",v))),
                                                       eval(parse(text=paste0("explain_",m,"_",v)))))
        }
        
        for(p in 1:length(pred_times)) {
          temp <- eval(parse(text=paste0("L_stability_",m,"_",v,"[[",p,"]]"))) 
          temp.ind <- which(names(eval(parse(text=paste0("explain_",m,"_",v))))==paste0("value_t_",pred_times[p]))
          eval(parse(text=paste0("L_stability_",m,"_",v,"[[",p,"]] <- rbind(temp, explain_",m,"_",v,"[-1,",temp.ind,"])")))
        }
      }
      
      temp_stability_SR <- NULL
      temp_stability_E <- NULL
      for(p in 1:length(pred_times)) {
        mat_SR <- philentropy::distance(eval(parse(text=paste0("L_stability_",m,"_",v,"[[",p,"]]"))), method="hassebrook", mute.message=TRUE)
        temp_stability_SR <- c(temp_stability_SR, mean(mat_SR[lower.tri(mat_SR)]))
        
        mat_E <- philentropy::distance(eval(parse(text=paste0("L_stability_",m,"_",v,"[[",p,"]]"))), method="euclidean", mute.message=TRUE)
        temp_stability_E <- c(temp_stability_E, mean(mat_E[lower.tri(mat_E)]))
      }
      
      eval(parse(text=paste0("stabilitySR_",m,"_",v," <- rbind(stabilitySR_",m,"_",v,", temp_stability_SR)")))
      eval(parse(text=paste0("stabilityE_",m,"_",v," <- rbind(stabilityE_",m,"_",v,", temp_stability_E)")))
      
    }
  }
}

# Individual summary explanation  -----------------------------------------
# compute an explanation over the individual's survival curve
# integrate over the survival curve and normalize it by the sum across all of their features 
# produces one metric per person per feature (averaged across time)
dat_expl <- NULL
dat_stability <- NULL
for(i in 1:nrow(dat_imp)) {
  temp_mod <- eval(parse(text=paste0("explanation_",dat_imp[i,1],"_",dat_imp[i,2])))
  
  #Integrate over an individual's survival curve
  avg_val <- by(temp_mod, temp_mod[,"case"], function(x) {
    apply(x, 1, function(z) {
      ys <- as.numeric(z[grep("value_t",names(temp_mod))])
      xs <- pred_times
      sum(diff(xs)*zoo::rollmean(ys,2))/(max(pred_times)-min(pred_times))
    })
  })
  temp_mod$summ_avg <- unlist(avg_val)
  
  absavg_val <- by(temp_mod, temp_mod[,"case"], function(x) {
    apply(x, 1, function(z) {
      ys <- as.numeric(z[grep("value_t",names(temp_mod))])
      ys2 <- abs(ys)
      xs <- pred_times
      sum(diff(xs)*zoo::rollmean(ys2,2))/(max(pred_times)-min(pred_times))
    })
  })
  temp_mod$summ_absavg <- unlist(absavg_val)
  
  sqavg_val <- by(temp_mod, temp_mod[,"case"], function(x) {
    apply(x, 1, function(z) {
      ys <- as.numeric(z[grep("value_t",names(temp_mod))])
      ys2 <- ys^2
      xs <- pred_times
      sum(diff(xs)*zoo::rollmean(ys2,2))/(max(pred_times)-min(pred_times))
    })
  })
  temp_mod$summ_sqavg <- unlist(sqavg_val)
  
  #measure of variability (similar to sd) to indicate whether there is evidence that some
  #variables have a lot of variation: are positive with high values and then
  #negative with high values 
  #VARIABILITY OF A FUNCTION (INTEGRAL OF THE SQUARE OF A FUNCTION)
  temp_mod$avgvalue_diff <- sqrt(temp_mod$summ_sqavg - temp_mod$summ_avg^2)
  
  #normalized values
  temp_mod2 <- temp_mod %>% 
    group_by(case) %>% 
    dplyr::summarize(summ_avg_norm = summ_avg/sum(summ_avg), 
                     summ_absavg_norm = summ_absavg/sum(summ_absavg), 
                     summ_sqavg_norm = summ_sqavg/sum(summ_sqavg), .groups="keep")
  
  temp_mod <- cbind(temp_mod, subset(temp_mod2, select=-case))
  
  # assign(paste0("explanation_",dat_imp[i,1],"_",dat_imp[i,2]), temp_mod)
  
  temp_mod$method <- dat_imp[i,"method"]
  temp_mod$vis <- dat_imp[i,"vis"]
  dat_expl <- rbind(dat_expl, temp_mod)
  
  #Combine stability info 
  temp_stabilitySR <- eval(parse(text=paste0("stabilitySR_",dat_imp[i,1],"_",dat_imp[i,2])))  
  temp_stabilityE <- eval(parse(text=paste0("stabilityE_",dat_imp[i,1],"_",dat_imp[i,2])))  
  temp_stabilitycomb <- rbind(colMeans(temp_stabilitySR), colMeans(temp_stabilityE))
  temp_stability <- data.frame(temp_stabilitycomb)
  temp_stability$stability <- c("SR","E")
  temp_stability$method <- dat_imp[i, "method"]
  temp_stability$vis <- dat_imp[i, "vis"]
  dat_stability <- rbind(dat_stability, temp_stability)
}
names(dat_stability)[1:5] <- paste0("stability_t_", pred_times)

# write.csv(table(dat$status)[1]/n*100, 
#           paste0(dir0, "/SimulationResults/S1/cens_",sim,".csv"), row.names=FALSE)
# write.csv(AUC,
#           paste0(dir0, "/SimulationResults/S1/AUC_",sim,".csv"), row.names=FALSE)
# write.csv(BS,
#           paste0(dir0, "/SimulationResults/S1/BS_",sim,".csv"), row.names=FALSE)
# write.csv(dat_expl,
#           paste0(dir0, "/SimulationResults/S1/expl_",sim,".csv"), row.names=FALSE)
# write.csv(dat_coef,
#           paste0(dir0, "/SimulationResults/S1/coef_",sim,".csv"), row.names=FALSE)
# write.csv(dat_stability,
#           paste0(dir0, "/SimulationResults/S1/stability_",sim,".csv"), row.names=FALSE)
