############################################################
# IMPROVED 3.0 - QUANTILE PIPELINE + MULTIVARIATE CV + PLOTS
############################################################

rm(list = ls())
set.seed(1234)

source("Original/Functions.R")

packages <- c("Hmisc","quantreg","fmsb","Amelia","sm")

invisible(lapply(packages, function(p){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)

############################################################
# SETTINGS
############################################################

TAUS_SHAPE <- seq(0.02, 0.98, 0.02)

TAUS_GROUPS <- list(
  Floor   = seq(0.02, 0.10, length.out = 20),
  Median  = seq(0.45, 0.55, length.out = 20),
  Ceiling = seq(0.90, 0.98, length.out = 20)
)

EPS <- 0.001
B_BOOT <- 50

############################################################
# SAFE FUNCTIONS
############################################################

safe_val <- function(x){
  if(is.null(x)) return(NA_real_)
  if(is.data.frame(x) || is.matrix(x)) x <- as.numeric(x[1,1])
  if(length(x) == 0) return(NA_real_)
  if(is.na(x)) return(NA_real_)
  as.numeric(x)
}

compute_cv <- function(vals){
  vals <- vals[!is.na(vals)]
  if(length(vals) < 2) return(NA_real_)
  sd(vals)/mean(vals)
}

fit_rq <- function(formula, data, tau = TAUS_SHAPE){
  rq(formula, tau = tau, method = "sfn", data = data)
}

pb_step <- function(pb, i, total, label = ""){
  setTxtProgressBar(pb, i)
  cat(sprintf("\r[%s] %d/%d", label, i, total))
}

############################################################
# DATA
############################################################

dati <- read.table("data/METRICS_2026_LB_OK.csv", header=TRUE, sep=",")

dati[,c(6,7)] <- rm_outlier_15iqr(dati[,c(6,7)])
dati <- na.omit(dati)

dati$SUB <- factor(dati$SUB)

GROUP <- dati$GROUP
SUB <- dati$SUB

Metrics <- dati[,10:ncol(dati)]
Variables_ST <- abs(dati[,c(6,7)])
Variables_ST$VEL[Variables_ST$VEL == 0] <- 0.001

fine_tass <- which(names(Metrics)=="Viviparidae")
Metrics[,1:fine_tass] <- log10(Metrics[,1:fine_tass] + 1)

logit <- function(x) log((x+EPS)/(1-x+EPS))

Metrics$logit_EPT_prop <- logit(Metrics$EPT_prop)
Metrics$logit_OCH_prop <- logit(Metrics$OCH_prop)
Metrics$logit_EPT_EPTOCH <- logit(Metrics$EPT_EPTOCH)

Metrics$EPT_abu <- log10(Metrics$EPT_abu + 1)
Metrics$OCH_abu <- log10(Metrics$OCH_abu + 1)
Metrics$ABUNDANCE <- log10(Metrics$ABUNDANCE + 1)

############################################################
# SHAPE SELECTION (WITH CV + PROGRESS)
############################################################

wiSHAPE <- vector("list", length(Metrics))
names(wiSHAPE) <- names(Metrics)

pb1 <- txtProgressBar(min = 0, max = length(Metrics), style = 3)
message("SHAPE SELECTION START")

for(i in seq_along(Metrics)){
  
  pb2 <- txtProgressBar(min = 0, max = ncol(Variables_ST), style = 3)
  message(sprintf("Metric %s (%d/%d)", names(Metrics)[i], i, length(Metrics)))
  
  VAR <- Metrics[,i]
  wiINVAR <- list()
  
  for(y in 1:ncol(Variables_ST)){
    
    pb_step(pb2, y, ncol(Variables_ST), "ShapeVar")
    
    INVARy <- Variables_ST[,y]
    df <- na.omit(data.frame(VAR, INVARy, GROUP, SUB))
    
    fits <- list(
      LIN = fit_rq(VAR ~ INVARy + GROUP + SUB, df),
      LOG = fit_rq(VAR ~ log10(INVARy) + GROUP + SUB, df),
      EXP = fit_rq(VAR ~ exp(INVARy) + GROUP + SUB, df),
      QUA = fit_rq(VAR ~ poly(INVARy,2) + GROUP + SUB, df)
    )
    
    wi_list <- lapply(fits, function(f){
      mean_wi(list(f), TAUS_SHAPE)
    })
    
    base <- sapply(wi_list, function(x) safe_val(x["full",1]))
    cv <- compute_cv(base)
    
    wiINVAR[[y]] <- data.frame(BASE = base, CV = cv)
  }
  
  names(wiINVAR) <- colnames(Variables_ST)
  wiSHAPE[[i]] <- wiINVAR
  
  close(pb2)
  setTxtProgressBar(pb1, i)
}

close(pb1)

############################################################
# VARIABLE SELECTION
############################################################

score_shape <- function(x){
  b <- safe_val(x$BASE)
  c <- safe_val(x$CV)
  if(is.na(b)) return(-Inf)
  if(is.na(c)) c <- 0
  0.7*b + 0.3*(1-c)
}

pick_var <- function(shape_list){
  s <- sapply(shape_list, score_shape)
  s[is.na(s)] <- -Inf
  names(s)[which.max(s)]
}

var_list <- lapply(seq_along(Metrics), function(i){
  data.frame(Variable = pick_var(wiSHAPE[[i]]), Shape = "LIN")
})

############################################################
# MULTIVARIATE MODELS
############################################################

build_model <- function(i){
  
  VAR <- Metrics[[i]]
  invar <- var_list[[i]]$Variable
  
  df <- data.frame(VAR, Variables_ST, GROUP, SUB)
  int <- paste0("(",invar,"):GROUP")
  
  formula <- as.formula(
    paste("VAR ~ 1 +", invar, "+", int, "+ GROUP + SUB")
  )
  
  full <- fit_rq(formula, df)
  res  <- fit_rq(VAR ~ GROUP + SUB, df)
  
  list(
    formula = formula,
    full = full,
    res = res,
    wi = mean_wi(list(full=full,res=res), TAUS_SHAPE)
  )
}

pb3 <- txtProgressBar(min=0, max=length(Metrics), style=3)
message("MULTIVARIATE FIT START")

var_final <- vector("list", length(Metrics))
names(var_final) <- names(Metrics)

for(i in seq_along(Metrics)){
  var_final[[i]] <- build_model(i)
  setTxtProgressBar(pb3, i)
}

close(pb3)

############################################################
# RESULTS TABLE
############################################################

results <- do.call(rbind, lapply(seq_along(var_final), function(i){
  
  x <- var_final[[i]]
  
  data.frame(
    Metrics = names(Metrics)[i],
    Variables = var_list[[i]]$Variable,
    
    Wifull = x$wi["full",1],
    
    Floor_Sig   = ifelse(mean_wi(list(x$full,x$res), TAUS_GROUPS$Floor)["full",1] > 0.999,"YES","NO"),
    Median_Sig  = ifelse(mean_wi(list(x$full,x$res), TAUS_GROUPS$Median)["full",1] > 0.999,"YES","NO"),
    Ceiling_Sig = ifelse(mean_wi(list(x$full,x$res), TAUS_GROUPS$Ceiling)["full",1] > 0.999,"YES","NO"),
    
    Formula = paste(deparse(x$formula), collapse="")
  )
}))

############################################################
# BOOTSTRAP CV
############################################################

pb4 <- txtProgressBar(min=0, max=length(Metrics), style=3)
message("BOOTSTRAP CV START")

model_cv <- numeric(length(Metrics))

for(i in seq_along(Metrics)){
  
  VAR <- Metrics[[i]]
  invar <- var_list[[i]]$Variable
  
  df <- data.frame(VAR, Variables_ST, GROUP, SUB)
  
  formula <- as.formula(
    paste("VAR ~ 1 +", invar, "+ (",invar,"):GROUP + GROUP + SUB")
  )
  
  boot <- replicate(B_BOOT, {
    
    idx <- sample(nrow(df), replace=TRUE)
    d <- df[idx,]
    
    m <- tryCatch(rq(formula, tau=TAUS_SHAPE, data=d),
                  error=function(e) NULL)
    
    if(is.null(m)) return(NA_real_)
    
    safe_val(mean_wi(list(m), TAUS_SHAPE)["full",1])
  })
  
  model_cv[i] <- compute_cv(boot)
  
  setTxtProgressBar(pb4, i)
}

close(pb4)

############################################################
# PLOTS (ONLY SIGNIFICANT MODELS)
############################################################

pb5 <- txtProgressBar(min=0, max=nrow(results), style=3)
message("PLOT GENERATION START")

for(i in seq_len(nrow(results))){
  
  pb_step(pb5, i, nrow(results), "Plots")
  
  if(!(results$Floor_Sig[i]=="YES" ||
       results$Median_Sig[i]=="YES" ||
       results$Ceiling_Sig[i]=="YES")) next
  
  df <- data.frame(VAR = Metrics[[i]], Variables_ST, GROUP, SUB)
  
  m <- tryCatch(
    fit_rq(var_final[[i]]$formula, df, tau=c(0.05,0.5,0.95)),
    error=function(e) NULL
  )
  
  if(is.null(m)) next
  
  png(paste0("results/plots/plot_", names(Metrics)[i], ".png"),
      width=1200, height=800)
  
  plot(m)
  
  title(main=paste("Metric:", names(Metrics)[i],
                   "\nVariable:", var_list[[i]]$Variable))
  
  dev.off()
}

close(pb5)

############################################################
# SAVE EVERYTHING
############################################################

stability <- data.frame(Metric=names(Metrics), Model_CV=model_cv)

write.csv(results, "results/improved3_results.csv", row.names=FALSE)
write.csv(stability, "results/improved3_stability.csv", row.names=FALSE)

saveRDS(var_final, "results/improved3_var_final.rds")
saveRDS(model_cv, "results/improved3_model_cv.rds")
