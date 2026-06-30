############################################################
# REFRACTOR 2.0 STABLE + PROGRESS BAR (2 LEVELS)
############################################################

rm(list = ls())
set.seed(1234)

############################################################
# FUNCTIONS
############################################################

source("Original/Functions.R")

progress_init <- function(total, label="Progress") {
  list(
    pb = txtProgressBar(min = 0, max = total, style = 3),
    i = 0,
    total = total,
    label = label
  )
}

progress_update <- function(p, msg=NULL) {
  p$i <- p$i + 1
  setTxtProgressBar(p$pb, p$i)
  
  if(!is.null(msg) && p$i %% 10 == 0){
    cat("\n", p$label, ":", msg, "\n")
  }
  
  p
}

progress_close <- function(p){
  close(p$pb)
}

############################################################
# PACKAGES
############################################################

PACKAGES <- c("Hmisc","quantreg","fmsb","Amelia","sm")

missing <- PACKAGES[!sapply(PACKAGES, require, character.only=TRUE)]
if(length(missing)>0) install.packages(missing)
lapply(PACKAGES, library, character.only=TRUE)

############################################################
# CONSTANTS
############################################################

DATA_FILE <- "data/METRICS_2026_LB_OK.csv"

EPS <- 0.001
VEL_MIN <- 0.001
TAUS <- seq(0.02,0.98,0.02)

TAUS_GROUPS <- list(
  Floor=seq(0.02, 0.10, length.out = 20),
  Median=seq(0.45, 0.55, length.out = 20),
  Ceiling=seq(0.90, 0.98, length.out = 20)
)

############################################################
# DATA
############################################################

dati <- read.table(DATA_FILE, header=TRUE, sep=",", na.strings="NA")

dati[,c(6,7)] <- rm_outlier_15iqr(dati[,c(6,7)])
dati <- dati[complete.cases(dati),]

dati$SUB <- factor(dati$SUB, levels=levels(factor(dati$SUB))[c(1,3,2)])
GROUP <- factor(dati$GROUP)
SUB <- dati$SUB

Metrics <- dati[,10:ncol(dati)]
Variables_ST <- abs(dati[,c(6,7)])
Variables_ST$VEL[Variables_ST$VEL==0] <- VEL_MIN

############################################################
# MODEL LIST
############################################################

Models <- list(
  
  null = function(df) rq(VAR ~ 1, tau=TAUS, method="sfn", data=df),
  
  SBS_GRP = function(df) rq(VAR ~ 1 + SUB + GROUP, tau=TAUS, method="sfn", data=df),
  
  LIN = function(df) rq(VAR ~ 1 + INVARy + GROUP + SUB, tau=TAUS, method="sfn", data=df),
  
  LOG = function(df) rq(VAR ~ 1 + log10(INVARy) + GROUP + SUB, tau=TAUS, method="sfn", data=df),
  
  EXP = function(df) rq(VAR ~ 1 + exp(INVARy) + GROUP + SUB, tau=TAUS, method="sfn", data=df),
  
  QUA = function(df) rq(VAR ~ 1 + poly(INVARy,2) + GROUP + SUB, tau=TAUS, method="sfn", data=df),
  
  LIN_INT = function(df) rq(VAR ~ 1 + INVARy*GROUP + SUB, tau=TAUS, method="sfn", data=df),
  
  LOG_INT = function(df) rq(VAR ~ 1 + log10(INVARy)*GROUP + SUB, tau=TAUS, method="sfn", data=df),
  
  EXP_INT = function(df) rq(VAR ~ 1 + exp(INVARy)*GROUP + SUB, tau=TAUS, method="sfn", data=df),
  
  QUA_INT = function(df) rq(VAR ~ 1 + poly(INVARy,2)*GROUP + SUB, tau=TAUS, method="sfn", data=df)
)

MODEL_NAMES <- names(Models)
NON_NULL <- 3:10

fit_models <- function(df){
  lapply(NON_NULL, function(i){
    tryCatch(Models[[i]](df), error=function(e) NULL)
  })
}

############################################################
# SHAPE SELECTION (PROGRESS LEVEL 1)
############################################################

wiSHAPE <- vector("list", length(Metrics))
names(wiSHAPE) <- names(Metrics)

total_shape <- length(Metrics) * ncol(Variables_ST)
p_shape <- progress_init(total_shape, "SHAPE SELECTION")

counter <- 0

for(i in seq_along(Metrics)){
  
  VAR <- Metrics[[i]]
  
  wiINVAR <- vector("list", ncol(Variables_ST))
  
  for(j in seq_along(Variables_ST)){
    
    INVARy <- Variables_ST[[j]]
    
    df <- na.omit(data.frame(VAR, INVARy, GROUP, SUB))
    
    models <- fit_models(df)
    
    wiINVAR[[j]] <- mean_wi(models, TAUS)
    
    counter <- counter + 1
    p_shape <- progress_update(
      p_shape,
      paste(names(Metrics)[i], colnames(Variables_ST)[j])
    )
  }
  
  wiSHAPE[[i]] <- wiINVAR
}

progress_close(p_shape)

############################################################
# VARIABLE SELECTION
############################################################

var_list <- vector("list", length(Metrics))
names(var_list) <- names(Metrics)

for(i in seq_along(Metrics)){
  
  sel <- data.frame(
    Variable = names(Variables_ST),
    Wi = NA_real_,
    Model = NA_character_
  )
  
  for(j in seq_along(Variables_ST)){
    sel$Wi[j] <- wiSHAPE[[i]][[j]][1,1]
    sel$Model[j] <- rownames(wiSHAPE[[i]][[j]])[1]
  }
  
  best <- sel[which.max(sel$Wi),]
  
  var_list[[i]] <- data.frame(
    Variable = best$Variable,
    Model = best$Model,
    Shape = substr(best$Model,1,4),
    stringsAsFactors=FALSE
  )
}

############################################################
# MULTIVARIATE LOOP (PROGRESS LEVEL 2)
############################################################

var_final <- vector("list", length(Metrics))
names(var_final) <- names(Metrics)

p_multi <- progress_init(length(Metrics), "MULTIVARIATE")

for(i in seq_along(Metrics)){
  
  VAR <- Metrics[[i]]
  
  if(is.null(var_list[[i]])){
    p_multi <- progress_update(p_multi)
    next
  }
  
  invar <- var_list[[i]]$Variable
  shape <- var_list[[i]]$Shape
  
  if(shape=="LOG") invar <- paste0("log10(",invar,")")
  if(shape=="EXP") invar <- paste0("exp(",invar,")")
  if(shape=="QUA") invar <- paste0("poly(",invar,",2)")
  
  invar_int <- paste0("(",invar,"):GROUP")
  
  df <- data.frame(VAR, Variables_ST, GROUP, SUB)
  
  full_formula <- as.formula(
    paste("VAR ~ 1 +", invar, "+", invar_int, "+ GROUP + SUB")
  )
  
  full <- rq(full_formula, tau=TAUS, method="sfn", data=df)
  
  res <- rq(VAR ~ 1 + GROUP + SUB, tau=TAUS, method="sfn", data=df)
  
  sig <- rownames(mean_wi(list(full,res), TAUS))[1]
  
  var_final[[i]] <- list(
    Formula = full_formula,
    Sig = sig
  )
  
  p_multi <- progress_update(p_multi, names(Metrics)[i])
}

progress_close(p_multi)

############################################################
# RESULTS
############################################################

results <- data.frame(
  Metrics = names(Metrics),
  Sig = NA,
  Formula = NA
)

for(i in seq_along(var_final)){
  if(is.null(var_final[[i]])) next
  results$Sig[i] <- var_final[[i]]$Sig
  results$Formula[i] <- paste(deparse(var_final[[i]]$Formula), collapse="")
}

results$Sig[is.na(results$Sig)] <- "NO"

############################################################
# SAVE
############################################################

dir.create("results", showWarnings=FALSE)

saveRDS(wiSHAPE, "results/refactor_wiSHAPE.rds")
saveRDS(var_final, "results/refactor_var_final.rds")
write.csv(results, "results/refactor_results.csv", row.names=FALSE)