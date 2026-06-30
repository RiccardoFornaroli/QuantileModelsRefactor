############################################################
# IMPROVED 3.1 - FINAL STABLE PIPELINE
############################################################

rm(list = ls())
set.seed(1234)

source("Original/Functions.R")

packages <- c("Hmisc","quantreg","fmsb","Amelia","sm")
invisible(lapply(packages, function(p){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

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

fit_rq <- function(formula, data, tau = TAUS_SHAPE){
  tryCatch(
    rq(formula, tau = tau, method = "sfn", data = data),
    error = function(e) NULL
  )
}

safe_mean_wi <- function(models, taus){
  models <- Filter(Negate(is.null), models)
  if(length(models) == 0){
    return(matrix(NA, 1, 1, dimnames = list("full","wi")))
  }
  tryCatch(
    mean_wi(models, taus),
    error = function(e)
      matrix(NA, 1, 1, dimnames = list("full","wi"))
  )
}

get_wi <- function(x){
  if(is.null(x)) return(NA_real_)
  if(is.matrix(x) || is.data.frame(x)) return(as.numeric(x[1,1]))
  as.numeric(x)
}

compute_cv <- function(x){
  x <- x[!is.na(x)]
  if(length(x) < 2) return(NA_real_)
  sd(x)/mean(x)
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
# PROGRESS BAR
############################################################

n_metrics <- ncol(Metrics)
n_vars <- ncol(Variables_ST)

pb_metric <- txtProgressBar(min = 0, max = n_metrics, style = 3)
cat("\n>>> SHAPE SELECTION START\n")

############################################################
# SHAPE SELECTION (FIXED + PROGRESS WORKING)
############################################################

wiSHAPE <- vector("list", n_metrics)

for(i in 1:n_metrics){
  
  cat("\nMetric:", i, "/", n_metrics, "-", names(Metrics)[i], "\n")
  
  VAR <- Metrics[, i]
  wiINVAR <- vector("list", n_vars)
  
  pb_var <- txtProgressBar(min = 0, max = n_vars, style = 3)
  
  for(j in 1:n_vars){
    
    INVAR <- Variables_ST[, j]
    
    df <- na.omit(data.frame(VAR, INVAR, GROUP, SUB))
    
    if(nrow(df) < 10){
      wiINVAR[[j]] <- list(BASE = NA, CV = NA)
      setTxtProgressBar(pb_var, j)
      next
    }
    
    fits <- list(
      LIN = fit_rq(VAR ~ INVAR + GROUP + SUB, df),
      LOG = fit_rq(VAR ~ log10(INVAR) + GROUP + SUB, df),
      EXP = fit_rq(VAR ~ exp(INVAR) + GROUP + SUB, df),
      QUA = fit_rq(VAR ~ poly(INVAR,2) + GROUP + SUB, df)
    )
    
    wi_list <- lapply(fits, function(m)
      safe_mean_wi(list(m), TAUS_SHAPE))
    
    base <- sapply(wi_list, get_wi)
    cv <- compute_cv(base)
    
    wiINVAR[[j]] <- list(BASE = base, CV = cv)
    
    setTxtProgressBar(pb_var, j)
  }
  
  close(pb_var)
  
  wiSHAPE[[i]] <- wiINVAR
  
  setTxtProgressBar(pb_metric, i)
}

close(pb_metric)
cat("\n>>> SHAPE DONE\n")

############################################################
# SCORE + VARIABLE SELECTION
############################################################

score_shape <- function(x){
  b <- x$BASE
  if(all(is.na(b))) return(-Inf)
  
  cv <- x$CV
  if(is.na(cv)) cv <- 0
  
  0.7*mean(b, na.rm=TRUE) + 0.3*(1-cv)
}

pick_var <- function(wi_list){
  scores <- sapply(wi_list, score_shape)
  names(scores)[which.max(scores)]
}

var_list <- lapply(seq_along(Metrics), function(i){
  data.frame(
    Variable = pick_var(wiSHAPE[[i]]),
    Shape = "LIN"
  )
})
names(var_list) <- names(Metrics)

############################################################
# MULTIVARIATE MODEL
############################################################

build_model <- function(i){
  
  VAR <- Metrics[, i]
  invar <- var_list[[i]]$Variable
  
  df <- data.frame(VAR, Variables_ST, GROUP, SUB)
  
  int <- paste0("(",invar,"):GROUP")
  
  f <- as.formula(
    paste("VAR ~ 1 +", invar, "+", int, "+ GROUP + SUB")
  )
  
  full <- fit_rq(f, df)
  res  <- fit_rq(VAR ~ GROUP + SUB, df)
  
  list(formula=f, full=full, res=res,
       wi=safe_mean_wi(list(full,res), TAUS_SHAPE))
}

var_final <- lapply(seq_along(Metrics), build_model)
names(var_final) <- names(Metrics)

############################################################
# RESULTS TABLE (FIXED)
############################################################

results <- do.call(rbind, lapply(seq_along(var_final), function(i){
  
  x <- var_final[[i]]
  
  floor  <- safe_mean_wi(list(x$full,x$res), TAUS_GROUPS$Floor)
  median <- safe_mean_wi(list(x$full,x$res), TAUS_GROUPS$Median)
  ceil   <- safe_mean_wi(list(x$full,x$res), TAUS_GROUPS$Ceiling)
  
  data.frame(
    Metrics = names(Metrics)[i],
    Variables = var_list[[i]]$Variable,
    
    Wifull = get_wi(x$wi),
    
    Floor_Sig   = ifelse(get_wi(floor) > 0.999, "YES", "NO"),
    Median_Sig  = ifelse(get_wi(median) > 0.999, "YES", "NO"),
    Ceiling_Sig = ifelse(get_wi(ceil) > 0.999, "YES", "NO"),
    
    Wi_Floor = get_wi(floor),
    Wi_Median = get_wi(median),
    Wi_Ceiling = get_wi(ceil),
    
    Formula = paste(deparse(x$formula), collapse="")
  )
}))

############################################################
# MODELS FOR PLOTS
############################################################

modelli_pronti <- lapply(seq_along(var_final), function(i){
  
  df <- data.frame(VAR = Metrics[,i], Variables_ST, GROUP, SUB)
  
  if(any(results[i, c("Floor_Sig","Median_Sig","Ceiling_Sig")] == "YES")){
    fit_rq(var_final[[i]]$formula, df, tau=c(0.05,0.5,0.95))
  } else NULL
})

names(modelli_pronti) <- names(Metrics)

############################################################
# SAVE
############################################################

dir.create("results", showWarnings = FALSE)

write.csv(results, "results/improved3_final_results.csv", row.names=FALSE)

saveRDS(var_final, "results/improved3_final_var_final.rds")
saveRDS(modelli_pronti, "results/improved3_final_models.rds")