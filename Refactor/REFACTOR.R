############################################################
# 0. GLOBAL SETTINGS
############################################################

rm(list=ls())

set.seed(1234)

source("Functions.R")

packages <- c("Hmisc","quantreg","fmsb","Amelia","sm")

tau_shape <- seq(0.02, 0.98, 0.02)

taus_floor   <- (2:10)/100
taus_median  <- (45:55)/100
taus_ceiling <- (90:98)/100

taus_groups <- list(
  Floor = taus_floor,
  Median = taus_median,
  Ceiling = taus_ceiling
)

NULL_MODELS <- c(1,2)

eps <- 0.001


############################################################
# 1. PACKAGES + DATA
############################################################

tmp <- which(lapply(packages, require, character.only = TRUE) == FALSE)
if(length(tmp) > 0) install.packages(packages[tmp])
lapply(packages, require, character.only = TRUE)

dati <- read.table("METRICS_2026_LB_OK.csv", header=TRUE, sep=",")

dati[,c(6,7)] <- rm_outlier_15iqr(dati[,c(6,7)])
dati <- na.omit(dati)

dati$SUB <- factor(dati$SUB)

GROUP <- dati$GROUP
SUB   <- dati$SUB


############################################################
# 2. METRICS + TRANSFORMATIONS
############################################################

Metrics <- dati[,10:ncol(dati)]

Variables_ST <- abs(dati[,c(6,7)])
Variables_ST$VEL[Variables_ST$VEL == 0] <- 0.001

fine_tass <- which(names(Metrics) == "Viviparidae")
Metrics[,1:fine_tass] <- log10(Metrics[,1:fine_tass] + 1)

logit_fun <- function(x) log((x + eps) / (1 - x + eps))

Metrics$logit_EPT_prop   <- logit_fun(Metrics$EPT_prop)
Metrics$logit_OCH_prop   <- logit_fun(Metrics$OCH_prop)
Metrics$logit_EPT_EPTOCH <- logit_fun(Metrics$EPT_EPTOCH)

Metrics$EPT_abu   <- log10(Metrics$EPT_abu + 1)
Metrics$OCH_abu   <- log10(Metrics$OCH_abu + 1)
Metrics$ABUNDANCE <- log10(Metrics$ABUNDANCE + 1)


############################################################
# 3. MODEL DEFINITIONS
############################################################

Models <- list(
  
  null = function(df) rq(VAR ~ 1, tau=tau_shape, method="sfn", data=df),
  
  SBS_GRP = function(df) rq(VAR ~ 1 + SUB + GROUP, tau=tau_shape, method="sfn", data=df),
  
  LIN = function(df) rq(VAR ~ INVARy + GROUP + SUB, tau=tau_shape, method="sfn", data=df),
  
  LOG = function(df) rq(VAR ~ log10(INVARy) + GROUP + SUB, tau=tau_shape, method="sfn", data=df),
  
  EXP = function(df) rq(VAR ~ exp(INVARy) + GROUP + SUB, tau=tau_shape, method="sfn", data=df),
  
  QUA = function(df) rq(VAR ~ poly(INVARy,2) + GROUP + SUB, tau=tau_shape, method="sfn", data=df),
  
  LIN_INT = function(df) rq(VAR ~ INVARy * GROUP + SUB, tau=tau_shape, method="sfn", data=df),
  
  LOG_INT = function(df) rq(VAR ~ log10(INVARy) * GROUP + SUB, tau=tau_shape, method="sfn", data=df),
  
  EXP_INT = function(df) rq(VAR ~ exp(INVARy) * GROUP + SUB, tau=tau_shape, method="sfn", data=df),
  
  QUA_INT = function(df) rq(VAR ~ poly(INVARy,2) * GROUP + SUB, tau=tau_shape, method="sfn", data=df)
)

MODEL_NAMES <- names(Models)


############################################################
# 4. SHAPE SELECTION (WI)
############################################################

wiSHAPE <- vector("list", length(Metrics))
names(wiSHAPE) <- names(Metrics)

for(i in seq_along(Metrics)) {
  
  VAR <- Metrics[,i]
  
  wiINVAR <- lapply(Variables_ST, function(INVARy) {
    
    df <- data.frame(VAR, INVARy, GROUP, SUB)
    df <- na.omit(df)
    
    fits <- lapply(Models, function(m) m(df))
    names(fits) <- MODEL_NAMES
    
    mean_wi(fits, tau_shape)
  })
  
  wiSHAPE[[i]] <- wiINVAR
}


############################################################
# 5. VARIABLE SELECTION
############################################################

wi_selected <- list()

for(v in seq_along(Metrics)) {
  
  sel <- data.frame(
    Index = 1:ncol(Variables_ST),
    Variable = colnames(Variables_ST),
    Selected = "YES",
    Model = NA,
    Wi = NA
  )
  
  for(i in 1:ncol(Variables_ST)) {
    sel$Model[i] <- rownames(wiSHAPE[[v]][[i]])[1]
    sel$Wi[i]    <- wiSHAPE[[v]][[i]][1,1]
  }
  
  wi_selected[[v]] <- sel
}

names(wi_selected) <- names(Metrics)


############################################################
# 6. VARIABLE LIST (SHAPE + MODEL)
############################################################

var_list <- list()

for(v in seq_along(Metrics)) {
  
  sel <- wi_selected[[v]]
  
  Variable <- sel$Variable[which.max(sel$Wi)]
  Model    <- sel$Model[which.max(sel$Wi)]
  
  Shape <- substr(Model, 1, 4)
  Shape <- gsub("[_M]", "", Shape)
  
  var_list[[v]] <- data.frame(
    Variable = Variable,
    Model = Model,
    Shape = Shape
  )
}

names(var_list) <- names(Metrics)


############################################################
# 7. MULTIVARIATE MODELS
############################################################

var_final <- list()

for(v in seq_along(Metrics)) {
  
  if(all(is.na(var_list[[v]]))) next
  
  VAR <- Metrics[[v]]
  
  invars <- var_list[[v]]$Variable
  
  fit_data <- data.frame(VAR, Variables_ST, GROUP, SUB)
  
  full_formula <- as.formula(
    paste("VAR ~ 1 +",
          paste(invars, collapse=" + "),
          "+ GROUP + SUB")
  )
  
  full <- rq(full_formula, tau=tau_shape, method="sfn", data=fit_data)
  
  res  <- rq(VAR ~ 1 + GROUP + SUB, tau=tau_shape, method="sfn", data=fit_data)
  
  meanWI <- mean_wi(list(full=full, res=res), tau_shape)
  
  var_final[[v]] <- list(
    Metric = names(Metrics)[v],
    Variables = invars,
    Formula = full_formula,
    Model = full,
    Wi = meanWI
  )
}

names(var_final) <- names(Metrics)


############################################################
# 8. RESULTS TABLE
############################################################

results <- data.frame(
  Metrics = names(Metrics),
  Variables = sapply(var_final, function(x) paste(x$Variables, collapse=" ")),
  Formula = sapply(var_final, function(x) paste(deparse(x$Formula), collapse="")),
  stringsAsFactors = FALSE
)


############################################################
# 9. MODELS FOR PLOTS
############################################################

modelli_pronti <- list()

for(i in seq_along(Metrics)) {
  
  m <- names(Metrics)[i]
  
  if(i %in% names(var_final)) {
    
    fit_data <- data.frame(VAR = Metrics[,m], Variables_ST, GROUP, SUB)
    
    modelli_pronti[[m]] <- tryCatch({
      
      rq(var_final[[m]]$Formula,
         tau = c(0.05, 0.5, 0.95),
         method = "sfn",
         data = fit_data)
      
    }, error = function(e) NULL)
  }
}


############################################################
# 10. SAVE OUTPUT
############################################################

save(results,
     var_final,
     modelli_pronti,
     file = "REFRACTOR_models.RData")