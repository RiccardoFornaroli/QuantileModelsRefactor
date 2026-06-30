################################################################################
# 0. CONFIGURAZIONE DELLE VARIABILI GLOBALI, SOGLIE E PARAMETRI
################################################################################
rm(list=ls())

# Sementi e riproducibilità
SET_SEED <- 1234
set.seed(SET_SEED)

# Parametri di Cross-Validazione
K_FOLDS <- 5

# Soglìe di Significatività Statistica, Costanti e Filtri CWI Rigidi
EPSILON_LOGIT  <- 0.001
VIF_THRESHOLD  <- 1.5
CWI_THRESHOLD  <- 0.75  # Unica doppia soglia applicata sia al modello globale che ai singoli quantili
MIN_NON_ZERO_PROP <- 0.10 # Filtro di ammissibilità metrica: almeno il 10% di valori diversi da zero

# Definizione dei Quantili (Taus)
TAUS_S <- seq(0.02, 0.98, 0.02)
TAUS_M <- list(
  Floor   = seq(0.02, 0.10, length.out = 20),
  Median  = seq(0.45, 0.55, length.out = 20),
  Ceiling = seq(0.90, 0.98, length.out = 20)
)

# Paths e File di Input/Output
INPUT_DATA_PATH         <- "data/METRICS_2026_LB_OK.csv"
FUNCTIONS_PATH          <- "Original/Functions.R"
OUTPUT_CSV_PATH         <- "results/optimized_results.csv"
OUTPUT_SUMMARY_CSV_PATH <- "results/tabellone_modelli_significativita.csv"
OUTPUT_HTML_PATH        <- "results/Report_Modelli_Quantilici.html"
OUTPUT_RDATA            <- "results/Tutto_Pronto_Per_Grafici.RData"
OUTPUT_PLOT_DIR         <- "results/plots/"

# Creazione cartelle di output se non existenti
if(!dir.exists("results")) dir.create("results")
if(!dir.exists(OUTPUT_PLOT_DIR)) dir.create(OUTPUT_PLOT_DIR)

################################################################################
# 1. CARICAMENTO PACCHETTI, PARALLELIZZAZIONE E FUNZIONI CUSTOM
################################################################################
source(FUNCTIONS_PATH)

package.list <- c("Hmisc", "quantreg", "fmsb", "Amelia", "sm", "caret", "ggplot2", "foreach", "doSNOW", "tidyr", "dplyr")
tmp.install  <- which(lapply(package.list, require, character.only = TRUE) == FALSE)
if(length(tmp.install) > 0) install.packages(package.list[tmp.install])
lapply(package.list, require, character.only = TRUE)

# Calcolo dei core disponibili della CPU
num_cores <- max(1, parallel::detectCores() - 1)

################################################################################
# 2. CARICAMENTO E PRE-PROCESSING DATI
################################################################################
dati <- read.table(INPUT_DATA_PATH, header = TRUE, sep = ",", na.strings = c("NA"))
dati[, c(6, 7)] <- rm_outlier_15iqr(dati[, c(6, 7)])
dati <- dati[complete.cases(dati), ]
dati$SUB <- factor(dati$SUB, levels = levels(factor(dati$SUB))[c(1, 3, 2)])

# Trasformazione Metriche biologiche
Metrics_or <- dati[c(10:length(dati))]
fine_tassonomia <- which(names(Metrics_or) == "Viviparidae")
Metrics_or[, 1:fine_tassonomia] <- log10(Metrics_or[, 1:fine_tassonomia] + 1)
Metrics_all <- Metrics_or

Metrics_all$logit_EPT_prop   <- log((Metrics_all$EPT_prop + EPSILON_LOGIT) / ((1 - Metrics_all$EPT_prop) + EPSILON_LOGIT))
Metrics_all$logit_OCH_prop   <- log((Metrics_all$OCH_prop + EPSILON_LOGIT) / ((1 - Metrics_all$OCH_prop) + EPSILON_LOGIT))
Metrics_all$logit_EPT_EPTOCH <- log((Metrics_all$EPT_EPTOCH + EPSILON_LOGIT) / ((1 - Metrics_all$EPT_EPTOCH) + EPSILON_LOGIT))

Metrics_all$EPT_abu   <- log10(Metrics_all$EPT_abu + 1)
Metrics_all$OCH_abu   <- log10(Metrics_all$OCH_abu + 1)
Metrics_all$ABUNDANCE <- log10(Metrics_all$ABUNDANCE + 1)

# --- APPLICAZIONE FILTRO DEL 10% VALORI NON-ZERO ---
non_zero_proportions <- colMeans(Metrics_all != 0, na.rm = TRUE)

valid_metrics_names  <- names(non_zero_proportions)[non_zero_proportions >= MIN_NON_ZERO_PROP]
excluded_metrics_names <- names(non_zero_proportions)[non_zero_proportions < MIN_NON_ZERO_PROP]

Metrics <- Metrics_all[, valid_metrics_names, drop = FALSE]

cat("\n========================================================================\n")
cat(" FILTRAGGIO METRICHE SPARSE (Soglia Minima Non-Zero:", MIN_NON_ZERO_PROP * 100, "%)\n")
cat(" Metriche Totali Caricate:", ncol(Metrics_all), "\n")
cat(" Metriche Ammesse al Calcolo:", ncol(Metrics), "\n")
cat(" Metriche Escluse (Troppi Zeri):", length(excluded_metrics_names), "\n")
if(length(excluded_metrics_names) > 0) {
  cat(" Lista Escluse:", paste(excluded_metrics_names, collapse = ", "), "\n")
}
cat("========================================================================\n\n")

# Sistemazione Variabili Idrauliche
Variables <- dati[c(6, 7)]
Variables <- abs(Variables)
Variables$VEL[Variables$VEL == 0] <- 0.01
Variables_ST <- as.data.frame(Variables)

GROUP <- factor(dati$GROUP)
SUB   <- dati$SUB

# Generazione dei fold per la Cross-Validazione bilanciata su GROUP
folds <- createFolds(GROUP, k = K_FOLDS, list = TRUE, returnTrain = FALSE)

################################################################################
# 3. DEFINIZIONE DEI MODELLI FUNZIONALI
################################################################################
Models <- list(
  null                = function(fit_data, taus_p) { rq(VAR ~ 1, tau = taus_p, method = "sfn", data = fit_data) },
  SBS_GRP             = function(fit_data, taus_p) { rq(VAR ~ 1 + SUB + GROUP, tau = taus_p, method = "sfn", data = fit_data) },
  LIN_SBS_GRP         = function(fit_data, taus_p) { rq(VAR ~ 1 + INVARy + GROUP + SUB, tau = taus_p, method = "sfn", data = fit_data) },
  LOG_SBS_GRP         = function(fit_data, taus_p) { rq(VAR ~ 1 + log10(INVARy) + GROUP + SUB, tau = taus_p, method = "sfn", data = fit_data) },
  EXP_SBS_GRP         = function(fit_data, taus_p) { rq(VAR ~ 1 + exp(INVARy) + GROUP + SUB, tau = taus_p, method = "sfn", data = fit_data) },
  QUA_SBS_GRP         = function(fit_data, taus_p) { rq(VAR ~ 1 + poly(INVARy, 2) + GROUP + SUB, tau = taus_p, method = "sfn", data = fit_data) },
  LIN_SBS_GRP_INT     = function(fit_data, taus_p) { rq(VAR ~ 1 + INVARy * GROUP + SUB, tau = taus_p, method = "sfn", data = fit_data) },
  LOG_SBS_GRP_INT     = function(fit_data, taus_p) { rq(VAR ~ 1 + log10(INVARy) * GROUP + SUB, tau = taus_p, method = "sfn", data = fit_data) },
  EXP_SBS_GRP_INT     = function(fit_data, taus_p) { rq(VAR ~ 1 + exp(INVARy) * GROUP + SUB, tau = taus_p, method = "sfn", data = fit_data) },
  QUA_SBS_GRP_INT     = function(fit_data, taus_p) { rq(VAR ~ 1 + poly(INVARy, 2) * GROUP + SUB, tau = fit_data, method = "sfn", data = fit_data) }
)
names(Models) <- c("null", "SBS_GRP", "LIN_SBS_GRP", "LOG_SBS_GRP", "EXP_SBS_GRP", 
                   "QUA_SBS_GRP", "LIN_SBS_GRP_INT", "LOG_SBS_GRP_INT", "EXP_SBS_GRP_INT", "QUA_SBS_GRP_INT")
NULL_Models     <- c(1:2)
NON_NULL_Models <- (1:length(Models))[(1:length(Models)) %nin% NULL_Models]

################################################################################
# 4. UNIVARIATE SHAPE SELECTION CON PARALLELIZZAZIONE E PROGRESS BAR
################################################################################
taus <- TAUS_S

cat(" FASE 1: Univariate Shape Selection (Calcolo parallelo su", num_cores, "core)\n")
cat(" Analisi delle forme funzionali ottimali con CV a", K_FOLDS, "fold per ogni metrica idonea\n")
cat("========================================================================\n\n")
flush.console()

cl <- parallel::makeCluster(num_cores)
doSNOW::registerDoSNOW(cl)

pb_shape <- txtProgressBar(max = length(Metrics), style = 3)
progress <- function(n) setTxtProgressBar(pb_shape, n)
opts <- list(progress = progress)

shape_results <- foreach(i = 1:length(Metrics), .packages = c("quantreg", "Hmisc"), .options.snow = opts) %dopar% {
  VAR <- Metrics[, i]
  metric_name <- names(Metrics)[i]
  wiINVAR <- list()
  metric_cv_values <- c()
  
  for (y in 1:ncol(Variables_ST)) {
    INVARy <- Variables_ST[, y]
    fit_data_all <- data.frame(VAR, INVARy, GROUP, SUB)
    fit_data_all <- fit_data_all[complete.cases(fit_data_all), ]
    
    cv_wi_list <- list()
    fold_top_wi_vals <- c()
    
    for(f in 1:K_FOLDS) {
      test_idx <- folds[[f]]
      train_data <- fit_data_all[-test_idx, ]
      
      fitted_model_cv <- list()
      for (m in NON_NULL_Models) {
        fit <- tryCatch({ Models[[m]](train_data, taus) }, error = function(e) { NULL })
        if(!is.null(fit)) {
          fitted_model_cv[[names(Models)[m]]] <- fit
        }
      }
      if(length(fitted_model_cv) > 0) {
        meanWI_cv <- mean_wi(fitted_model_cv, taus)
        cv_wi_list[[f]] <- meanWI_cv
        fold_top_wi_vals <- c(fold_top_wi_vals, meanWI_cv[1, 1])
      }
    }
    
    if(length(cv_wi_list) > 0) {
      wiINVAR[[colnames(Variables_ST)[y]]] <- cv_wi_list[[1]]
      metric_cv_values <- c(metric_cv_values, sd(fold_top_wi_vals))
    } else {
      wiINVAR[[colnames(Variables_ST)[y]]] <- NULL
      metric_cv_values <- c(metric_cv_values, NA)
    }
  }
  list(wiINVAR = wiINVAR, cv_stability = mean(metric_cv_values, na.rm = TRUE))
}
close(pb_shape)
parallel::stopCluster(cl)

wiSHAPE <- list()
cv_stability_list <- list()
for(i in 1:length(Metrics)) {
  m_name <- names(Metrics)[i]
  wiSHAPE[[m_name]] <- shape_results[[i]]$wiINVAR
  cv_stability_list[[m_name]] <- shape_results[[i]]$cv_stability
}

wi_selection <- list()
wi_selected  <- list()
for (v in 1:ncol(Metrics)) {
  selection <- as.data.frame(matrix(ncol = 5, nrow = ncol(Variables_ST), NA))
  names(selection) <- c("Index", "Variable", "Selected", "Model", "Wi")
  selection$Variable <- colnames(Variables_ST)
  selection$Index    <- 1:ncol(Variables_ST)
  
  for (i in 1:ncol(Variables_ST)) {
    selection$Selected[i] <- "YES"
    if(!is.null(wiSHAPE[[v]][[i]])) {
      selection$Model[i]    <- rownames(wiSHAPE[[v]][[i]])[1]
      selection$Wi[i]       <- wiSHAPE[[v]][[i]][1, 1]
    }
    selected <- selection[selection$Selected == "YES", ]
  }
  wi_selection[[names(Metrics)[v]]] <- selection
  wi_selected[[names(Metrics)[v]]]  <- selected
}

var_list <- list()
for (v in 1:ncol(Metrics)) {
  if (nrow(wi_selected[[v]]) > 1) {
    Variable <- try(vif_func(in_frame = Variables_ST[, wi_selected[[v]]$Index], thresh = VIF_THRESHOLD, trace = FALSE), silent = TRUE)
    if (inherits(Variable, "try-error") || length(Variable) == 0 || Variable[1] == "Error in if (vif_max < thresh) break : argument is of length zero\n") {
      Variable <- wi_selected[[v]]$Variable[which.max(wi_selected[[v]]$Wi)]
    }
    Model <- wi_selected[[v]]$Model[wi_selected[[v]]$Variable %in% Variable]
    Shape <- gsub('[M_]', '', substr(Model, start = 1, stop = 4))
    sel   <- cbind(Variable, Model, Shape)
    var_list[[v]] <- sel
  } else if (nrow(wi_selected[[v]]) > 0) {
    Variable <- colnames(Variables_ST)[wi_selected[[v]]$Index]
    Model    <- wi_selected[[v]]$Model[wi_selected[[v]]$Variable %in% Variable]
    Shape    <- gsub('[M_]', '', substr(Model, start = 1, stop = 4))
    sel      <- cbind(Variable, Model, Shape)
    var_list[[v]] <- sel
  } else {
    var_list[[v]] <- NA
  }
}
names(var_list) <- names(Metrics)

################################################################################
# 5. MULTIVARIATE FITTING & STEPWISE IN PARALLELO CON PROGRESS BAR
################################################################################
build_robust_formula <- function(response, main_effects, interactions, basic_covariates) {
  components <- c()
  if (length(main_effects) > 0 && any(main_effects != "")) components <- c(components, main_effects)
  if (length(interactions) > 0 && any(interactions != "")) components <- c(components, interactions)
  if (length(basic_covariates) > 0 && any(basic_covariates != "")) components <- c(components, basic_covariates)
  
  if (length(components) == 0) {
    return(as.formula(paste(response, "~ 1")))
  } else {
    return(as.formula(paste(response, "~ 1 +", paste(components, collapse = " + "))))
  }
}

cat("\n\n========================================================================\n")
cat(" FASE 2: Multivariate Stepwise Selection (Calcolo parallelo su", num_cores, "core)\n")
cat(" Selezione robusta all'indietro e verifica delle interazioni con i gruppi\n")
cat("========================================================================\n\n")
flush.console()

cl <- parallel::makeCluster(num_cores)
doSNOW::registerDoSNOW(cl)

pb_multi <- txtProgressBar(max = length(var_list), style = 3)
progress_multi <- function(n) setTxtProgressBar(pb_multi, n)
opts_multi <- list(progress = progress_multi)

var_final <- foreach(v = 1:length(var_list), .packages = c("quantreg", "Hmisc"), .options.snow = opts_multi) %dopar% {
  
  var_met_final <- list()
  var_met_final[[1]] <- names(var_list)[v]
  
  if (all(!is.na(var_list[[v]]))) {
    VAR <- Metrics[[v]]
    invar <- c()
    for (i in 1:nrow(var_list[[v]])) {
      if (var_list[[v]][i, 3] == "LIN") invar <- append(invar, as.character(var_list[[v]][i, 1]))
      if (var_list[[v]][i, 3] == "LOG") invar <- append(invar, paste0("log10(", var_list[[v]][i, 1], ")"))
      if (var_list[[v]][i, 3] == "EXP") invar <- append(invar, paste0("exp(", var_list[[v]][i, 1], ")"))
      if (var_list[[v]][i, 3] == "QUA") invar <- append(invar, paste0("poly(", var_list[[v]][i, 1], ",2)"))
    }
    
    var_met_final[[2]] <- wi_selected[[v]]
    var_met_final[[3]] <- invar
    fit_data <- data.frame(VAR, Variables_ST, GROUP, SUB)
    fit_data <- fit_data[complete.cases(fit_data), ]
    
    invar_int      <- paste0("(", invar, "):GROUP")
    invar_int_iniz <- invar_int
    
    full_formula <- build_robust_formula("VAR", invar, invar_int, c("GROUP", "SUB"))
    full <- rq(full_formula, tau = taus, method = "sfn", data = fit_data)
    mod  <- list(full)
    names(mod) <- c("start")
    meanWI     <- mean_wi(mod, taus)
    
    while (rownames(meanWI)[1] != "full" && length(invar_int) > 1) {
      mod <- list(full)
      for (b in 1:length(invar_int)) {
        res_formula <- build_robust_formula("VAR", invar, invar_int[-b], c("GROUP", "SUB"))
        res <- rq(res_formula, tau = taus, method = "sfn", data = fit_data)
        mod <- lappend(mod, res)
      }
      names(mod) <- c("full", invar_int)
      meanWI <- mean_wi(mod, taus)
      invar_int <- invar_int[invar_int != rownames(meanWI)[1]]
      full_formula <- build_robust_formula("VAR", invar, invar_int, c("GROUP", "SUB"))
      full <- rq(full_formula, tau = taus, method = "sfn", data = fit_data)
    }
    
    if (length(invar_int) == 1) {
      res_formula <- build_robust_formula("VAR", invar, c(), c("GROUP", "SUB"))
      last_res    <- rq(res_formula, tau = taus, method = "sfn", data = fit_data)
      mod         <- list(full, last_res)
      names(mod)  <- c("full", "last_res")
      meanWI      <- mean_wi(mod, taus)
      if(rownames(meanWI)[1] != "full") invar_int <- invar_int[-1]
    }
    
    all_int <- if(length(invar_int) == length(invar_int_iniz)) TRUE else FALSE
    invar_to_mantain <- which(invar_int_iniz %in% invar_int)
    
    full_formula <- build_robust_formula("VAR", invar, invar_int, c("GROUP", "SUB"))
    full <- rq(full_formula, tau = taus, method = "sfn", data = fit_data)
    mod  <- list(full)
    names(mod) <- c("start")
    meanWI     <- mean_wi(mod, taus)
    
    while (rownames(meanWI)[1] != "full" && (length(invar) - length(invar_to_mantain)) > 1 && all_int == FALSE) {
      mod = list(full)
      for (b in 1:length(invar)) {
        res_formula <- build_robust_formula("VAR", invar[-b], invar_int, c("GROUP", "SUB"))
        res <- rq(res_formula, tau = taus, method = "sfn", data = fit_data)
        mod <- lappend(mod, res)
      }
      names(mod) <- c("full", invar)
      meanWI     <- mean_wi(mod, taus)
      invar      <- invar[invar != rownames(meanWI)[1]]
      full_formula <- build_robust_formula("VAR", invar, invar_int, c("GROUP", "SUB"))
      full <- rq(full_formula, tau = taus, method = "sfn", data = fit_data)
    }
    
    if ((length(invar) - length(invar_to_mantain)) == 1) {
      res_formula <- build_robust_formula("VAR", invar[invar_to_mantain], invar_int, c("GROUP", "SUB"))
      last_res    <- rq(res_formula, tau = taus, method = "sfn", data = fit_data)
      mod         <- list(full, last_res)
      names(mod)  <- c("full", "last_res")
      meanWI      <- mean_wi(mod, taus)
      if(rownames(meanWI)[1] != "full") invar <- invar[invar_to_mantain]
    }
    
    full_formula    <- build_robust_formula("VAR", invar, invar_int, c("GROUP", "SUB"))
    full            <- rq(full_formula, tau = taus, method = "sfn", data = fit_data)
    res_grp_formula <- build_robust_formula("VAR", invar, invar_int, "SUB")
    res_group       <- rq(res_grp_formula, tau = taus, method = "sfn", data = fit_data)
    res_sub_formula <- build_robust_formula("VAR", invar, invar_int, "GROUP")
    res_sub         <- rq(res_sub_formula, tau = taus, method = "sfn", data = fit_data)
    
    mod        <- list(full, res_group, res_sub)
    names(mod) <- c("full", "res_group", "res_sub")
    meanWI     <- mean_wi(mod, taus)
    
    if (rownames(meanWI)[1] == "full") {
      sig <- "full"
    } else if (rownames(meanWI)[1] == "res_group") {
      full_formula    <- res_grp_formula
      full            <- rq(full_formula, tau = taus, method = "sfn", data = fit_data)
      res_sub_formula <- build_robust_formula("VAR", invar, invar_int, c())
      res_sub         <- rq(res_sub_formula, tau = taus, method = "sfn", data = fit_data)
      mod             <- list(full, res_sub)
      names(mod)      <- c("full", "res_sub")
      meanWI          <- mean_wi(mod, taus)
      sig             <- if(rownames(meanWI)[1] == "full") "no_group" else "no"
      if(sig == "no") { full <- res_sub; full_formula <- res_sub_formula }
    } else if (rownames(meanWI)[1] == "res_sub") {
      full_formula    <- res_sub_formula
      full            <- rq(full_formula, tau = taus, method = "sfn", data = fit_data)
      res_grp_formula <- build_robust_formula("VAR", invar, invar_int, c())
      res_group       <- rq(res_grp_formula, tau = taus, method = "sfn", data = fit_data)
      mod             <- list(full, res_group)
      names(mod)      <- c("full", "res_group")
      meanWI          <- mean_wi(mod, taus)
      sig             <- if(rownames(meanWI)[1] == "full") "no_sub" else "no"
      if(sig == "no") { full <- res_group; full_formula <- res_grp_formula }
    }
    
    var_met_final[[4]] <- invar
    var_met_final[[5]] <- invar_int
    var_met_final[[6]] <- sig
    var_met_final[[7]] <- full_formula
    
    res_formula <- if (sig == "full") as.formula("VAR ~ 1 + GROUP + SUB") else 
      if (sig == "no_group") as.formula("VAR ~ 1 + SUB") else 
        if (sig == "no_sub") as.formula("VAR ~ 1 + GROUP") else as.formula("VAR ~ 1")
    
    res        <- rq(res_formula, tau = taus, method = "sfn", data = fit_data)
    mod        <- list(full, res)
    names(mod) <- c("full", "res")
    meanWI     <- mean_wi(mod, taus)
    var_met_final[[8]] <- meanWI
    var_met_final[[9]] <- full
    
    for (i in 1:length(TAUS_M)) {
      taus_i <- TAUS_M[[i]]
      full_m <- rq(full_formula, tau = taus_i, method = "sfn", data = fit_data)
      res_m  <- rq(res_formula, tau = taus_i, method = "sfn", data = fit_data)
      mod_m  <- list(full_m, res_m)
      names(mod_m) <- c("full", "res")
      var_met_final <- lappend(var_met_final, mean_wi(mod_m, taus_i))
    }
    names(var_met_final) <- c("Metric", "Selected Variables", "Invar_Inizial", "Invar_Selected", "Invar_int_Selected", 
                              "Intercept_grouping", "Formula", "Significance", "Model", "Wi_floor", "Wi_median", "Wi_ceiling")
    return(var_met_final)
  } else {
    return(NULL)
  }
}
close(pb_multi)
parallel::stopCluster(cl)
names(var_final) <- names(Metrics)

cat("\n\n------------------------------------------------------------------------\n")
cat(" Elaborazione parallela completata con successo.\n")
cat(paste(" Generazione dei grafici basati sulla doppia soglia CWI >=", CWI_THRESHOLD, "...\n"))
cat("------------------------------------------------------------------------\n\n")
flush.console()

################################################################################
# 6. APPLICAZIONE CRITERIO CWI SU MODELLO GLOBALE E SINGOLI QUANTILI
################################################################################
results <- as.data.frame(matrix(ncol = 18, nrow = length(Metrics_all), NA))
names(results) <- c("Metrics", "Variables", "Interaction", "Intercept_grouping", 
                    "Wifull", "CWI_Global", "Model_Approved", "Floor_Sig", "Wi_Floor", "CWI_Floor", 
                    "Median_Sig", "Wi_Median", "CWI_Median", "Ceiling_Sig", "Wi_Ceiling", "CWI_Ceiling", 
                    "CV_Uncertainty", "Formula")

tabellone_export <- as.data.frame(matrix(ncol = 16, nrow = length(Metrics_all), NA))
names(tabellone_export) <- c("Metrica", "Variabili_Idrauliche", "Forma_Funzionale", "Interazioni_Attive", 
                             "Fattori_Intercetta", "Modello_Globale_Approvato",
                             "Floor_Sig_0.05", "CWI_Floor", "Median_Sig_0.50", "CWI_Median", 
                             "Ceiling_Sig_0.95", "CWI_Ceiling", "Numero_Quantili_Plottabili", "Formula_Riferimento")

for (i in 1:length(Metrics_all)) {
  m_name <- names(Metrics_all)[i]
  results$Metrics[i]    <- m_name
  tabellone_export$Metrica[i] <- m_name
  
  if (m_name %in% excluded_metrics_names) {
    results$Model_Approved[i] <- "EXCLUDED_SPARSE"
    tabellone_export$Modello_Globale_Approvato[i] <- "Esclusa (Troppi Zeri)"
    tabellone_export$Variabili_Idrauliche[i]       <- "N/A"
    tabellone_export$Forma_Funzionale[i]           <- "N/A"
    tabellone_export$Interazioni_Attive[i]         <- "N/A"
    tabellone_export$Fattori_Intercetta[i]         <- "N/A"
    tabellone_export$Numero_Quantili_Plottabili[i] <- 0
    next
  }
  
  if(!is.null(var_final[[m_name]])) {
    results$Variables[i]          <- paste(var_final[[m_name]]$Invar_Selected, collapse = " ")
    results$Interaction[i]        <- paste(var_final[[m_name]]$Invar_int_Selected, collapse = " ")
    results$Intercept_grouping[i] <- var_final[[m_name]]$Intercept_grouping
    results$Wifull[i]             <- var_final[[m_name]]$Significance["full", ][1]
    results$CV_Uncertainty[i]     <- cv_stability_list[[m_name]]
    
    results$CWI_Global[i]         <- results$Wifull[i] * (1 - results$CV_Uncertainty[i])
    results$Model_Approved[i]     <- ifelse(results$CWI_Global[i] >= CWI_THRESHOLD, "YES", "NO")
    tabellone_export$Modello_Globale_Approvato[i] <- results$Model_Approved[i]
    
    results$Wi_Floor[i]           <- var_final[[m_name]]$Wi_floor["full", ][1]
    results$CWI_Floor[i]          <- results$Wi_Floor[i] * (1 - results$CV_Uncertainty[i])
    results$Floor_Sig[i]          <- ifelse(results$Model_Approved[i] == "YES" & results$CWI_Floor[i] >= CWI_THRESHOLD, "YES", "NO")
    
    results$Wi_Median[i]          <- var_final[[m_name]]$Wi_median["full", ][1]
    results$CWI_Median[i]         <- results$Wi_Median[i] * (1 - results$CV_Uncertainty[i])
    results$Median_Sig[i]         <- ifelse(results$Model_Approved[i] == "YES" & results$CWI_Median[i] >= CWI_THRESHOLD, "YES", "NO")
    
    results$Wi_Ceiling[i]         <- var_final[[m_name]]$Wi_ceiling["full", ][1]
    results$CWI_Ceiling[i]        <- results$Wi_Ceiling[i] * (1 - results$CV_Uncertainty[i])
    results$Ceiling_Sig[i]        <- ifelse(results$Model_Approved[i] == "YES" & results$CWI_Ceiling[i] >= CWI_THRESHOLD, "YES", "NO")
    
    results$Formula[i]            <- paste(deparse(var_final[[m_name]]$Formula), collapse = "")
    
    raw_shape_info <- var_final[[m_name]]$`Selected Variables`
    tabellone_export$Variabili_Idrauliche[i] <- paste(raw_shape_info$Variable, collapse = "; ")
    tabellone_export$Forma_Funzionale[i]     <- paste(raw_shape_info$Shape, collapse = "; ")
    
    interazioni_attive_testo <- var_final[[m_name]]$Invar_int_Selected
    tabellone_export$Interazioni_Attive[i]  <- if(length(interazioni_attive_testo) > 0 && any(interazioni_attive_testo != "")) {
      paste(gsub("[:\\(\\)]", "", interazioni_attive_testo), collapse = "; ")
    } else { "Nessuna" }
    
    tabellone_export$Fattori_Intercetta[i]  <- var_final[[m_name]]$Intercept_grouping
    
    tabellone_export$Floor_Sig_0.05[i]      <- results$Floor_Sig[i]
    tabellone_export$CWI_Floor[i]           <- round(results$CWI_Floor[i], 4)
    
    tabellone_export$Median_Sig_0.50[i]     <- results$Median_Sig[i]
    tabellone_export$CWI_Median[i]          <- round(results$CWI_Median[i], 4)
    
    tabellone_export$Ceiling_Sig_0.95[i]    <- results$Ceiling_Sig[i]
    tabellone_export$CWI_Ceiling[i]         <- round(results$CWI_Ceiling[i], 4)
    
    num_plottabili <- sum(c(results$Floor_Sig[i], results$Median_Sig[i], results$Ceiling_Sig[i]) == "YES")
    tabellone_export$Numero_Quantili_Plottabili[i] <- num_plottabili
    tabellone_export$Formula_Riferimento[i]        <- results$Formula[i]
  }
}

write.csv(results, OUTPUT_CSV_PATH, row.names = FALSE)
write.csv(tabellone_export, OUTPUT_SUMMARY_CSV_PATH, row.names = FALSE)

################################################################################
# 7. GENERAZIONE PLOT CON CORREZIONE DELLA LEGENDA E DEL BUG DI POLY()
################################################################################
modelli_pronti <- list()
saved_plots <- list() 

for (i in 1:length(Metrics)) {
  m <- names(Metrics)[i]
  if (!is.null(var_final[[m]])) {
    
    if (results$Model_Approved[results$Metrics == m] == "NO") next
    
    taus_da_disegnare <- c()
    if(results$Floor_Sig[results$Metrics == m] == "YES")   taus_da_disegnare <- c(taus_da_disegnare, 0.05)
    if(results$Median_Sig[results$Metrics == m] == "YES")  taus_da_disegnare <- c(taus_da_disegnare, 0.50)
    if(results$Ceiling_Sig[results$Metrics == m] == "YES") taus_da_disegnare <- c(taus_da_disegnare, 0.95)
    
    if(length(taus_da_disegnare) > 0) {
      
      fit_data_all_vars <- as.data.frame(cbind(VAR = Metrics[, m], Variables_ST, GROUP, SUB))
      vera_formula      <- as.formula(results$Formula[results$Metrics == m])
      
      modelli_pronti[[m]] <- tryCatch({
        rq(vera_formula, tau = taus_da_disegnare, method = "sfn", data = fit_data_all_vars)
      }, error = function(e) { NULL })
      
      if(!is.null(modelli_pronti[[m]])) {
        
        var_idrauliche_coinvolte <- c()
        for (col_check in colnames(Variables_ST)) {
          if (any(grepl(col_check, var_final[[m]]$Invar_Selected))) {
            var_idrauliche_coinvolte <- c(var_idrauliche_coinvolte, col_check)
          }
        }
        
        has_interaction <- length(var_final[[m]]$Invar_int_Selected) > 0 && any(var_final[[m]]$Invar_int_Selected != "")
        
        if(length(var_idrauliche_coinvolte) > 0) {
          metric_saved_paths <- c()
          
          for(var_idraulica_attiva in var_idrauliche_coinvolte) {
            
            x_seq <- seq(min(Variables_ST[[var_idraulica_attiva]], na.rm=TRUE), 
                         max(Variables_ST[[var_idraulica_attiva]], na.rm=TRUE), 
                         length.out = 200)
            
            if(has_interaction) {
              grid_df <- expand.grid(X_Var = x_seq, GROUP = levels(GROUP))
              colnames(grid_df)[1] <- var_idraulica_attiva
              grid_df$SUB <- factor(levels(SUB)[1], levels = levels(SUB))
              
              altre_vars <- colnames(Variables_ST)[colnames(Variables_ST) != var_idraulica_attiva]
              if(length(altre_vars) > 0) {
                for(av in altre_vars) grid_df[[av]] <- median(Variables_ST[[av]], na.rm=TRUE)
              }
              
              preds_grid <- as.data.frame(predict(modelli_pronti[[m]], newdata = grid_df))
              colnames(preds_grid) <- paste0("Tau_", taus_da_disegnare)
              plot_df <- cbind(grid_df, preds_grid)
              
              plot_long <- plot_df %>%
                pivot_longer(cols = starts_with("Tau_"), names_to = "Quantile", values_to = "Prediction") %>%
                mutate(Quantile = gsub("Tau_", "Quantile ", Quantile))
              
            } else {
              grid_df <- data.frame(x_seq)
              colnames(grid_df) <- var_idraulica_attiva
              
              altre_vars <- colnames(Variables_ST)[colnames(Variables_ST) != var_idraulica_attiva]
              if(length(altre_vars) > 0) {
                for(av in altre_vars) grid_df[[av]] <- median(Variables_ST[[av]], na.rm=TRUE)
              }
              
              preds_grid <- matrix(0, nrow = nrow(grid_df), ncol = length(taus_da_disegnare))
              colnames(preds_grid) <- paste0("Tau_", taus_da_disegnare)
              
              for(k in 1:nrow(grid_df)) {
                synthetic_row <- fit_data_all_vars[1, , drop=FALSE]
                synthetic_row[[var_idraulica_attiva]] <- grid_df[[var_idraulica_attiva]][k]
                if(length(altre_vars) > 0) {
                  for(av in altre_vars) synthetic_row[[av]] <- median(Variables_ST[[av]], na.rm=TRUE)
                }
                
                p_single <- predict(modelli_pronti[[m]], newdata = synthetic_row)
                if(length(taus_da_disegnare) == 1) {
                  preds_grid[k, 1] <- p_single
                } else {
                  preds_grid[k, ] <- p_single
                }
              }
              
              plot_df <- cbind(grid_df, as.data.frame(preds_grid))
              plot_long <- plot_df %>%
                pivot_longer(cols = starts_with("Tau_"), names_to = "Quantile", values_to = "Prediction") %>%
                mutate(Quantile = gsub("Tau_", "Quantile ", Quantile))
            }
            
            is_logit_metric <- grepl("^logit_", m)
            if(is_logit_metric) {
              plot_long$Prediction <- (exp(plot_long$Prediction) * (1 + EPSILON_LOGIT) - EPSILON_LOGIT) / (1 + exp(plot_long$Prediction))
              plot_long$Prediction[plot_long$Prediction < 0] <- 0
              plot_long$Prediction[plot_long$Prediction > 1] <- 1
              y_label <- "Proporzione Attesa Metrica (Scala Originale 0-1)"
            } else {
              y_label <- "Valore Atteso Metrica"
            }
            
            if(has_interaction) {
              p <- ggplot(plot_long, aes(x = .data[[var_idraulica_attiva]], y = Prediction, 
                                         color = GROUP, linetype = Quantile, group = interaction(GROUP, Quantile))) +
                geom_line(size = 1.1) +
                labs(title = paste("Effetto Interazione Modellizzato (CWI >=", CWI_THRESHOLD, ") - Metrica:", m),
                     subtitle = paste("Curve back-transformed condizionate per gruppo su:", var_idraulica_attiva),
                     x = paste("Variabile Idraulica:", var_idraulica_attiva),
                     y = y_label) +
                scale_color_discrete(name = "Gruppo Idromorfologico") +
                scale_linetype_manual(values = c("Quantile 0.05" = "dotted", "Quantile 0.5" = "solid", "Quantile 0.95" = "dashed"), name = "Fascia Quantile") +
                guides(color = guide_legend(ncol = 4, title.position = "top"),
                       linetype = guide_legend(ncol = 1, title.position = "top")) +
                theme_minimal()
            } else {
              p <- ggplot(plot_long, aes(x = .data[[var_idraulica_attiva]], y = Prediction, color = Quantile, linetype = Quantile)) +
                geom_line(size = 1.25) +
                labs(title = paste("Andamento Quantili Convalida (CWI >=", CWI_THRESHOLD, ") - Metrica:", m),
                     subtitle = paste("Modello Globale Marginale (Non-Group Specific) su:", var_idraulica_attiva),
                     x = paste("Variabile Idraulica:", var_idraulica_attiva),
                     y = y_label) +
                scale_color_manual(values = c("Quantile 0.05" = "#377eb8", "Quantile 0.5" = "#4daf4a", "Quantile 0.95" = "#e41a1c"), name = "Fascia Quantile") +
                scale_linetype_manual(values = c("Quantile 0.05" = "dashed", "Quantile 0.5" = "solid", "Quantile 0.95" = "dashed"), name = "Fascia Quantile") +
                theme_minimal()
            }
            
            if(is_logit_metric) {
              p <- p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
            }
            
            p <- p + theme(legend.position = "bottom", legend.box = "vertical",
                           panel.grid.minor = element_blank(),
                           plot.title = element_text(face = "bold", size = 11))
            
            fig_path <- paste0(OUTPUT_PLOT_DIR, m, "_", var_idraulica_attiva, "_curve_plot.png")
            ggsave(filename = fig_path, plot = p, width = 7.5, height = 5.5)
            
            metric_saved_paths <- c(metric_saved_paths, fig_path)
          }
          saved_plots[[m]] <- metric_saved_paths
        }
      }
    }
  }
}

################################################################################
# 8. GENERAZIONE REPORT HTML CON SEZIONE DEDICATA ALLE METRICHE ESCLUSE
################################################################################
html_file <- file(OUTPUT_HTML_PATH, "w", encoding = "UTF-8")
writeLines("<!DOCTYPE html><html><head><meta charset='UTF-8'><title>Report Modelli Quantilici Rigorosi</title>", html_file)
writeLines("<style>
  body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 30px; background-color: #f8f9fa; color: #333; }
  h1, h2 { color: #2c3e50; }
  .metric-card { background: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); margin-bottom: 30px; }
  .metric-card-excluded { background: #fdf3f3; padding: 15px; border-radius: 8px; border-left: 5px solid #dc3545; margin-bottom: 15px; }
  table { width: 100%; border-collapse: collapse; margin-top: 15px; background: #fff; }
  th, td { border: 1px solid #dee2e6; padding: 12px; text-align: left; }
  th { background-color: #e9ecef; color: #495057; }
  tr:nth-child(even) { background-color: #f1f3f5; }
  .badge { padding: 5px 10px; border-radius: 4px; font-weight: bold; font-size: 12px; }
  .badge-yes { background-color: #d4edda; color: #155724; }
  .badge-no { background-color: #f8d7da; color: #721c24; }
  .badge-ex { background-color: #e2e3e5; color: #383d41; }
  .flex-container { display: flex; flex-direction: column; gap: 20px; }
  .table-container { width: 100%; }
  .plots-wrapper { display: flex; flex-wrap: wrap; gap: 15px; justify-content: flex-start; margin-top: 15px; }
  .image-container { flex: 1; min-width: 380px; max-width: 49%; text-align: center; background: #fafafa; border: 1px solid #eee; padding: 10px; border-radius: 6px; }
  img { max-width: 100%; height: auto; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.15); }
</style></head><body>", html_file)

writeLines("<h1>Report Quantile Regression & CWI Strict Double-Filtering</h1>", html_file)
writeLines(paste0("<p>Generato in data: ", Sys.time(), " | <strong>Filtro CWI &ge; ", CWI_THRESHOLD, "</strong> | <strong>Soglia Presenza Metrica &ge; ", MIN_NON_ZERO_PROP*100, "%</strong></p><hr>"), html_file)

writeLines("<h2>1. Metriche Idonee Analizzate</h2>", html_file)
for(i in 1:nrow(results)) {
  m <- results$Metrics[i]
  if(results$Model_Approved[i] == "EXCLUDED_SPARSE") next
  if(is.na(results$Wifull[i])) next
  
  has_plots <- results$Model_Approved[i] == "YES" && (results$Floor_Sig[i] == "YES" || results$Median_Sig[i] == "YES" || results$Ceiling_Sig[i] == "YES")
  
  writeLines("<div class='metric-card'>", html_file)
  writeLines(paste0("<h3>Metrica: ", m, "</h3>"), html_file)
  writeLines("<div class='flex-container'><div class='table-container'>", html_file)
  
  writeLines("<table>", html_file)
  writeLines(paste0("<tr><th>Variabili Selezionate</th><td>", results$Variables[i], "</td></tr>"), html_file)
  writeLines(paste0("<tr><th>Interazioni Attive</th><td>", results$Interaction[i], "</td></tr>"), html_file)
  writeLines(paste0("<tr><th>Modello Intercetta</th><td><span class='badge' style='background:#e2e3e5;'>", results$Intercept_grouping[i], "</span></td></tr>"), html_file)
  writeLines(paste0("<tr><th>Formula di Riferimento</th><td><code>", results$Formula[i], "</code></td></tr>"), html_file)
  writeLines(paste0("<tr><th>Incertezza Totale CV (&sigma;)</th><td>", round(results$CV_Uncertainty[i], 4), "</td></tr>"), html_file)
  
  badge_global = paste0("<span class='badge ", ifelse(results$Model_Approved[i] == "YES", "badge-yes", "badge-no"), "'>", results$Model_Approved[i], "</span>")
  writeLines(paste0("<tr><th><strong>Approvazione Modello Globale</strong></th><td>", badge_global, " (CWI Globale: ", round(results$CWI_Global[i], 4), ")</td></tr>"), html_file)
  
  badge_floor  = paste0("<span class='badge ", ifelse(results$Floor_Sig[i] == "YES", "badge-yes", "badge-no"), "'>", results$Floor_Sig[i], "</span>")
  badge_median = paste0("<span class='badge ", ifelse(results$Median_Sig[i] == "YES", "badge-yes", "badge-no"), "'>", results$Median_Sig[i], "</span>")
  badge_ceil   = paste0("<span class='badge ", ifelse(results$Ceiling_Sig[i] == "YES", "badge-yes", "badge-no"), "'>", results$Ceiling_Sig[i], "</span>")
  
  writeLines(paste0("<tr><th>Ammissibilità Floor (0.05)</th><td>", badge_floor, " (CWI Quantile: ", round(results$CWI_Floor[i],4), ")</td></tr>"), html_file)
  writeLines(paste0("<tr><th>Ammissibilità Mediana (0.50)</th><td>", badge_median, " (CWI Quantile: ", round(results$CWI_Median[i],4), ")</td></tr>"), html_file)
  writeLines(paste0("<tr><th>Ammissibilità Ceiling (0.95)</th><td>", badge_ceil, " (CWI Quantile: ", round(results$CWI_Ceiling[i],4), ")</td></tr>"), html_file)
  writeLines("</table></div>", html_file)
  
  writeLines("<div class='plots-wrapper'>", html_file)
  if(has_plots && m %in% names(saved_plots)) {
    for(plot_path in saved_plots[[m]]) {
      writeLines("<div class='image-container'>", html_file)
      writeLines(paste0("<img src='plots/", basename(plot_path), "' alt='Grafico Quantile Reale'>"), html_file)
      writeLines("</div>", html_file)
    }
  } else {
    # CORREZIONE DINAMICA RICHIESTA: Il testo ora stampa il valore esatto di CWI_THRESHOLD impostato all'inizio (0.75)
    writeLines(paste0("<p style='color:#dc3545; padding: 15px 0; font-weight:500;'>Grafici non generati: il modello globale o i singoli quantili non hanno raggiunto la stabilità critica di CWI &ge; ", CWI_THRESHOLD, ".</p>"), html_file)
  }
  writeLines("</div></div></div>", html_file)
}

if(length(excluded_metrics_names) > 0) {
  writeLines("<hr><h2>2. Metriche Escluse Automaticamente</h2>", html_file)
  writeLines(paste0("<p>Le seguenti metriche presentano valori diversi da zero in una percentuale di campioni inferiore al limite richiesto del ", MIN_NON_ZERO_PROP*100, "% e non sono state modellate per evitare bias strutturali:</p>"), html_file)
  
  for(ex_m in excluded_metrics_names) {
    prop_effettiva <- round(non_zero_proportions[ex_m] * 100, 2)
    writeLines("<div class='metric-card-excluded'>", html_file)
    writeLines(paste0("<strong>Metrica:</strong> <code>", ex_m, "</code> | <strong>Percentuale campioni presenti:</strong> ", prop_effettiva, "% (Soglia minima richiesta: ", MIN_NON_ZERO_PROP*100, "%)"), html_file)
    writeLines("</div>", html_file)
  }
}

writeLines("</body></html>", html_file)
close(html_file)

save(results, modelli_pronti, file = OUTPUT_RDATA)
cat("\nProcesso terminato. Script configurato con soglia 0.75 e testi dinamici aggiornati.\n")
