################################################################################
# 0. CONFIGURAZIONE DELLE VARIABILI GLOBALI, SOGLIE E PARAMETRI
################################################################################
rm(list=ls())

# Sementi e riproducibilità
SET_SEED <- 1234
set.seed(SET_SEED)

# Parametri di Cross-Validazione
K_FOLDS <- 5

# Soglie di Significatività Statistica e Costanti
EPSILON_LOGIT <- 0.001
VIF_THRESHOLD <- 1.5
SIG_THRESHOLD <- 0.999

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

# Creazione cartelle di output se non esistenti
if(!dir.exists("results")) dir.create("results")
if(!dir.exists(OUTPUT_PLOT_DIR)) dir.create(OUTPUT_PLOT_DIR)

################################################################################
# 1. CARICAMENTO PACCHETTI, PARALLELIZZAZIONE E FUNZIONI CUSTOM
################################################################################
source(FUNCTIONS_PATH)

package.list <- c("Hmisc", "quantreg", "fmsb", "Amelia", "sm", "caret", "ggplot2", "foreach", "doSNOW")
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
Metrics <- Metrics_or

Metrics$logit_EPT_prop   <- log((Metrics$EPT_prop + EPSILON_LOGIT) / ((1 - Metrics$EPT_prop) + EPSILON_LOGIT))
Metrics$logit_OCH_prop   <- log((Metrics$OCH_prop + EPSILON_LOGIT) / ((1 - Metrics$OCH_prop) + EPSILON_LOGIT))
Metrics$logit_EPT_EPTOCH <- log((Metrics$EPT_EPTOCH + EPSILON_LOGIT) / ((1 - Metrics$EPT_EPTOCH) + EPSILON_LOGIT))

Metrics$EPT_abu   <- log10(Metrics$EPT_abu + 1)
Metrics$OCH_abu   <- log10(Metrics$OCH_abu + 1)
Metrics$ABUNDANCE <- log10(Metrics$ABUNDANCE + 1)

# Sistemazione Variabili Idrauliche
Variables <- dati[c(6, 7)]
Variables <- abs(Variables)
Variables$VEL[Variables$VEL == 0] <- 0.001
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
  QUA_SBS_GRP_INT     = function(fit_data, taus_p) { rq(VAR ~ 1 + poly(INVARy, 2) * GROUP + SUB, tau = sati, method = "sfn", data = fit_data) }
)
names(Models) <- c("null", "SBS_GRP", "LIN_SBS_GRP", "LOG_SBS_GRP", "EXP_SBS_GRP", 
                   "QUA_SBS_GRP", "LIN_SBS_GRP_INT", "LOG_SBS_GRP_INT", "EXP_SBS_GRP_INT", "QUA_SBS_GRP_INT")
NULL_Models     <- c(1:2)
NON_NULL_Models <- (1:length(Models))[(1:length(Models)) %nin% NULL_Models]

################################################################################
# 4. UNIVARIATE SHAPE SELECTION CON PARALLELIZZAZIONE E PROGRESS BAR
################################################################################
taus <- TAUS_S

# Stampa forzata a console del testo descrittivo PRIMA di inizializzare il cluster
cat("\n\n========================================================================\n")
cat(" FASE 1: Univariate Shape Selection (Calcolo parallelo su", num_cores, "core)\n")
cat(" Analisi delle forme funzionali ottimali con CV a", K_FOLDS, "fold per ogni metrica\n")
cat("========================================================================\n\n")
flush.console()

# Avvio cluster specifico per doSNOW
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

# Ricostruzione liste dai risultati paralleli
wiSHAPE <- list()
cv_stability_list <- list()
for(i in 1:length(Metrics)) {
  m_name <- names(Metrics)[i]
  wiSHAPE[[m_name]] <- shape_results[[i]]$wiINVAR
  cv_stability_list[[m_name]] <- shape_results[[i]]$cv_stability
}

# Estrazione Forme Selezionate
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

# Stampa forzata a console del testo descrittivo PRIMA del secondo blocco parallelo
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
    fit_data <- as.data.frame(cbind(VAR, Variables_ST, GROUP, SUB))
    
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
cat(" Generazione dei grafici e scrittura del report HTML...\n")
cat("------------------------------------------------------------------------\n\n")
flush.console()

################################################################################
# 6. COSTRUZIONE DATASET RISULTATI E ESPORTAZIONE TABELLONE RIASSUNTIVO (CSV)
################################################################################
results <- as.data.frame(matrix(ncol = 14, nrow = length(Metrics), NA))
names(results) <- c("Metrics", "Variables", "Interaction", "Intercept_grouping", 
                    "Wifull", "Floor_Sig", "Wi_Floor", "Median_Sig", "Wi_Median", "Ceiling_Sig", "Wi_Ceiling", 
                    "CV_Uncertainty", "CWI_Index", "Formula")

# Tabellone strutturato ad hoc richiesto dall'utente
tabellone_export <- as.data.frame(matrix(ncol = 13, nrow = length(Metrics), NA))
names(tabellone_export) <- c("Metrica", "Variabili_Idrauliche", "Forma_Funzionale", "Interazioni_Attive", 
                             "Fattori_Intercetta", "Significativo_In_Almeno_Un_Quantile",
                             "Floor_Sig_0.05", "Wi_Floor", "Median_Sig_0.50", "Wi_Median", 
                             "Ceiling_Sig_0.95", "Wi_Ceiling", "Formula_Riferimento")

for (i in 1:length(Metrics)) {
  m_name <- names(Metrics)[i]
  results$Metrics[i]    <- m_name
  tabellone_export$Metrica[i] <- m_name
  
  if(!is.null(var_final[[m_name]])) {
    # Dataset di calcolo interno standard
    results$Variables[i]          <- paste(var_final[[m_name]]$Invar_Selected, collapse = " ")
    results$Interaction[i]        <- paste(var_final[[m_name]]$Invar_int_Selected, collapse = " ")
    results$Intercept_grouping[i] <- var_final[[m_name]]$Intercept_grouping
    results$Wifull[i]             <- var_final[[m_name]]$Significance["full", ][1]
    
    results$Wi_Floor[i]           <- var_final[[m_name]]$Wi_floor["full", ][1]
    results$Floor_Sig[i]          <- ifelse(results$Wi_Floor[i] > SIG_THRESHOLD, "YES", "NO")
    
    results$Wi_Median[i]          <- var_final[[m_name]]$Wi_median["full", ][1]
    results$Median_Sig[i]         <- ifelse(results$Wi_Median[i] > SIG_THRESHOLD, "YES", "NO")
    
    results$Wi_Ceiling[i]         <- var_final[[m_name]]$Wi_ceiling["full", ][1]
    results$Ceiling_Sig[i]        <- ifelse(results$Wi_Ceiling[i] > SIG_THRESHOLD, "YES", "NO")
    
    results$CV_Uncertainty[i]     <- cv_stability_list[[m_name]]
    results$CWI_Index[i]          <- results$Wifull[i] * (1 - results$CV_Uncertainty[i])
    results$Formula[i]            <- paste(deparse(var_final[[m_name]]$Formula), collapse = "")
    
    # Riempimento del Tabellone Esteso pulito per Excel/CSV
    raw_shape_info <- var_final[[m_name]]$`Selected Variables`
    tabellone_export$Variabili_Idrauliche[i] <- paste(raw_shape_info$Variable, collapse = "; ")
    tabellone_export$Forma_Funzionale[i]     <- paste(raw_shape_info$Shape, collapse = "; ")
    
    interazioni_attive_testo <- var_final[[m_name]]$Invar_int_Selected
    tabellone_export$Interazioni_Attive[i]  <- if(length(interazioni_attive_testo) > 0 && any(interazioni_attive_testo != "")) {
      paste(gsub("[:\\(\\)]", "", interazioni_attive_testo), collapse = "; ")
    } else { "Nessuna" }
    
    tabellone_export$Fattori_Intercetta[i]  <- var_final[[m_name]]$Intercept_grouping
    
    tabellone_export$Floor_Sig_0.05[i]      <- results$Floor_Sig[i]
    tabellone_export$Wi_Floor[i]            <- round(results$Wi_Floor[i], 4)
    
    tabellone_export$Median_Sig_0.50[i]     <- results$Median_Sig[i]
    tabellone_export$Wi_Median[i]           <- round(results$Wi_Median[i], 4)
    
    tabellone_export$Ceiling_Sig_0.95[i]    <- results$Ceiling_Sig[i]
    tabellone_export$Wi_Ceiling[i]          <- round(results$Wi_Ceiling[i], 4)
    
    any_sig <- (results$Floor_Sig[i] == "YES" || results$Median_Sig[i] == "YES" || results$Ceiling_Sig[i] == "YES")
    tabellone_export$Significativo_In_Almeno_Un_Quantile[i] <- ifelse(any_sig, "SI", "NO")
    tabellone_export$Formula_Riferimento[i]                 <- results$Formula[i]
  }
}

# Esportazione di entrambi i file CSV nelle cartelle dei risultati
write.csv(results, OUTPUT_CSV_PATH, row.names = FALSE)
write.csv(tabellone_export, OUTPUT_SUMMARY_CSV_PATH, row.names = FALSE)

################################################################################
# 7. GENERAZIONE GRAFICI RIGIDI SU COVARIATE ORDINATE ED EQUISPAZIATE
################################################################################
modelli_pronti <- list()
saved_plots <- c()

for (i in 1:length(Metrics)) {
  m <- names(Metrics)[i]
  if (!is.null(var_final[[m]])) {
    
    taus_da_disegnare <- c()
    if(results$Floor_Sig[i] == "YES")   taus_da_disegnare <- c(taus_da_disegnare, 0.05)
    if(results$Median_Sig[i] == "YES")  taus_da_disegnare <- c(taus_da_disegnare, 0.50)
    if(results$Ceiling_Sig[i] == "YES") taus_da_disegnare <- c(taus_da_disegnare, 0.95)
    
    if(length(taus_da_disegnare) > 0) {
      
      # Ricostruzione e stima del modello finale globale sui dati completi
      fit_data_all_vars <- as.data.frame(cbind(VAR = Metrics[, m], Variables_ST, GROUP, SUB))
      vera_formula      <- as.formula(results$Formula[i])
      
      modelli_pronti[[m]] <- tryCatch({
        rq(vera_formula, tau = taus_da_disegnare, method = "sfn", data = fit_data_all_vars)
      }, error = function(e) { NULL })
      
      if(!is.null(modelli_pronti[[m]])) {
        
        # Individuazione della variabile idraulica attiva pulita (rimuovendo poly, log10, etc.)
        raw_string <- var_final[[m]]$Invar_Selected[1]
        var_idraulica_attiva <- NA
        for (col_check in colnames(Variables_ST)) {
          if (grepl(col_check, raw_string)) {
            var_idraulica_attiva <- col_check
            break
          }
        }
        
        if(!is.na(var_idraulica_attiva)) {
          # COSTRUZIONE GRID DI PREDIZIONE ORDINATA ED EQUISPAZIATA (200 punti)
          x_seq <- seq(min(Variables_ST[[var_idraulica_attiva]], na.rm=TRUE), 
                       max(Variables_ST[[var_idraulica_attiva]], na.rm=TRUE), 
                       length.out = 200)
          
          grid_df <- as.data.frame(matrix(ncol = ncol(Variables_ST), nrow = 200, NA))
          colnames(grid_df) <- colnames(Variables_ST)
          grid_df[[var_idraulica_attiva]] <- x_seq
          
          # Fissaggio delle altre variabili idrauliche sulla loro mediana campionaria
          altre_vars <- colnames(Variables_ST)[colnames(Variables_ST) != var_idraulica_attiva]
          if(length(altre_vars) > 0) {
            for(av in altre_vars) {
              grid_df[[av]] <- median(Variables_ST[[av]], na.rm=TRUE)
            }
          }
          
          # Fissaggio dei fattori sui livelli di riferimento più frequenti (o mediani)
          grid_df$GROUP <- factor(levels(GROUP)[1], levels = levels(GROUP))
          grid_df$SUB   <- factor(levels(SUB)[1], levels = levels(SUB))
          
          # Calcolo delle predizioni pulite lungo la sequenza ordinata
          preds_grid <- as.data.frame(predict(modelli_pronti[[m]], newdata = grid_df))
          colnames(preds_grid) <- paste0("Tau_", taus_da_disegnare)
          
          plot_df <- data.frame(X_Var = x_seq)
          plot_df <- cbind(plot_df, preds_grid)
          
          # Plotting geometrico puro (Senza geom_point)
          p <- ggplot(plot_df, aes(x = X_Var)) +
            labs(title = paste("Andamento Quantili Significativi - Metrica:", m),
                 x = paste("Variabile Idraulica:", var_idraulica_attiva),
                 y = "Valore Atteso Metrica") +
            theme_minimal() +
            theme(legend.position = "bottom")
          
          for(tau_curr in taus_da_disegnare) {
            col_name <- paste0("Tau_", tau_curr)
            p <- p + geom_line(aes(y = .data[[col_name]], color = col_name), size = 1.2)
          }
          
          p <- p + scale_color_brewer(palette = "Set1", name = "Quantili Sgn.")
          fig_path <- paste0(OUTPUT_PLOT_DIR, m, "_curve_plot.png")
          ggsave(filename = fig_path, plot = p, width = 7, height = 5)
          saved_plots[m] <- fig_path
        }
      }
    }
  }
}

################################################################################
# 8. GENERAZIONE REPORT HTML AUTOMATICO NATIVO
################################################################################
html_file <- file(OUTPUT_HTML_PATH, "w", encoding = "UTF-8")
writeLines("<!DOCTYPE html><html><head><meta charset='UTF-8'><title>Report Modelli Quantilici Ottimizzati</title>", html_file)
writeLines("<style>
  body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 30px; background-color: #f8f9fa; color: #333; }
  h1, h2 { color: #2c3e50; }
  .metric-card { background: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); margin-bottom: 30px; }
  table { width: 100%; border-collapse: collapse; margin-top: 15px; background: #fff; }
  th, td { border: 1px solid #dee2e6; padding: 12px; text-align: left; }
  th { background-color: #e9ecef; color: #495057; }
  tr:nth-child(even) { background-color: #f1f3f5; }
  .badge { padding: 5px 10px; border-radius: 4px; font-weight: bold; font-size: 12px; }
  .badge-yes { background-color: #d4edda; color: #155724; }
  .badge-no { background-color: #f8d7da; color: #721c24; }
  .flex-container { display: flex; flex-wrap: wrap; gap: 20px; }
  .table-container { flex: 1; min-width: 500px; }
  .image-container { flex: 1; min-width: 400px; text-align: center; }
  img { max-width: 100%; height: auto; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.15); }
</style></head><body>", html_file)

writeLines("<h1>Report Quantile Regression & Cross-Validation Stability (High-Speed Mode)</h1>", html_file)
writeLines(paste("<p>Generato in data:", Sys.time(), " | Fold Cross-Validazione k =", K_FOLDS, "</p><hr>"), html_file)

for(i in 1:nrow(results)) {
  m <- results$Metrics[i]
  if(is.na(results$Wifull[i])) next
  
  is_sig <- (results$Floor_Sig[i] == "YES" || results$Median_Sig[i] == "YES" || results$Ceiling_Sig[i] == "YES")
  
  writeLines("<div class='metric-card'>", html_file)
  writeLines(paste0("<h2>Metrica: ", m, "</h2>"), html_file)
  writeLines("<div class='flex-container'><div class='table-container'>", html_file)
  
  writeLines("<table>", html_file)
  writeLines(paste0("<tr><th>Variabili Selezionate</th><td>", results$Variables[i], "</td></tr>"), html_file)
  writeLines(paste0("<tr><th>Interazioni Attive</th><td>", results$Interaction[i], "</td></tr>"), html_file)
  writeLines(paste0("<tr><th>Modello Intercetta</th><td>", results$Intercept_grouping[i], "</td></tr>"), html_file)
  writeLines(paste0("<tr><th>Formula di Riferimento</th><td><code>", results$Formula[i], "</code></td></tr>"), html_file)
  writeLines(paste0("<tr><th>Peso Informativo (Wifull)</th><td>", round(results$Wifull[i], 4), "</td></tr>"), html_file)
  writeLines(paste0("<tr><th>Incertezza CV (&sigma;)</th><td>", round(results$CV_Uncertainty[i], 4), "</td></tr>"), html_file)
  writeLines(paste0("<tr><th><strong>Indice Composto CWI</strong></th><td><strong>", round(results$CWI_Index[i], 4), "</strong></td></tr>"), html_file)
  
  badge_floor  = paste0("<span class='badge ", ifelse(results$Floor_Sig[i] == "YES", "badge-yes", "badge-no"), "'>", results$Floor_Sig[i], "</span>")
  badge_median = paste0("<span class='badge ", ifelse(results$Median_Sig[i] == "YES", "badge-yes", "badge-no"), "'>", results$Median_Sig[i], "</span>")
  badge_ceil   = paste0("<span class='badge ", ifelse(results$Ceiling_Sig[i] == "YES", "badge-yes", "badge-no"), "'>", results$Ceiling_Sig[i], "</span>")
  
  writeLines(paste0("<tr><th>Significatività Floor (0.05)</th><td>", badge_floor, " (Wi: ", round(results$Wi_Floor[i],4), ")</td></tr>"), html_file)
  writeLines(paste0("<tr><th>Significatività Mediana (0.50)</th><td>", badge_median, " (Wi: ", round(results$Wi_Median[i],4), ")</td></tr>"), html_file)
  writeLines(paste0("<tr><th>Significatività Ceiling (0.95)</th><td>", badge_ceil, " (Wi: ", round(results$Wi_Ceiling[i],4), ")</td></tr>"), html_file)
  writeLines("</table></div>", html_file)
  
  writeLines("<div class='image-container'>", html_file)
  if(is_sig && m %in% names(saved_plots)) {
    writeLines(paste0("<img src='plots/", basename(saved_plots[m]), "' alt='Grafico Quantili Significativi'>"), html_file)
  } else {
    writeLines("<p style='color:#777; padding-top:40px;'>Nessun quantile significativo individuato. Grafico non generato.</p>", html_file)
  }
  writeLines("</div></div></div>", html_file)
}

writeLines("</body></html>", html_file)
close(html_file)

save(results, modelli_pronti, file = OUTPUT_RDATA)
cat("\nProcesso completato con successo.\n")
cat("Tabellone CSV esportato in:", OUTPUT_SUMMARY_CSV_PATH, "\n")
