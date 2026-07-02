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
INPUT_DATA_PATH  <- "data/METRICS_2026_LB_OK.csv"
FUNCTIONS_PATH   <- "Original/Functions.R"
OUTPUT_CSV_PATH  <- "results/optimized_results.csv"
OUTPUT_RDATA     <- "results/Tutto_Pronto_Per_Grafici.RData"
OUTPUT_PLOT_DIR  <- "results/plots/"

# Creazione cartelle di output se non esistenti
if(!dir.exists("results")) dir.create("results")
if(!dir.exists(OUTPUT_PLOT_DIR)) dir.create(OUTPUT_PLOT_DIR)

################################################################################
# 1. CARICAMENTO PACCHETTI E FUNZIONI CUSTOM
################################################################################
source(FUNCTIONS_PATH)

package.list <- c("Hmisc", "quantreg", "fmsb", "Amelia", "sm", "caret", "ggplot2")
tmp.install  <- which(lapply(package.list, require, character.only = TRUE) == FALSE)
if(length(tmp.install) > 0) install.packages(package.list[tmp.install])
lapply(package.list, require, character.only = TRUE)

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
  QUA_SBS_GRP_INT     = function(fit_data, taus_p) { rq(VAR ~ 1 + poly(INVARy, 2) * GROUP + SUB, tau = taus_p, method = "sfn", data = fit_data) }
)

names(Models) <- c("null", "SBS_GRP", "LIN_SBS_GRP", "LOG_SBS_GRP", "EXP_SBS_GRP", 
                   "QUA_SBS_GRP", "LIN_SBS_GRP_INT", "LOG_SBS_GRP_INT", "EXP_SBS_GRP_INT", "QUA_SBS_GRP_INT")
NULL_Models     <- c(1:2)
NON_NULL_Models <- (1:length(Models))[(1:length(Models)) %nin% NULL_Models]

################################################################################
# 4. UNIVARIATE SHAPE SELECTION CON CROSS-VALIDAZIONE
################################################################################
taus <- TAUS_S
wiSHAPE <- list()

for (i in 1:length(Metrics)) {
  VAR <- Metrics[, i]
  metric <- names(Metrics[i])
  wiINVAR <- list()
  
  for (y in 1:ncol(Variables_ST)) {
    perc <- ((((i - 1) * ncol(Variables_ST)) + y) / (length(Metrics) * ncol(Variables_ST))) * 100
    cat(paste0("\rShape CV - Metric: ", i, "/", length(Metrics), " | Var: ", y, "/", ncol(Variables_ST), " | Progresso: ", round(perc, 1), "%"))
    
    INVARy <- Variables_ST[, y]
    fit_data_all <- data.frame(VAR, INVARy, GROUP, SUB)
    fit_data_all <- fit_data_all[complete.cases(fit_data_all), ]
    
    cv_wi_list <- list()
    
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
      }
    }
    
    if(length(cv_wi_list) > 0) {
      combined_wi <- cv_wi_list[[1]]
      wiINVAR[[colnames(Variables_ST)[y]]] <- combined_wi
    } else {
      wiINVAR[[colnames(Variables_ST)[y]]] <- NULL
    }
  }
  wiSHAPE[[metric]] <- wiINVAR
}
cat("\n")

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
# 5. MULTIVARIATE FITTING & STEPWISE CON FUNZIONE DI COSTRUZIONE FORMULA ROBUSTA
################################################################################
# Funzione helper interna per evitare formule malformate (stringhe vuote o + isolati)
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

var_final <- list()
taus <- TAUS_S

for (v in 1:length(var_list)) {
  perc <- (v / length(Metrics)) * 100
  cat(paste0("\rMultivariate Stepwise CV - Metric: ", v, "/", length(Metrics), " | Progresso: ", round(perc, 1), "%"))
  
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
    
    # Formula iniziale Completa
    full_formula <- build_robust_formula("VAR", invar, invar_int, c("GROUP", "SUB"))
    full <- rq(full_formula, tau = taus, method = "sfn", data = fit_data)
    mod  <- list(full)
    names(mod) <- c("start")
    meanWI     <- mean_wi(mod, taus)
    
    # Backward Stepwise Interazioni
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
    
    # Selezione Variabili Principali
    full_formula <- build_robust_formula("VAR", invar, invar_int, c("GROUP", "SUB"))
    full <- rq(full_formula, tau = taus, method = "sfn", data = fit_data)
    mod  <- list(full)
    names(mod) <- c("start")
    meanWI     <- mean_wi(mod, taus)
    
    while (rownames(meanWI)[1] != "full" && (length(invar) - length(invar_to_mantain)) > 1 && all_int == FALSE) {
      mod <- list(full)
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
    
    # Selezione dei Fattori di Raggruppamento (Intercept)
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
    
    # Modello Ridotto di riferimento
    res_formula <- if (sig == "full") as.formula("VAR ~ 1 + GROUP + SUB") else 
      if (sig == "no_group") as.formula("VAR ~ 1 + SUB") else 
        if (sig == "no_sub") as.formula("VAR ~ 1 + GROUP") else as.formula("VAR ~ 1")
    
    res        <- rq(res_formula, tau = taus, method = "sfn", data = fit_data)
    mod        <- list(full, res)
    names(mod) <- c("full", "res")
    meanWI     <- mean_wi(mod, taus)
    var_met_final[[8]] <- meanWI
    var_met_final[[9]] <- full
    
    # Validazione sui tre gruppi di quantili principali
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
    var_final[[names(Metrics)[v]]] <- var_met_final
  }
}
cat("\n")

################################################################################
# 6. GENERAZIONE REPORT DI SIGNIFICATIVITÀ
################################################################################
results <- as.data.frame(matrix(ncol = 12, nrow = length(Metrics), NA))
names(results) <- c("Metrics", "Variables", "Interaction", "Intercept_grouping", 
                    "Wifull", "Floor_Sig", "Wi_Floor", "Median_Sig", "Wi_Median", "Ceiling_Sig", "Wi_Ceiling", "Formula")

for (i in 1:length(Metrics)) {
  m_name <- names(Metrics)[i]
  results$Metrics[i] <- m_name
  
  if(!is.null(var_final[[m_name]])) {
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
    
    results$Formula[i]            <- paste(deparse(var_final[[m_name]]$Formula), collapse = "")
  }
}

write.csv(results, OUTPUT_CSV_PATH, row.names = FALSE)
print("--- REPORT SIGNIFICATIVITÀ GENERATO ---")
print(head(results))

################################################################################
# 7. GRAFICI DEI MODELLI SIGNIFICATIVI E SALVATAGGIO
################################################################################
modelli_pronti <- list()

for (i in 1:length(Metrics)) {
  m <- names(Metrics)[i]
  if (!is.null(var_final[[m]]) && (results$Floor_Sig[i] == "YES" || results$Median_Sig[i] == "YES" || results$Ceiling_Sig[i] == "YES")) {
    
    fit_data$VAR <- Metrics[, m]
    vera_formula <- as.formula(results$Formula[i])
    
    modelli_pronti[[m]] <- tryCatch({
      rq(vera_formula, tau = c(0.05, 0.50, 0.95), method = "sfn", data = fit_data)
    }, error = function(e) { NULL })
    
    if(!is.null(modelli_pronti[[m]])) {
      tryCatch({
        preds <- predict(modelli_pronti[[m]])
        plot_df <- data.frame(Observed = fit_data$VAR, Predicted_Med = preds[, 2], GROUP = fit_data$GROUP)
        
        p <- ggplot(plot_df, aes(x = Predicted_Med, y = Observed, color = GROUP)) +
          geom_point(alpha = 0.6) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
          labs(title = paste("Modello Quantilico Significativo:", m),
               x = "Valori Predetti (Mediana tau=0.50)", y = "Valori Osservati") +
          theme_minimal()
        
        ggsave(filename = paste0(OUTPUT_PLOT_DIR, m, "_significant_plot.png"), plot = p, width = 7, height = 5)
      }, error = function(e) { message(paste("Impossibile generare il grafico per la metrica:", m)) })
    }
  }
}

save(results, modelli_pronti, file = OUTPUT_RDATA)
print(paste("Processo completato senza errori. Ambiente salvato in:", OUTPUT_RDATA))