
################################################################################
# 0. CONFIGURAZIONE DELLE VARIABILI GLOBALI, SOGLIE E PARAMETRI
################################################################################
rm(list=ls())

# Sementi e riproducibilità
SET_SEED <- 1234
set.seed(SET_SEED)

# Parametri di Cross-Validazione
K_FOLDS <- 10

# Soglie di Significatività Statistica e Costanti
EPSILON_LOGIT <- 0.001
VIF_THRESHOLD <- 1.5
SIG_THRESHOLD <- 0.750
MIN_NON_ZERO_PROP <- 0.20 # Soglia per esclusione

# Definizione dei Quantili (Taus)
# TAUS_S <- seq(0.02, 0.98, length.out = 5)
# TAUS_M <- list(
#   Floor   = seq(0.02, 0.10, length.out = 5),
#   Median  = seq(0.45, 0.55, length.out = 5),
#   Ceiling = seq(0.90, 0.98, length.out = 5)
# )
TAUS_S <- seq(0.02, 0.98, length.out = 97)
TAUS_M <- list(
  Floor   = seq(0.02, 0.10, length.out = 100),
  Median  = seq(0.45, 0.55, length.out = 100),
  Ceiling = seq(0.90, 0.98, length.out = 100)
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

################################################################################
# 2.1 CONTROLLO SPARSITÀ METRICHE
################################################################################

# Calcolo della proporzione di valori > 0 per ogni colonna in Metrics
non_zero_proportions <- sapply(Metrics, function(x) {
  sum(x > 0, na.rm = TRUE) / length(na.omit(x))
})

# Identificazione delle metriche da escludere
excluded_metrics_names <- names(non_zero_proportions[non_zero_proportions < MIN_NON_ZERO_PROP])

# Filtro effettivo del set di dati per la modellazione
metrics_to_model <- Metrics[, non_zero_proportions >= MIN_NON_ZERO_PROP]

cat("Metriche escluse per eccessiva sparsità (< 10%):", paste(excluded_metrics_names, collapse = ", "), "\n")
cat("Metriche procedenti alla modellazione:", ncol(metrics_to_model), "\n")

# Aggiornamento dell'oggetto Metrics per il resto dello script
Metrics <- metrics_to_model

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

    if(length(cv_wi_list) > 0){

        all_models <- unique(unlist(lapply(cv_wi_list, rownames)))

        wi_matrix <- matrix(NA,
                            nrow = length(all_models),
                            ncol = length(cv_wi_list),
                            dimnames = list(all_models,
                                            paste0("Fold",1:length(cv_wi_list))))

        for(f in seq_along(cv_wi_list)){
          wi_matrix[rownames(cv_wi_list[[f]]), f] <- cv_wi_list[[f]]$Wi
        }

        combined_wi <- data.frame(
          MeanWi = rowMeans(wi_matrix, na.rm = TRUE),
          SDWi   = apply(wi_matrix, 1, sd, na.rm = TRUE),
          SEWi   = apply(wi_matrix, 1, sd, na.rm = TRUE) /
            sqrt(colSums(!is.na(wi_matrix))[1]),
          BestFrequency = apply(wi_matrix, 1, function(x)
            mean(x == max(x, na.rm = TRUE), na.rm = TRUE))
        )

        combined_wi$StableScore <- combined_wi$MeanWi * (1 - combined_wi$SDWi)

        combined_wi <- combined_wi[order(combined_wi$StableScore, decreasing = TRUE), ]
        wiINVAR[[colnames(Variables_ST)[y]]] <- combined_wi

      }
      else {
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

  selection <- data.frame(
    Index = 1:ncol(Variables_ST),
    Variable = colnames(Variables_ST),
    Selected = NA,
    Model = NA,
    MeanWi = NA,
    StableScore = NA
  )

  for (i in 1:ncol(Variables_ST)) {

    selection$Selected[i] <- "YES"

    if (!is.null(wiSHAPE[[v]][[i]])) {
      selection$Model[i]       <- rownames(wiSHAPE[[v]][[i]])[1]
      selection$MeanWi[i]      <- wiSHAPE[[v]][[i]]$MeanWi[1]
      selection$StableScore[i] <- wiSHAPE[[v]][[i]]$StableScore[1]
    }

  }

  selected <- selection[selection$Selected == "YES", ]

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
###############################################################################
# 5. MULTIVARIATE FITTING & STEPWISE (LOGICA INTEGRALE)
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

cv_mean_wi <- function(models, fit_data, taus, folds){
  cv_wi_list <- list()
  for(f in seq_along(folds)){
    test_idx <- folds[[f]]
    train_data <- fit_data[-test_idx, ]
    fitted_models <- list()
    for(m in seq_along(models)){
      fit <- tryCatch(rq(models[[m]], tau = taus, method = "sfn", data = train_data), error = function(e) NULL)
      if(!is.null(fit)) fitted_models[[names(models)[m]]] <- fit
    }
    if(length(fitted_models) > 1) cv_wi_list[[f]] <- mean_wi(fitted_models, taus)
  }
  if(length(cv_wi_list) == 0) return(NULL)

  all_models <- unique(unlist(lapply(cv_wi_list, rownames)))
  wi_matrix <- matrix(NA, nrow = length(all_models), ncol = length(cv_wi_list),
                      dimnames = list(all_models, paste0("Fold_", seq_along(cv_wi_list))))
  for(f in seq_along(cv_wi_list)){
    if(!is.null(cv_wi_list[[f]])) wi_matrix[rownames(cv_wi_list[[f]]), f] <- cv_wi_list[[f]]$Wi
  }
  return(list(wi_matrix = wi_matrix))
}

aggregate_cv_wi <- function(cv_list_obj){
  wi_matrix <- cv_list_obj$wi_matrix
  data.frame(
    MeanWi = rowMeans(wi_matrix, na.rm = TRUE),
    SDWi   = apply(wi_matrix, 1, sd, na.rm = TRUE),
    StableScore = rowMeans(wi_matrix, na.rm = TRUE) * (1 - apply(wi_matrix, 1, sd, na.rm = TRUE))
  )
}

multivariate_cv_stepwise <- function(VAR, var_list_row, Metrics_row, Variables_ST, GROUP, SUB, TAUS_S, TAUS_M, folds, wi_selected, debug = TRUE){

  fit_data <- as.data.frame(cbind(VAR, Variables_ST, GROUP, SUB))
  var_met_final <- list()
  var_met_final[[1]] <- names(var_list_row)

  invar <- c()
  for(i in 1:nrow(var_list_row)){
    if(var_list_row[i,3]=="LIN") invar <- c(invar, as.character(var_list_row[i,1]))
    if(var_list_row[i,3]=="LOG") invar <- c(invar, paste0("log10(", var_list_row[i,1], ")"))
    if(var_list_row[i,3]=="EXP") invar <- c(invar, paste0("exp(", var_list_row[i,1], ")"))
    if(var_list_row[i,3]=="QUA") invar <- c(invar, paste0("poly(", var_list_row[i,1], ",2)"))
  }

  var_met_final[[2]] <- wi_selected
  var_met_final[[3]] <- invar
  invar_int <- paste0("(", invar, "):GROUP")

  # --- STEP 1: INTERAZIONI ---
  repeat{
    base <- build_robust_formula("VAR", invar, invar_int, c("GROUP","SUB"))
    formulas <- list(full = base)
    for(i in seq_along(invar_int)){
      formulas[[invar_int[i]]] <- build_robust_formula("VAR", invar, invar_int[-i], c("GROUP","SUB"))
    }
    score <- aggregate_cv_wi(cv_mean_wi(formulas, fit_data, TAUS_S, folds))
    best <- rownames(score)[which.max(score$StableScore)]
    if(best=="full" || length(invar_int)<=1) break
    invar_int <- setdiff(invar_int, best)
  }

  # --- STEP 2: VARIABILI ---
  repeat{
    base <- build_robust_formula("VAR", invar, invar_int, c("GROUP","SUB"))
    formulas <- list(full = base)
    for(i in seq_along(invar)){
      formulas[[invar[i]]] <- build_robust_formula("VAR", invar[-i], invar_int, c("GROUP","SUB"))
    }
    score <- aggregate_cv_wi(cv_mean_wi(formulas, fit_data, TAUS_S, folds))
    best <- rownames(score)[which.max(score$StableScore)]
    if(best=="full" || length(invar)<=1) break
    invar <- setdiff(invar, best)
  }

  # --- STEP 3: GROUP / SUB ---
  formulas <- list(
    full = build_robust_formula("VAR", invar, invar_int, c("GROUP","SUB")),
    no_group = build_robust_formula("VAR", invar, invar_int, "SUB"),
    no_sub = build_robust_formula("VAR", invar, invar_int, "GROUP")
  )
  score <- aggregate_cv_wi(cv_mean_wi(formulas, fit_data, TAUS_S, folds))
  sig <- rownames(score)[1]

  # --- MODELLI FINALI ---
  full_formula <- build_robust_formula("VAR", invar, invar_int, c("GROUP","SUB"))
  res_formula <- if(sig=="full") as.formula("VAR ~ 1 + GROUP + SUB") else if(sig=="no_group") as.formula("VAR ~ 1 + SUB") else if(sig=="no_sub") as.formula("VAR ~ 1 + GROUP") else as.formula("VAR ~ 1")

  # Calcolo statistiche per il modello ridotto
  cv_res_stats <- cv_mean_wi(list(res = res_formula), fit_data, TAUS_S, folds)
  stats_res <- if(!is.null(cv_res_stats)) aggregate_cv_wi(cv_res_stats) else data.frame(MeanWi=NA, SDWi=NA, StableScore=NA)

  # Fit dei modelli
  res_model <- rq(res_formula, tau = TAUS_S, method = "sfn", data = fit_data)
  full_model <- rq(full_formula, tau = TAUS_S, method = "sfn", data = fit_data)

  # Aggiungiamo le statistiche come attributi del modello (non rompe la compatibilità con rq)
  attr(res_model, "stats") <- stats_res

  var_met_final[[4]] <- invar
  var_met_final[[5]] <- invar_int
  var_met_final[[6]] <- sig
  var_met_final[[7]] <- full_formula
  var_met_final[[8]] <- res_model  # Questo è l'oggetto rq
  var_met_final[[9]] <- full_model

  # --- VALIDAZIONE QUANTILI ---
  # Creiamo una lista temporanea per i risultati dei quantili
  results_quantiles <- list()

  for(i in seq_along(TAUS_M)){
    taus_i <- TAUS_M[[i]]
    cv_q <- cv_mean_wi(list(full = full_formula, res = res_formula), fit_data, taus_i, folds)

    if(!is.null(cv_q)) {
      wi_q <- aggregate_cv_wi(cv_q)
    } else {
      wi_q <- data.frame(MeanWi = NA, SDWi = NA, StableScore = NA)
    }
    results_quantiles[[i]] <- wi_q
  }

  var_met_final <- c(var_met_final, results_quantiles)
  names(var_met_final) <- c("Metric", "Selected Variables", "Invar_Inizial", "Invar_Selected", "Invar_int_Selected",
                            "Intercept_grouping", "Formula", "Significance", "Model", "Wi_floor", "Wi_median", "Wi_ceiling")

    if(debug) print(score)
  return(var_met_final)
}

# --- ESECUZIONE ---
var_final <- list()
for(v in 1:length(var_list)){
  cat(paste0("\rMultivariate CV Stepwise: ", v,"/",length(var_list)))
  if(all(!is.na(var_list[[v]]))){
    var_final[[names(Metrics)[v]]] <- multivariate_cv_stepwise(
      VAR = Metrics[[v]], var_list_row = var_list[[v]], Metrics_row = Metrics[[v]],
      Variables_ST = Variables_ST, GROUP = GROUP, SUB = SUB,
      TAUS_S = TAUS_S, TAUS_M = TAUS_M, folds = folds, wi_selected = wi_selected[[v]],
      debug=F
    )
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

    # Estraiamo le statistiche salvate nell'attributo
    stats <- attr(var_final[[m_name]]$Significance, "stats")

    results$Wifull[i]       <- stats$StableScore[1]

    # Per i quantili, assicurati che i nomi siano quelli giusti ("full" non esiste in Wi_floor)
    results$Wi_Floor[i]     <- var_final[[m_name]]$Wi_floor[1, 1]
    results$Floor_Sig[i]    <- ifelse(results$Wi_Floor[i] > SIG_THRESHOLD, "YES", "NO")

    results$Wi_Median[i]    <- var_final[[m_name]]$Wi_median[1, 1]
    results$Median_Sig[i]   <- ifelse(results$Wi_Median[i] > SIG_THRESHOLD, "YES", "NO")

    results$Wi_Ceiling[i]   <- var_final[[m_name]]$Wi_ceiling[1, 1]
    results$Ceiling_Sig[i]  <- ifelse(results$Wi_Ceiling[i] > SIG_THRESHOLD, "YES", "NO")
    results$Formula[i]            <- paste(deparse(var_final[[m_name]]$Formula), collapse = "")
  }
}

write.csv(results, OUTPUT_CSV_PATH, row.names = FALSE)
print("--- REPORT SIGNIFICATIVITÀ GENERATO ---")
print(head(results))


# --- SEZIONE 7: Visualizzazione Corretta (versione finale) ---

library(ggplot2)
library(dplyr)
library(quantreg)
library(tidyr)

if (!dir.exists("results/plots")) {
  dir.create("results/plots", recursive = TRUE)
}

plot_models_final <- function(var_final, Variables_ST, GROUP, SUB, Metrics, metric_name) {
  
  base_df <- cbind(Variables_ST, GROUP = GROUP, SUB = SUB)
  
  for (m_name in names(var_final)) {
    
    if (!m_name %in% colnames(Metrics)) next
    
    df_plot <- base_df
    df_plot$VAL <- Metrics[[m_name]]
    
    f_orig <- var_final[[m_name]]$Formula
    if (is.null(f_orig)) next
    
    f_char <- as.character(f_orig)
    f_new <- as.formula(paste("VAL", f_char[3], sep = " ~ "))
    
    message("  > Processando: ", m_name)
    
    fit_rq <- tryCatch(
      rq(f_new, data = df_plot, tau = c(0.05, 0.5, 0.95), method = "fn"),
      error = function(e) return(NULL)
    )
    
    if (is.null(fit_rq)) next
    
    # SELEZIONE VARIABILE X
    f_txt <- paste(deparse(f_orig), collapse = " ")
    numeric_vars <- names(df_plot)[sapply(df_plot, is.numeric)]
    numeric_vars <- setdiff(numeric_vars, "VAL")
    
    present_vars <- numeric_vars[sapply(numeric_vars, function(v) grepl(paste0("\\b", v, "\\b"), f_txt))]
    
    if (length(present_vars) == 0) next
    
    priority <- intersect(c("VEL", "DEPTH"), present_vars)
    clean_x <- if (length(priority) > 0) priority[1] else present_vars[1]
    
    # CONTROLLO INTERAZIONE REALE
    has_interaction <- grepl(paste0(clean_x, ".*:.*GROUP|GROUP.*:.*", clean_x), f_txt) | 
      grepl(paste0(clean_x, ".*\\*GROUP|GROUP.*\\*", clean_x), f_txt)
    
    # COSTRUZIONE NEWDATA
    pred_df <- df_plot[rep(1, 100), , drop = FALSE]
    all_vars <- setdiff(names(pred_df), "VAL")
    
    for (v in all_vars) {
      if (is.numeric(df_plot[[v]])) {
        pred_df[[v]] <- median(df_plot[[v]], na.rm = TRUE)
      } else if (is.factor(df_plot[[v]])) {
        pred_df[[v]] <- factor(levels(df_plot[[v]])[1], levels = levels(df_plot[[v]]))
      } else {
        pred_df[[v]] <- df_plot[[v]][1]
      }
    }
    
    pred_df[[clean_x]] <- seq(min(df_plot[[clean_x]], na.rm = TRUE), max(df_plot[[clean_x]], na.rm = TRUE), length.out = 100)
    
    # INTERAZIONE GROUP
    if (has_interaction) {
      pred_df <- bind_rows(
        lapply(levels(df_plot$GROUP), function(g) {
          tmp <- pred_df
          tmp$GROUP <- factor(g, levels = levels(df_plot$GROUP))
          tmp
        })
      )
    } else {
      pred_df$GROUP <- factor(levels(df_plot$GROUP)[1], levels = levels(df_plot$GROUP))
    }
    
    # PREDIZIONE
    preds <- predict(fit_rq, newdata = pred_df, type = "matrix")
    plot_data <- cbind(pred_df, preds)
    colnames(plot_data)[(ncol(pred_df)+1):ncol(plot_data)] <- c("Q05", "Q50", "Q95")
    
    # FILTRO QUANTILI SIGNIFICATIVI
    sig_quantiles <- c("Q05", "Q50", "Q95")[c(!is.na(var_final[[m_name]]$Wi_floor["full", "StableScore"]) && var_final[[m_name]]$Wi_floor["full", "StableScore"] >= SIG_THRESHOLD,
                                              !is.na(var_final[[m_name]]$Wi_median["full", "StableScore"]) && var_final[[m_name]]$Wi_median["full", "StableScore"] >= SIG_THRESHOLD,
                                              !is.na(var_final[[m_name]]$Wi_ceiling["full", "StableScore"]) && var_final[[m_name]]$Wi_ceiling["full", "StableScore"] >= SIG_THRESHOLD)]
    
    if (length(sig_quantiles) == 0) next
    
    plot_data_long <- plot_data %>% pivot_longer(cols = all_of(sig_quantiles), names_to = "Quantile", values_to = "Pred")
    
    # PLOT
    p <- ggplot(plot_data_long, aes(x = .data[[clean_x]], y = Pred, color = Quantile)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c(Q05 = "red", Q50 = "black", Q95 = "green4")[sig_quantiles],
                         labels = c(Q05 = "0.05", Q50 = "0.50", Q95 = "0.95")[sig_quantiles], drop = FALSE) +
      labs(title = paste("Modello Quantile:", m_name), x = clean_x, y = metric_name, color = "Tau") +
      theme_minimal()
    
    if (has_interaction) p <- p + facet_wrap(~GROUP)
    
    ggsave(file.path("results/plots", paste0("Plot_", m_name, ".png")), p, width = 8, height = 6)
    message("  > Salvato: ", m_name)
  }
}

plot_models_final(var_final, Variables_ST, GROUP, SUB, Metrics, "Valore_Metrica")






# SEZIONE 8: Esportazione Dati e Reportistica (Versione a prova di errore)
library(rmarkdown)
library(knitr)
library(dplyr)

# 1. Estrazione dati
write.csv(results, "results/results_summary.csv", row.names = FALSE)

# 2. Creazione file Rmd (Senza caratteri di escape complessi)
rmd_path <- "results/Report_Analisi.Rmd"

sink(rmd_path)
cat("---\n")
cat("title: 'Report Analisi Quantilica'\n")
cat("output: html_document\n")
cat("---\n\n")

cat("## Riepilogo Significatività\n")
cat("```{r echo=FALSE}\n")
cat("library(knitr)\n")
cat("df <- read.csv('results_summary.csv')\n")
cat("kable(df)\n")
cat("```\n\n")

cat("## Grafici Modelli\n")
cat("```{r echo=FALSE, results='asis'}\n")
# Usiamo ../plots per uscire da 'results' e andare in 'plots'
cat("plots <- list.files('../results/plots', pattern = '.png', full.names = TRUE)\n")
cat("if(length(plots) > 0) {\n")
cat("  for (p in plots) {\n")
cat("    nome_grafico <- gsub('Plot_', '', basename(p))\n")
cat("    nome_grafico <- gsub('.png', '', nome_grafico)\n")
cat("    cat('\\n### ', nome_grafico, '\\n')\n")
# Qui scriviamo il percorso relativo corretto per il file HTML finale
cat("    cat('![](', p, '){width=80%}\\n\\n')\n")
cat("  }\n")
cat("} else {\n")
cat("  cat('Nessun grafico trovato. Percorso cercato:', list.files('../results/plots', full.names=TRUE))\n")
cat("}\n")
cat("```\n")
sink()

# 3. Compilazione
render(rmd_path, output_dir = "results")
message("  > Report generato: results/Report_Analisi.html")
