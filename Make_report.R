############################################################
# COMPARE RESULTS - ORIGINAL vs REFACTOR vs IMPROVED vs IMPROVED2
############################################################

rm(list = ls())

library(dplyr)
library(rmarkdown)

############################################################
# 1. LOAD RESULTS
############################################################

orig <- read.csv("results/original_results.csv")
# ref  <- read.csv("results/refactor_results.csv")
imp  <- read.csv("results/improved_results.csv")

############################################################
# 2. STANDARDIZE COLUMN NAMES
############################################################

standardize <- function(df, name) {
  
  df %>%
    mutate(Model = name) %>%
    select(
      Metrics,
      Variables,
      Formula,
      Model,
      everything()
    )
}

orig <- standardize(orig, "ORIGINAL")
# ref  <- standardize(ref, "REFACTOR")
imp  <- standardize(imp, "IMPROVED")

############################################################
# 3. MERGE ALL RESULTS
############################################################

all_results <- bind_rows(orig, imp)

############################################################
# 4. METRIC-LEVEL SUMMARY
############################################################

summary_table <- all_results %>%
  group_by(Metrics, Model) %>%
  summarise(
    Variables = first(Variables),
    Formula   = first(Formula),
    .groups = "drop"
  )

############################################################
# 5. CV COMPARISON (ONLY IMPROVED HAS CV)
############################################################

if("CV_full" %in% colnames(imp)) {
  
  cv_table <- imp %>%
    mutate(
      CV_gain = CV_null - CV_full,
      CV_relative = (CV_null - CV_full) / abs(CV_null + 1e-8)
    ) %>%
    select(Metrics, CV_full, CV_null, CV_gain, CV_relative)
  
} else {
  cv_table <- NULL
}

############################################################
# 6. BEST MODEL PER METRIC (BASED ON SIMPLE RULE)
############################################################

# regola: preferisci CV se disponibile, altrimenti Wi
best_models <- all_results %>%
  mutate(
    Score = case_when(
      Model == "IMPROVED" ~ -CV_full,   # min CV = best
      TRUE ~ -1  # placeholder per altri
    )
  ) %>%
  group_by(Metrics) %>%
  slice_min(order_by = Score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Metrics, Model, Variables, Formula)

############################################################
# 7. MODEL COUNTS (robustezza confronto)
############################################################

model_counts <- all_results %>%
  group_by(Model) %>%
  summarise(
    n_metrics = n(),
    n_unique_vars = n_distinct(Variables),
    .groups = "drop"
  )

############################################################
# 8. FINAL MASTER TABLE
############################################################

final_table <- list(
  all_results = all_results,
  summary = summary_table,
  best_models = best_models,
  model_counts = model_counts,
  cv_table = cv_table
)

############################################################
# 9. SAVE OUTPUTS (LIGHT + GITHUB SAFE)
############################################################

dir.create("results", showWarnings = FALSE)

write.csv(all_results, "results/COMPARE_all_results.csv", row.names = FALSE)
write.csv(best_models, "results/COMPARE_best_models.csv", row.names = FALSE)
write.csv(model_counts, "results/COMPARE_model_counts.csv", row.names = FALSE)

if(!is.null(cv_table)) {
  write.csv(cv_table, "results/COMPARE_cv_table.csv", row.names = FALSE)
}

saveRDS(final_table, "results/COMPARE_full_object.rds")

############################################################
# 10. QUICK PRINT
############################################################

print("=== MODEL COUNTS ===")
print(model_counts)

print("=== BEST MODELS (TOP 10) ===")
print(head(best_models, 10))

library(rmarkdown)

rmarkdown::render(
  input = "report_template.Rmd",
  output_file = "results/FINAL_REPORT.html"
)
