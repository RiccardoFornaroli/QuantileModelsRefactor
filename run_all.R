set.seed(1234)

dir.create("results", showWarnings = FALSE)

source("Original/Script_2026.R")
write.csv(results, "results/original_results.csv", row.names = FALSE)

source("Refactor/REFACTOR.R")
write.csv(results, "results/refactor_results.csv", row.names = FALSE)

source("Improved/IMPROVED.R")
write.csv(results, "results/improved_results.csv", row.names = FALSE)