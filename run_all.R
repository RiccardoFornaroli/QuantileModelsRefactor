
set.seed(1234)

source("Original/Script_2026.R")
save(results, file="results/original.RData")

source("Refactor/REFACTOR.R")
save(results, file="results/refactor.RData")

source("Improved/IMPROVED.R")
save(results, file="results/improved.RData")

