############################################################
# COMPARE ORIGINAL vs REFACTOR vs IMPROVED
############################################################

rm(list=ls())

library(ggplot2)


############################################################
# 1. LOAD RESULTS
############################################################

load("results/original.RData")
res_orig <- results

load("results/refactor.RData")
res_ref <- results

load("results/improved.RData")
res_imp <- results


############################################################
# 2. BASIC COMPARISON TABLE
############################################################

comparison <- data.frame(
  Metric = res_orig$Metrics,
  
  Orig_Variable = res_orig$Variables,
  Ref_Variable  = res_ref$Variables,
  Imp_Variable  = res_imp$Variable,
  
  stringsAsFactors = FALSE
)

############################################################
# 3. VARIABLE STABILITY
############################################################

comparison$Stable_Refactor <- comparison$Orig_Variable == comparison$Ref_Variable
comparison$Stable_Improved <- comparison$Orig_Variable == comparison$Imp_Variable


############################################################
# 4. CV COMPARISON (ONLY IMPROVED)
############################################################

if("CV_error" %in% names(res_imp)) {
  
  cv_table <- data.frame(
    Metric = res_imp$Metrics,
    CV_error = res_imp$CV_error
  )
  
} else {
  cv_table <- NULL
}


############################################################
# 5. SUMMARY STATISTICS
############################################################

summary_table <- data.frame(
  Version = c("Original","Refactor","Improved"),
  
  Stable_vs_Original = c(
    NA,
    mean(comparison$Stable_Refactor, na.rm=TRUE),
    mean(comparison$Stable_Improved, na.rm=TRUE)
  ),
  
  Mean_CV_Error = c(
    NA,
    NA,
    if(!is.null(cv_table)) mean(cv_table$CV_error, na.rm=TRUE) else NA
  )
)


############################################################
# 6. VISUALISATION
############################################################

# --- Stability plot
stability_plot <- data.frame(
  Version = c("Refactor","Improved"),
  Stability = c(
    mean(comparison$Stable_Refactor, na.rm=TRUE),
    mean(comparison$Stable_Improved, na.rm=TRUE)
  )
)

p1 <- ggplot(stability_plot, aes(x=Version, y=Stability)) +
  geom_bar(stat="identity") +
  ylim(0,1) +
  ggtitle("Variable Selection Stability vs Original")


# --- CV plot
if(!is.null(cv_table)) {
  
  p2 <- ggplot(cv_table, aes(x=Metric, y=CV_error)) +
    geom_point() +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    ggtitle("Cross-Validation Error (Improved only)")
}


############################################################
# 7. OUTPUT PRINT
############################################################

print("===== VARIABLE COMPARISON =====")
print(head(comparison, 10))

print("===== SUMMARY =====")
print(summary_table)

if(!is.null(cv_table)) {
  print("===== CV TABLE =====")
  print(head(cv_table, 10))
}

print(p1)

if(exists("p2")) print(p2)


############################################################
# 8. SAVE OUTPUTS
############################################################

write.csv(comparison, "results/comparison_table.csv", row.names=FALSE)
write.csv(summary_table, "results/summary_table.csv", row.names=FALSE)

if(!is.null(cv_table)) {
  write.csv(cv_table, "results/cv_table.csv", row.names=FALSE)
}

ggsave("results/stability_plot.png", p1, width=7, height=4)

if(exists("p2")) {
  ggsave("results/cv_plot.png", p2, width=10, height=5)
}
