# Stratified by sex and check the association
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(dplyr)
library(pheatmap)
library(tidyr)

load("input/protine_time_points.rdata")
duplicate_all$time_point = substr(duplicate_all$Assay, 1,2)

both = data.frame(rbind(duplicate_all, no_duplicate)) # 459

# We only focus on the T1 data 
both = both %>% filter(time_point == "T1") # 378

both2 = both[,2:65]
rownames(both2) = both[,1]
sum(is.na(both2))
both2 <- both2[apply(both2, 1, function(x) !all(is.na(x))), ]

both3 = as.data.frame(lapply(both2, as.numeric))
column_means <- colMeans(both3, na.rm = TRUE)
both3_filled <- apply(both3, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
rownames(both3_filled) = rownames(both2)

dt_3 = read_excel("input/PRIME_studydata.xlsx")

anno = dt_3 %>% filter(Record_Id %in% substr(rownames(both2),3, 10))

rownames(both2) = substr(rownames(both2),3, 10)
both = both %>% filter(ID %in% rownames(both2))
pheno = both %>% select("Assay", "ID", "Plate.ID", "QC.Warning", "time_point")

pheno2 = merge(pheno, anno, by.x = "ID", by.y = "Record_Id")
pheno2 = pheno2[match(rownames(both2), pheno2$ID), ]

rownames(both3_filled) = substr(rownames(both3_filled),3, 10)
sum(pheno2$ID == rownames(both3_filled))

EAA = read.csv("input/BCG-PRIME.csv")
EAA_filter = EAA %>% filter(SID %in% rownames(both3_filled))

sum(rownames(both3_filled) == EAA_filter$SID)
EAA_filter_match = EAA_filter[match(rownames(both3_filled), EAA_filter$SID),]

sum(rownames(both3_filled) == EAA_filter_match$SID)
sum(rownames(both3_filled) == pheno2$ID)

pheno2$frailty_factor <- factor(pheno2$frailty_scale, levels = c("1 - Very fit", "2 - Well", "3 - Managing well", "4 - Vulnerable", "5 - Mildly frail", "6 - Moderately frail", "7 - Severely frail"))
pheno2$frailty_numeric <- as.numeric(pheno2$frailty_factor)

Category = c(
  "AgeAccelPheno",
  "AgeAccelGrim",
  "AgeAccelerationResidualHannum",
  "AgeAccelerationResidual",
  "DNAmAge",
  "DNAmAgeHannum",
  "DNAmPhenoAge",
  "DNAmGrimAgeBasedOnPredictedAge",
  "Age"
)

sex_group = "Male" # Male or Females

all_results <- list()

for (j in Category) {
  
  result <- data.frame()
  
  for (i in 1:64){
    
    data = data.frame(protein = both3_filled[,i], EAA = EAA_filter_match[[j]], sex = pheno2$sex, bmi = pheno2$bmi_calc, frailty = pheno2$frailty_numeric) %>% filter(sex == sex_group)
    data$nor <- qnorm((rank(data$protein, na.last="keep") - 0.5) / sum(!is.na(data$protein)))  
    mod <- lm(formula = nor ~ EAA + bmi, data = data)
    cf <- summary(mod)$coefficients
    ci <- confint(mod, level = 0.95)
    res <- data.frame(
      protein = colnames(both3_filled)[i],
      Estimate = cf[2, "Estimate"],
      P_value = cf[2, "Pr(>|t|)"],
      CI_lower = ci[2, "2.5 %"],
      CI_upper = ci[2, "97.5 %"])
    result <- rbind(result, res)
  }
  
  rownames(result) <- colnames(both3_filled)[1:64]
  colnames(result) <- c("protein","estimate", "p", "CI_lower", "CI_upper")
  result$padj <- p.adjust(result$p, method = "BH")
  result$protein = rownames(result)
  result$Category = j
  
  all_results[[j]] <- result
  
}

final_result <- do.call(rbind, all_results)

final_result$stars <- cut(final_result$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

all_result = final_result 

all_result <- all_result %>%
  mutate(Category = recode(Category,
                           "Age" = "Age",
                           "DNAmAge" = "Horvath",
                           "DNAmAgeHannum" = "Hannum",
                           "DNAmGrimAgeBasedOnPredictedAge" = "GrimAge",
                           "DNAmPhenoAge" = "PhenoAge",
                           "AgeAccelerationResidual" = "EAA_Horvath",
                           "AgeAccelerationResidualHannum" = "EAA_Hannum",
                           "AgeAccelGrim" = "EAA_Grim",
                           "AgeAccelPheno" = "EAA_Pheno"))

male_result = all_result


sex_group = "Female" 

all_results <- list()

for (j in Category) {
  
  result <- data.frame()
  
  for (i in 1:64){
    
    data = data.frame(protein = both3_filled[,i], EAA = EAA_filter_match[[j]], sex = pheno2$sex, bmi = pheno2$bmi_calc, frailty = pheno2$frailty_numeric) %>% filter(sex == sex_group)
    data$nor <- qnorm((rank(data$protein, na.last="keep") - 0.5) / sum(!is.na(data$protein)))  
    mod <- lm(formula = nor ~ EAA + bmi, data = data)
    cf <- summary(mod)$coefficients
    ci <- confint(mod, level = 0.95)
    res <- data.frame(
      protein = colnames(both3_filled)[i],
      Estimate = cf[2, "Estimate"],
      P_value = cf[2, "Pr(>|t|)"],
      CI_lower = ci[2, "2.5 %"],
      CI_upper = ci[2, "97.5 %"])
    result <- rbind(result, res)
  }
  
  rownames(result) <- colnames(both3_filled)[1:64]
  colnames(result) <- c("protein","estimate", "p", "CI_lower", "CI_upper")
  result$padj <- p.adjust(result$p, method = "BH")
  result$protein = rownames(result)
  result$Category = j
  
  all_results[[j]] <- result
  
}

final_result <- do.call(rbind, all_results)

final_result$stars <- cut(final_result$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

all_result = final_result 

all_result <- all_result %>%
  mutate(Category = recode(Category,
                           "Age" = "Age",
                           "DNAmAge" = "Horvath",
                           "DNAmAgeHannum" = "Hannum",
                           "DNAmGrimAgeBasedOnPredictedAge" = "GrimAge",
                           "DNAmPhenoAge" = "PhenoAge",
                           "AgeAccelerationResidual" = "EAA_Horvath",
                           "AgeAccelerationResidualHannum" = "EAA_Hannum",
                           "AgeAccelGrim" = "EAA_Grim",
                           "AgeAccelPheno" = "EAA_Pheno"))

female_result = all_result


all_results <- list()
for (j in Category) {
  
  result <- data.frame()
  
  for (i in 1:64){
    
    data = data.frame(protein = both3_filled[,i], EAA = EAA_filter_match[[j]], sex = pheno2$sex, bmi = pheno2$bmi_calc, frailty = pheno2$frailty_numeric)
    data$nor <- qnorm((rank(data$protein, na.last="keep") - 0.5) / sum(!is.na(data$protein)))  
    mod <- lm(formula = nor ~ EAA + bmi + sex, data = data)
    cf <- summary(mod)$coefficients
    ci <- confint(mod, level = 0.95)
    res <- data.frame(
      protein = colnames(both3_filled)[i],
      Estimate = cf[2, "Estimate"],
      P_value = cf[2, "Pr(>|t|)"],
      CI_lower = ci[2, "2.5 %"],
      CI_upper = ci[2, "97.5 %"])
    result <- rbind(result, res)
  }
  
  rownames(result) <- colnames(both3_filled)[1:64]
  colnames(result) <- c("protein","estimate", "p", "CI_lower", "CI_upper")
  result$padj <- p.adjust(result$p, method = "BH")
  result$protein = rownames(result)
  result$Category = j
  
  all_results[[j]] <- result
  
}

final_result <- do.call(rbind, all_results)

final_result$stars <- cut(final_result$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

all_result = final_result 

all_result <- all_result %>%
  mutate(Category = recode(Category,
                           "Age" = "Age",
                           "DNAmAge" = "Horvath",
                           "DNAmAgeHannum" = "Hannum",
                           "DNAmGrimAgeBasedOnPredictedAge" = "GrimAge",
                           "DNAmPhenoAge" = "PhenoAge",
                           "AgeAccelerationResidual" = "EAA_Horvath",
                           "AgeAccelerationResidualHannum" = "EAA_Hannum",
                           "AgeAccelGrim" = "EAA_Grim",
                           "AgeAccelPheno" = "EAA_Pheno"))

male_result$Group <- "Male"
female_result$Group <- "Female" 
all_result$Group <- "All"

combined_data <- rbind(all_result, male_result, female_result)

significant_proteins <- all_result %>%
  filter(padj < 0.05) %>%
  pull(protein) %>%
  unique()

cat("Found", length(significant_proteins), "significant proteins in All group:\n")
print(significant_proteins)

output_dir <- "/output/03_protein/forest_plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#significant_proteins = c("CXCL9", "CXCL10", "CCL11", "CCL4", "IL18")
if (length(significant_proteins) > 0) {
  significant_data <- combined_data %>% 
    filter(protein %in% significant_proteins & Category %in% c("EAA_Pheno", "PhenoAge")) %>%
    mutate(
      Group = factor(Group, levels = c("All", "Male", "Female")),
      protein = factor(protein)
    )
  
  p_combined <- ggplot(significant_data, aes(x = estimate, y = protein, color = Group, shape = Group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_point(position = position_dodge(width = 0.7), size = 2.5) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), 
                   height = 0.15, position = position_dodge(width = 0.7), size = 0.6) +
    scale_color_manual(values = c("All" = "#2E86AB", "Male" = "#A23B72", "Female" = "#F18F01")) +
    scale_shape_manual(values = c("All" = 16, "Male" = 17, "Female" = 15)) +
    facet_wrap(~ Category, scales = "free_y", ncol = 2) +
    labs(
      title = " ",
      x = "Effect Size (Beta Coefficient)",
      y = "Protein",
      color = "Group",
      shape = "Group"
    ) +
    theme_classic() +
    theme(
      #panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(face = "bold")
    )
  
  n_proteins <- length(significant_proteins)
  plot_height <- 10
  plot_width <- 10
  
  ggsave(file.path(output_dir, "combined_forest_plot_all_significantpheno.png"), 
         p_combined, width = plot_width, height = plot_height, dpi = 300, bg = "white")
  
  cat("Created combined forest plot with", n_proteins, "significant proteins\n")
}

if (length(significant_proteins) > 0) {
  significant_data <- combined_data %>% 
    filter(protein %in% significant_proteins & Category %in% c("EAA_Grim", "GrimAge")) %>%
    mutate(
      Group = factor(Group, levels = c("All", "Male", "Female")),
      protein = factor(protein)
    )
  
  p_combined <- ggplot(significant_data, aes(x = estimate, y = protein, color = Group, shape = Group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_point(position = position_dodge(width = 0.7), size = 2.5) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), 
                   height = 0.15, position = position_dodge(width = 0.7), size = 0.6) +
    scale_color_manual(values = c("All" = "#2E86AB", "Male" = "#A23B72", "Female" = "#F18F01")) +
    scale_shape_manual(values = c("All" = 16, "Male" = 17, "Female" = 15)) +
    facet_wrap(~ Category, scales = "free_y", ncol = 2) +
    labs(
      title = " ",
      x = "Effect Size (Beta Coefficient)",
      y = "Protein",
      color = "Group",
      shape = "Group"
    ) +
    theme_classic() +
    theme(
      #panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(face = "bold")
    )
  
  n_proteins <- length(significant_proteins)
  plot_height <- 10
  plot_width <- 10
  
  ggsave(file.path(output_dir, "combined_forest_plot_all_significant_grim.png"), 
         p_combined, width = plot_width, height = plot_height, dpi = 300, bg = "white")
  
  cat("Created combined forest plot with", n_proteins, "significant proteins\n")
}

