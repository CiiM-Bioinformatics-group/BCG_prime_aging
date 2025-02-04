rm(list = ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

EAA = read.csv("input/BCG-PRIME.csv")

dt_3 = read_excel("input/PRIME_studydata.xlsx")
pheno_filter = dt_3 %>% filter(Record_Id %in% EAA$SID)
sum(pheno_filter$Record_Id == EAA$SID) # 384

# diabetes mellitus, Chronic Obstructive Pulmonary Disease,Chronic kidney disease (CKD), Malignancy, dementia
comorbidy = c("comorb_hypertension", "comorb_cvd", "comorb_stroke", "comorb_dm", "comorb_copd", "comorb_asthma", 
              "comorb_pulm_real", "comorb_ckd", "comorb_malign", "comorb_demen")

pheno = pheno_filter[, c("Record_Id", comorbidy)]
pheno$yes_count <- rowSums(pheno[, 2:ncol(pheno)] == "Yes")

sum(pheno$Record_Id == EAA$SID) # 384


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

shapiro.test(pheno$yes_count)

result <- data.frame()

for (j in Category) {
  
  data = data.frame(comorbidy = pheno$yes_count, sex = pheno_filter$sex, EAA = EAA[[j]]) %>% na.omit()
  data$sex = as.factor(data$sex)
  #data$comorbidy <- qnorm((rank(data$comorbidy, na.last="keep") - 0.5) / sum(!is.na(data$comorbidy)))  
  
  mod <- lm(formula = comorbidy ~ EAA + sex, data = data)
  cf <- summary(mod)$coefficients
  res <- cf[2, c("Estimate", "Pr(>|t|)")]
  result <- rbind(result, res)
  
}

rownames(result) <- Category
colnames(result) <- c("estimate", "p")
result$padj <- p.adjust(result$p, method = "BH")
result$value = -log10(result$p) * result$estimate

result$stars <- cut(result$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))
final_result = result

all_result = final_result %>% filter(!Category %in% c("IEAA", "EEAA"))

all_result$estimate <- as.numeric(all_result$estimate)
all_result$sig = ifelse(all_result$padj < 0.05, "sig", "no")

all_result$Category = rownames(all_result)
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

category_order <- c("Age", "Horvath", "Hannum", "GrimAge", "PhenoAge", "EAA_Horvath", "EAA_Hannum", "EAA_Grim", "EAA_Pheno")
all_result$Category <- factor(all_result$Category, levels = category_order)

p1 = ggplot(all_result, aes(x = estimate, y = Category, color = sig)) +
  geom_bar(stat = "identity", fill = "transparent") +
  #facet_wrap(~ Category, scales = "free_y", nrow = 1) +  # Facet by Category
  labs(x = "Regression Coefficient", y = "Multimorbidity (number)") +  # Labels
  #geom_text(aes(label = stars, x = estimate + 0.05 * diff(range(estimate, na.rm = TRUE)), y = protein), size = 6, vjust = 0.7, hjust = 0) +  # Add text labels for significance at the end of each bar
  
  geom_text(aes(label = stars, x = estimate + 0.03 * diff(range(estimate, na.rm = TRUE)), y = Category), size = 6, vjust = 0.7, hjust = 0) +  # Add text labels for significance at the end of each bar
  theme_classic() +
  scale_color_manual(values = c("black", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + # Rotate x-axis labels
  xlim(0,0.05)


### heatmap ###
rm(list = ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

EAA = read.csv("input/BCG-PRIME.csv")

dt_3 = read_excel("input/PRIME_studydata.xlsx")
pheno_filter = dt_3 %>% filter(Record_Id %in% EAA$SID)
sum(pheno_filter$Record_Id == EAA$SID) # 384

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

comorbidy = c("comorb_hypertension", "comorb_cvd", "comorb_stroke", "comorb_dm", "comorb_copd", "comorb_asthma", 
              "comorb_pulm_real", "comorb_ckd", "comorb_malign")

all_results <- list()

# with smoking_household and bmi as covariable
for (j in Category){
  
  result = data.frame()
  for(i in comorbidy){
    data = data.frame(Comorbidy = pheno_filter[[i]], sex = pheno_filter$sex, EAA = EAA[[j]]) %>% na.omit()
    data$Comorbidy = as.factor(data$Comorbidy)
    data$sex = as.factor(data$sex)
    model <- glm(Comorbidy ~ EAA + sex, data = data, family = "binomial")
    p_value <- summary(model)$coefficients["EAA", "Pr(>|z|)"]
    coefficient <- coef(model)["EAA"]
    res = data.frame(Comorbidity = i, pvalue = p_value, coefficient = coefficient)
    result <- rbind(result, res)
  }
  
  result$padj <- p.adjust(result$pvalue, method = "BH")
  result$Category = j
  all_results[[j]] <- result
}

final_result <- do.call(rbind, all_results)

final_result$stars <- cut(final_result$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

final_result <- final_result %>%
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

final_result <- final_result %>%
  mutate(Comorbidity = recode(Comorbidity,
                              "comorb_hypertension" = "hypertension",
                              "comorb_cvd" = "cvd",
                              "comorb_stroke" = "stroke",
                              "comorb_dm" = "diabetes mellitus",
                              "comorb_copd" = "copd",
                              "comorb_asthma" = "asthma",
                              "comorb_pulm_real" = "pulm_real",
                              "comorb_ckd" = "ckd",
                              "comorb_malign" = "malign"))

category_order <- c("Age", "Horvath", "Hannum", "GrimAge", "PhenoAge", "EAA_Horvath", "EAA_Hannum", "EAA_Grim", "EAA_Pheno")
final_result$Category <- factor(final_result$Category, levels = category_order)

final_result$value = -log10(final_result$pvalue) * final_result$coefficient

p2 = ggplot(final_result, aes(x = Comorbidity, y = Category, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = muted("blue"), mid = "white",
                       high = muted("red"), midpoint = 0, name = "-log10(p)*estimate")+
  geom_text(aes(label=stars), color="black", size=5)+
  theme_linedraw() +
  #ylab(label = "Covariates") +
  #xlab(label = "PCs") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Comorbidity") + ylab(" ")

ppi = 300
png("output/ms_figure2c_cormorbidity_age_DNAmage_EAA_protein.png", width = 3.2 * ppi, height = 3.8 * ppi, res = ppi)
p1         
dev.off()

ppi = 300
png("output/ms_figure2d_cormorbidity_age_DNAmage_EAA_protein.png", width = 5.8 * ppi, height = 4 * ppi, res = ppi)
p2         
dev.off()

### figure 2e ###
for (j in Category){
  
  result = data.frame()
  for(i in comorbidy){
    data = data.frame(Comorbidy = pheno_filter[[i]], sex = pheno_filter$sex, EAA = EAA[[j]], smoking = pheno_filter$base_smoking_household) %>% na.omit()
    data$Comorbidy = as.factor(data$Comorbidy)
    data$sex = as.factor(data$sex)
    data$smoking = as.factor(data$smoking)
    model <- glm(Comorbidy ~ EAA + sex + smoking, data = data, family = "binomial")
    p_value <- summary(model)$coefficients["EAA", "Pr(>|z|)"]
    coefficient <- coef(model)["EAA"]
    res = data.frame(Comorbidity = i, pvalue = p_value, coefficient = coefficient)
    result <- rbind(result, res)
  }
  
  result$padj <- p.adjust(result$pvalue, method = "BH")
  result$Category = j
  all_results[[j]] <- result
}

final_result <- do.call(rbind, all_results)

final_result$stars <- cut(final_result$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

final_result <- final_result %>%
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

final_result <- final_result %>%
  mutate(Comorbidity = recode(Comorbidity,
                              "comorb_hypertension" = "hypertension",
                              "comorb_cvd" = "cvd",
                              "comorb_stroke" = "stroke",
                              "comorb_dm" = "diabetes mellitus",
                              "comorb_copd" = "copd",
                              "comorb_asthma" = "asthma",
                              "comorb_pulm_real" = "pulm_real",
                              "comorb_ckd" = "ckd",
                              "comorb_malign" = "malign"))

category_order <- c("Age", "Horvath", "Hannum", "GrimAge", "PhenoAge", "EAA_Horvath", "EAA_Hannum", "EAA_Grim", "EAA_Pheno")
final_result$Category <- factor(final_result$Category, levels = category_order)

final_result$value = -log10(final_result$pvalue) * final_result$coefficient

p3 = ggplot(final_result, aes(x = Comorbidity, y = Category, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = muted("blue"), mid = "white",
                       high = muted("red"), midpoint = 0, name = "-log10(p)*estimate")+
  geom_text(aes(label=stars), color="black", size=5)+
  theme_linedraw() +
  #ylab(label = "Covariates") +
  #xlab(label = "PCs") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Comorbidity") + ylab(" ")

ppi = 300
png("output/ms_figure2e_cormorbidity_age_DNAmage_EAA_protein_smoking.png", width = 5.8 * ppi, height = 4 * ppi, res = ppi)
p3         
dev.off()



