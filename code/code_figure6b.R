rm(list = ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(reshape2)
library(tibble)
library(tidyr)

cytokine = read.csv("input/BCGprime_cytokine_nor.csv", row.names = 1)

EAA = read.csv("input/BCG-PRIME.csv")

idx = intersect(rownames(cytokine), EAA$SID)

cytokine_filter = cytokine[which(rownames(cytokine) %in% idx),]
EAA_filter = EAA %>% filter(SID %in% idx)

sum(rownames(cytokine_filter) == EAA_filter$SID)

EAA_filter_match = EAA_filter[match(rownames(cytokine_filter), EAA_filter$SID),]

sum(rownames(cytokine_filter) == EAA_filter_match$SID)

cytokine_filter <- cytokine_filter[, !grepl("rpmi", colnames(cytokine_filter))]

colnames(cytokine_filter)

dt_3 = read_excel("input/PRIME_studydata.xlsx")

anno = dt_3 %>% filter(Record_Id %in% rownames(cytokine_filter))

anno_match = anno[match(rownames(cytokine_filter), anno$Record_Id),]

sum(anno_match$Record_Id == EAA_filter_match$SID)
sum(tolower(anno_match$sex) == tolower(EAA_filter_match$Sex))
sum(rownames(cytokine_filter) == EAA_filter_match$SID)

Category = c(
  "AgeAccelPheno",
  "AgeAccelGrim",
  "AgeAccelerationResidualHannum",
  "AgeAccelerationResidual",
  "DNAmAge",
  "DNAmAgeHannum",
  "DNAmPhenoAge",
  "DNAmGrimAgeBasedOnPredictedAge",
  "Age")


all_results <- list()

for (j in Category) {
  
  result <- data.frame()
  
  for (i in 1:37){
    
    data = data.frame(protein = cytokine_filter[,i], EAA = EAA_filter_match[[j]], sex = EAA_filter_match$Sex, bmi = anno_match$bmi_calc)
    mod <- lm(formula = protein ~ EAA + sex + bmi, data = data)
    cf <- summary(mod)$coefficients
    res <- cf[2, c("Estimate", "Pr(>|t|)")]
    result <- rbind(result, res)
  }
  
  rownames(result) <- colnames(cytokine_filter)[1:37]
  colnames(result) <- c("estimate", "p")
  result$padj <- p.adjust(result$p, method = "BH")
  result$value = -log10(result$p) * result$estimate
  result$protein = rownames(result)
  result$Category = j
  
  all_results[[j]] <- result
  
}


final_result <- do.call(rbind, all_results)

final_result$stars <- cut(final_result$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

all_result <- final_result
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

all_result$sig = ifelse(all_result$padj < 0.05, "sig", "no")

stimuli = c("r848", "rpmi", "lps", "c_albicans", "delta", "omicron", "s_aureus", "wuhan", "influenza")
cyto = c("il_1ra", "il_6", "il_10", "il_1ra", "il_1b", "if_na", "tn_fa", "if_ng")

all_result <- all_result %>%
  rowwise() %>%
  mutate(
    stimuli = stringr::str_extract(protein, paste(stimuli, collapse = "|")),
    cyto = stringr::str_extract(protein, paste(cyto, collapse = "|"))
  ) %>%
  ungroup()

stimuli = c("R848", "rpmi", "LPS", "C.albicans", "Delta", "Omicron", "S.aureus", "Wuhan", "Influenza")
cyto = c("IL1RA", "IL6", "IL10",  "IL1b", "IFNA", "TNFA", "IFNy")

head(all_result$stimuli)
head(all_result$cyto)

stimuli_map <- c("r848" = "R848", "rpmi" = "rpmi", "lps" = "LPS", 
                 "c_albicans" = "C.albicans", "delta" = "Delta", 
                 "omicron" = "Omicron", "s_aureus" = "S.aureus", 
                 "wuhan" = "Wuhan", "influenza" = "Influenza")

cyto_map <- c("il_1ra" = "IL1RA", "il_6" = "IL6", "il_10" = "IL10", 
              "il_1b" = "IL1b", "if_na" = "IFNA", "tn_fa" = "TNFA", 
              "if_ng" = "IFNy")

all_result <- all_result %>%
  mutate(
    stimuli_id = recode(stimuli, !!!stimuli_map),
    cyto_id = recode(cyto, !!!cyto_map)
  )

all_result$cytokine = paste0(all_result$cyto_id, "_", all_result$stimuli_id)

ppi <- 300
png("output/figure5_heatmap_age_DNAmage_EAA_cytokine_no_rpmi_order_by_cytokine.png", width = 10 * ppi, height = 3 * ppi, res = ppi)
ggplot(all_result, aes(x = cytokine, y = Category, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = muted("blue"), mid = "white",
                       high = muted("red"), midpoint = 0, name = "-log10(p)*estimate")+
  geom_text(aes(label=stars), color="black", size=5)+
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("BCG-PRIME") +
  xlab("Cytokine")
dev.off()

