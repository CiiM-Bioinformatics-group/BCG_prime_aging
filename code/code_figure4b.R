# correlation analysis between inflammatory proteins and clinical phenotype.
rm(list = ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

# we load in the clinical phenotype data
dt_3 = read_excel("input/PRIME_studydata.xlsx")

# we load in the protein data
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

rownames(both3_filled) = substr(rownames(both3_filled), 3, 10)
pheno_filter = dt_3 %>% filter(Record_Id %in% rownames(both3_filled))

sum(pheno_filter$Record_Id == rownames(both3_filled)) 

pheno_filter_match = pheno_filter[match(rownames(both3_filled), pheno_filter$Record_Id),]

sum(pheno_filter_match$Record_Id == rownames(both3_filled)) # 368 samples 

# the number of cormorbidity
result <- data.frame()
for (j in 1:64){
  data = data.frame(comorbidy = pheno_filter_match$comorb_num, sex = pheno_filter_match$sex, protein = both3_filled[,j]) %>% na.omit()
  data$sex = as.factor(data$sex)
  model <- lm(comorbidy ~ protein + sex, data = data)
  p_value <- summary(model)$coefficients["protein", "Pr(>|t|)"]
  coefficient <- coef(model)["protein"]
  res = data.frame(pvalue = p_value, coefficient = coefficient, protein = colnames(both3_filled)[j])
  result <- rbind(result, res)
}
result$padj <- p.adjust(result$pvalue, method = "BH")

result$stars <- cut(result$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " "))

result$sig = ifelse(result$padj < 0.05, "sig", "no")

result <- result %>%
  arrange(desc(coefficient))

ppi = 300
png("output/ms_figure4_barplot_protein_number_of_comorbidity.png", width = 5*ppi, height = 10*ppi, res = ppi) 
ggplot(result, aes(x = coefficient, y = reorder(protein, coefficient), color = sig)) +
  geom_bar(stat = "identity", fill = "transparent") +
  labs(x = "Regression Coefficient", y = "Protein") +  # Labels
  geom_text(aes(label = stars, x = coefficient + 0.05 * diff(range(coefficient, na.rm = TRUE)), y = protein), size = 6, vjust = 0.7, hjust = 0) +  # Add text labels for significance at the end of each bar
  
  #geom_text(aes(label = stars, x = estimate, y = protein), size = 6, vjust = 0.5, hjust = 1) +  # Add text labels for significance
  theme_classic() +
  scale_color_manual(values = c("black", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")  # Rotate x-axis labels
dev.off()





