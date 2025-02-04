rm(list = ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(reshape2)
library(tibble)
library(tidyr)
library(stringr)

cytokine = read.csv("input/500fg_pheno_91cytokines_4Raul.csv")

split_data <- str_split(cytokine$X, "_", simplify = TRUE)

cytokine <- cbind(cytokine, as.data.frame(split_data))
colnames(cytokine)[491:494] <- c("cytokine", "Stimulus", "Cell_Type", "Time")  # 重命名新列

cytokine = cytokine %>% filter(Cell_Type == "PBMC")

# Step 1: Check the distribution of the cytokine data 
cytokine_long <- gather(cytokine[,1:490], key = "sampleid", value = "value", -X)
colnames(cytokine_long) = c("cytokine", "sampleid", "value")

shapiro_test <- function(data) {
  shapiro.test(data)$p.value
}

results <- cytokine_long %>%
  group_by(cytokine) %>%
  summarize(p_value = shapiro_test(value))

non_normal_cytokines <- results %>%
  filter(p_value < 0.05)

print(results)
print(non_normal_cytokines)

# Step 2: check the missing values 
cytokine = cytokine[,1:490]
rownames(cytokine) = cytokine[,1]
cytokine.t = t(cytokine)
cytokine.t = cytokine.t[-1,]
cytokine.t = data.frame(cytokine.t)

missing_values <- cytokine.t %>%
  summarise(across(everything(), ~mean(is.na(.)), .names = "missing_{col}"))

cytokine = cytokine.t


# Step 3: Check same value
check_same_value <- function(x) {
  value_counts <- table(x)
  max(value_counts) / length(x) > 0.5
}

same_value_columns <- cytokine %>%
  summarise(across(everything(), check_same_value, .names = "same_value_{col}"))

selected_columns <- same_value_columns %>%
  select(starts_with("same_value_")) %>%
  pivot_longer(cols = everything(), names_to = "cytokine", values_to = "same_value") %>%
  filter(!same_value) %>%
  mutate(cytokine = sub("same_value_", "", cytokine)) %>%
  pull(cytokine)
selected_columns
setdiff(colnames(cytokine), selected_columns)


set.seed(123)
rank_normalize <- function(value) {
  qnorm((rank(value, na.last = "keep", ties.method = "random") - 0.5) / sum(!is.na(value)))
}

results <- cytokine %>%
  mutate(across(everything(), rank_normalize)) %>%
  summarise(across(everything(), ~shapiro.test(.x)$p.value, .names = "p_value_{col}"))

non_normal_cytokines <- results %>%
  pivot_longer(cols = everything(), names_to = "cytokine", values_to = "p_value") %>%
  filter(p_value < 0.05)


cytokine = cytokine %>%
  mutate(across(everything(), rank_normalize)) %>% 
  data.frame()

EAA = read.csv("input/500FG.csv")

idx = intersect(rownames(cytokine), EAA$SID)

cytokine_filter = cytokine[which(rownames(cytokine) %in% idx),]
EAA_filter = EAA %>% filter(SID %in% idx)

sum(rownames(cytokine_filter) == EAA_filter$SID)

pheno = read.csv("input/500FG_basicPhenos.csv")
pheno_filter = pheno %>% filter(ID_500fg %in% idx)

sum(rownames(cytokine_filter) == pheno_filter$ID_500fg)

anno_match = pheno_filter[match(rownames(cytokine_filter), pheno_filter$ID_500fg),]

sum(anno_match$ID_500fg == EAA_filter$SID)
sum(anno_match$ID_500fg == rownames(cytokine_filter))
sum(anno_match$Gender == EAA_filter$Sex)

anno_match$Gender = as.factor(anno_match$Gender)

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

all_results <- list()

for (j in Category) {
  
  result <- data.frame()
  
  for (i in 1:67){
    
    data = data.frame(protein = cytokine_filter[,i], EAA = EAA_filter[[j]], sex = anno_match$Gender, bmi = anno_match$BMI)
    mod <- lm(formula = protein ~ EAA + sex + bmi, data = data)
    cf <- summary(mod)$coefficients
    res <- cf[2, c("Estimate", "Pr(>|t|)")]
    result <- rbind(result, res)
  }
  
  rownames(result) <- colnames(cytokine_filter)[1:67]
  colnames(result) <- c("estimate", "p")
  result$padj <- p.adjust(result$p, method = "BH")
  #table(result$p < 0.05)
  #table(result$padj < 0.05)
  #table(result$estimate < 0)
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
all_result$protein <- sub("^(([^_]*_[^_]*))_.*$", "\\1", all_result$protein)

ppi <- 300
png("output/figure5_500fg.png", width = 15 * ppi, height = 4 * ppi, res = ppi)
ggplot(all_result, aes(x = protein, y = Category, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = muted("blue"), mid = "white",
                       high = muted("red"), midpoint = 0, name = "-log10(p)*estimate")+
  geom_text(aes(label=stars), color="black", size=5)+
  theme_linedraw() +
  #ylab(label = "Covariates") +
  #xlab(label = "PCs") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Cytokine")+
  ggtitle("500FG") 
dev.off()






