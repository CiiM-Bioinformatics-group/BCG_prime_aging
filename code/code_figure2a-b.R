rm(list = ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(reshape2)
library(ggpubr)

EAA = read.csv("input/BCG-PRIME.csv")

dt_3 = read_excel("input/PRIME_studydata.xlsx")
pheno_filter = dt_3 %>% filter(Record_Id %in% EAA$SID)
sum(pheno_filter$Record_Id == EAA$SID) # 384
pheno_filter$frailty_score = substr(pheno_filter$frailty_scale, 1, 2)

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

result <- NULL
for(j in Category){
    data = data.frame(comorbidy = pheno_filter$frailty_score, sex = pheno_filter$sex, EAA = EAA[[j]]) %>% na.omit()
    data$comorbidy = as.numeric(data$comorbidy)
    data$sex = as.factor(data$sex)
    model <- lm(comorbidy ~ EAA + sex, data = data)
    p_value <- summary(model)$coefficients["EAA", "Pr(>|t|)"]
    coefficient <- coef(model)["EAA"]
    res = data.frame(pvalue = p_value, coefficient = coefficient, Category = j)
    result <- rbind(result, res)
}

final_result <- result
final_result$padj = p.adjust(final_result$pvalue, method = "BH")
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


category_order <- c("Age", "Horvath", "Hannum", "GrimAge", "PhenoAge", "EAA_Horvath", "EAA_Hannum", "EAA_Grim", "EAA_Pheno")
final_result$Category <- factor(final_result$Category, levels = category_order)

final_result$sig = ifelse(final_result$padj < 0.05, "sig", "no")

sum(pheno_filter$Record_Id == EAA$SID) # 384

EAA_select = EAA[,c("SID", Category)]
EAA_select$frailty = pheno_filter$frailty_scale

EAA_long <- melt(EAA_select, 
                 id.vars = c("SID", "frailty"),  
                 value.name = "Value")  

EAA_long <- EAA_long %>%
  mutate(Category = recode(variable,
                           "Age" = "Age",
                           "DNAmAge" = "Horvath",
                           "DNAmAgeHannum" = "Hannum",
                           "DNAmGrimAgeBasedOnPredictedAge" = "GrimAge",
                           "DNAmPhenoAge" = "PhenoAge",
                           "AgeAccelerationResidual" = "EAA_Horvath",
                           "AgeAccelerationResidualHannum" = "EAA_Hannum",
                           "AgeAccelGrim" = "EAA_Grim",
                           "AgeAccelPheno" = "EAA_Pheno"))

p1 = ggplot(final_result, aes(x = coefficient, y = Category, color = sig)) +
  geom_bar(stat = "identity", fill = "transparent") +
  labs(x = "Regression Coefficient") +   
  geom_text(aes(label = stars, x = coefficient + 0.03 * diff(range(coefficient, na.rm = TRUE)), y = Category), size = 6, vjust = 0.7, hjust = 0) + 
  theme_classic() +
  scale_color_manual(values = c("black", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + 
  xlim(0,0.05)+
  theme(plot.margin = margin(b = 65, t=5))  


p2 = ggplot(data = EAA_long, aes(x=frailty, y=Value)) + geom_boxplot() + 
  facet_wrap(~Category,scales = "free_y", ncol = 5) +
  theme_classic()+
  geom_smooth(method = "lm", aes(group = 1), color = "red", se = FALSE) +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
  ylab("Age/EpiAge/EAA") +
  xlab("Frailty")

ppi <- 300
png("output/figure2ab_age_DNAmage_EAA_frailty.png", width = 9 * ppi, height = 4.5 * ppi, res = ppi)
ggarrange(p1,p2,
          widths = c(1, 2.1),
          heights = c(1, 2.1),
          ncol = 2, 
          nrow = 1)
dev.off()



