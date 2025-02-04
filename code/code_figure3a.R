# All four cohorts, EAA, real age, protein, heatmap
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(scales)

load("input/BCGprime_final_result_all_age_eaa.rdata")

bcg_prime_all_result = final_result
rm(final_result)

bcg_prime_all_result$protein <- gsub("\\.", "_", bcg_prime_all_result$protein)
bcg_prime_all_result$protein <- gsub("X4E_BP1", "4E_BP1", bcg_prime_all_result$protein)
bcg_prime_all_result = bcg_prime_all_result[,-1]

protein = bcg_prime_all_result %>% filter(padj < 0.05) %>% pull(protein) %>% unique()
bcg_prime_all_result_sig = bcg_prime_all_result[which(bcg_prime_all_result$protein %in% protein),]

fg500 = read_excel("input/500FG_final_result_all_age_eaa.xlsx")

length(intersect(fg500$protein, bcg_prime_all_result_sig$protein))

bcg300 = read_excel("input/300BCG_final_result_all_age_eaa.xlsx")

bcg300$protein <- gsub("X4E_BP1", "4E_BP1", bcg300$protein)

imed = read_excel("input/iMED_final_result_all_age_eaa.xlsx")
length(intersect(imed$protein, bcg_prime_all_result_sig$protein))

imed$protein <- gsub("IL17C", "IL_17C", imed$protein)
imed$protein <- gsub("TGFA", "TGF_alpha", imed$protein)
imed$protein <- gsub("MMP1", "MMP_1", imed$protein)
imed$protein <- gsub("LIFR", "LIF_R", imed$protein)
imed$protein <- gsub("IL10RB", "IL_10RB", imed$protein)
imed$protein <- gsub("IL18R1", "IL_18R1", imed$protein)
imed$protein <- gsub("IL12B", "IL_12B", imed$protein)
imed$protein <- gsub("MMP10", "MMP_10", imed$protein)
imed$protein <- gsub("FLT3LG", "Flt3L", imed$protein)
imed$protein <- gsub("IFNG", "IFN_gamma", imed$protein)
imed$protein <- gsub("CSF1", "CSF_1", imed$protein)
length(intersect(imed$protein, bcg_prime_all_result_sig$protein))

idx =  intersect(fg500$protein ,intersect(bcg300$protein, intersect(bcg_prime_all_result_sig$protein, imed$protein)))

bcg_prime_2 = bcg_prime_all_result_sig %>% filter(!Category %in% c("IEAA", "EEAA"))
imed_2 = imed %>% filter(!Category %in% c("IEAA", "EEAA"))
fg500_2 = fg500 %>% filter(!Category %in% c("IEAA", "EEAA"))
bcg300_2 = bcg300 %>% filter(!Category %in% c("IEAA", "EEAA"))

reference_proteins <- unique(bcg_prime_all_result_sig$protein)

prepare_data <- function(df, reference_proteins) {
  df %>%
    filter(protein %in% reference_proteins)
}

bcg_prime_2 <- prepare_data(bcg_prime_2, reference_proteins)
imed_2 <- prepare_data(imed_2, reference_proteins)
fg500_2 <- prepare_data(fg500_2, reference_proteins)
bcg300_2 <- prepare_data(bcg300_2, reference_proteins)


# Combine all protein lists and count occurrences
all_proteins <- data.frame(protein = unique(c(reference_proteins, fg500_2$protein, bcg300_2$protein, imed_2$protein)))
all_proteins <- all_proteins %>%
  mutate(fg500 = protein %in% fg500_2$protein,
         bcg300 = protein %in% bcg300_2$protein,
         bcg_prime_2 = protein %in% bcg_prime_2$protein,
         imed = protein %in% imed_2$protein,
         count = fg500 + bcg300 + bcg_prime_2 + imed) %>%
  arrange(desc(count), protein)

# Step 2: Ensure consistent protein order
ordered_proteins <- all_proteins$protein

prepare_data <- function(df, ordered_proteins) {
  df %>%
    right_join(data.frame(protein = ordered_proteins), by = "protein") %>%
    arrange(match(protein, ordered_proteins)) %>%
    mutate(protein = factor(protein, levels = ordered_proteins))
}

bcg_prime_2 <- prepare_data(bcg_prime_2, ordered_proteins)
imed_2 <- prepare_data(imed_2, ordered_proteins)
fg500_2 <- prepare_data(fg500_2, ordered_proteins)
bcg300_2 <- prepare_data(bcg300_2, ordered_proteins)


add_stars <- function(df) {
  df %>%
    mutate(padj = p.adjust(p, method = "BH"),
           stars = cut(padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", " ")))
}

#bcg_prime_2 <- add_stars(bcg_prime_2)
imed_2 <- add_stars(imed_2)
fg500_2 <- add_stars(fg500_2)
bcg300_2 <- add_stars(bcg300_2)

recode_category <- function(df) {
  df %>%
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
}

bcg_prime_2 <- recode_category(bcg_prime_2)
imed_2 <- recode_category(imed_2)
fg500_2 <- recode_category(fg500_2)
bcg300_2 <- recode_category(bcg300_2)

category_order <- c("Age", "Horvath", "Hannum", "GrimAge", "PhenoAge", "EAA_Horvath", "EAA_Hannum", "EAA_Grim", "EAA_Pheno")
bcg_prime_2$Category <- factor(bcg_prime_2$Category, levels = category_order)
imed_2$Category <- factor(imed_2$Category, levels = category_order)
fg500_2$Category <- factor(fg500_2$Category, levels = category_order)
bcg300_2$Category <- factor(bcg300_2$Category, levels = category_order)

plot_heatmap <- function(df, title) {
  ggplot(df, aes(x = protein, y = Category, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = muted("blue"), mid = "white",
                         high = muted("red"), midpoint = 0, name = "-log10(p)*sign(estimate)"
                         #limits = c(-7, 15),
                         #breaks = seq(-7, 15, by = 2)
                         ) +
    geom_text(aes(label = stars), color = "black", size = 5) +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title)
}

p1 <- plot_heatmap(bcg_prime_2, "Old (BCG-PRIME:n=368)")
p2 <- plot_heatmap(imed_2, "Old (iMed:n=165)")
p3 <- plot_heatmap(bcg300_2, "Young (300BCG:n=283)")
p4 <- plot_heatmap(fg500_2, "Young (500FG:n=212)")

ppi = 300
png("output/fig3_heatmap_eaa_age_protein_sig_in_bcgprime_all_cohorts.png", width = 14*ppi, height = 14*ppi, res = ppi) 
ggarrange(p1, p2, p3, p4, ncol = 1)
dev.off()


