# Figures for the manuscript #
library(dplyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggsci)
library(cowplot)

# Figure 1B
datasets <- list(
  `BCG-PRIME` = "input/BCG-PRIME.csv",
  iMed = "input/iMed.csv",
  `300BCG` = "input/300BCG.csv",  
  `500FG` = "input/500FG.csv"
)


generate_plot <- function(file_path, cohort_name) {
  dt <- read.csv(file_path)
  
  horvath <- dt %>% dplyr::select(SID, Age, DNAmAge) %>% mutate(Age_type = "Horvath") %>% na.omit()
  hannum <- dt %>% dplyr::select(SID, Age, DNAmAgeHannum) %>% mutate(Age_type = "Hannum") %>% na.omit()
  phenoage <- dt %>% dplyr::select(SID, Age, DNAmPhenoAge) %>% mutate(Age_type = "PhenoAge") %>% na.omit()
  grimage <- dt %>% dplyr::select(SID, Age, DNAmGrimAgeBasedOnRealAge) %>% mutate(Age_type = "GrimAge") %>% na.omit()
  
  colnames(horvath)[3] <- "DNAm_age"
  colnames(hannum)[3] <- "DNAm_age"
  colnames(phenoage)[3] <- "DNAm_age"
  colnames(grimage)[3] <- "DNAm_age"
  
  plot_data <- rbind(horvath, hannum, phenoage, grimage)
  
  plot_data$Age = as.numeric(plot_data$Age)
  plot_data$DNAm_age = as.numeric(plot_data$DNAm_age)
  
  cor_results <- plot_data %>%
    group_by(Age_type) %>%
    summarise(
      correlation = cor(Age, DNAm_age, method = "spearman"),
      p_value = cor.test(Age, DNAm_age, method = "spearman")$p.value
    )
  
  annotation_table <- cor_results %>%
    mutate(
      label = sprintf("%s: r = %.2f, p = %.3g", Age_type, correlation, p_value)
    )
  
  table_text <- paste(annotation_table$label, collapse = "\n")
  
  p <- ggplot(data = plot_data, aes(x = Age, y = DNAm_age, color = Age_type)) + 
    geom_point(alpha = 0.8) + 
    theme_classic() + 
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
    ggtitle(cohort_name) +
    xlab("Chronological Age") +
    ylab("DNA Methylation Age") +
    scale_color_npg() +
    theme(
      legend.position = "none", 
      plot.margin = unit(c(1, 1, 1, 1), "lines")
    )
  
  annotation_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = table_text, hjust = 0.5, vjust = 0.5, size = 4) +
    theme_void() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  return(list(main_plot = p, annotation_plot = annotation_plot))
}

plots <- lapply(names(datasets), function(name) {
  generate_plot(datasets[[name]], name)
})

main_plots <- lapply(plots, function(p) p$main_plot)
annotation_plots <- lapply(plots, function(p) p$annotation_plot)

shared_legend <- get_legend(
  main_plots[[1]] +
    theme(legend.position = "right", legend.title = element_blank())
)

combined_plots <- mapply(function(main, annotation) {
  grid.arrange(main, annotation, nrow = 2, heights = c(3, 1))
}, main_plots, annotation_plots, SIMPLIFY = FALSE)

ppi <- 300
png("output/figure1b_scatter_plots_realage_dnamage.png", 
    width = 12*ppi, height = 4*ppi, res = ppi)
grid.arrange(
  grobs = c(combined_plots, list(shared_legend)), 
  ncol = 5,  
  widths = c(4, 4, 4, 4, 1.5)  
)
dev.off()



