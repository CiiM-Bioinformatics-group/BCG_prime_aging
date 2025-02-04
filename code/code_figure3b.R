# Make the exposure list
library(data.table)
library(dplyr)
library(ieugwasr)
library(TwoSampleMR)
library(dplyr)

# get all folders with inflammation, 736 proteins
base_path <- "UKB-PPP/European/unzip/"
folders <- list.files(base_path, full.names = TRUE)
inflammation_folders <- folders[grepl("Inflammation", folders)]

# define p threshold
p_threshold <- 5e-08
log_p_threshold <- -log10(p_threshold)

read_and_filter <- function(folder_path) {

  files <- list.files(folder_path, full.names = TRUE)
  
  combined_data <- data.frame()
  
  for (file in files) {
    data <- fread(file)
    
    filtered_data <- data %>% filter(LOG10P > log_p_threshold)
    
    combined_data <- rbind(combined_data, filtered_data)
  }
  
  combined_data
}

output_base_path <- "input/MR_input/"

for (folder in inflammation_folders) {
  cat("Processing folder:", folder, "\n")
   
  filtered_snps <- read_and_filter(folder)
   
  output_file <- file.path(output_base_path, paste0(basename(folder), "_significant_snps.csv"))
  
  if (nrow(filtered_snps) > 0) {
    fwrite(filtered_snps, output_file)
    cat("Results saved to:", output_file, "\n")
  } else {
    cat("No significant SNPs found in:", folder, "\n")
  }
}

# sig. inflammatory names, 
protein_sig = read.csv("BCG-PRIME.csv") %>% filter(padj < 0.05) %>% pull(protein) %>% unique()

protein_sig <- gsub("\\.", "_", protein_sig)

file_path <- "input/MR_input"
files <- list.files(file_path, full.names = TRUE, pattern = "_significant_snps.csv$")

selected_files <- files[sapply(files, function(file) {
  base_name <- strsplit(basename(file), "_")[[1]][1]
  any(base_name %in% protein_sig)
})]

print(selected_files)
selected_files <- c(selected_files, "input/MR_input/TGFA_P01135_OID20600_v1_Inflammation_significant_snps.csv",
	"input/MR_input/LIFR_P42702_OID20606_v1_Inflammation_significant_snps.csv",
	"input/MR_input/IL10RB_Q08334_OID20696_v1_Inflammation_significant_snps.csv",
	"input/MR_input/IL18R1_Q13478_OID20652_v1_Inflammation_significant_snps.csv",
	"input/MR_input/IFNG_P01579_OID20495_v1_Inflammation_significant_snps.csv",
	"input/MR_input/CSF1_P09603_OID20719_v1_Inflammation_significant_snps.csv",
	"input/MR_input/MMP1_P03956_OID20672_v1_Inflammation_significant_snps.csv",
	"input/MR_input/MMP10_P09238_OID20687_v1_Inflammation_significant_snps.csv")
# Now we should use these select_files to combine and do clumping

# 23 proteins have been selected. 

setwd("input/MR_input/rsid_maps")
files <- list.files(pattern = "\\.tsv\\.gz$")

data_list <- list()

# Loop through each file, read it, and store it in the list
for (file in files) {
  data <- fread(file, header = TRUE, sep = "\t")
  data_list[[file]] <- data
}

mapping_combined_data <- rbindlist(data_list)

format_exposure_data <- function(file) {
  data <- fread(file)
  
  protein_name <- strsplit(basename(file), "_")[[1]][1]
  
  combined_data_filter <- mapping_combined_data %>% filter(ID %in% data$ID)
  
  matched_data <- data[match(combined_data_filter$ID, data$ID),]
  
  matched_data$rsid <- combined_data_filter$rsid
  
  matched_data$p_value <- 10^(-matched_data$LOG10P)
  
  matched_data$SNP <- matched_data$rsid
  matched_data$beta <- matched_data$BETA
  matched_data$se <- matched_data$SE
  matched_data$effect_allele <- matched_data$ALLELE1
  matched_data$other_allele <- matched_data$ALLELE0
  matched_data$eaf <- matched_data$A1FREQ
  matched_data$pval <- matched_data$p_value
  matched_data$samplesize <- matched_data$N
  matched_data$chr_name <- matched_data$CHROM
  matched_data$chrom_start <- matched_data$GENPOS
  matched_data$protein <- protein_name
  
  exp_dat <- matched_data %>% select(SNP, chr_name, chrom_start, beta, se, effect_allele, other_allele, eaf, pval, samplesize, protein)
  
  exp_dat <- as.data.frame(exp_dat)

  format_data(exp_dat, type = "exposure", chr_col = "chr_name", pos_col = "chrom_start")
}

exposure_list <- list()
combined_exposure_data <- data.table()

for (file in selected_files) {
  exposure_data <- format_exposure_data(file)
  
  protein_name <- strsplit(basename(file), "_")[[1]][1]
  
  exposure_list[[protein_name]] <- exposure_data
  
  combined_exposure_data <- rbind(combined_exposure_data, exposure_data, fill = TRUE)
}

## Do clumping for each exposure list and then combind with outcome data, run MR 

out_data = fread("input/MR_input/PhenoAge_EUR_summary_statistics.txt")

out_data$A1 <- toupper(out_data$A1)
out_data$A2 <- toupper(out_data$A2)


res_list = list()
heterogeneity_res_list = list()

for (i in 1:length(exposure_list)) {
  exposure_data <- exposure_list[[i]]
  protein_name <- names(exposure_list)[i]

  ld_reflookup(exposure_data$SNP, pop = "EUR", opengwas_jwt = get_opengwas_jwt())

  ld_clump_res <- ld_clump(
    dplyr::tibble(rsid=exposure_data$SNP, pval=exposure_data$pval.exposure, id="exposure"),
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = "input/MR_input/ref/EUR",
    clump_kb = 500) 

  exp_dat_clump = exposure_data %>% filter(SNP %in% ld_clump_res$rsid)

  out_data_filter = out_data %>% filter(rsID %in% ld_clump_res$rsid)


  out_data_filter$SNP = out_data_filter$rsID
  out_data_filter$beta = out_data_filter$Effect
  out_data_filter$se = out_data_filter$SE
  out_data_filter$effect_allele = out_data_filter$A1
  out_data_filter$other_allele = out_data_filter$A2
  out_data_filter$eaf = out_data_filter$Freq1
  out_data_filter$pval = out_data_filter$P
  out_data_filter$samplesize = out_data_filter$N
  out_data_filter$chr_name = out_data_filter$chr
  out_data_filter$chrom_start = out_data_filter$bp
  out_data_filter_new = out_data_filter %>% select(SNP, beta, se, effect_allele, other_allele, eaf, pval, samplesize, chr_name, chrom_start) %>% data.frame()


  out_dat <- format_data(out_data_filter_new, type = "outcome",
    chr_col = "chr_name",
    pos_col = "chrom_start")

  exp_dat_clump = exp_dat_clump %>% filter(SNP %in% out_dat$SNP)

  sum(exp_dat_clump$SNP == out_dat$SNP)
  out_dat_match = out_dat[match(exp_dat_clump$SNP, out_dat$SNP),]
  sum(exp_dat_clump$SNP == out_dat_match$SNP)

  hamor_data = harmonise_data(exp_dat_clump, out_dat_match, action = 2)
  res_list[[protein_name]] = mr(hamor_data)

  heterogeneity_res_list[[protein_name]] = mr_heterogeneity(hamor_data)
  
}

save(res_list, heterogeneity_res_list, file = "output/EAA_protein_MR_res.rdata")






