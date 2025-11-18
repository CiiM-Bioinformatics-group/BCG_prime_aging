# trans_MR.R

# Configuration 
config <- list(
  # Input / output
  pqtl_dir = "../Downloads/MR_input",            # folder containing pQTL _significant_snps.csv files
  rsid_map_dir = "../Downloads/rsid_map",       # folder containing chr*_rsid.tsv.gz files
  outcome_file = "../GrimAge_EUR_summary_statistics.txt", # outcome GWAS file (change as needed)
  mapping_output = "../Downloads/all_trans_pqtl_rsid_mapping.csv.gz",
  outdir = "../Downloads/MR_results_trans",
  
  # Clumping settings
  clump_kb = 1000,
  clump_r2 = 0.001,
  F_threshold = 10,
  
  # LD reference (plink prefix: .bed/.bim/.fam expected)
  plink_bfile = "../Downloads/ref/EUR",
  
  # MR-PRESSO
  mr_presso_nb = 1000,
  
  # Runtime options
  verbose = TRUE
)

# Load libraries 
library(TwoSampleMR)
library(dplyr)
library(data.table)
library(ieugwasr)
library(biomaRt)
library(MRPRESSO)

# log message
log_msg <- function(...) { if (isTRUE(config$verbose)) message(paste0("[", Sys.time(), "] ", paste0(...))) }

safe_fread <- function(fp, ...) { if (!file.exists(fp)) stop(sprintf("File not found: %s", fp)); data.table::fread(fp, ...) }

if (!dir.exists(config$outdir)) dir.create(config$outdir, recursive = TRUE)

# 1) Load significant proteins list and map names to HGNC symbols 
protein_sig_path <- "../BCG_prime/output/03_protein/EAA_protein_correlation.csv"
if (!file.exists(protein_sig_path)) stop("protein_sig input CSV not found; update path in config section.")

protein_sig_raw <- data.table::fread(protein_sig_path)
protein_sig <- unique(protein_sig_raw[padj < 0.05, protein])
protein_sig <- gsub("\\.", "_", protein_sig)

# map non-standard names to HGNC symbols 
protein_map <- c(
  "TGF_alpha" = "TGFA",
  "LIF_R" = "LIFR",
  "IL_10RB" = "IL10RB",
  "MMP_1" = "MMP1",
  "MMP_10" = "MMP10",
  "IFN_gamma" = "IFNG",
  "IL_18R1" = "IL18R1",
  "FGF_21" = "FGF21",
  "CSF_1" = "CSF1",
  "PD_L1" = "CD274",
  "Flt3L" = "FLT3LG",
  "X4E_BP1" = "EIF4EBP1",
  "EN_RAGE" = "S100A12",
  "MCP_2" = "CCL8",
  "LAP_TGF_beta_1" = "TGFB1"
)
protein_sig_rename <- ifelse(protein_sig %in% names(protein_map), protein_map[protein_sig], protein_sig)

# 2) Discover pqtl files to analyze
pqtl_files_all <- list.files(config$pqtl_dir, full.names = TRUE, pattern = "_significant_snps\\.csv$")
if (length(pqtl_files_all) == 0) stop("No pqtl files found in pqtl_dir")

pqtl_files_selected <- pqtl_files_all[sapply(pqtl_files_all, function(fp) {
  base_name <- strsplit(basename(fp), "_")[[1]][1]
  base_name %in% protein_sig_rename
})]
log_msg("Found ", length(pqtl_files_selected), " protein files for trans-MR analysis")

# 3) Collect all trans SNP IDs and map to rsIDs
all_trans_snp_ids <- character()
for (pqtl_fp in pqtl_files_selected) {
  pqtl <- data.table::fread(pqtl_fp)
  # For trans analysis we use all SNPs (no cis filtering)
  all_trans_snp_ids <- unique(c(all_trans_snp_ids, pqtl$ID))
}
log_msg("Total unique trans SNP IDs to map: ", length(all_trans_snp_ids))

# read mapping files from rsid_map_dir (assume tab-delimited gz files with columns: ID, rsid, chrom, pos, ref, alt)
rsid_map_files <- list.files(config$rsid_map_dir, full.names = TRUE, pattern = "\\.tsv\\.gz$")
if (length(rsid_map_files) == 0) stop("No rsid mapping files found in rsid_map_dir")

# iterate and collect matches 
all_rsid_mapping <- data.table()
for (rf in rsid_map_files) {
  log_msg("Reading rsid map: ", basename(rf))
  cols_keep <- c("ID", "rsid", "chrom", "pos", "ref", "alt")
  tmp <- data.table::fread(rf, select = intersect(cols_keep, names(data.table::fread(rf, nrows = 0))))
  if (nrow(tmp) == 0) next
  matched <- tmp[ID %in% all_trans_snp_ids]
  if (nrow(matched) > 0) all_rsid_mapping <- data.table::rbindlist(list(all_rsid_mapping, matched), use.names = TRUE, fill = TRUE)
}
all_rsid_mapping <- unique(all_rsid_mapping)
log_msg("Successfully mapped ", nrow(all_rsid_mapping), " variants (", round(nrow(all_rsid_mapping) / length(all_trans_snp_ids) * 100, 2), "%)")

# save mapping for provenance
data.table::fwrite(all_rsid_mapping, config$mapping_output)

# 4) Main loop: for each protein, prepare exposure (rsid-mapped), clump, harmonize with outcome and run MR
all_mr_results <- list()
all_heterogeneity <- list()
all_pleiotropy <- list()
all_presso <- list()

for (pqtl_fp in pqtl_files_selected) {
  protein_name <- strsplit(basename(pqtl_fp), "_")[[1]][1]
  log_msg("Processing protein: ", protein_name)
  
  tryCatch({
    pqtl <- data.table::fread(pqtl_fp)
    pqtl_rsmatch <- all_rsid_mapping[ID %in% pqtl$ID]
    if (nrow(pqtl_rsmatch) == 0) { log_msg("  No rsid mapping for this protein - skipping"); next }
    
    matched <- pqtl[match(pqtl_rsmatch$ID, pqtl$ID), ]
    matched$rsid <- pqtl_rsmatch$rsid
    matched$SNP <- matched$rsid
    
    # standardize columns for TwoSampleMR
    matched$p_value <- 10^(-matched$LOG10P)
    matched$beta <- matched$BETA
    matched$se <- matched$SE
    matched$effect_allele <- matched$ALLELE1
    matched$other_allele <- matched$ALLELE0
    matched$eaf <- matched$A1FREQ
    matched$pval <- matched$p_value
    matched$samplesize <- matched$N
    matched$chr_name <- matched$CHROM
    matched$chrom_start <- matched$GENPOS
    
    exp_df <- as.data.frame(matched %>% dplyr::select(SNP, chr_name, chrom_start, beta, se, effect_allele, other_allele, eaf, pval, samplesize))
    exposure_data <- TwoSampleMR::format_data(exp_df, type = "exposure", chr_col = "chr_name", pos_col = "chrom_start")
    
    # LD clumping
    log_msg("  Clumping exposure SNPs (kb=", config$clump_kb, ", r2=", config$clump_r2, ")")
    ld_clump_res <- TwoSampleMR::ld_clump(
      dplyr::tibble(rsid = exposure_data$SNP, pval = exposure_data$pval.exposure, id = "exposure"),
      plink_bin = 'plink',
      bfile = config$plink_bfile,
      clump_kb = config$clump_kb,
      clump_r2 = config$clump_r2
    )
    if (nrow(ld_clump_res) == 0) { log_msg("  No SNPs after clumping - skipping"); next }
    
    exposure_clumped <- exposure_data %>% filter(SNP %in% ld_clump_res$rsid)
    
    # Read outcome GWAS and subset to clumped SNPs
    outcome_dt <- data.table::fread(config$outcome_file)
    # standardize allele case if present
    if ("A1" %in% names(outcome_dt)) outcome_dt$A1 <- toupper(outcome_dt$A1)
    if ("A2" %in% names(outcome_dt)) outcome_dt$A2 <- toupper(outcome_dt$A2)
    if (!"rsID" %in% names(outcome_dt)) stop("Outcome file must contain 'rsID' column; adjust script if different.")
    
    outcome_sub <- outcome_dt[rsID %in% ld_clump_res$rsid]
    if (nrow(outcome_sub) == 0) { log_msg("  No overlapping SNPs with outcome - skipping"); next }
    
    outcome_sub$SNP <- outcome_sub$rsID
    outcome_sub$beta <- outcome_sub$Effect
    outcome_sub$se <- outcome_sub$SE
    outcome_sub$effect_allele <- outcome_sub$A1
    outcome_sub$other_allele <- outcome_sub$A2
    outcome_sub$eaf <- outcome_sub$Freq1
    outcome_sub$pval <- outcome_sub$P
    outcome_sub$samplesize <- outcome_sub$N
    outcome_sub$chr_name <- outcome_sub$chr
    outcome_sub$chrom_start <- outcome_sub$bp
    
    out_df <- as.data.frame(outcome_sub %>% dplyr::select(SNP, beta, se, effect_allele, other_allele, eaf, pval, samplesize, chr_name, chrom_start))
    out_dat <- TwoSampleMR::format_data(out_df, type = "outcome", chr_col = "chr_name", pos_col = "chrom_start")
    
    # Keep overlapping SNPs
    exposure_clumped <- exposure_clumped %>% filter(SNP %in% out_dat$SNP)
    if (nrow(exposure_clumped) == 0) { log_msg("  No overlapping SNPs after formatting - skipping"); next }
    
    # Filter weak instruments by F-stat
    exposure_clumped$F_stat <- (exposure_clumped$beta.exposure^2) / (exposure_clumped$se.exposure^2)
    exposure_clumped <- exposure_clumped %>% filter(F_stat > config$F_threshold)
    if (nrow(exposure_clumped) == 0) { log_msg("  No instruments exceed F threshold - skipping"); next }
    
    out_dat_match <- out_dat[match(exposure_clumped$SNP, out_dat$SNP), ]
    harmonized <- TwoSampleMR::harmonise_data(exposure_clumped, out_dat_match, action = 2)
    
    # Run MR
    mr_res <- TwoSampleMR::mr(harmonized)
    mr_res$protein <- protein_name
    mr_res$ci_lower <- mr_res$b - 1.96 * mr_res$se
    mr_res$ci_upper <- mr_res$b + 1.96 * mr_res$se
    all_mr_results[[protein_name]] <- mr_res
    
    # Heterogeneity & pleiotropy when applicable
    if (nrow(harmonized) >= 2) {
      all_heterogeneity[[protein_name]] <- tryCatch(TwoSampleMR::mr_heterogeneity(harmonized), error = function(e) NULL)
      all_pleiotropy[[protein_name]] <- tryCatch(TwoSampleMR::mr_pleiotropy_test(harmonized), error = function(e) NULL)
    }
    
    # MR-PRESSO
    mr_input <- data.frame(BetaOutcome = harmonized$beta.outcome,
                           BetaExposure = harmonized$beta.exposure,
                           SdOutcome = harmonized$se.outcome,
                           SdExposure = harmonized$se.exposure)
    presso_res <- tryCatch({
      MRPRESSO::mr_presso(BetaOutcome = "BetaOutcome", BetaExposure = "BetaExposure",
                          SdOutcome = "SdOutcome", SdExposure = "SdExposure",
                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                          NbDistribution = config$mr_presso_nb, SignifThreshold = 0.05,
                          data = mr_input)
    }, error = function(e) { log_msg("  MR-PRESSO error: ", e$message); NULL })
    all_presso[[protein_name]] <- presso_res
    
    # Save per-protein outputs
    write.csv(harmonized, file = file.path(config$outdir, paste0(protein_name, "_harmonized.csv")), row.names = FALSE)
    write.csv(mr_res, file = file.path(config$outdir, paste0(protein_name, "_mr_results.csv")), row.names = FALSE)
    if (!is.null(all_heterogeneity[[protein_name]])) write.csv(all_heterogeneity[[protein_name]], file = file.path(config$outdir, paste0(protein_name, "_heterogeneity.csv")), row.names = FALSE)
    if (!is.null(all_pleiotropy[[protein_name]])) write.csv(all_pleiotropy[[protein_name]], file = file.path(config$outdir, paste0(protein_name, "_pleiotropy.csv")), row.names = FALSE)
    
  }, error = function(e) { log_msg("  Error processing ", protein_name, ": ", e$message) })
}

# 5) Combine & save summary results
if (length(all_mr_results) > 0) {
  mr_combined <- dplyr::bind_rows(all_mr_results, .id = "protein")
  write.csv(mr_combined, file = file.path(config$outdir, "trans_mr_combined_results.csv"), row.names = FALSE)
  log_msg("Saved combined MR results to: ", file.path(config$outdir, "trans_mr_combined_results.csv"))
}

if (length(all_presso) > 0) {
  saveRDS(all_presso, file = file.path(config$outdir, "trans_mr_presso_results.rds"))
  log_msg("Saved MR-PRESSO results to: ", file.path(config$outdir, "trans_mr_presso_results.rds"))
}

# Save session info
capture.output(sessionInfo(), file = file.path(config$outdir, "sessionInfo_trans_mr.txt"))
log_msg("trans-MR pipeline finished.")
