# cis-MR
# clump_kb threshold = 10,000
# clump_r2 threshold = 0.01
# cis_window = +-500000 bp
# MR methods:  one SNP wald ratio, at least two SNPs, inverse weighted, MR-PRESSO, 
# output_file_path = "../Downloads/GrimAge_EUR_summary_statistics.txt" or "../Downloads/PhenoAge_EUR_summary_statistics.txt"
# clump_kb_threshold = 10000
# clump_r2_threshold = 0.01
# cis_window = 500000
# F_stat_threshold = 10

config <- list(
  # Input / output
  pqtl_dir = "../Downloads/MR_input",            # folder containing pQTL _significant_snps.csv files
  rsid_map_dir = "../Downloads/rsid_map",        # folder containing chr*_rsid.tsv.gz files
  outcome_file = "../Downloads/GrimAge_EUR_summary_statistics.txt", # outcome GWAS file (change as needed)
  outdir = "../Downloads/MR_results",            # where to save results
  
  # Clumping / cis settings
  clump_kb = 10000,
  clump_r2 = 0.01,
  cis_window = 500000,       # +- cis window around TSS in bp
  F_threshold = 10,
  
  # LD reference (plink prefix; .bed/.bim/.fam expected)
  plink_bfile = "../Downloads/ref/EUR",
  
  # MR-PRESSO
  mr_presso_nb = 1000,
  
  # Runtime options
  verbose = TRUE
)

# Load library
library(TwoSampleMR)
library(dplyr)
library(data.table)
library(ieugwasr)
library(biomaRt)
library(VariantAnnotation)
library(GenomicRanges)
library(data.table)
library(MRPRESSO)

# log message
log_msg <- function(...) {
  if (isTRUE(config$verbose)) message(paste0("[", Sys.time(), "] ", paste0(...)))
}

safe_fread <- function(fp, ...) {
  if (!file.exists(fp)) stop(sprintf("File not found: %s", fp))
  data.table::fread(fp, ...)
}

# Ensure output dir
if (!dir.exists(config$outdir)) dir.create(config$outdir, recursive = TRUE)

# 1) Get transcripts and TSS from Ensembl (GRCh38)
log_msg("Fetching transcript info from Ensembl (GRCh38)...")
ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

transcripts_info <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "chromosome_name",
    "transcript_start",
    "transcript_end",
    "strand",
    "transcription_start_site",
    "gene_biotype",
    "transcript_biotype"
  ),
  filters = "transcript_biotype",
  values = "protein_coding",
  mart = ensembl
)

# Compute a single TSS per transcript row: use transcript_start for positive strand,
# transcript_end for negative strand. We'll later select one row per gene.
transcripts_info$TSS <- ifelse(transcripts_info$strand == 1,
                               transcripts_info$transcript_start,
                               transcripts_info$transcript_end)

# 2) Load list of significant proteins and harmonize gene names
protein_sig_path <- "/BCG_prime/output/03_protein/EAA_protein_correlation.csv"
if (!file.exists(protein_sig_path)) stop("protein_sig input CSV not found; update `protein_sig_path` in config section.")

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

# 3) Filter transcripts to only the significant proteins 
transcripts_info_sig <- transcripts_info %>%
  filter(external_gene_name %in% protein_sig_rename)

# remove non-numeric chromosomes 
transcripts_info_sig_clean <- transcripts_info_sig[!grepl("[A-Za-z]", transcripts_info_sig$chromosome_name) |
                                                     transcripts_info_sig$chromosome_name %in% c("X", "Y"), ]
transcripts_info_sig_unique <- transcripts_info_sig_clean[!duplicated(transcripts_info_sig_clean$ensembl_gene_id), ]

log_msg("Selected ", nrow(transcripts_info_sig_unique), " genes for MR analysis.")

# 4) Find pQTL files that match selected proteins
pqtl_files_all <- list.files(config$pqtl_dir, full.names = TRUE, pattern = "_significant_snps\\.csv$")
get_base <- function(path) tools::file_path_sans_ext(basename(path))

pqtl_files_selected <- pqtl_files_all[sapply(pqtl_files_all, function(fp) {
  base_name <- strsplit(basename(fp), "_")[[1]][1]
  base_name %in% protein_sig_rename
})]

log_msg("Found ", length(pqtl_files_selected), " pQTL files matching proteins of interest.")

# 5) rsID mapping files downloaded from pQTL literature
rsid_files <- list.files(config$rsid_map_dir, full.names = TRUE, pattern = "\\.tsv\\.gz$")
if (length(rsid_files) == 0) stop("No rsid mapping files found in rsid_map_dir")

load_rsid_map <- function(fp, ids_to_keep) {
  # fread can read gz compressed tsv
  dt <- data.table::fread(fp)
  if (!"ID" %in% colnames(dt)) stop("rsid mapping file must contain column 'ID'")
  dt[ID %in% ids_to_keep]
}

# 6) Main loop: iterate over each protein's pQTL file and run cis-MR
all_mr_results <- list()
all_presso_results <- list()

for (pqtl_fp in pqtl_files_selected) {
  protein_name <- strsplit(basename(pqtl_fp), "_")[[1]][1]
  log_msg("Processing protein: ", protein_name)
  
  tryCatch({
    pqtl <- data.table::fread(pqtl_fp)
    
    gene_info <- transcripts_info_sig_unique %>% filter(external_gene_name == protein_name)
    if (nrow(gene_info) == 0) {
      log_msg("  No gene info found for ", protein_name, " - skipping")
      next
    }
    
    chr_location <- as.character(gene_info$chromosome_name[1])
    tss <- as.integer(gene_info$TSS[1])
    
    # filter cis pQTLs
    cis_pqtl <- pqtl[CHROM == chr_location &
                       GENPOS >= (tss - config$cis_window) &
                       GENPOS <= (tss + config$cis_window)]
    
    if (nrow(cis_pqtl) == 0) {
      log_msg("  No cis pQTLs in window for ", protein_name)
      next
    }
    
    target_pattern <- paste0("chr", chr_location, "_")
    matching_rsid_files <- rsid_files[grepl(target_pattern, basename(rsid_files))]
    
    if (length(matching_rsid_files) == 0) {
      log_msg("  No rsid mapping files for chromosome ", chr_location, " - skipping")
      next
    }
    
    # collect mapped rsids
    mapping_list <- list()
    ids_to_find <- cis_pqtl$ID
    for (rf in matching_rsid_files) {
      mm <- load_rsid_map(rf, ids_to_find)
      if (nrow(mm) > 0) mapping_list[[rf]] <- mm
    }
    
    if (length(mapping_list) == 0) {
      log_msg("  No rsid mappings found for cis variants - skipping")
      next
    }
    
    mapping_combined <- data.table::rbindlist(mapping_list, fill = TRUE)
    
    # match and prepare exposure file for TwoSampleMR formatting
    matched <- cis_pqtl[match(mapping_combined$ID, cis_pqtl$ID), ]
    matched$rsid <- mapping_combined$rsid
    matched$SNP <- matched$rsid
    
    # convert values and name columns consistently
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
    
    exp_dat <- as.data.frame(matched %>%
                               dplyr::select(SNP, chr_name, chrom_start, beta, se, effect_allele, other_allele, eaf, pval, samplesize))
    
    exposure_data <- TwoSampleMR::format_data(exp_dat, type = "exposure", chr_col = "chr_name", pos_col = "chrom_start")
    
    # 6a) LD clumping using TwoSampleMR::ld_clump
    log_msg("  Performing LD clumping (plinker) ...")
    ld_clump_res <- TwoSampleMR::ld_clump(
      dplyr::tibble(rsid = exposure_data$SNP, pval = exposure_data$pval.exposure, id = "exposure"),
      plink_bin = 'plink',
      bfile = config$plink_bfile,
      clump_kb = config$clump_kb,
      clump_r2 = config$clump_r2
    )
    
    if (nrow(ld_clump_res) == 0) {
      log_msg("  No SNPs passed LD clumping for ", protein_name)
      next
    }
    
    exposure_data_clumped <- exposure_data %>% filter(SNP %in% ld_clump_res$rsid)
    
    # 6b) Read outcome GWAS and filter to clumped SNPs
    log_msg("  Reading outcome GWAS and subsetting to clumped SNPs...")
    outcome_dt <- data.table::fread(config$outcome_file)
    # normalize allele case
    if ("A1" %in% names(outcome_dt)) outcome_dt$A1 <- toupper(outcome_dt$A1)
    if ("A2" %in% names(outcome_dt)) outcome_dt$A2 <- toupper(outcome_dt$A2)
    
    if (!"rsID" %in% names(outcome_dt)) stop("Outcome file must contain 'rsID' column. Edit script if different.")
    
    outcome_filter <- outcome_dt[rsID %in% ld_clump_res$rsid]
    if (nrow(outcome_filter) == 0) {
      log_msg("  No overlapping SNPs between exposure and outcome after clumping - skipping")
      next
    }
    
    # format outcome for TwoSampleMR
    outcome_filter$SNP <- outcome_filter$rsID
    outcome_filter$beta <- outcome_filter$Effect
    outcome_filter$se <- outcome_filter$SE
    outcome_filter$effect_allele <- outcome_filter$A1
    outcome_filter$other_allele <- outcome_filter$A2
    outcome_filter$eaf <- outcome_filter$Freq1
    outcome_filter$pval <- outcome_filter$P
    outcome_filter$samplesize <- outcome_filter$N
    outcome_filter$chr_name <- outcome_filter$chr
    outcome_filter$chrom_start <- outcome_filter$bp
    
    out_dat <- TwoSampleMR::format_data(
      as.data.frame(outcome_filter %>% dplyr::select(SNP, beta, se, effect_allele, other_allele, eaf, pval, samplesize, chr_name, chrom_start)),
      type = "outcome",
      chr_col = "chr_name",
      pos_col = "chrom_start"
    )
    
    # keep only SNPs that are present in both datasets
    exposure_data_clumped <- exposure_data_clumped %>% filter(SNP %in% out_dat$SNP)
    if (nrow(exposure_data_clumped) == 0) {
      log_msg("  No overlapping SNPs after formatting - skipping")
      next
    }
    
    # compute F-statistic and filter weak instruments
    exposure_data_clumped$F_stat <- (exposure_data_clumped$beta.exposure^2) / (exposure_data_clumped$se.exposure^2)
    exposure_data_clumped <- exposure_data_clumped %>% filter(F_stat > config$F_threshold)
    if (nrow(exposure_data_clumped) == 0) {
      log_msg("  No instruments exceed F threshold for ", protein_name)
      next
    }
    
    # subset outcome to matching SNPs in the same order
    out_dat_match <- out_dat[match(exposure_data_clumped$SNP, out_dat$SNP), ]
    
    # harmonize
    harmonized <- TwoSampleMR::harmonise_data(exposure_data_clumped, out_dat_match, action = 2)
    
    # run MR analysis (TwoSampleMR::mr implements many methods)
    mr_res <- TwoSampleMR::mr(harmonized)
    mr_res$protein <- protein_name
    mr_res$ci_lower <- mr_res$b - 1.96 * mr_res$se
    mr_res$ci_upper <- mr_res$b + 1.96 * mr_res$se
    
    # record results
    all_mr_results[[protein_name]] <- mr_res
    
    # If at least two instruments, run heterogeneity and pleiotropy tests
    if (nrow(harmonized) >= 2) {
      het <- tryCatch(TwoSampleMR::mr_heterogeneity(harmonized), error = function(e) NULL)
      pleio <- tryCatch(TwoSampleMR::mr_pleiotropy_test(harmonized), error = function(e) NULL)
    } else {
      het <- NULL; pleio <- NULL
    }
    
    # MR-PRESSO (requires a data.frame with specific column names)
    mr_input <- data.frame(
      BetaOutcome = harmonized$beta.outcome,
      BetaExposure = harmonized$beta.exposure,
      SdOutcome = harmonized$se.outcome,
      SdExposure = harmonized$se.exposure
    )
    
    presso_res <- tryCatch({
      MRPRESSO::mr_presso(
        BetaOutcome = "BetaOutcome",
        BetaExposure = "BetaExposure",
        SdOutcome = "SdOutcome",
        SdExposure = "SdExposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        NbDistribution = config$mr_presso_nb,
        SignifThreshold = 0.05,
        data = mr_input
      )
    }, error = function(e) {
      log_msg("  MR-PRESSO failed for ", protein_name, ": ", e$message)
      NULL
    })
    
    all_presso_results[[protein_name]] <- presso_res
    
    # Save intermediate harmonized & MR results
    write.csv(harmonized, file = file.path(config$outdir, paste0(protein_name, "_harmonized.csv")), row.names = FALSE)
    write.csv(mr_res, file = file.path(config$outdir, paste0(protein_name, "_mr_results.csv")), row.names = FALSE)
    if (!is.null(het)) write.csv(het, file = file.path(config$outdir, paste0(protein_name, "_heterogeneity.csv")), row.names = FALSE)
    if (!is.null(pleio)) write.csv(pleio, file = file.path(config$outdir, paste0(protein_name, "_pleiotropy.csv")), row.names = FALSE)
    
  }, error = function(e) {
    log_msg("  Error processing ", protein_name, ": ", e$message)
  })
}

# 7) Combine and save summary results
if (length(all_mr_results) > 0) {
  mr_combined <- dplyr::bind_rows(all_mr_results, .id = "protein")
  write.csv(mr_combined, file = file.path(config$outdir, "mr_combined_results.csv"), row.names = FALSE)
  log_msg("Wrote MR combined results to: ", file.path(config$outdir, "mr_combined_results.csv"))
}

if (length(all_presso_results) > 0) {
  # mr_presso returns nested objects; we save the list as RDS
  saveRDS(all_presso_results, file = file.path(config$outdir, "mr_presso_results.rds"))
  log_msg("Saved MR-PRESSO results to: ", file.path(config$outdir, "mr_presso_results.rds"))
}

# 8) Save provenance information
session_file <- file.path(config$outdir, "sessionInfo.txt")
capture.output(sessionInfo(), file = session_file)
log_msg("Saved sessionInfo to: ", session_file)
log_msg("cis-MR pipeline finished.")

