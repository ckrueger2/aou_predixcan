#if needed, install packages
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("coloc", quietly = TRUE)) install.packages("coloc")
if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse")
#load packages
library(data.table)
library(dplyr)
library(coloc)
library(argparse)
#set up argparse
parser <- ArgumentParser()
parser$add_argument("--phecode", help="all of us phenotype ID")
parser$add_argument("--pop", help="all of us population ID")
args <- parser$parse_args()
#find bucket
my_bucket <- Sys.getenv('WORKSPACE_BUCKET')
#read in MESA pQTL table
name_of_pqtl_file <- paste0(args$pop, "_mesa_pqtls.txt")
pqtl_command <- paste0("gsutil cp ", my_bucket, "/data/", name_of_pqtl_file, " .")
system(pqtl_command, intern=TRUE)
qtl_data <- fread(name_of_pqtl_file, header=TRUE)
#read in AoU GWAS table
name_of_gwas_file <- paste0(args$pop, "_formatted_mesa_", args$phecode, ".tsv")
gwas_command <- paste0("gsutil cp ", my_bucket, "/data/", name_of_gwas_file, " .")
system(gwas_command, intern=TRUE)
gwas_data <- fread(name_of_gwas_file, header=TRUE)
#extract unique phenotypes
unique_phenotypes <- unique(qtl_data$phenotype_id)
for (phenotype in unique_phenotypes) {
  cat("Processing phenotype:", phenotype, "\n")
  
  #filter QTL data for current phenotype
  qtl_subset <- qtl_data[qtl_data$phenotype_id == phenotype, ]
  
  #find common SNPs for each phenotype
  common_variants <- intersect(qtl_subset$variant_id, gwas_data$ID)
  
  if (length(common_variants) > 0) {
    # Filter to common SNPs
    qtl_coloc <- qtl_subset[qtl_subset$variant_id %in% common_variants, ]
    gwas_coloc <- gwas_data[gwas_data$ID %in% common_variants, ]
    
    #merge tables
    merged_data <- inner_join(gwas_coloc, qtl_coloc, by = c("ID" = "variant_id"))
    head(merged_data)
    
    pre_filter <- (nrow(merged_data))
    
    #remove duplicate SNPs
    duplicate_snps <- merged_data$ID[duplicated(merged_data$ID)]
    if (length(duplicate_snps) > 0) {
      cat("Found", length(duplicate_snps), "duplicate SNPs. Removing duplicates...\n")
      #keep most significant
      merged_data <- merged_data %>%
        group_by(ID) %>%
        slice_min(pval_nominal, n = 1, with_ties = FALSE) %>%
        ungroup()
    }
    post_filter <-(nrow(merged_data))
    cat("Pre-filter SNP count: ", pre_filter, " Post-filter SNP count: ", post_filter, "\n")
    
    #prepare datasets
    dataset1 <- list(
      beta = merged_data$slope,
      varbeta = merged_data$slope_se^2,
      snp = merged_data$ID,
      type = "quant",
      N = if (args$pop == "EUR") 1270 else 2953, #META: 2953, EUR: 1270
      MAF = merged_data$af
    )
    
    dataset2 <- list(
      beta = merged_data$BETA,
      varbeta = merged_data$SE^2,
      snp = merged_data$ID,
      type = "quant",
      N = 277593,
      sdY = 1
    )
    
    #run coloc
    result <- coloc.abf(dataset1, dataset2)
    
    #view summary
    cat("Results for", phenotype, ":\n")
    print(result$summary)
    
    #show top SNPs if colocalization found
    #show top SNPs if colocalization found
    if (result$summary["PP.H4.abf"] > 0.8) {
      cat("Colocalizing SNPs found:\n")
      top_snps <- result$results[order(-result$results$SNP.PP.H4), ]
      print(head(top_snps[, c("snp", "SNP.PP.H4")], 10))
      
      #install locuscomparer if needed
      if (!requireNamespace("locuscomparer", quietly = TRUE)) {
        if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
        devtools::install_github("boxiangliu/locuscomparer")
      }
      library(locuscomparer)
      
      #prepare data for locuscomparer
      pqtl_data <- merged_data %>%
        select(rsid = ID, 
               chr = CHR,
               pos = POS, 
               pval = pval_nominal) %>%
        write.table("/tmp/pqtl.tsv", sep="\t", row.names=FALSE, quote=FALSE)
      
      gwas_data_formatted <- merged_data %>%
        select(rsid = ID,
               chr = CHR, 
               pos = POS,
               pval = P) %>%
        write.table("/tmp/gwas.tsv", sep="\t", row.names=FALSE, quote=FALSE)
      
      #get lead SNP (most probable colocalization SNP)
      lead_snp <- top_snps$SNP.PP.H4[which.min(top_snps$SNP.PP.H4)]
      
      #retrieve LD information
      ld_data <- retrieve_LD(chr = merged_data$CHR[1], 
                             start = min(merged_data$POS), 
                             end = max(merged_data$POS),
                             lead_snp = lead_snp,
                             population = "META")
      
      #create the three-panel plot
      plot_filename <- paste0(phenotype, "_", args$phecode, "_locuscompare.png")
      png(plot_filename, width = 1200, height = 400)
      locuscompare(in_fn1 = "/tmp/gwas.tsv", 
                   in_fn2 = "/tmp/pqtl.tsv",
                   title1 = paste("GWAS:", args$phecode),
                   title2 = paste("pQTL:", phenotype),
                   ld = ld_data,
                   color_scheme = c("blue", "green", "orange", "red"))
      dev.off()
      
      #clean up tmp files
      file.remove("/tmp/pqtl.tsv", "/tmp/gwas.tsv")
      
      cat("Locuscompare plot saved as:", plot_filename, "\n")
    }
  } else {
    cat("No common variants found for", phenotype, "\n\n")
  }
}
