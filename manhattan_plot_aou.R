library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)

filtered_table <- fread("~/META_cis_region_CV_404.tsv", header=TRUE)

filtered_table$alleles_formatted <- gsub('\\["', "", filtered_table$alleles)
filtered_table$alleles_formatted <- gsub('"\\]', "", filtered_table$alleles_formatted) 
filtered_table$alleles_formatted <- gsub('","', ":", filtered_table$alleles_formatted)
filtered_table <- filtered_table %>% separate(alleles_formatted, into = c("REF", "ALT"), sep = ":", remove=F)
filtered_table$ID <- paste0(filtered_table$locus, ":", filtered_table$REF, ":", filtered_table$ALT)

filtered_table <- filtered_table %>% select(ID, Pvalue, Pvalue_log10, BETA, SE, CHR, POS)


#gsutil plink --bfile gs://fc-secure-915a27e9-c961-4aa0-b53f-09825278579c/data/plink --ld chr1:109274968:G:T chr1:109274570:A:G

write.table(filtered_table, "~/SORT1_cis_region.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

add_ld <- function(df) {
  df$ld <- NA
  
  for(i in 1:nrow(df)) {
    #plink command
    cmd <- paste("plink --bfile ~/plink --ld chr1:109274968:G:T", df$ID[i], "--noweb")
    output <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    
    #extract r2 value
    rsq_line <- grep("R-sq", output, value = TRUE)
    print(rsq_line)
    if(length(rsq_line) > 0) {
      rsq_value <- as.numeric(sub(".*R-sq = ([0-9.e+-]+).*", "\\1", rsq_line))
      df$ld[i] <- rsq_value
      cat(df$ld[i], "\n")
    }
    
    #progress
    if(i %% 10 == 0) cat("Done:", i, "\n")
  }
  
  return(df)
}

coloc_results <- add_ld(filtered_table)

write.table(coloc_results, "~/SORT1_ld.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#cis region: 108309568	110397918

coloc_results$ld <- as.numeric(coloc_results$ld)
coloc_results$ld[coloc_results$ld > 1] <- 0.1
coloc_results$ld_category <- cut(coloc_results$ld, 
                                 breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                                 labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"),
                                 include.lowest = TRUE)

p <- ggplot(coloc_results, aes(x = POS, y = Pvalue_log10)) +
  geom_point(aes(fill = ld_category), color = "black", shape = 21, size = 2.5, stroke = 0.2, alpha = 0.8) +
  scale_fill_manual(
    values = c("0-0.2" = "royalblue3",
               "0.2-0.4" = "lightblue2",
               "0.4-0.6" = "springgreen3",
               "0.6-0.8" = "orange1",
               "0.8-1.0" = "red2"
               ),
    name = expression(r^2)
  ) +
  scale_color_identity() +
  labs(
    x = "Chromosome 1 SNP Position",
    y = "AoU SORT1 CDS ± 1 Mb -log10(P)"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Add this line
    axis.title = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black", size = 0.8)
  ) +
  xlim(108309567, 110397919) +
  ylim(0, max(coloc_results$Pvalue_log10, na.rm = TRUE) * 1.05)

#highlight reference SNP
top_snp <- coloc_results[ID == "chr1:109274968:G:T", ]
p <- p +
  geom_point(data = top_snp,
             aes(x = POS, y = Pvalue_log10),
             color = "purple", size = 3.5, shape = 18) +
  geom_text(data = top_snp,
            aes(x = POS, y = Pvalue_log10),
            label = "chr1:109274968:G:T (rs12740374)",
            vjust = -0.8, hjust = 1, size = 4) #hjust > 1 moves label left

print(p)
ggsave("~/locuszoomplot_SORT1_ld.png", p, width = 8, height = 6, dpi = 300)


# gsutil -m cp -r -v /home/wheelerlab3/Data/predictdb_models/scPrediXcan_models/l-ctPred_models_for_immune_cell_types_from_OneK1K_dataset/ gs://fc-secure-d160515e-6f0a-4e0d-bd0a-71ee788a7b0c/data/
#   gsutil -m cp -r -v /home/wheelerlab3/Data/predictdb_models/scPrediXcan_models/l-ctPred_models_for_islet_cell_types_from_T2D_dataset/ gs://fc-secure-d160515e-6f0a-4e0d-bd0a-71ee788a7b0c/data/
