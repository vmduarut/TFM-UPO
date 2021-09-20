library(biomaRt)
setwd("F:/TFM/popGenome_analysis_october/30000_500/outliers")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
results_negative <- vector()
results_positive <- vector()
for (j in 1:22) {
    
  chr_negative <- read.table(paste("negative_outliers_chr", j,"_30000_500_IBS_updated.tsv", sep = ""), header = TRUE, sep = "\t")
  plot(chr_negative$Tajima.D ~ chr_negative$Window_start, main = paste("negative_chromosome_", j, sep = ""), xlab = "Window start", ylab = "Tajima")
  regions <- vector()
  
  for (i in 1:nrow(chr_negative)) {
    regions[i] <- paste(j, chr_negative$Window_start[i], chr_negative$Window_end[i], sep = ":")
    ## Aqu? en el bucle se introducir?a un condicional if por ejemplo para coger solo unas venatanas espec?ficas
  }
  genes <- getBM(attributes = c("hgnc_symbol"), filters = c("chromosomal_region"), values = list(regions), mart = ensembl)
  results_negative <- c(results_negative, genes$hgnc_symbol)

  chr_positive <- read.table(paste("positive_outliers_chr", j,"_30000_500_IBS_updated.tsv", sep = ""), header = TRUE, sep = "\t")
  plot(chr_positive$Tajima.D ~ chr_positive$Window_start, main = paste("positive_chromosome_", j, sep = ""), xlab = "Window start", ylab = "Tajima")
  regions <- vector()
  
  for (i in 1:nrow(chr_positive)) {
    regions[i] <- paste(j, chr_positive$Window_start[i], chr_positive$Window_end[i], sep = ":")
  }
  genes <- getBM(attributes = c("hgnc_symbol"), filters = c("chromosomal_region"), values = list(regions), mart = ensembl)
  results_positive <- c(results_positive, genes$hgnc_symbol)
}

write.table(results_negative, file = "negative_genes_hgnc_symbol.panel", row.names = FALSE, col.names = FALSE)
write.table(results_positive, file = "positive_genes_hgnc_symbol.panel", row.names = FALSE, col.names = FALSE)