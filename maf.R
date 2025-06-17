
# Minor Allele Frequency (MAF) < 1%
calculate_maf <- function(snp) {
  # Count the number of each genotype
  n0 <- sum(snp == -1, na.rm = TRUE)  # Homozygous allele
  n1 <- sum(snp == 0, na.rm = TRUE)  # Heterozygous
  n2 <- sum(snp == 1, na.rm = TRUE)  # Homozygous allele
  
  # Calculate the allele counts
  total_alleles <- 2 * (n0 + n1 + n2)
  ref_alleles <- 2*n0+n1
  alt_alleles <- 2*n2 + n1
  ref_freq <- ref_alleles/total_alleles
  alt_freq <- alt_alleles/total_alleles
  maf <- min(ref_freq, alt_freq)
  return(maf)
}

# for all chr

maf_values <- apply(ara_gen, 2, calculate_maf)
snps_index <- which(maf_values > 0.01) 
ara_gen_filter <- ara_gen[,snps_index]
ara_map_filter <- ara_map[snps_index, ]
write.csv(ara_map_filter, "ara_LES_map_maffilter.csv",row.names = FALSE)
