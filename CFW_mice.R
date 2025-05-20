pheno <- read.csv("pheno.csv")
geno <- read.csv(gzfile("geno.txt.gz"), sep=" ")
map <- read.csv("map.txt", sep=" ")

### phenotype selection
pheno_var <- pheno[(pheno$discard=="no"),c("id","TA", "EDL","bw1", "taillength", "tibia", "testisweight")]
pheno_var_noNa <- na.omit(pheno_var)

p = ncol(pheno_var_noNa)-1 # number of traits

# plot the variables
library("PerformanceAnalytics")
chart.Correlation(pheno_var_noNa[,2:ncol(pheno_var_noNa)], histogram=TRUE, pch=19)


### genotype 
geno_nodiscard = geno[(geno$discard=="no"),]

k = ncol(geno_nodiscard) -2# number of mar

PandG <- merge(geno_nodiscard,pheno_var_noNa, by = "id") # inner join (only include common row)

P <- scale(PandG[,(3+k):ncol(PandG)])
write.table(P, "cfw_gemma_pheno.txt",sep="\t", col.names = FALSE, row.names = FALSE)
write.csv(P, "cfw_pheno.csv", row.names = FALSE)

G <- PandG[,3:(2+k)]
G_gemma <- cbind(map[, c("id","alt","ref")], t(G))
write.table(G_gemma, "cfw_gemma_geno.txt", sep=",",col.names = FALSE, row.names = FALSE)

map_gemma<- map[,c("id","pos","chr")]
write.table(map_gemma, "cfw_gemma_annot.txt",sep="\t", col.names = FALSE, row.names = FALSE)

G <- round(G)-1 # the original datast is 0, 1,2

maf_values <- apply(G, 2, calculate_maf) 

which(is.na(maf_values)) # no NA
snps_index <- which(maf_values > 0.05) 
cfw_gen_filter <- G[,snps_index]
cfw_map_filter <- map[snps_index, ]
write.csv(cfw_gen_filter, "cfw_gen_maffilter.csv",row.names = FALSE)
write.csv(cfw_map_filter, "cfw_map_maffilter.csv",row.names = FALSE)


cfw_map_filter <- read.csv("cfw_map_maffilter.csv")
cfw_gen_filter <-read.csv("cfw_gen_maffilter.csv")

## kinship matrix
Z = as.matrix(cfw_gen_filter)
K = Z%*%t(Z)
p = ncol(Z)
K = K/ p 
write.csv(K, "K_cfw.csv", row.names = FALSE)


## profiling gwas
cfw_gen_filter <- read.csv("cfw_gen_maffilter.csv")
P <- read.csv("cfw_pheno.csv")

Z = as.matrix(cfw_gen_filter)
K = as.matrix(read.csv("K_cfw.csv"))
source("C:/Users/faf04/Dropbox/jason/Abroad/UCR/Research in UCR/Dr. Xu/project/transformed_method_BLUP/R_function/profiling_gwas.R")
source("C:/Users/jason/Dropbox/jason/Abroad/UCR/Research in UCR/Dr. Xu/project/transformed_method_BLUP/R_function/profiling_gwas.R")
ID = PandG[,1]
result_cfw <- profiling_gwas(as.matrix(P), t(Z), K, ID, "CFW")

# manhatton plot
library(qqman)
gwasinfo <- data.frame(
  SNP = cfw_map_filter$id,
  CHR =  cfw_map_filter$chr,
  BP = cfw_map_filter$pos,
  P =result_cfw$log10p
)
max_p = max(result_cfw$log10p)

genomewideline =-log10(0.05/k)
manhattan(gwasinfo,col = c("blue4","cyan3"),main="CFW",ylim=c(0,max_p+2),cex=2,
          logp =FALSE,
          genomewideline=genomewideline ,#50.17457,
          ylab="-log10p",
          suggestiveline=F)
