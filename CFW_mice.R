pheno <- read.csv("pheno.csv")
geno <- read.csv(gzfile("geno.txt.gz"), sep=" ")
map <- read.csv("map.txt", sep=" ")


pheno_var <- pheno[(pheno$discard=="no"),c("id","TA", "EDL","bw1", "taillength", "tibia", "testisweight")]
pheno_var_noNa <- na.omit(pheno_var)

# plot the variables
library("PerformanceAnalytics")
chart.Correlation(pheno_var_noNa[,2:ncol(pheno_var_noNa)], histogram=TRUE, pch=19)


G = geno[(geno$discard=="no"),3:ncol(geno)]
P_var = pheno_var_noNa[,c("id", "BMD", "TA", "EDL", "taillength", "p120b1")]#,"PPIweight"

P_nodiscard_noNa <- na.omit(P_var)

PandG <- merge(P_nodiscard_noNa,geno, by = "id") # should combine in the end

P <- scale(PandG[,2:6])
G <- PandG[,8:ncol(PandG)]
G <- round(G)-1 # the original datast is 0, 1,2

maf_values <- apply(G, 2, calculate_maf) 

which(is.na(maf_values)) # no NA
snps_index <- which(maf_values > 0.01) 
cfw_gen_filter <- G[,snps_index]
cfw_map_filter <- map[snps_index, ]
write.csv(cfw_gen_filter, "cfw_gen_maffilter.csv",row.names = FALSE)
write.csv(ara_map_filter, "cfw_map_maffilter.csv",row.names = FALSE)


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
Z = as.matrix(cfw_gen_filter)
K = as.matrix(read.csv("K_cfw.csv"))
source("C:/Users/faf04/Dropbox/jason/Abroad/UCR/Research in UCR/Dr. Xu/project/transformed_method_BLUP/R_function/profiling_gwas.R")
source("C:/Users/jason/Dropbox/jason/Abroad/UCR/Research in UCR/Dr. Xu/project/transformed_method_BLUP/R_function/profiling_gwas.R")
ID = PandG[,1]
result_cfw <- profiling_gwas(as.matrix(P), t(Z), K, ID, "CFW")



## other eda