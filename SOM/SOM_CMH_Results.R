setwd("~/Dropbox/Yeast_Time_Series/Revisions_1/SOM/")

#list of sigs
library(poolSeq)
#read in SNP table
snps2018_v2 <- read.table("~/Dropbox/Yeast_Time_Series/Revisions_1/SNP_table.txt", header = TRUE)
snps2018_v2=subset(snps2018_v2,snps2018_v2$ancmaf_0_0>0.025 & snps2018_v2$ancmaf_0_0 <0.975)

#Filter out "bad" populations (1,2,4,5 and 6)
X <- rep(c(3,7:12),each =16)
Y <- c(1:15,18)
maf <- paste("maf_", X, "_", Y, sep="")
cov <- paste("cov_", X, "_", Y, sep="")
good <- c(rbind(maf,cov))
snps2018_v2 <- snps2018_v2[c("chr","pos","ref","alt","ancmaf_0_0","anccov_0_0",good)]

#create inputs
A0 <- snps2018_v2[,seq(7,ncol(snps2018_v2),by=32)] * snps2018_v2[,seq(8,ncol(snps2018_v2),by=32)]
a0 <- snps2018_v2[,seq(8,ncol(snps2018_v2),by=32)] - A0
A0 <- t(A0)
a0 <- t(a0)

At <- snps2018_v2[,seq(37,ncol(snps2018_v2),by=32)] * snps2018_v2[,seq(38,ncol(snps2018_v2),by=32)]
at <- snps2018_v2[,seq(38,ncol(snps2018_v2),by=32)] - At
At <- t(At)
at <- t(at)

pval <- cmh.test(A0=A0, a0=a0, At=At, at=at)
cmh_results <- na.omit(cbind(snps2018_v2[,1:2], pval))
cmh_results$pval <-  -log10(cmh_results$pval)
sig_snps <- subset(cmh_results, pval >= 13)
sig_list <- data.frame(sig_snps$chr, sig_snps$pos)
colnames(sig_list) <- c("chr","pos")

#create trans table
snps2018_v2 <- read.table("~/Dropbox/Yeast_Time_Series/Revisions_1/SNP_table.txt", header = TRUE)
snps2018_v2=subset(snps2018_v2,snps2018_v2$ancmaf_0_0>0.025 & snps2018_v2$ancmaf_0_0 <0.975)
snps2018_v2 <- subset(snps2018_v2, chr != "chrmito")
#Filter out "bad" populations (1,2,4,5 and 6)
X <- rep(c(3,7:12),each =16)
Y <- c(1:15,18)
maf <- paste("maf_", X, "_", Y, sep="")
cov <- paste("cov_", X, "_", Y, sep="")
good <- c(rbind(maf,cov))
snps2018_v2 <- snps2018_v2[c("chr","pos","ref","alt","ancmaf_0_0","anccov_0_0",good)]
snps2018_v2 <- merge(snps2018_v2, sig_list)
maf <- snps2018_v2[,c(seq(7,ncol(snps2018_v2),by = 2))]

#find snps where mean final freq lower than starting
test <- cbind(snps2018_v2[,c(1:2,5)],maf[,c(seq(16,ncol(maf), by = 16))])
mean_alt <- rowMeans(test[,4:ncol(test)])
test2 <- cbind(test, mean_alt)
test3 <- subset(test2, mean_alt < ancmaf_0_0) #lower than starting
test4 <-  subset(test2, mean_alt > ancmaf_0_0) #higher than starting
lower <- test3[,1:2]
higher <- test4[,1:2]

#subset lower ending than start and transforms
snps2018_v3 <- merge(snps2018_v2, lower)
snps2018_v3[,c(seq(5, ncol(snps2018_v3), by = 2))] <- 1- snps2018_v3[,c(seq(5, ncol(snps2018_v3), by = 2))]
snps2018_v3[,3:4] <- snps2018_v3[,4:3]

#subset higher ending, combing with lower ending
snps2018_v4 <- merge(snps2018_v2, higher)

#recombine
snps2018_v5 <- rbind(snps2018_v3,snps2018_v4)


#make new traj table
maf <- snps2018_v5[,c(seq(7,ncol(snps2018_v5),by = 2))]
Sig_Traj_table_uni_dir <- data.frame(matrix(ncol = 1, nrow = 17))

for (i in 1:nrow(maf)){
  
  chr <- snps2018_v5[i,1]
  position <- snps2018_v5[i,2]
  pop <- c(3,7:12)
  col_names <- paste("pop",pop,chr,position, sep = "_")
  maf_snp <- as.numeric(maf[i,])
  matrix_maf <- data.frame(matrix(maf_snp, nrow = 16, ncol = 7))
  ancestor_maf <- rep(snps2018_v5$ancmaf_0_0[i],16)
  maf_final <- rbind(ancestor_maf, matrix_maf)
  colnames(maf_final) <- col_names
  Sig_Traj_table_uni_dir <- cbind(Sig_Traj_table_uni_dir, maf_final)
  
  print(i)
}

Sig_Traj_table_uni_dir <- Sig_Traj_table_uni_dir[,2:ncol(Sig_Traj_table_uni_dir)]


Trans_traj_sig_uni_dir <- t(Sig_Traj_table_uni_dir)

#save(Trans_traj_sig_uni_dir , file="SOM_table_sig_all.rda")
library("som", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

Full_sig_uni_traj <- normalize(Trans_traj_sig_uni_dir)
Full_sig_uni <- som(Full_sig_uni_traj, xdim=3, ydim=1)
plot(Full_sig_uni, main = "Candidate SNPs (CMH)", ylim=c(-2, 2), cex.axis=2,cex.main =2)

