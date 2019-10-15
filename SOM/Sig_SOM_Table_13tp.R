setwd("~/Dropbox/Yeast_Time_Series/Revisions_1/SOM/")
#list of sigs
Full <- read.table("~/Dropbox/Yeast_Time_Series/Revisions_1/Regression/Pval_table_13tp_ModelA.txt", header = TRUE)

Full$p.week <- -log10(Full$p.week)
sig_snps <- subset(Full, p.week >= 7)
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
#downsample to oringinal 
X <- rep(c(3,7:12),each =12)
Y <- c(1,2,3,5,7,9,10,11,12,13,15,18)
maf <- paste("maf_", X, "_", Y, sep="")
cov <- paste("cov_", X, "_", Y, sep="")
the_rest <- c(rbind(maf,cov))
snps2018_v5 <- snps2018_v5[c("chr","pos","ref","alt","ancmaf_0_0","anccov_0_0",the_rest)]
snps2018_v5 <- merge(snps2018_v5, sig_list)
maf <- snps2018_v5[,c(seq(7,ncol(snps2018_v5),by = 2))]

Sig_Traj_table_uni_dir_ds13 <- data.frame(matrix(ncol = 1, nrow = 13))


for (i in 1:nrow(maf)){
  
  chr <- snps2018_v5[i,1]
  position <- snps2018_v5[i,2]
  pop <- c(3,7:12)
  col_names <- paste("pop",pop,chr,position, sep = "_")
  maf_snp <- as.numeric(maf[i,])
  matrix_maf <- data.frame(matrix(maf_snp, nrow = 12, ncol = 7))
  ancestor_maf <- rep(snps2018_v5$ancmaf_0_0[i],12)
  maf_final <- rbind(ancestor_maf, matrix_maf)
  colnames(maf_final) <- col_names
  Sig_Traj_table_uni_dir_ds13 <- cbind(Sig_Traj_table_uni_dir_ds13, maf_final)
  
  print(i)
}

Sig_Traj_table_uni_dir_ds13 <- Sig_Traj_table_uni_dir_ds13[,2:ncol(Sig_Traj_table_uni_dir_ds13)]



Trans_traj_sig_uni_dir_ds13 <- t(Sig_Traj_table_uni_dir_ds13)
save(Trans_traj_sig_uni_dir_ds13 , file="SOM_table_sig_13tp.rda")


