#change file paths as needed
setwd("~/Dropbox/Yeast_Time_Series/Revisions_1/Regression")

#read in SNP table
snps2018_v2 <- read.table("~/Dropbox/Yeast_Time_Series/Revisions_1/SNP_table.txt", header = TRUE)
snps2018_v2=subset(snps2018_v2,snps2018_v2$ancmaf_0_0>0.025 & snps2018_v2$ancmaf_0_0 <0.975)

#Filter out "bad" populations (1,2,4,5 and 6) and downsample to 13 timepoints
X <- rep(c(3,7:12),each =12)
Y <- c(1,2,3,5,7,9,10,11,12,13,15,18)
maf <- paste("maf_", X, "_", Y, sep="")
cov <- paste("cov_", X, "_", Y, sep="")
the_rest <- c(rbind(maf,cov))
snps2018_v2 <- snps2018_v2[c("chr","pos","ref","alt","ancmaf_0_0","anccov_0_0",the_rest)]

#prep for function
index=seq(5,ncol(snps2018_v2),2)
names=colnames(snps2018_v2)[index] 
col2 <- sapply(names,function(x) unlist(strsplit(x,"_"))[2]) 
col3 <- sapply(names,function(x) unlist(strsplit(x,"_"))[3]) 
format <- data.frame(id=index,line=as.factor(col2),week=as.numeric(col3)) 




#fitting function
fiterms <- function(xx,format){
  
  freq <- cbind(asin(sqrt(as.numeric(xx[format$id]))),as.numeric(xx[format$id+1]),format)
  colnames(freq)[1:2] <- c("Y","cov")
  
  #modelA: time as a continuous variable
  time.reg <- lm(freq$Y~freq$week,weights=sqrt(freq$cov))  						
  p.reg <- anova(time.reg)$"Pr(>F)"[1]
  
  data.frame("p.week"=p.reg)
}

#apply function
models <- apply(snps2018_v2,1,function(x) fiterms(x,format))
models <-  do.call(rbind.data.frame,models)



#output table
p_table <- cbind(snps2018_v2[,1:4],models)

write.table(p_table,"Pval_table_13tp_ModelA.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")
