#change file paths as needed
setwd("~/Dropbox/Yeast_Time_Series/Revisions_1/Regression")

#read in SNP table
snps2018_v2 <- read.table("~/Dropbox/Yeast_Time_Series/Revisions_1/SNP_table.txt", header = TRUE)
snps2018_v2=subset(snps2018_v2,snps2018_v2$ancmaf_0_0>0.025 & snps2018_v2$ancmaf_0_0 <0.975)

#Filter out "bad" populations (1,2,4,5 and 6)
X <- rep(c(3,7:12),each =12)
Y <- c(1:15,18)
maf <- paste("maf_", X, "_", Y, sep="")
cov <- paste("cov_", X, "_", Y, sep="")
good <- c(rbind(maf,cov))
snps2018_v2 <- snps2018_v2[c("chr","pos","ref","alt","ancmaf_0_0","anccov_0_0",good)]

#run regression analysis
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
  
  
  #modelB:time plus quadratic term
  timeQ <- lm(freq$Y~freq$week+I(freq$week^2),weights=sqrt(freq$cov))					
  pQ<- anova(timeQ)$"Pr(>F)"[1]
  
  
  #modelC:time and pop 
  ANC1=lm(freq$Y~freq$week+freq$line,weights=sqrt(freq$cov))											#model4:time and pop 
  ANC2=lm(freq$Y~freq$week+freq$line+freq$line:freq$week,weights=sqrt(freq$cov))						#model4:same as above, plus all interaction terms
  p.ancova=anova(ANC1,ANC2)$"Pr(>F)"[2]																# if significant, the extended model fits better (ANC2)
  
  #modelD:time, pop, and the quadratic term
  ANC3=lm(freq$Y~freq$week+I(freq$week^2)+freq$line,weights=sqrt(freq$cov))														
  ANC4=lm(freq$Y~freq$week+freq$line+I(freq$week^2)+freq$line:I(freq$week^2)+freq$line:freq$week,weights=sqrt(freq$cov))  		#model5:same as above, plus all interaction terms
  p.ancovaQ=anova(ANC3,ANC4)$"Pr(>F)"[2]																							# if significant, the extended model including a quadratic term fits better (ANC4)

  
  
  data.frame("p.week"=p.reg,"pQ"=pQ,"p.ancova"=p.ancova,"p.ancovaQ"=p.ancovaQ)
  
}

#apply function
models <- apply(snps2018_v2,1,function(x) fiterms(x,format))
models <-  do.call(rbind.data.frame,models)



#output table
p_table <- cbind(snps2018_v2[,1:4],models)
write.table(p_table,"Pval_table_17tp_All_Models.txt", quote= FALSE, row.names= TRUE,col.names= TRUE, sep="\t")


