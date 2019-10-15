setwd("~/Dropbox/Yeast_Time_Series/Revisions_1/SOM/")
library("som", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

par(mfrow=c(4,2))
par(mfrow=c(4,2),mar=c(2,1,3,2),mgp=c(3,1,0), oma=c(1,1,1,1))





#####DS4 SIG#########
load("SOM_table_sig_4tp.rda")
DS4_sig_uni_traj <- normalize(Trans_traj_sig_uni_dir_ds4)
DS4_sig_uni <- som(DS4_sig_uni_traj, xdim=3, ydim=1)
plot(DS4_sig_uni, main = "A. Candidate SNPs (T=4)", ylim=c(-2, 2), cex.axis=2, cex.main =2)

#####DS4 Rand#########
load("SOM_table_random_4tp.rda")
DS4_rand_uni_traj <- normalize(Trans_traj_rand_uni_dir_ds4)
DS4_rand_uni <- som(DS4_rand_uni_traj, xdim=3, ydim=1)
plot(DS4_rand_uni, main = "F. Random SNPs (T=4)", ylim=c(-2, 2), cex.axis=2,cex.main =2)



#####DS8 SIG#########
load("SOM_table_sig_8tp.rda")
DS8_sig_uni_traj <- normalize(Trans_traj_sig_uni_dir_ds8)
DS8_sig_uni <- som(DS8_sig_uni_traj, xdim=3, ydim=1)
plot(DS8_sig_uni, main = "B. Candidate SNPs (T=8)", ylim=c(-2, 2), cex.axis=2,cex.main =2)

#####DS8 Rand#########
load("SOM_table_random_8tp.rda")
DS8_rand_uni_traj <- normalize(Trans_traj_rand_uni_dir_ds8)
DS8_rand_uni <- som(DS8_rand_uni_traj, xdim=3, ydim=1)
plot(DS8_rand_uni, main = "G. Random SNPs (T=8)", ylim=c(-2, 2), cex.axis=2,cex.main =2)


#####DS13 SIG#########
load("SOM_table_sig_13tp.rda")
DS13_sig_uni_traj <- normalize(Trans_traj_sig_uni_dir_ds13)
DS13_sig_uni <- som(DS13_sig_uni_traj, xdim=3, ydim=1)
plot(DS13_sig_uni, main = "C. Candidate SNPs (T=13)", ylim=c(-2, 2), cex.axis=2,cex.main =2)

#####DS13 Rand#########
load("SOM_table_random_13tp.rda")
DS13_rand_uni_traj <- normalize(Trans_traj_rand_uni_dir_ds13)
DS13_rand_uni <- som(DS13_rand_uni_traj, xdim=3, ydim=1)
plot(DS13_rand_uni, main = "H. Random SNPs (T=13)", ylim=c(-2, 2), cex.axis=2,cex.main =2)


######## Sig SNPs Full########
load("SOM_table_sig_all.rda")
Full_sig_uni_traj <- normalize(Trans_traj_sig_uni_dir)
Full_sig_uni <- som(Full_sig_uni_traj, xdim=3, ydim=1)
plot(Full_sig_uni, main = "E. Candidate SNPs (T=17)", ylim=c(-2, 2), cex.axis=2,cex.main =2)


######## Random SNPs Full########
load("SOM_table_ramdom_all.rda")
Full_rand_uni_traj <- normalize(Trans_traj_rand_uni_dir)
Full_rand_uni <- som(Full_rand_uni_traj, xdim=3, ydim=1)
plot(Full_rand_uni, main = "I. Random SNPs (T=17)", ylim=c(-2, 2), cex.axis=2,cex.main =2)