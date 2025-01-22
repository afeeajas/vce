#/usr/bin/Rscript


pedVCE <- function(pheno="/scratch/afeesa/GSE_redo/phenoGS",yc=5){
library(data.table)

pheno_dat <- read.table(pheno,header=T,stringsAsFactor=F)
pheno_dat <- pheno_dat[pheno_dat$trait_no ==yc,]
write.table(pheno_dat,"phen.dat",col.names=F,row.names=F,quote=F)
system("wombat")

}

