#/usr/bin/Rscript


VCEest <- function(pheno="/scratch/afeesa/GSE_redo/phenoGS",yc=5,ycname="YC2016A",accpredpath="/scratch/afeesa/GSE_redo/within_pop/YC2016A/accPred"){
library(data.table)
pheno_dat<- read.table(pheno,header=T,stringsAsFactor=F)
pheno_dat <- pheno_dat[pheno_dat$trait_no==yc,]
write.table(pheno_dat[,c(2,2)],paste0(ycname),col.names=F,row.names=F,quote=F)


system(paste0("plink --bfile /scratch/afeesa/YC2017/GENO_EATFISH --keep ", ycname ," --nonfounders --chr-set 30 --recodeA --out ", ycname)) #This is being edited with selected SNP 
genofile = paste0(ycname,".raw")
geno <- read.table(genofile,header=T,stringsAsFactor=F)
pheno_dat <- pheno_dat[pheno_dat$animal %in% geno$IID,]
write.table(pheno_dat,"phen.dat",col.names=F,row.names=F,quote=F)
geno <- geno[geno$IID %in% pheno_dat$animal,]
pedfile <- pheno_dat[,c(2:4)]
pedfile[,c(2,3)] <- 0
pedfile$gtype <- 2
write.table(pedfile,paste0(ycname,"_ped"),col.names=F,row.names=F,quote=F)
write.table(geno[,-c(1,3,4:6)],"MarkerCounts.dat",col.names=F,row.names=F,quote=F)

nmrk <- ncol(geno)-6
grmpar <- c("RUNOP --hinv -v","COM Build G-inverse with default settings",
paste0("PED ",ycname,"_ped"),"MRK MarkerCounts.dat","SPECIAL",paste0("  HINVERSE SNP ",nmrk),"END")
write.table(grmpar,"grm.par",col.names=F,row.names=F,quote=F)

system("wombat grm.par")
system("mv Hinverse.gin animal.gin")
system("mv Hinverse.codes animal.codes")
system("wombat")
cat("Variance component estimated with selected SNPS")
#system(paste0("mv animal.gin ",accpredpath,"/animal.gin"))
#system(paste0("mv animal.codes ",accpredpath,"/animal.codes"))
system("rm *.raw *.bin *.gin MarkerCounts.dat")
}
