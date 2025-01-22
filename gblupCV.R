#/usr/bin/Rscript


gblupCV <- function(pheno="/scratch/afeesa/GSE_redo/phenoGS",yc=5,ycname="YC2016A",
              h2PBLUP=0.21,
              validGRP="/scratch/afeesa/GSE_redo/within_pop/YC2016A/CV/adg5_CVtypefamily_rep"){
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

#file.copy(from='../VCE/PBLUP/BestPoint',to='.')

######## accuracy and bias #########

accuracy <- NULL
bias <- NULL

for(i in 1:50){
system("rm *.bin")
test_id <- read.table(paste0(validGRP,i,".dat"),header=F,stringsAsFactor=F)
test_dat<- pheno_dat[pheno_dat$animal %in% test_id$V1,]
train_dat <- pheno_dat[!(pheno_dat$animal %in% test_id$V1),]


## write.table(test_dat,"test_data.txt",col.names = F, row.names = F, quote = F)
write.table(train_dat,"train_dat.txt",col.names = F, row.names = F, quote = F)

cat("... starting",i,"cross validation ...","\n")
system("wombat CV.par")
solutions <- fread("RnSoln_animal.dat",header=F,skip=1)

validation <- merge(test_dat,solutions,by.x="animal",by.y="V2")
accuracy[i]<- with(validation,cor(V4,AGD))/sqrt(h2PBLUP)

#bias
reg <- lm(validation$AGD ~ validation$V4)
bias[i]<- reg$coefficients[2]

cat("... Done with",i,"cross validation ...","\n")
}
result <- cbind(accuracy,bias)
cat("The mean accuracy is", mean(result[,1]),"with standard error", sd(result[,1]),"\n")
cat("The mean bias is ", mean(result[,2]), "with standard error", sd(result[,2]),"\n")

write.table(result,"result_accuracy",col.names=T,row.names=F)
}
