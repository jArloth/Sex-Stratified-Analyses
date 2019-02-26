#update pheno
load("../ESet_dex_n578.rda")
exp=exprs(ESet_dex)
edata_sort=exp[,order(colnames(exp))]

d1= read.delim("../../Expression/SVA/SV.txt", head=T)
pheno <- ESet_dex@phenoData@data[ sampleNames(ESet_dex)%in% sampleNames(ESet_dex),]
pheno_sv <- merge(pheno, d1, by.x="RNA_ID", by.y="row.names")
#ids = read.delim("../../SampleSheet/sample_sheet_exp_geno_meth.txt",head=T)[,c(1,32)]
#pheno_sv2 <- merge(pheno_sv, ids, by="RNA_ID")
rownames(pheno_sv)=pheno_sv$RNA_ID
pheno=pheno_sv[order(rownames(pheno_sv)),]

if(all(rownames(pheno_sv)==colnames(edata_sort))){
  phenoD=new("AnnotatedDataFrame", data= pheno)#extend to AnnotatedDataFrame (requierd for ESet)
  ESet = new("ExpressionSet", exprs= as.matrix(edata_sort), phenoData=phenoD)
}
ESet_dex=ESet
save(ESet_dex,file="ESet_dex_sva.rda")

########################
test = pData(ESet_dex)
M=cor(test[,c(20:30)],test[,c(4:16)],use="complete" ) #NA

for(i in 1:11){
  tmp=which.max(abs(M[i,]))
  write.table(paste("SV",i,":",names(tmp),":",round(max(abs(M[i,])),3),sep=""),file="05_corr_m_sv.txt",row.names=F,quote=F,col.names=F,sep="\t",append=T)
}

library(corrplot)
pdf("05_cor_M_sv.pdf")
corrplot(M,method = "square")
dev.off()

############################

#calc new residuals add the 2 SNP chips as batch
for(i in 1: dim(ESet_dex)[1]){
  tmp=as.data.frame(cbind(t(exprs(ESet_dex[i,])),ESet_dex@phenoData@data$Dex,as.character(ESet_dex@phenoData@data$DNA_ID)))
  #dex
  dex0=tmp[ tmp$V2==0,]
  colnames(dex0)=c("dex0","base","NID")
  
  dex1=tmp[ tmp$V2==1,]
  colnames(dex1)=c("dex1","stim","NID")
  
  data=merge(dex0,dex1,by.x=c("NID"),by.y=c("NID"))
  delta_new=(as.numeric(as.character(data$dex1))-as.numeric(as.character(data$dex0)))/as.numeric(as.character(data$dex0))
  
  data$delta_new=delta_new
  data =merge(data, ESet_dex@phenoData@data[ ESet_dex@phenoData@data$Dex==0,], by.y="DNA_ID", by.x="NID")
  
  # data$res_new=lm(delta_new~ as.numeric(as.character(age))+as.numeric(as.character(BMI))+Status +as.numeric(as.character(V1))+as.numeric(as.character(V2))+as.numeric(as.character(V3))+
  #                 as.numeric(as.character(V4))+as.numeric(as.character(V5))+as.numeric(as.character(V6))+ as.numeric(as.character(V7))+as.numeric(as.character(V8))+as.numeric(as.character(V9))+
  #                 as.numeric(as.character(V10))+as.numeric(as.character(V11))+as.numeric(as.character(V12))+as.numeric(as.character(V13))+as.numeric(as.character(V14))+as.numeric(as.character(V15))+
  #                 as.numeric(as.character(V16))+as.numeric(as.character(V17))+as.numeric(as.character(V18))+as.numeric(as.character(V19))+as.numeric(as.character(V20))+as.numeric(as.character(V21))+
  #                 as.numeric(as.character(V22))+as.numeric(as.character(V23))+as.numeric(as.character(V24)),
  #                 data=data)$residuals

  data$res_new=lm(delta_new~ as.numeric(as.character(Age))+as.numeric(as.character(BMI_D1))+Status +Sex +as.numeric(as.character(V1))+as.numeric(as.character(V2))+as.numeric(as.character(V3)),
                  data=data)$residuals
  
  sub_data=subset(data,select=c("NID","NID","res_new"))
  colnames(sub_data) <- c("FID","IID",paste(featureNames(ESet_dex)[i],"_residuals",sep=""))
  write.table(sub_data, file=paste("Residuals_SV/",colnames(sub_data)[3],".txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}  

#run 05_02_permutation.sh





#### just for sarah
load("ESet_dex_sva.rda")

#calc new residuals add the 2 SNP chips as batch
for(i in 1: dim(ESet_dex)[1]){
  tmp=as.data.frame(cbind(t(exprs(ESet_dex[i,])),ESet_dex@phenoData@data$Dex,as.character(ESet_dex@phenoData@data$DNA_ID)))
  #dex
  dex0=tmp[ tmp$V2==0,]
  colnames(dex0)=c("dex0","base","NID")
  
  dex1=tmp[ tmp$V2==1,]
  colnames(dex1)=c("dex1","stim","NID")
  
  data=merge(dex0,dex1,by.x=c("NID"),by.y=c("NID"))
  delta_new=(as.numeric(as.character(data$dex1))-as.numeric(as.character(data$dex0)))/as.numeric(as.character(data$dex0))
  
  data$delta_new=delta_new
  data =merge(data, ESet_dex@phenoData@data[ ESet_dex@phenoData@data$Dex==0,], by.y="DNA_ID", by.x="NID")
  
  # data$res_new=lm(delta_new~ as.numeric(as.character(age))+as.numeric(as.character(BMI))+Status +as.numeric(as.character(V1))+as.numeric(as.character(V2))+as.numeric(as.character(V3))+
  #                 as.numeric(as.character(V4))+as.numeric(as.character(V5))+as.numeric(as.character(V6))+ as.numeric(as.character(V7))+as.numeric(as.character(V8))+as.numeric(as.character(V9))+
  #                 as.numeric(as.character(V10))+as.numeric(as.character(V11))+as.numeric(as.character(V12))+as.numeric(as.character(V13))+as.numeric(as.character(V14))+as.numeric(as.character(V15))+
  #                 as.numeric(as.character(V16))+as.numeric(as.character(V17))+as.numeric(as.character(V18))+as.numeric(as.character(V19))+as.numeric(as.character(V20))+as.numeric(as.character(V21))+
  #                 as.numeric(as.character(V22))+as.numeric(as.character(V23))+as.numeric(as.character(V24)),
  #                 data=data)$residuals
  
  data$res_new=lm(delta_new~ as.numeric(as.character(Age))+as.numeric(as.character(BMI_D1))+Status +as.numeric(as.character(V1))+as.numeric(as.character(V2))+as.numeric(as.character(V3)),
                  data=data)$residuals
  
  sub_data=subset(data,select=c("NID","NID","res_new"))
  colnames(sub_data) <- c("FID","IID",paste(featureNames(ESet_dex)[i],"_residuals",sep=""))
  write.table(sub_data, file=paste("Residuals_SV_noSexAdjustment/",colnames(sub_data)[3],".txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}  

#run 05_02_permutation.sh

####
#check
#load("/Users/Janine/Downloads/delta_vals_all_april20.RData")
#test.probe= read.delim("Residuals_SV/ILMN_1755321_residuals.txt", head=T) #no
#test.probe= read.delim("Residuals_SV_noSexAdjustment/ILMN_1755321_residuals.txt", head=T)
