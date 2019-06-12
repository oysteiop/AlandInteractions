##########################################################################
######################### - Åland Interactions - #########################
##########################################################################

###########################################
#### - Single-species col/ext models - ####
###########################################
rm(list=ls())
setwd("Z:/data/aland/")

library(reshape2)
library(mvtnorm)  
library(lme4)
library(MuMIn)

#Read datafiles
list.files()
newdat=read.csv("colextdat.csv")

#Insert patch mean for missing PL data
pmeanPL=data.frame(tapply(newdat$PLM2,newdat$Patch,mean,na.rm=T))
pmeanPL$Patch=row.names(pmeanPL)
names(pmeanPL)=c("PL","Patch")
newPL=subset(newdat,select=c("Patch","PLM2"))
for(i in 1:nrow(newPL)){
  if(is.na(newPL$PLM2[i])){
    newPL$PLM2[i]=pmeanPL$PL[which(pmeanPL$Patch==newPL$Patch[i])]
  }
}
newdat$PLM2=newPL$PLM2

#Insert patch mean for missing VS data
newdat$VS[which(newdat$VS>3)]=3

pmeanVS=data.frame(tapply(newdat$VS,newdat$Patch,mean,na.rm=T))
pmeanVS$Patch=row.names(pmeanVS)
names(pmeanVS)=c("VS","Patch")
newVS=subset(newdat,select=c("Patch","VS"))
for(i in 1:nrow(newVS)){
  if(is.na(newVS$VS[i])){
    newVS$VS[i]=pmeanVS$VS[which(pmeanVS$Patch==newVS$Patch[i])]
  }
}
newdat$VS=newVS$VS

patches=read.csv("patch_data.csv")
nhdata=read.csv("NEWneighbourhooddata.csv")
Cnhdata=read.csv("NEWCotesia_neighbourhooddata.csv")
Mnhdata=read.csv("Mildew_neighbourhooddata.csv")
#Mnhdata=read.csv("Mildew_neighbourhooddata_PatchArea.csv")
TotalNH=read.csv("TotalNH.csv")

#### - Format reduced datafile for model fitting - ####
mdat=subset(newdat,select=c("Patch","Year","PatchYear","PLM2","VSM2","PL","VS","PLVS","MelCc","MelCe","PhoPc","PhoPe","CotPc","CotPe","MelClast","PhoPlast","CotPlast","Statelast"))
mdat$Patch=as.factor(mdat$Patch)
r=rowSums(mdat[,9:14]>-1,na.rm=T)==0
mdat=mdat[-which(r),] #Remove all-NA rows

#Environmental data
#list.files()
pdat=read.csv("patch_precipitation.csv")

#May prec
mayprec=pdat[,c(1,25,35,45,55,65,75,85,95,105,133,144,154,164,174,184,194,204,214)]
colnames(mayprec)=c("Patch",paste(2001:2018))
mayprec=melt(mayprec,id.var="Patch")
mayprec$PatchYear=paste(mayprec$Patch,mayprec$variable,sep="_")
colnames(mayprec)=c("Patch","Year","MayPrec","PatchYear")
mayprec=mayprec[mayprec$MayPrec>0,]

#June prec
juneprec=pdat[,c(1,26,36,46,56,66,76,86,96,106,134,145,155,165,175,185,195,205,215)]
colnames(juneprec)=c("Patch",paste(2001:2018))
juneprec=melt(juneprec,id.var="Patch")
juneprec$PatchYear=paste(juneprec$Patch,juneprec$variable,sep="_")
colnames(juneprec)=c("Patch","Year","junePrec","PatchYear")
juneprec=juneprec[juneprec$junePrec>0,]

#July prec
julyprec=pdat[,c(1,27,37,47,57,67,77,87,97,107,135,146,156,166,176,186,196,206,216)]
colnames(julyprec)=c("Patch",paste(2001:2018))
julyprec=melt(julyprec,id.var="Patch")
julyprec$PatchYear=paste(julyprec$Patch,julyprec$variable,sep="_")
colnames(julyprec)=c("Patch","Year","julyPrec","PatchYear")
julyprec=julyprec[julyprec$julyPrec>0,]

#August prec
augustprec=pdat[,c(1,28,38,48,58,68,78,88,98,108,137,147,157,167,177,187,197,207,217)]
colnames(augustprec)=c("Patch",paste(2001:2018))
augustprec=melt(augustprec,id.var="Patch")
augustprec$PatchYear=paste(augustprec$Patch,augustprec$variable,sep="_")
colnames(augustprec)=c("Patch","Year","augustPrec","PatchYear")
augustprec=augustprec[augustprec$augustPrec>0,]

#May temp #MISSING MAY 2015 DATA
#list.files()
#tdat=read.csv("patch_temp.csv")
#cbind(names(tdat))
#maytemp=tdat[,c(1,25,37,49,61,72,84,96,108,120,131,143,155,167,179,NA, 202,214]
#head(maytemp)
#colnames(maytemp)=c("Patch",paste(2001:2017))
#maytemp=melt(maytemp,id.var="Patch")
#maytemp$PatchYear=paste(maytemp$Patch,maytemp$variable,sep="_")
#colnames(maytemp)=c("Patch","Year","MayTemp","PatchYear")
#head(maytemp)
#maytemp[which(maytemp$MayTemp<0),3]=NA

#Combine
comb=merge(mdat,mayprec[,3:4],by="PatchYear",all.x=T)
comb=merge(comb,juneprec[,3:4],by="PatchYear",all.x=T)
comb=merge(comb,julyprec[,3:4],by="PatchYear",all.x=T)
comb=merge(comb,augustprec[,3:4],by="PatchYear",all.x=T)

head(patches)
colnames(patches)=c("Patch","Lat","Lon","Area","Road")
comb=merge(comb,patches[,c(1,4,5)],by="Patch",all.x=T)

#Connectivity data
colnames(nhdata)=c("Patch",paste(2001:2019)) #Offset 1 year because these values were computed with Nt rather than Nt-i
nhood=melt(nhdata,id.var="Patch")
nhood$PatchYear=paste(nhood$Patch,nhood$variable,sep="_")
colnames(nhood)=c("Patch","Year","NH","PatchYear")
comb=merge(comb,nhood[,3:4],by="PatchYear",all.x=T)

colnames(Cnhdata)=c("Patch",paste(2001:2019)) #Same
Cnhood=melt(Cnhdata,id.var="Patch")
Cnhood$PatchYear=paste(Cnhood$Patch,Cnhood$variable,sep="_")
colnames(Cnhood)=c("Patch","Year","NH_Cot","PatchYear")
comb=merge(comb,Cnhood[,3:4],by="PatchYear",all.x=T)

colnames(Mnhdata)=c("Patch",paste(2001:2019)) #Same
Mnhood=melt(Mnhdata,id.var="Patch")
Mnhood$PatchYear=paste(Mnhood$Patch,Mnhood$variable,sep="_")
colnames(Mnhood)=c("Patch","Year","NH_Mil","PatchYear")
comb=merge(comb,Mnhood[,3:4],by="PatchYear",all.x=T)

colnames(TotalNH)
comb=merge(comb,TotalNH,by="Patch",all.x=T)

#Variable transformations
comb$logMayPrec=log(comb$MayPrec)
comb$logJunePrec=log(comb$junePrec)
comb$logJulyPrec=log(comb$julyPrec)
comb$logAugPrec=log(comb$augustPrec)

comb$logPLM2=log(comb$PLM2+1)
comb$logArea=log(comb$Area)
comb$logNH=log(comb$NH)
comb$logNH_Mil=log(comb$NH_Mil)
comb$logNH_Cot=sqrt(comb$NH_Cot)
comb$logTotalNH=log(comb$TotalNH)

comb$Patch=as.factor(comb$Patch)
comb$Year=as.factor(comb$Year)

comb$MelClast=as.numeric(comb$MelClast>0)
comb$PhoPlast=as.numeric(comb$PhoPlast>0)
comb$CotPlast=as.numeric(comb$CotPlast>0)
comb2=comb[which(comb$Statelast %in% c("0_NA_0", "1_0_0", "1_1_0", "1_1_1","1_0_1","0_NA_1")),]
comb2$Statelast=factor(comb2$Statelast)
str(comb)
summary(comb)
#write.csv(comb,"combdat.csv",row.names=F)

names(comb)
apply(comb[,9:14],2, sum, na.rm=T)
apply(comb[,9:14]>-Inf,2, sum, na.rm=T)
signif(apply(comb[,9:14],2, sum, na.rm=T)/apply(comb[,9:14]>-Inf,2, sum, na.rm=T),2)

d=comb[which(comb$CotPe>0),]
table(d$MelCe)

d=comb[which(comb$CotPc>0),]
table(d$MelCc)

##########################################
#### - Fitting the candidate models - ####
##########################################

#### - Butterfly colonisation models - ####
MelCc_dat=na.omit(subset(comb,select=c("MelCc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","VS","logNH","Road","Year","Patch")))

mMelCc=glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc2=glmer(MelCc~logPLM2+VS+logNH+Road+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc3=glmer(MelCc~logPLM2+VS+logNH+Road+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc4=glmer(MelCc~logPLM2+VS+logNH+Road+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc5=glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc6=glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc7=glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc8=glmer(MelCc~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc9=glmer(MelCc~logPLM2+VS+logNH+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc10=glmer(MelCc~logPLM2+VS+logNH+Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc11=glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc12=glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc13=glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc14=glmer(MelCc~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)
mMelCc15=glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCc_dat)

MelCc_modList=list(mMelCc,mMelCc2,mMelCc3,mMelCc4,mMelCc5,mMelCc6,mMelCc7,mMelCc8,mMelCc9,mMelCc10,mMelCc11,mMelCc12,mMelCc13,mMelCc14,mMelCc15)
names(MelCc_modList)=c("mMelCc",paste0("mMelCc",2:15))
save(MelCc_modList,file="fitted_mods/MelCc_modList.RData")

#### - Butterfly extinction models - ####
MelCe_dat=na.omit(subset(comb,select=c("MelCe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","VS","logNH","Road","Year","Patch")))

mMelCe=glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe2=glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe3=glmer(MelCe~logPLM2+VS+logNH+Road+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe4=glmer(MelCe~logPLM2+VS+logNH+Road+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe5=glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe6=glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe7=glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe8=glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe9=glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe10=glmer(MelCe~logPLM2+VS+logNH+Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe11=glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe12=glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe13=glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe14=glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe15=glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)

MelCe_modList=list(mMelCe,mMelCe2,mMelCe3,mMelCe4,mMelCe5,mMelCe6,mMelCe7,mMelCe8,mMelCe9,mMelCe10,mMelCe11,mMelCe12,mMelCe13,mMelCe14,mMelCe15)
names(MelCe_modList)=c("mMelCe",paste0("mMelCe",2:15))
save(MelCe_modList,file="fitted_mods/MelCe_modList.RData")

#### - Mildew colonisation models - ####
PhoPc_dat=na.omit(subset(comb,select=c("PhoPc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","logNH_Mil","logTotalNH","Road","Year","Patch")))

mPhoPc=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logMayPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc2=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc3=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc4=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc5=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc6=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc7=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc8=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc9=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc10=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc11=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc12=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc13=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc14=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc15=glmer(PhoPc~logPLM2+logNH_Mil + logTotalNH +Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)

PhoPc_modList=list(mPhoPc,mPhoPc2,mPhoPc3,mPhoPc4,mPhoPc5,mPhoPc6,mPhoPc7,mPhoPc8,mPhoPc9,mPhoPc10,mPhoPc11,mPhoPc12,mPhoPc13,mPhoPc14,mPhoPc15)
names(PhoPc_modList)=c("mPhoPc",paste0("mPhoPc",2:15))
save(PhoPc_modList,file="fitted_mods/PhoPc_modList.RData")

#### - Mildew extinction models - ####
PhoPe_dat=na.omit(subset(comb,select=c("PhoPe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","logNH_Mil","logTotalNH","Road","Year","Patch")))

mPhoPe=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logMayPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe2=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe3=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH +Road+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe4=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe5=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe6=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe7=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe8=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe9=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe10=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe11=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe12=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe13=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe14=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe15=glmer(PhoPe~logPLM2+logNH_Mil + logTotalNH+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)

PhoPe_modList=list(mPhoPe,mPhoPe2,mPhoPe3,mPhoPe4,mPhoPe5,mPhoPe6,mPhoPe7,mPhoPe8,mPhoPe9,mPhoPe10,mPhoPe11,mPhoPe12,mPhoPe13,mPhoPe14,mPhoPe15)
names(PhoPe_modList)=c("mPhoPe",paste0("mPhoPe",2:15))
save(PhoPe_modList,file="fitted_mods/PhoPe_modList.RData")

#### - Parasitoid colonisation models - ####
CotPc_dat=na.omit(subset(comb,select=c("CotPc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPc_dat=CotPc_dat[CotPc_dat$Year!="2018",]

mCotPc=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc2=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc3=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc4=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc5=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc6=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc7=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc8=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc9=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc10=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc11=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc12=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc13=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc14=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc15=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPc_dat)

CotPc_modList=list(mCotPc,mCotPc2,mCotPc3,mCotPc4,mCotPc5,mCotPc6,mCotPc7,mCotPc8,mCotPc9,mCotPc10,mCotPc11,mCotPc12,mCotPc13,mCotPc14,mCotPc15)
names(CotPc_modList)=c("mCotPc",paste0("mCotPc",2:15))
save(CotPc_modList,file="fitted_mods/CotPc_modList.RData")

#### - Parasitoid extinction models - ####
CotPe_dat=na.omit(subset(comb,select=c("CotPe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPe_dat=CotPe_dat[CotPe_dat$Year!="2018",]

mCotPe=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe2=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe3=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe4=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe5=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe6=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe7=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe8=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe9=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe10=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe11=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe12=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe13=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe14=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe15=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)

CotPe_modList=list(mCotPe,mCotPe2,mCotPe3,mCotPe4,mCotPe5,mCotPe6,mCotPe7,mCotPe8,mCotPe9,mCotPe10,mCotPe11,mCotPe12,mCotPe13,mCotPe14,mCotPe15)
names(CotPe_modList)=c("mCotPe",paste0("mCotPe",2:15))
save(CotPe_modList,file="fitted_mods/CotPe_modList.RData")

############################################################
#### - Parameter estimates from highest ranked models - ####
############################################################

#Butterfly col
MelCc_dat=na.omit(subset(comb,select=c("MelCc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH","logTotalNH","Road","Year","Patch")))
load("fitted_mods/MelCc_modList.RData")

AICtab=data.frame(AIC=unlist(lapply(MelCc_modList,AIC)))
AICtab$delta=AICtab$AIC-min(AICtab$AIC)
AICtab=AICtab[order(AICtab$delta),]
AICtab #Mod1 best, 11 good

mMelCc=MelCc_modList[["mMelCc"]]
summary(mMelCc)
r2=r.squaredGLMM(mMelCc)
r2

#Varcomps
mod=mMelCc
dat=MelCc_dat

xdat=as.matrix(data.frame(Intercept=rep(1,nrow(dat)), subset(dat, select = rownames(summary(mod)$coef)[-1])))
fixed=t(summary(mod)$coef[,1])%*%cov(xdat)%*%summary(mod)$coef[,1]
yr=as.numeric(VarCorr(mod)$Year)
pa=as.numeric(VarCorr(mod)$Patch)
total=fixed+pa+yr
df=cbind(fixed,yr,pa)
round(df/sum(df)*100,2)

pl=(summary(mod)$coef[2,1]^2)*var(dat$logPLM2)
vs=(summary(mod)$coef[3,1]^2)*var(dat$VS)
round(pl/total*100, 2)
round(vs/total*100, 2)

#Butterfly ext
MelCe_dat=na.omit(subset(comb,select=c("MelCe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH","logTotalNH","Road","Year","Patch")))
load("fitted_mods/MelCe_modList.RData")

AICtab=data.frame(AIC=unlist(lapply(MelCe_modList,AIC)))
AICtab$delta=AICtab$AIC-min(AICtab$AIC)
AICtab=AICtab[order(AICtab$delta),]
AICtab #Mod10 best but poor convergence, use mod14

mMelCe14=MelCe_modList[["mMelCe14"]]
summary(mMelCe14)
r2=r.squaredGLMM(mMelCe14)
r2

#Mildew col
PhoPc_dat=na.omit(subset(comb,select=c("PhoPc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Mil","logTotalNH","Road","Year","Patch")))
load("fitted_mods/PhoPc_modList.RData")

AICtab=data.frame(AIC=unlist(lapply(PhoPc_modList,AIC)))
AICtab$delta=AICtab$AIC-min(AICtab$AIC)
AICtab=AICtab[order(AICtab$delta),]
AICtab #Mod9 best

mPhoPc9=PhoPc_modList[["mPhoPc9"]]
summary(mPhoPc9)
r2=r.squaredGLMM(mPhoPc9)
r2

#Varcomps
mod=mPhoPc9
dat=PhoPc_dat

xdat=as.matrix(data.frame(Intercept=rep(1,nrow(dat)), subset(dat, select = rownames(summary(mod)$coef)[-1])))
fixed=t(summary(mod)$coef[,1])%*%cov(xdat)%*%summary(mod)$coef[,1]
yr=as.numeric(VarCorr(mod)$Year)
pa=as.numeric(VarCorr(mod)$Patch)
total=fixed+pa+yr
df=cbind(fixed,yr,pa)
round(df/sum(df)*100,2)

pl=(summary(mod)$coef[2,1]^2)*var(dat$logPLM2)
nhMil=(summary(mod)$coef[3,1]^2)*var(dat$logNH_Mil)
nh=(summary(mod)$coef[4,1]^2)*var(dat$logTotalNH)
june=(summary(mod)$coef[6,1]^2)*var(dat$logJunePrec)
prec=t((summary(mod)$coef[6:7,1])) %*% cov(subset(dat, select=c("logJunePrec", "logAugPrec"))) %*% (summary(mod)$coef[6:7,1])

round(pl/total*100, 2)
round(nhMil/total*100, 2)
round(nh/total*100, 2)
round(prec/total*100, 2)

#Mildew extinction
PhoPe_dat=na.omit(subset(comb,select=c("PhoPe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Mil","Road","Year","Patch")))
load("fitted_mods/PhoPe_modList.RData")

AICtab=data.frame(AIC=unlist(lapply(PhoPe_modList,AIC)))
AICtab$delta=AICtab$AIC-min(AICtab$AIC)
AICtab=AICtab[order(AICtab$delta),]
AICtab #Mod3 best

mPhoPe3=PhoPe_modList[["mPhoPe3"]]
summary(mPhoPe3)
r2=r.squaredGLMM(mPhoPe3)
r2

#Parasitoid col
CotPc_dat=na.omit(subset(comb,select=c("CotPc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPc_dat=CotPc_dat[CotPc_dat$Year!="2018",]
load("fitted_mods/CotPc_modList.RData")

AICtab=data.frame(AIC=unlist(lapply(CotPc_modList,AIC)))
AICtab$delta=AICtab$AIC-min(AICtab$AIC)
AICtab=AICtab[order(AICtab$delta),]
AICtab #Mod 2 best, Mod 1 and Mod 4 equal

mCotPc2=CotPc_modList[["mCotPc2"]]
summary(mCotPc2)
r2=r.squaredGLMM(mCotPc2)
r2

#Varcomps
mod=mCotPc2
dat=CotPc_dat

xdat=as.matrix(data.frame(Intercept=rep(1,nrow(dat)), subset(dat, select = rownames(summary(mod)$coef)[-1])))
fixed=t(summary(mod)$coef[,1])%*%cov(xdat)%*%summary(mod)$coef[,1]
yr=as.numeric(VarCorr(mod)$Year)
pa=as.numeric(VarCorr(mod)$Patch)
total=fixed+pa+yr
df=cbind(fixed,yr,pa)
round(df/sum(df)*100,2)

pl=(summary(mod)$coef[2,1]^2)*var(dat$logPLM2)
vs=(summary(mod)$coef[3,1]^2)*var(dat$VS)
round(pl/total*100, 1)
round(vs/total*100, 1)

#Parasitoid ext
CotPe_dat=na.omit(subset(comb,select=c("CotPe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPe_dat=CotPe_dat[CotPe_dat$Year!="2018",]
load("fitted_mods/CotPe_modList.RData")

AICtab=data.frame(AIC=unlist(lapply(CotPe_modList,AIC)))
AICtab$delta=AICtab$AIC-min(AICtab$AIC)
AICtab=AICtab[order(AICtab$delta),]
AICtab #Mod2 best

mCotPe2=CotPe_modList[["mCotPe2"]]
summary(mCotPe2)
r2=r.squaredGLMM(mCotPe2)
r2

#############################################################
#### - Adding patch state to the highest ranked models - ####
#############################################################

#Butterfly col
MelCc_dat=na.omit(subset(comb2,select=c("MelCc","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH","logTotalNH","Road","Year","Patch")))

mMelCca=glmer(MelCc~logPLM2 + VS + logNH + Road + logMayPrec + (1|Patch) + (1|Year), family="binomial", data=MelCc_dat)
mMelCcb=glmer(MelCc~logPLM2 + VS + logNH + Road + logMayPrec + Statelast + (1|Patch) + (1|Year), family="binomial", data=MelCc_dat)
AIC(mMelCca, mMelCcb)
AIC(mMelCca, mMelCcb)$AIC[2]-AIC(mMelCca, mMelCcb)$AIC[1]

mMelCcc=glmer(MelCc~scale(logPLM2, scale=F) + scale(VS, scale=F) + scale(logNH, scale=F) + scale(Road, scale=F) +
                scale(logMayPrec, scale=F) + Statelast - 1 + (1|Patch) + (1|Year), family="binomial", data = MelCc_dat)

summary(mMelCcc)

#000
s=summary(mMelCcc)$coef[6,1]
sse=summary(mMelCcc)$coef[6,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#001
s=summary(mMelCcc)$coef[7,1]
sse=summary(mMelCcc)$coef[7,2]
round(1/(1+exp(-s))*100,2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#Butterfly ext
MelCe_dat=na.omit(subset(comb2,select=c("MelCe","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH","logTotalNH","Road","Year","Patch")))

mMelCe14a=glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe14b=glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+logAugPrec+Statelast + (1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
AIC(mMelCe14a,mMelCe14b)
AIC(mMelCe14a,mMelCe14b)$AIC[2]-AIC(mMelCe14a,mMelCe14b)$AIC[1]

mMelCe14c=glmer(MelCe~scale(logPLM2, scale=F) + scale(VS, scale=F) + scale(logNH, scale=F)+
                 scale(Road, scale=F) + scale(logJunePrec, scale=F) +scale(logJulyPrec, scale=F) + 
                 scale(logAugPrec, scale=F) + Statelast -1 + (1|Patch)+(1|Year),family="binomial",data=MelCe_dat)

summary(mMelCe14c)

#100
s=summary(mMelCe14c)$coef[8,1]
sse=summary(mMelCe14c)$coef[8,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#101
s=summary(mMelCe14c)$coef[9,1]
sse=summary(mMelCe14c)$coef[9,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#110
s=summary(mMelCe14c)$coef[10,1]
sse=summary(mMelCe14c)$coef[10,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#111
s=summary(mMelCe14c)$coef[11,1]
sse=summary(mMelCe14c)$coef[11,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#Mildew col
PhoPc_dat=na.omit(subset(comb2,select=c("PhoPc","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Mil","logTotalNH","Road","Year","Patch")))

mPhoPc9a=glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc9b=glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+logAugPrec+Statelast + (1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
AIC(mPhoPc9a,mPhoPc9b)
AIC(mPhoPc9a,mPhoPc9b)$AIC[2]-AIC(mPhoPc9a,mPhoPc9b)$AIC[1]

mPhoPc9c=glmer(PhoPc~scale(logPLM2, scale=F) + scale(logNH_Mil, scale=F) + scale(logTotalNH, scale=F) + scale(Road, scale=F)+
                  scale(logJunePrec, scale=F) + scale(logAugPrec, scale=F) + Statelast -1 + (1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)

summary(mPhoPc9c)

#000
s=summary(mPhoPc9c)$coef[7,1]
sse=summary(mPhoPc9c)$coef[7,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#100
s=summary(mPhoPc9c)$coef[8,1]
sse=summary(mPhoPc9c)$coef[8,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#110
s=summary(mPhoPc9c)$coef[9,1]
sse=summary(mPhoPc9c)$coef[9,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#Mildew ext
PhoPe_dat=na.omit(subset(comb2,select=c("PhoPe","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Mil","logTotalNH","Road","Year","Patch")))

mPhoPe3a=glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logJulyPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
mPhoPe3b=glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logJulyPrec+ Statelast + (1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)
AIC(mPhoPe3a,mPhoPe3b)
AIC(mPhoPe3a,mPhoPe3b)$AIC[2]-AIC(mPhoPe3a,mPhoPe3b)$AIC[1]

mPhoPe3c=glmer(PhoPe~scale(logPLM2, scale=F) + scale(logNH_Mil, scale=F) + scale(logTotalNH, scale=F) + scale(Road, scale=F) + 
                 scale(logJulyPrec, scale=F) + Statelast -1 + (1|Patch)+(1|Year),family="binomial",data=PhoPe_dat)

summary(mPhoPe3c)

#001
s=summary(mPhoPe3c)$coef[6,1]
sse=summary(mPhoPe3c)$coef[6,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#101
s=summary(mPhoPe3c)$coef[7,1]
sse=summary(mPhoPe3c)$coef[7,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#111
s=summary(mPhoPe3c)$coef[8,1]
sse=summary(mPhoPe3c)$coef[8,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#Parasitoid col
CotPc_dat=na.omit(subset(comb2,select=c("CotPc","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPc_dat=CotPc_dat[CotPc_dat$Year!="2018",]

mCotPc2a=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec + (1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
mCotPc2b=glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec+Statelast + (1|Patch)+(1|Year),family="binomial",data=CotPc_dat)
AIC(mCotPc2a, mCotPc2b)
AIC(mCotPc2a, mCotPc2b)$AIC[2]-AIC(mCotPc2a, mCotPc2b)$AIC[1]

mCotPc2c=glmer(CotPc~scale(logPLM2, scale=F) + scale(VS, scale=F) + scale(logNH_Cot, scale=F) + scale(Road, scale=F) + 
                scale(logJunePrec, scale=F) + Statelast - 1 + (1|Patch)+(1|Year),family="binomial",data=CotPc_dat)

summary(mCotPc2c)

#000
s=summary(mCotPc2c)$coef[6,1]
sse=summary(mCotPc2c)$coef[6,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#001
s=summary(mCotPc2c)$coef[7,1]
sse=summary(mCotPc2c)$coef[7,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#100
s=summary(mCotPc2c)$coef[8,1]
sse=summary(mCotPc2c)$coef[8,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#101
s=summary(mCotPc2c)$coef[9,1]
sse=summary(mCotPc2c)$coef[9,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#Parasitoid ext
CotPe_dat=na.omit(subset(comb2,select=c("CotPe","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPe_dat=CotPe_dat[CotPe_dat$Year!="2018",]

mCotPe2a=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
mCotPe2b=glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+Statelast+(1|Patch)+(1|Year),family="binomial",data=CotPe_dat)
AIC(mCotPe2a,mCotPe2b)
AIC(mCotPe2a,mCotPe2b)$AIC[2]-AIC(mCotPe2a,mCotPe2b)$AIC[1]

mCotPe2c=glmer(CotPe~scale(logPLM2, scale=F) + scale(VS, scale=F) + scale(logNH_Cot, scale=F)+
                 scale(Road, scale=F) + scale(logJunePrec, scale=F) +Statelast -1 + (1|Patch)+(1|Year),family="binomial",data=CotPe_dat)

summary(mCotPe2c)

#110
s=summary(mCotPe2c)$coef[6,1]
sse=summary(mCotPe2c)$coef[6,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#111
s=summary(mCotPe2c)$coef[7,1]
sse=summary(mCotPe2c)$coef[7,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

######################################################
#### - Generating predicted population dynamics - ####
######################################################

invlogit=function(x){1/(1+exp(-x))}

#Butterfly
load("fitted_mods/MelCc_modList.RData")
mMelCc=MelCc_modList[["mMelCc"]]
load("fitted_mods/MelCe_modList.RData")
mMelCe14=MelCe_modList[["mMelCe14"]]

TotalPatchDat=read.csv("TotalPatchDat.csv")

br=quantile(newdat$PLM2,c(seq(0,1,length.out=15)),na.rm=T)
test=cut(newdat$PLM2,breaks=br)
obs=tapply(newdat$MelC>0,test,sum,na.rm=T)/tapply(newdat$MelC>-1,test,sum,na.rm=T)

#Setup
years=unique(newdat$Year)
predN=matrix(NA,nrow=length(years)-1,ncol=100)
predP=matrix(NA,nrow=length(years)-1,ncol=100)
alldat=list()
cmod=mMelCc
emod=mMelCe14

for(y in 1:(length(years)-1)){
  year=years[y+1]
  print(paste("Year",year))
  
  lastyear=subset(newdat[newdat$Year==year-1,],select=c("Patch","MelC"))
  lastyear$MelC=as.numeric(lastyear$MelC>0)
  thisyear=comb[comb$Year==year,]
  
  preddat=merge(lastyear,thisyear,by="Patch",all.x=T)
  preddat$PLcat=cut(preddat$PLM2,breaks=br) 
  outdat=data.frame(Year=rep(paste(year),nrow(preddat)),Patch=preddat$Patch,PLcat=preddat$PLcat,Pred=rep(NA,nrow(preddat)))
  
  for(j in 1:100){
    
    #Col
    rr=MASS::mvrnorm(1,mu=summary(cmod)$coef[,1],Sigma=vcov(cmod))
    rYear=ranef(cmod)$Year[which(row.names(ranef(cmod)$Year)==year),1]
    #rYear=0
    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(rr[-1])))
    rr[1]=rr[1]+rYear
    predc=invlogit(as.matrix(pd)%*%as.matrix(rr))
  
    #Ext
    rr=MASS::mvrnorm(1,mu=summary(emod)$coef[,1],Sigma=vcov(emod))
    rYear=ranef(emod)$Year[which(row.names(ranef(emod)$Year)==year),1]
    #rYear=0
    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(rr[-1])))
    rr[1]=rr[1]+rYear
    prede=invlogit(as.matrix(pd)%*%as.matrix(rr))
  
  preddat$pred=rep(NA,nrow(preddat))
  preddat$pred[which(preddat$MelC==1)]=(1-prede[which(preddat$MelC==1)])
  preddat$pred[which(preddat$MelC==0)]=predc[which(preddat$MelC==0)]
  
  preddat$predBin=rbinom(nrow(preddat),1,p=preddat$pred)
  predN[y,j]=sum(preddat$predBin,na.rm=T)
  predP[y,j]=sum(preddat$predBin,na.rm=T)/(sum(preddat$predBin>-1,na.rm=T))
  
  outdat[,j+3]=preddat$predBin
  }
  alldat[[y]]=outdat
  }

allCinxiaDat=alldat

df=data.frame(Year=years[1:18]+1,predN=apply(predN,1,mean),
              predNupper=apply(predN,1,quantile,c(.025,.975),na.rm=T)[2,],
              predNlower=apply(predN,1,quantile,c(.025,.975),na.rm=T)[1,],
              predP=apply(predP,1,mean),
              predPupper=apply(predP,1,quantile,c(.025,.975),na.rm=T)[2,],
              predPlower=apply(predP,1,quantile,c(.025,.975),na.rm=T)[1,])

colr=rgb(0,0,1, alpha=.5)

yearlyDat=read.csv("yearlyDat.csv")
x11()
par(mfrow=c(1,3))
plot(years,yearlyDat$MelC,type="b",pch=16,ylim=c(0,.5),xlab="",las=1,ylab="Proportion of patches occupied",main="Melitaea cinxia")
xx=2001:2018
polygon(c(xx,rev(xx)),c(df$predPlower[1:18],rev(df$predPupper[1:18])),col = colr, border = FALSE)
points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years,yearlyDat$MelC,type="b",pch=16)

#Mildew
test=cut(newdat$PLM2,breaks=br)
obs=tapply(newdat$PhoP>0,test,sum,na.rm=T)/tapply(newdat$PhoP>-1,test,sum,na.rm=T)

load("fitted_mods/PhoPc_modList.RData")
mPhoPc9=PhoPc_modList[["mPhoPc9"]]
load("fitted_mods/PhoPe_modList.RData")
mPhoPe3=PhoPe_modList[["mPhoPe3"]]

#Parameter uncertainty
years=unique(newdat$Year)
predN=matrix(NA,nrow=length(years)-1,ncol=100)
predP=matrix(NA,nrow=length(years)-1,ncol=100)
alldat=list()
cmod=mPhoPc9
emod=mPhoPe3

a=Sys.time()
for(y in 1:(length(years)-1)){
  year=years[y+1]
  print(paste("Year",year))
  lastyear=subset(newdat[newdat$Year==year-1,],select=c("Patch","PhoP"))
  lastyear$PhoP=as.numeric(lastyear$PhoP>0)
  thisyear=comb[comb$Year==year,]
  
  preddat=merge(lastyear,thisyear,by="Patch",all.x=T)
  preddat$PLcat=cut(preddat$PLM2,breaks=br) 
  outdat=data.frame(Year=rep(paste(year),nrow(preddat)),Patch=preddat$Patch,PLcat=preddat$PLcat,Pred=rep(NA,nrow(preddat)))
  
  
  for(j in 1:100){
    
    #Col
    rr=MASS::mvrnorm(1,mu=summary(cmod)$coef[,1],Sigma=vcov(cmod))
    rYear=ranef(cmod)$Year[which(row.names(ranef(cmod)$Year)==year),1]
    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(rr[-1])))
    rr[1]=rr[1]+rYear
    predc=invlogit(as.matrix(pd)%*%as.matrix(rr))
    
    #Ext
    rr=MASS::mvrnorm(1,mu=summary(emod)$coef[,1],Sigma=vcov(emod))
    rYear=ranef(emod)$Year[which(row.names(ranef(emod)$Year)==year),1]
    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(rr[-1])))
    rr[1]=rr[1]+rYear
    prede=invlogit(as.matrix(pd)%*%as.matrix(rr))
    
    preddat$pred=rep(NA,nrow(preddat))
    preddat$pred[which(preddat$PhoP==1)]=(1-prede[which(preddat$PhoP==1)])
    preddat$pred[which(preddat$PhoP==0)]=predc[which(preddat$PhoP==0)]
    
    preddat$predBin=rbinom(nrow(preddat),1,p=preddat$pred)
    predN[y,j]=sum(preddat$predBin,na.rm=T)
    predP[y,j]=sum(preddat$predBin,na.rm=T)/(sum(preddat$predBin>-1,na.rm=T))
    
    outdat[,j+3]=preddat$predBin
    }
  alldat[[y]]=outdat
}

allMildewDat=alldat

df=data.frame(Year=years[1:18]+1,predN=apply(predN,1,mean),
              predNupper=apply(predN,1,quantile,c(.025,.975),na.rm=T)[2,],
              predNlower=apply(predN,1,quantile,c(.025,.975),na.rm=T)[1,],
              predP=apply(predP,1,mean),
              predPupper=apply(predP,1,quantile,c(.025,.975),na.rm=T)[2,],
              predPlower=apply(predP,1,quantile,c(.025,.975),na.rm=T)[1,])

plot(years,yearlyDat$PhoPcomb,type="b",pch=16,ylim=c(0,.4),xlab="",las=1,ylab="",main="Podosphaera plantaginis")
xx=2001:2018
polygon(c(xx,rev(xx)),c(df$predPlower[1:18],rev(df$predPupper[1:18])),col = colr, border = FALSE)
points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years,yearlyDat$PhoPcomb,type="b",pch=16)

#Cotesia
test=cut(newdat$PLM2,breaks=br)
obs=tapply(newdat$CotP>0,test,sum,na.rm=T)/tapply(newdat$CotP>-1,test,sum,na.rm=T)

load("fitted_mods/CotPc_modList.RData")
mCotPc2=CotPc_modList[["mCotPc2"]]
load("fitted_mods/CotPe_modList.RData")
mCotPe2=CotPe_modList[["mCotPe2"]]

#Setup
years=unique(newdat$Year)
predN=matrix(NA,nrow=length(years)-1,ncol=100)
predP=matrix(NA,nrow=length(years)-1,ncol=100)
alldat=list()
coldat=list()
extdat=list()
cmod=mCotPc2
emod=mCotPe2

for(y in 1:(length(years)-1)){
  year=years[y+1]
  print(paste("Year",year))
  lastyear=subset(newdat[newdat$Year==year-1,],select=c("Patch","CotP"))
  lastyear$CotP=as.numeric(lastyear$CotP>0)
  thisyear=comb[comb$Year==year,]
  
  preddat=merge(lastyear,thisyear,by="Patch",all.x=T)
  preddat$PLcat=cut(preddat$PLM2,breaks=br)
  
  outdat=data.frame(Year=rep(paste(year),nrow(preddat)),Patch=preddat$Patch,PLcat=preddat$PLcat,Pred=rep(NA,nrow(preddat)))
  cdat=data.frame(Year=rep(paste(year),nrow(preddat)),Patch=preddat$Patch,PLcat=preddat$PLcat,Pred=rep(NA,nrow(preddat)))
  edat=data.frame(Year=rep(paste(year),nrow(preddat)),Patch=preddat$Patch,PLcat=preddat$PLcat,Pred=rep(NA,nrow(preddat)))
  
  
  for(j in 1:100){

    #Col
    rr=MASS::mvrnorm(1,mu=summary(cmod)$coef[,1],Sigma=vcov(cmod))
    rYear=ranef(cmod)$Year[which(row.names(ranef(cmod)$Year)==year),1]
    if(length(rYear)<1){rYear=0} #Fix to sampling
    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(rr[-1])))
    rr[1]=rr[1]+rYear
    predc=invlogit(as.matrix(pd)%*%as.matrix(rr))
    
    #Ext
    rr=MASS::mvrnorm(1,mu=summary(emod)$coef[,1],Sigma=vcov(emod))
    rYear=ranef(emod)$Year[which(row.names(ranef(emod)$Year)==year),1]
    if(length(rYear)<1){rYear=0} #Fix to sampling
    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(rr[-1])))
    rr[1]=rr[1]+rYear
    prede=invlogit(as.matrix(pd)%*%as.matrix(rr))
    
    preddat$pred=rep(NA,nrow(preddat))
    preddat$pred[which(preddat$CotP==1)]=(1-prede[which(preddat$CotP==1)])
    preddat$pred[which(preddat$CotP==0)]=predc[which(preddat$CotP==0)]
    
    preddat$predBin=rbinom(nrow(preddat),1,p=preddat$pred)
    predN[y,j]=sum(preddat$predBin,na.rm=T)
    predP[y,j]=sum(preddat$predBin,na.rm=T)/length(preddat$predBin) #Scale by total patches surveyed
    
    outdat[,j+3]=preddat$predBin
    cdat[,j+3]=predc
    edat[,j+3]=prede
    
    #predN[y,j]=sum(preddat$pred,na.rm=T)
    #predP[y,j]=sum(preddat$pred,na.rm=T)/(sum(preddat$pred>-1,na.rm=T))
    #predP[y,j]=sum(preddat$pred,na.rm=T)/nrow(preddat) #Total Patches surveyes
      }
  alldat[[y]]=outdat
  coldat[[y]]=cdat
  extdat[[y]]=edat
  
}

allCotesiaDat=alldat

df=data.frame(Year=years[1:18]+1,predN=apply(predN,1,mean),
              predNupper=apply(predN,1,quantile,c(.025,.975),na.rm=T)[2,],
              predNlower=apply(predN,1,quantile,c(.025,.975),na.rm=T)[1,],
              predP=apply(predP,1,mean),
              predPupper=apply(predP,1,quantile,c(.025,.975),na.rm=T)[2,],
              predPlower=apply(predP,1,quantile,c(.025,.975),na.rm=T)[1,])

plot(years,yearlyDat$CotPAdj,type="b",pch=16,ylim=c(0,.06),xlab="",ylab="",las=1,main="Cotesia melitaearum")
xx=2001:2018
polygon(c(xx,rev(xx)),c(df$predPlower[1:18],rev(df$predPupper[1:18])),col = colr, border = FALSE)
points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years,yearlyDat$CotPAdj,type="b",pch=16)

######################################
#### - Plot population dynamics - ####
######################################
yearlyDat=read.csv("yearlyDat.csv")

cols=c(rgb(0,0,1,.5), rgb(0,1,0,.5), rgb(1,0,0,.5))

#Cinxia
predP=matrix(NA,nrow=length(allCinxiaDat),ncol=100)
for(i in 1:(length(allCinxiaDat))){
  sums=apply(allCinxiaDat[[i]][,4:103], 2, sum, na.rm=T)
  ns=apply(allCinxiaDat[[i]][,4:103]>-1, 2, sum, na.rm=T)
  predP[i,]=sums/ns
}

df=data.frame(Year=years[1:18]+1,
              predP=apply(predP,1,mean,na.rm=T),
              predPupper=apply(predP,1,quantile,c(.025,.975),na.rm=T)[2,],
              predPlower=apply(predP,1,quantile,c(.025,.975),na.rm=T)[1,])

x11()
par(mfrow=c(1,1),mar=c(4,4,2,4))

plot(years,yearlyDat$MelC,type="b",pch=16,ylim=c(0,.4),xlab="",las=1,ylab="",main="")
xx=2001:2018
polygon(c(xx,rev(xx)),c(df$predPlower[1:18],rev(df$predPupper[1:18])),col = cols[1], border = FALSE)
#points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years,yearlyDat$MelC,type="b",pch=16)

mtext("Proportion of patches occupied", 2, line=2.5)

#Mildew
nc=ncol(allMildewDat[[1]])
predP=matrix(NA,nrow=length(allMildewDat),ncol=100)
for(i in 1:(length(allMildewDat))){
  sums=apply(allMildewDat[[i]][,4:nc], 2, sum, na.rm=T)
  ns=apply(allMildewDat[[i]][,4:nc]>-1, 2, sum, na.rm=T)
  predP[i,]=sums/ns
}

df=data.frame(Year=years[1:18]+1,
              predP=apply(predP,1,mean,na.rm=T),
              predPupper=apply(predP,1,quantile,c(.025,.975),na.rm=T)[2,],
              predPlower=apply(predP,1,quantile,c(.025,.975),na.rm=T)[1,])

xx=2001:2018
polygon(c(xx,rev(xx)),c(df$predPlower[1:18],rev(df$predPupper[1:18])),col = cols[2], border = FALSE)
#points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years,yearlyDat$PhoPcomb,type="b",pch=16)

#Cotesia
predP=matrix(NA,nrow=length(allCotesiaDat),ncol=100)
for(i in 1:(length(allCotesiaDat))){
  sums=apply(allCotesiaDat[[i]][,4:103], 2, sum, na.rm=T)
  ns=apply(allCotesiaDat[[i]][,4:103]>-1, 2, length)
  predP[i,]=sums/ns
}

df=data.frame(Year=years[1:18]+1,
              predP=apply(predP,1,mean,na.rm=T),
              predPupper=apply(predP,1,quantile,c(.025,.975),na.rm=T)[2,],
              predPlower=apply(predP,1,quantile,c(.025,.975),na.rm=T)[1,])

xx=2001:2018
polygon(c(xx,rev(xx)),c(df$predPlower[1:18],rev(df$predPupper[1:18])),col = cols[3], border = FALSE)
#points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years,yearlyDat$CotPAdj,type="b",pch=16)

#Plantago
par(new=T)
plot(yearlyDat$Year,yearlyDat$mPLM2,type="l",lwd=3,col="darkgreen",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4,c(5,10,15,20,25),c(5,10,15,20,25),las=1)
mtext(4,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)

legend("topleft",col=c("darkgreen", cols),lwd=3, legend=c("Plantago lanceolata", "Butterfly", "Mildew", "Parasitoid"),bty="n")

plot(yearlyDat$Year,yearlyDat$mVS,type="l",lwd=3,col="darkgreen")

###################################
#### - Community predictions - ####
###################################

allCinxiaDat[[1]][1:10,1:10]
allCotesiaDat[[1]][1:10,1:10]
allMildewDat[[1]][1:10,1:10]

dist=matrix(NA,ncol=6,nrow=100)
medians=matrix(NA,ncol=6,nrow=18)
upper=matrix(NA,ncol=6,nrow=18)
lower=matrix(NA,ncol=6,nrow=18)
w=list()

for(y in 1:(length(years)-1)){
  for(i in 1:100){
mc=allCinxiaDat[[y]][,c(1:2,i+3)]
names(mc)=c("Year","Patch","MelC")
cm=allCotesiaDat[[y]][,c(1:2,i+3)]
names(cm)=c("Year","Patch","CotP")
pp=allMildewDat[[y]][,c(1:2,i+3)]
names(pp)=c("Year","Patch","PhoP")

rep=cbind(mc,cm,pp)
rep=rep[,c(1:3,6,9)]
head(rep,10)

rep$CotP[which(rep$MelC==0)]=NA #'Butterfly extinction'
rep$CotP[which(is.na(rep$MelC))]=NA #Missing data

Mel_col_patches=which(rep$MelC==1 & is.na(rep$CotP))
rep$CotP[Mel_col_patches] = rbinom(length(Mel_col_patches),1,coldat[[y]][Mel_col_patches,c(i+3)]) #Joint colonisations
head(rep,20)

rep$State=paste(as.numeric(rep$MelC>0),as.numeric(rep$CotP>0),as.numeric(rep$PhoP>0),sep="_")

rep$State=factor(rep$State,levels=c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1"))

dist[i,]=table(rep$State)
}
  a=apply(dist,2,quantile,c(.025,.5,.975))
  medians[y,]=a[2,]
  upper[y,]=a[3,]
  lower[y,]=a[1,]
}

colnames(medians)=colnames(upper)=colnames(lower)=levels(rep$State)
rownames(medians)=rownames(upper)=rownames(lower)=years[-1]
wide=read.csv("wide.csv")
rownames(wide)=2000:2017
colnames(wide)=sub("X","",colnames(wide))

mm=match(colnames(wide),colnames(medians))
medians=medians[,mm]
upper=upper[,mm]
lower=lower[,mm]
dim(wide)
dim(medians)

head(wide,20) #Observed data, 2009=NA
head(medians,20)

rs=rowSums(medians)
wide=t(apply(wide,1,function(x){x/sum(x)}))
medians=t(apply(medians,1,function(x){x/sum(x)}))

upper2=upper
for(i in 1:nrow(upper2)){upper2[i,]=upper[i,]/rs[i]}
upper=upper2

lower2=lower
for(i in 1:nrow(lower2)){lower2[i,]=lower[i,]/rs[i]}
lower=lower2

xx=2000:2017
colr=rgb(0,0,1, alpha=.5)

x11()
par(mfrow=c(2,3))
for(i in 1:6){
  plot(xx+1, medians[,i], type = "n",ylab="Proportion of patches",las=1,xlab="",xlim=c(2000,2018),
       ylim=c(min(c(lower[,i],wide[,i]),na.rm=T),max(c(upper[,i],wide[,i]),na.rm=T)),main=colnames(wide)[i])
  polygon(c(xx+1,rev(xx+1)),c(lower[,i],rev(upper[,i])),col = colr, border = FALSE)
  
  lines(rownames(wide),wide[,i],lwd=2)
  lines(xx+1,medians[,i],lwd=2,col="darkblue")
  }

medians
ord=order(medians[,1])
barplot(t(medians[,]),legend=F,col=topo.colors(6))

#### - Communities with Mildew lag - ####
dist=matrix(NA,ncol=6,nrow=100)
medians=matrix(NA,ncol=6,nrow=18)
upper=matrix(NA,ncol=6,nrow=18)
lower=matrix(NA,ncol=6,nrow=18)
w=list()

y=3
j=i=1
for(y in 2:(length(years)-1)){
  for(i in 1:100){
    mc=allCinxiaDat[[y]][,c(1:2,i+3)]
    names(mc)=c("Year","Patch","MelC")
    cm=allCotesiaDat[[y]][,c(1:2,i+2)]
    names(cm)=c("Year","Patch","CotP")
    pp=allMildewDat[[y-1]][,c(1:2,i+2)]
    names(pp)=c("Year","Patch","PhoP")
    
    rep=cbind(mc,cm)
    rep=rep[,c(1:3,6)]
    
    rep=merge(rep, pp, by="Patch", all=T)
    rep=rep[,c(1:4,6)]
    colnames(rep)[2]="Year"
    rep=rep[which(rep$Patch%in%coldat[[y]]$Patch),]
    head(rep)
    
    rep$CotP[which(rep$MelC==0)]=NA #'Butterfly extinction'
    rep$CotP[which(is.na(rep$MelC))]=NA #Missing data
    
    Mel_col_patches=which(rep$MelC==1 & is.na(rep$CotP))
    rep$CotP[Mel_col_patches] = rbinom(length(Mel_col_patches),1,coldat[[y]][Mel_col_patches,c(i+2)]) #Joint colonisations
    head(rep,20)
    
    rep$State=paste(as.numeric(rep$MelC>0),as.numeric(rep$CotP>0),as.numeric(rep$PhoP>0),sep="_")
    
    rep$State=factor(rep$State,levels=c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1"))
    
    dist[i,]=table(rep$State)
  }
  a=apply(dist,2,quantile,c(.025,.5,.975))
  medians[y,]=a[2,]
  upper[y,]=a[3,]
  lower[y,]=a[1,]
}

colnames(medians)=colnames(upper)=colnames(lower)=levels(rep$State)
rownames(medians)=rownames(upper)=rownames(lower)=years[-1]
wideLag=read.csv("wideLag.csv")
rownames(wideLag)=2000:2017
colnames(wideLag)=sub("X","",colnames(wide))

mm=match(colnames(wideLag),colnames(medians))
medians=medians[,mm]
upper=upper[,mm]
lower=lower[,mm]
dim(wide)
dim(medians)

head(wideLag,20) #Observed data, 2009=NA, 2000=NULL
head(medians,20) #Predicted data, 2001=NA, 2010 = ?
head(upper,20)
head(lower,20)

rs=rowSums(medians)
wideLag=t(apply(wideLag,1,function(x){x/sum(x)}))
medians=t(apply(medians,1,function(x){x/sum(x)}))

upper2=upper
for(i in 1:nrow(upper2)){upper2[i,]=upper[i,]/rs[i]}
upper=upper2

lower2=lower
for(i in 1:nrow(lower2)){lower2[i,]=lower[i,]/rs[i]}
lower=lower2

xx=2000:2017
x11()
par(mfrow=c(2,3))
for(i in 1:6){
  #polygon(c(xx,rev(xx)),c(lower[1:9,i],rev(upper[1:9,i])),col = "grey75", border = FALSE)
  plot(xx+1, medians[,i], type = "n",ylab="Proportion of patches",las=1,xlab="",xlim=c(2000,2018),
       ylim=c(min(c(lower[,i],wide[,i]),na.rm=T),max(c(upper[,i],wide[,i]),na.rm=T)),main=colnames(wide)[i])
  polygon(c(xx+1,rev(xx+1)),c(lower[,i],rev(upper[,i])),col = colr, border = FALSE)
  
  #lines(xx,medians[1:9,i], lwd = 2)
  #lines(xx+1, upper[,i], col="red",lty=2)
  #lines(xx+1, lower[,i], col="red",lty=2)
  
  lines(rownames(wideLag),wideLag[,i],lwd=2)
  lines(xx+1,medians[,i],lwd=2,col="darkblue")
  
  
}

plot(rownames(wideLag),wideLag[,2],pch=16)
points(rownames(wideLag),wide[,2])

############################################################
#### - Community predictions for the Plantago gradient- ####
############################################################

allCinxiaDat[[1]][1:10,1:10]
allCotesiaDat[[1]][1:10,1:10]
allMildewDat[[1]][1:10,1:10]
dim(allCinxiaDat[[1]])
dim(allMildewDat[[1]])
dim(allCotesiaDat[[1]])

dist=matrix(NA,ncol=6,nrow=100)
medians=matrix(NA,ncol=6,nrow=18)
upper=matrix(NA,ncol=6,nrow=18)
lower=matrix(NA,ncol=6,nrow=18)
w=list()

medianlist=list()
lowerlist=list()
upperlist=list()

y=j=i=1
y=3
for(y in 1:(length(years)-1)){
  for(i in 1:100){
    mc=allCinxiaDat[[y]][,c(1:3,i+3)]
    names(mc)=c("Year","Patch","PLcat","MelC")
    cm=allCotesiaDat[[y]][,c(1:3,i+3)]
    names(cm)=c("Year","Patch","PLcat","CotP")
    pp=allMildewDat[[y]][,c(1:3,i+3)]
    names(pp)=c("Year","Patch","PLcat","PhoP")
    
    rep=cbind(mc,cm,pp)
    rep=rep[,c(1:4,8,12)]
    head(rep,10)
    
    rep$CotP[which(rep$MelC==0)]=NA #'Butterfly extinction'
    rep$CotP[which(is.na(rep$MelC))]=NA #Missing data
    
    Mel_col_patches=which(rep$MelC==1 & is.na(rep$CotP))
    rep$CotP[Mel_col_patches] = rbinom(length(Mel_col_patches),1,coldat[[y]][Mel_col_patches,c(i+3)]) #Joint colonisations
    head(rep,20)
    
    rep$State=paste(as.numeric(rep$MelC>0),as.numeric(rep$CotP>0),as.numeric(rep$PhoP>0),sep="_")
    
    rep$State=factor(rep$State,levels=c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1"))
    
    l=lapply(tapply(rep$State, rep$PLcat, table), function(x) as.data.frame(x))
    
    for(ii in 1:length(l)){
      if(dim(l[[ii]])[1]<1){
          l[[ii]]=data.frame(Var1=factor(c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1")), Freq=rep(NA,6))}}
    
    w[[i]]=sapply(l, function(x) x[,2])
    rownames(w[[i]])=c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1")
    
    dist[i,]=table(rep$State)
  }
  a=apply(dist,2,quantile,c(.025,.5,.975))
  medians[y,]=a[2,]
  upper[y,]=a[3,]
  lower[y,]=a[1,]
  medianlist[[y]]=apply(simplify2array(w), 1:2, mean, na.rm=T)
  upperlist[[y]]=apply(simplify2array(w), 1:2, quantile, .975, na.rm=T)
  lowerlist[[y]]=apply(simplify2array(w), 1:2, quantile, .025, na.rm=T)
  }


temp=medianlist
temp=lapply(medianlist, function(x) apply(x, 2, function(xx){xx/sum(xx)}))

tempmean=apply(simplify2array(temp), 1:2, mean, na.rm=T)
upper=apply(simplify2array(temp), 1:2, quantile, .975, na.rm=T)
lower=apply(simplify2array(temp), 1:2, quantile, .025, na.rm=T)

#Observed
newdat$PLcat=cut(newdat$PLM2,breaks=br) #br from TotalPatchDat
newdat$State=factor(newdat$State,levels=c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1"))
tapply(newdat$State, newdat$PLcat, table)
l=lapply(tapply(newdat$State, newdat$PLcat, table), function(x) as.data.frame(x))
obs=sapply(l, function(x) x[,2])
obs=apply(obs, 2, function(x){x/sum(x)})

obs=obs[c(1,5,2,6,3,4),]
tempmean=tempmean[c(1,5,2,6,3,4),]
upper=upper[c(1,5,2,6,3,4),]
lower=lower[c(1,5,2,6,3,4),]

titles=rownames(tempmean)
titles=c("No species", "Mildew only", "Butterfly only", "Butterfly & Mildew", "Butterfly & Parasitoid", "All species")

br2=NULL #Midpoints
for(i in 1:(length(br)-1)){
  br2[i]=(br[i]+br[i+1])/2
}

xx=br2
xx=log(br2)

x11()
par(mfrow=c(2,3),mar=c(2,2,2,2), oma=c(3,3,0,0))
for(i in 1:6){
  plot(xx, tempmean[i,], type = "n",ylab="",las=1, xlab="",xaxt="n",
       ylim=c(min(lower[i,]), max(upper[i,])), main=as.expression(titles[i]))
  axis(1,seq(-4,8,2),signif(exp(seq(-4,8,2)),1))
  arrows(xx,upper[i,],xx,lower[i,],length=0,col="darkgrey")
  
  #polygon(c(xx+1,rev(xx+1)),c(lower[,i],rev(upper[,i])),col = colr, border = FALSE)
  points(xx, tempmean[i,], pch=16, col="darkblue")
  #points(xx, obs[i,],pch=16, col="black")
 }
mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=1, outer=T)
mtext(2,text="Proportion of patches",line=1, outer=T)

x11()
par(mar=c(4,4,2,9))
barplot(tempmean,col=topo.colors(6),legend.text = titles, axisnames = F,las=1,ylab="",
        args.legend = list(x=22.25, bty="n", xpd=T), ylim=c(.3,1),xpd=F)
axis(1,seq(.75,16.25, length=14),labels=signif(exp(xx),2))
mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5,cex=1.1)
mtext(2,text="Proporton of patches",line=2.5,cex=1.1)


##############################################################
##### - Visualising the effects of host-plant abundance - ####
##############################################################
easyPredCI <- function(model,newdata=NULL,alpha=0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form=NA,newdata=newdata)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## inverse-link function
  linkinv <- family(model)$linkinv
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(conf.low=pred0-crit*pred.se,
                conf.high=pred0+crit*pred.se))
}

#Datasets
#MelCc_dat=na.omit(subset(comb,select=c("MelCc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","VS","logNH","Road","Year","Patch")))
#MelCe_dat=na.omit(subset(comb,select=c("MelCe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","VS","logNH","Road","Year","Patch")))
#PhoPc_dat=na.omit(subset(comb,select=c("PhoPc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","logNH_Mil","logTotalNH","Road","Year","Patch")))
#PhoPe_dat=na.omit(subset(comb,select=c("PhoPe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","logNH_Mil","logTotalNH","Road","Year","Patch")))
#CotPc_dat=na.omit(subset(comb,select=c("CotPc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
#CotPc_dat=CotPc_dat[CotPc_dat$Year!="2018",]
#CotPe_dat=na.omit(subset(comb,select=c("CotPe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
#CotPe_dat=CotPe_dat[CotPe_dat$Year!="2018",]

#Models
#load("fitted_mods/MelCc_modList.RData")
#mMelCc=MelCc_modList[["mMelCc"]]
#load("fitted_mods/MelCe_modList.RData")
#mMelCe14=MelCe_modList[["mMelCe14"]]
#load("fitted_mods/PhoPc_modList.RData")
#mPhoPc9=PhoPc_modList[["mPhoPc9"]]
#load("fitted_mods/PhoPe_modList.RData")
#mPhoPe3=PhoPe_modList[["mPhoPe3"]]
#load("fitted_mods/CotPc_modList.RData")
#mCotPc2=CotPc_modList[["mCotPc2"]]
#load("fitted_mods/CotPe_modList.RData")
#mCotPe2=CotPe_modList[["mCotPe2"]]


#### - Colonisation - ####
x11()
#par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,0,0))
par(mfrow=c(2,4),mar=c(2,2,2,2),oma=c(2,2,0,0))

#Butterly colonisation
head(MelCc_dat)
cmod=mMelCcb
summary(cmod)

newdata=data.frame(logPLM2=seq(min(MelCc_dat$logPLM2,na.rm=T),max(MelCc_dat$logPLM2,na.rm=T),length.out=1000),
                   VS=rep(mean(MelCc_dat$VS,na.rm=T),1000),
                   logNH=rep(mean(MelCc_dat$logNH,na.rm=T),1000),
                   Road=rep(mean(MelCc_dat$Road,na.rm=T),1000),
                   logMayPrec=rep(mean(MelCc_dat$logMayPrec,na.rm=T),1000),
                   Statelast=rep("0_NA_0", 1000))
levels(newdata$Statelast)=levels(factor(MelCc_dat$Statelast))

preds=predict(cmod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(cmod, newdata=newdata)

xx=newdata$logPLM2
plot(xx,preds,type="l",lwd=2,ylim=c(0,1),las=1,col="white",xaxt="n",
     xlab="", ylab="",main="Butterfly")
polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,0,1,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1,c(0,2,4,6,8),signif(exp(c(0,2,4,6,8)),1))

newdata$Statelast=rep("0_NA_1", 1000)
newdata$Statelast=factor(newdata$Statelast)
levels(newdata$Statelast)=rev(levels(factor(MelCc_dat$Statelast)))

preds=predict(cmod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(cmod, newdata=newdata)

polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,1,0,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

legend("topleft", c("Mildew present", "Mildew absent"),col=c(rgb(0,1,0), rgb(0,0,1)),lty=1,lwd=2,bty="n")

mtext(2,text="Colonisation probability",line=2.5, outer=F, xpd=T)

#Cotesia colonisation
head(CotPc_dat)
cmod=mCotPc2b
summary(cmod)


newdata=data.frame(logPLM2=seq(min(CotPc_dat$logPLM2,na.rm=T),max(CotPc_dat$logPLM2,na.rm=T),length.out=1000),
                   VS=rep(mean(CotPc_dat$VS,na.rm=T),1000),
                   logNH_Cot=rep(mean(CotPc_dat$logNH_Cot,na.rm=T),1000),
                   Road=rep(mean(CotPc_dat$Road,na.rm=T),1000),
                   logJunePrec=rep(mean(CotPc_dat$logJunePrec,na.rm=T),1000),
                   Statelast=rep("0_NA_0", 1000))
levels(newdata$Statelast)=levels(factor(CotPc_dat$Statelast))
preds=predict(cmod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(cmod,newdata=newdata)

xx=newdata$logPLM2
plot(xx,preds,type="l",lwd=2,ylim=c(0,.1),las=1,col="white",xaxt="n",
     xlab="", ylab="", main="Parasitoid")
polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,0,1,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1,c(0,2,4,6,8),signif(exp(c(0,2,4,6,8)),1))

newdata$Statelast=rep("1_0_0", 1000)
newdata$Statelast=factor(newdata$Statelast)
levels(newdata$Statelast)=levels(factor(CotPc_dat$Statelast))[c(3,1,2,4)]

preds=predict(cmod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(cmod, newdata=newdata)

polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,1,0,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

legend("topleft", c("Butterfly present", "Butterfly absent"),col=c(rgb(0,1,0), rgb(0,0,1)),lty=1,lwd=2,bty="n")

#Cotesia colonisation, mildew effect
cmod=mCotPc2b
summary(cmod)

newdata=data.frame(logPLM2=seq(min(CotPc_dat$logPLM2,na.rm=T),max(CotPc_dat$logPLM2,na.rm=T),length.out=1000),
                   VS=rep(mean(CotPc_dat$VS,na.rm=T),1000),
                   logNH_Cot=rep(mean(CotPc_dat$logNH_Cot,na.rm=T),1000),
                   Road=rep(mean(CotPc_dat$Road,na.rm=T),1000),
                   logJunePrec=rep(mean(CotPc_dat$logJunePrec,na.rm=T),1000),
                   Statelast=rep("0_NA_0", 1000))
levels(newdata$Statelast)=levels(factor(CotPc_dat$Statelast))
preds=predict(cmod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(cmod, newdata=newdata)

xx=newdata$logPLM2
plot(xx,preds,type="l",lwd=2,ylim=c(0,.1),las=1,col="white",xaxt="n",
     xlab="", ylab="", main="Parasitoid")
polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,0,1,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1,c(0,2,4,6,8),signif(exp(c(0,2,4,6,8)),1))

newdata$Statelast=rep("0_NA_1", 1000)
newdata$Statelast=factor(newdata$Statelast)
levels(newdata$Statelast)=levels(factor(CotPc_dat$Statelast))[c(2,1,3,4)]

preds=predict(cmod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(cmod, newdata=newdata)

polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,1,0,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

legend("topleft", c("Mildew present", "Mildew absent"),col=c(rgb(0,1,0), rgb(0,0,1)),lty=1,lwd=2,bty="n")

#Mildew colonisation
head(PhoPc_dat)
cmod=mPhoPc9b
summary(cmod)

newdata=data.frame(logPLM2=seq(min(PhoPc_dat$logPLM2,na.rm=T),max(PhoPc_dat$logPLM2,na.rm=T),length.out=1000),
                   logNH_Mil=rep(mean(PhoPc_dat$logNH_Mil,na.rm=T),1000),
                   logTotalNH=rep(mean(PhoPc_dat$logTotalNH,na.rm=T),1000),
                   Road=rep(mean(PhoPc_dat$Road,na.rm=T),1000),
                   logJunePrec=rep(mean(PhoPc_dat$logJunePrec,na.rm=T),1000),
                   logAugPrec=rep(mean(PhoPc_dat$logAugPrec,na.rm=T),1000),
                   Statelast=rep("0_NA_0", 1000))
levels(newdata$Statelast)=levels(factor(PhoPc_dat$Statelast))
head(newdata)

preds=predict(cmod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(cmod, newdata=newdata)

xx=newdata$logPLM2
plot(xx,preds,type="l",lwd=2,ylim=c(0,1),las=1,col="white",xaxt="n",
     xlab="", ylab="",main="Mildew")
polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,0,1,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1,c(0,2,4,6,8),signif(exp(c(0,2,4,6,8)),1))

newdata$Statelast=rep("1_0_0", 1000)
newdata$Statelast=factor(newdata$Statelast)
levels(newdata$Statelast)=levels(factor(PhoPc_dat$Statelast))[c(2,1,3)]

preds=predict(cmod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(cmod, newdata=newdata)

polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,1,0,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

legend("topleft", c("Butterfly present", "Butterfly absent"),col=c(rgb(0,1,0), rgb(0,0,1)),lty=1,lwd=2,bty="n")

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=0.5, outer=T)
#mtext(2,text="Colonisation probability",line=0.5, outer=T)



#### - Extinction - ####
#x11()
#par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,0,0))

#Butterly extinction: Mildew effects
head(MelCe_dat)
emod=mMelCe14b
summary(emod)

newdata=data.frame(logPLM2=seq(min(MelCe_dat$logPLM2,na.rm=T),max(MelCe_dat$logPLM2,na.rm=T),length.out=1000),
                   VS=rep(mean(MelCe_dat$VS,na.rm=T),1000),
                   logNH=rep(mean(MelCe_dat$logNH,na.rm=T),1000),
                   Road=rep(mean(MelCe_dat$Road,na.rm=T),1000),
                   logJunePrec=rep(mean(MelCe_dat$logJunePrec,na.rm=T),1000),
                   logJulyPrec=rep(mean(MelCe_dat$logJulyPrec,na.rm=T),1000),
                   logAugPrec=rep(mean(MelCe_dat$logAugPrec,na.rm=T),1000),
                   Statelast=rep("1_0_0", 1000))
levels(newdata$Statelast)=levels(factor(MelCe_dat$Statelast))
head(newdata)

preds=predict(emod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(emod, newdata=newdata)

xx=newdata$logPLM2
plot(xx,preds,type="l",lwd=2,ylim=c(0,1),las=1,col="white",xaxt="n",
     xlab="", ylab="",main="Butterfly")
polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,0,1,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1,c(0,2,4,6,8),signif(exp(c(0,2,4,6,8)),1))

newdata$Statelast=rep("1_0_1", 1000)
newdata$Statelast=factor(newdata$Statelast)
levels(newdata$Statelast)=levels(factor(MelCe_dat$Statelast))[c(2,1,3,4)]

preds=predict(emod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(emod, newdata=newdata)

polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,1,0,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

legend("topleft", c("Mildew present", "Mildew absent"),col=c(rgb(0,1,0), rgb(0,0,1)),lty=1,lwd=2,bty="n")

mtext(2,text="Extinction probability",line=2.5, outer=F, xpd=T)

#Butterly extinction: Cotesia effects
emod=mMelCe14b
summary(emod)

newdata=data.frame(logPLM2=seq(min(MelCe_dat$logPLM2,na.rm=T),max(MelCe_dat$logPLM2,na.rm=T),length.out=1000),
                   VS=rep(mean(MelCe_dat$VS,na.rm=T),1000),
                   logNH=rep(mean(MelCe_dat$logNH,na.rm=T),1000),
                   Road=rep(mean(MelCe_dat$Road,na.rm=T),1000),
                   logJunePrec=rep(mean(MelCe_dat$logJunePrec,na.rm=T),1000),
                   logJulyPrec=rep(mean(MelCe_dat$logJulyPrec,na.rm=T),1000),
                   logAugPrec=rep(mean(MelCe_dat$logAugPrec,na.rm=T),1000),
                   Statelast=rep("1_0_0", 1000))
levels(newdata$Statelast)=levels(factor(MelCe_dat$Statelast))
head(newdata)

preds=predict(emod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(emod, newdata=newdata)

xx=newdata$logPLM2
plot(xx,preds,type="l",lwd=2,ylim=c(0,1),las=1,col="white",xaxt="n",
     xlab="", ylab="",main="Butterfly")
polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,0,1,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1,c(0,2,4,6,8),signif(exp(c(0,2,4,6,8)),1))

newdata$Statelast=rep("1_1_0", 1000)
newdata$Statelast=factor(newdata$Statelast)
levels(newdata$Statelast)=levels(factor(MelCe_dat$Statelast))[c(3,1,2,4)]

preds=predict(emod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(emod, newdata=newdata)

polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,1,0,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

legend("topleft", c("Parasitoid present", "Parasitoid absent"),col=c(rgb(0,1,0), rgb(0,0,1)),lty=1,lwd=2,bty="n")

#Cotesia extinction: Mildew effect
head(CotPe_dat)
emod=mCotPe2b
summary(emod)


newdata=data.frame(logPLM2=seq(min(CotPe_dat$logPLM2,na.rm=T),max(CotPe_dat$logPLM2,na.rm=T),length.out=1000),
                   VS=rep(mean(CotPe_dat$VS,na.rm=T),1000),
                   logNH_Cot=rep(mean(CotPe_dat$logNH_Cot,na.rm=T),1000),
                   Road=rep(mean(CotPe_dat$Road,na.rm=T),1000),
                   logJunePrec=rep(mean(CotPe_dat$logJunePrec,na.rm=T),1000),
                   Statelast=rep("1_1_0", 1000))
levels(newdata$Statelast)=levels(factor(CotPe_dat$Statelast))
preds=predict(emod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(emod, newdata=newdata)

xx=newdata$logPLM2
plot(xx,preds,type="l",lwd=2,ylim=c(0,1),las=1,col="white",xaxt="n",main="Parasitoid",
     xlab="", ylab="")
polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,0,1,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1,c(0,2,4,6,8),signif(exp(c(0,2,4,6,8)),1))

newdata$Statelast=rep("1_1_1", 1000)
newdata$Statelast=factor(newdata$Statelast)
levels(newdata$Statelast)=levels(factor(CotPe_dat$Statelast))[c(2,1)]

preds=predict(emod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(emod, newdata=newdata)

polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,1,0,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

legend("topleft", c("Mildew present", "Mildew absent"),col=c(rgb(0,1,0), rgb(0,0,1)),lty=1,lwd=2,bty="n")

#Mildew extinction: Butterfly effect
head(PhoPe_dat)
emod=mPhoPe3b
summary(emod)

newdata=data.frame(logPLM2=seq(min(PhoPe_dat$logPLM2,na.rm=T),max(PhoPe_dat$logPLM2,na.rm=T),length.out=1000),
                   logNH_Mil=rep(mean(PhoPe_dat$logNH_Mil,na.rm=T),1000),
                   logTotalNH=rep(mean(PhoPe_dat$logTotalNH,na.rm=T),1000),
                   Road=rep(mean(PhoPe_dat$Road,na.rm=T),1000),
                   logJulyPrec=rep(mean(PhoPe_dat$logJulyPrec,na.rm=T),1000),
                   Statelast=rep("0_NA_1", 1000))
levels(newdata$Statelast)=levels(factor(PhoPe_dat$Statelast))
head(newdata)

preds=predict(emod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(emod, newdata=newdata)

xx=newdata$logPLM2
plot(xx,preds,type="l",lwd=2,ylim=c(0,1),las=1,col="white",xaxt="n",
     xlab="", ylab="",main="Mildew")
polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,0,1,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1,c(0,2,4,6,8),signif(exp(c(0,2,4,6,8)),1))

newdata$Statelast=rep("1_0_1", 1000)
newdata$Statelast=factor(newdata$Statelast)
levels(newdata$Statelast)=levels(factor(PhoPe_dat$Statelast))[c(2,1,3)]

preds=predict(emod, newdata=newdata,re.form=NA, type="response")
predCI=easyPredCI(emod, newdata=newdata)

polygon(c(xx,rev(xx)),c(predCI[,1],rev(predCI[,2])),col = rgb(0,1,0,.5), border = FALSE)
lines(xx,preds,type="l",lwd=2)

legend("topleft", c("Butterfly present", "Butterfly absent"),col=c(rgb(0,1,0), rgb(0,0,1)),lty=1,lwd=2,bty="n")


mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=1, outer=T)
#mtext(2,text="Extinction probability",line=0.5, outer=T)

###################################################################
###################################################################
###################################################################

##############
invlogit=function(x){1/(1+exp(-x))}

####Compile transition matrices####

summary(mMelCcb)
summary(mMelCe14b)

summary(mPhoPc9b)
summary(mPhoPe3b)

summary(mCotPc2b)
summary(mCotPe2b)

#Overall support for inclusion of patch states
full=AIC(mMelCcb)+AIC(mMelCe14b)+AIC(mPhoPc9b)+AIC(mPhoPe3b)+AIC(mCotPc2b)+AIC(mCotPe2b)
null=AIC(mMelCca)+AIC(mMelCe14a)+AIC(mPhoPc9a)+AIC(mPhoPe3a)+AIC(mCotPc2a)+AIC(mCotPe2a)
full
null
null-full

#000
c1=invlogit(summary(mMelCcc)$coef[6,1])
c2=invlogit(summary(mCotPc2c)$coef[6,1])
c3=invlogit(summary(mPhoPc9c)$coef[6,1])
  
p000_000=(1-c1)*(1-c3) #No Cot
p000_001=(1-c1)*(c3) #No Cot
p000_100=c1*(1-c2)*(1-c3)
p000_110=c1*c2*(1-c3)
p000_101=c1*(1-c2)*c3
p000_111=c1*c2*c3

#001
c1=invlogit(summary(mMelCcc)$coef[7,1])
c2=invlogit(summary(mCotPc2c)$coef[7,1])
e3=invlogit(summary(mPhoPe3c)$coef[5,1])

p001_000=(1-c1)*e3
p001_001=(1-c1)*(1-e3)
p001_100=c1*(1-c2)*e3
p001_110=c1*c2*e3
p001_101=c1*(1-c2)*(1-e3)
p001_111=c1*c2*(1-e3)

#100
c2=invlogit(summary(mCotPc2c)$coef[8,1])
c3=invlogit(summary(mPhoPc9c)$coef[7,1])
e1=invlogit(summary(mMelCe14c)$coef[8,1])
  
p100_000=e1*(1-c3)
p100_001=e1*c3
p100_100=(1-e1)*(1-c2)*(1-c3)
p100_110=(1-e1)*c2*(1-c3)
p100_101=(1-e1)*(1-c2)*c3
p100_111=(1-e1)*c2*c3

#110
c3=invlogit(summary(mPhoPc9c)$coef[8,1])
e1=invlogit(summary(mMelCe14c)$coef[10,1])
e2=invlogit(summary(mCotPe2c)$coef[6,1])

p110_000=e1*(1-c3)
p110_001=e1*c3
p110_100=(1-e1)*e2*(1-c3)
p110_110=(1-e1)*(1-e2)*(1-c3)
p110_101=(1-e1)*e2*c3
p110_111=(1-e1)*(1-e2)*c3

#101
c2=invlogit(summary(mCotPc2c)$coef[9,1])
e1=invlogit(summary(mMelCe14c)$coef[9,1])
e3=invlogit(summary(mPhoPe3c)$coef[6,1])

p101_000=e1*e3
p101_001=e1*(1-e3)
p101_100=(1-e1)*(1-c2)*e3
p101_110=(1-e1)*c2*e3
p101_101=(1-e1)*(1-c2)*(1-e3)
p101_111=(1-e1)*c2*(1-e3)

#111
e1=invlogit(summary(mMelCe14c)$coef[11,1])
e2=invlogit(summary(mCotPe2c)$coef[7,1])
e3=invlogit(summary(mPhoPe3c)$coef[7,1])
  
p111_000=e1*e3
p111_001=e1*(1-e3)
p111_100=(1-e1)*e2*e3
p111_110=(1-e1)*(1-e2)*e3
p111_101=(1-e1)*e2*(1-e3)
p111_111=(1-e1)*(1-e2)*(1-e3)

vals=c(p000_000,p000_001,p000_100,p000_110,p000_101,p000_111,
       p001_000,p001_001,p001_100,p001_110,p001_101,p001_111,
       p100_000,p100_001,p100_100,p100_110,p100_101,p100_111,
       p110_000,p110_001,p110_100,p110_110,p110_101,p110_111,
       p101_000,p101_001,p101_100,p101_110,p101_101,p101_111,
       p111_000,p111_001,p111_100,p111_110,p111_101,p111_111)

mat=matrix(vals,nrow=6,byrow=T)
mat
rowSums(mat)
colnames(mat)=rownames(mat)=c("000","001","100","110","101","111")

library(corrplot)
x11()
corrplot(sqrt(mat),is.corr=F,method="color",tl.col="black",title="",mar=c(0,0,0,4), cl.pos="n",col=colorRampPalette(c("Red","White","Blue"))(200))
par(xpd=T)
colorlegend(xlim=c(6.7,7.7),ylim=c(0.5,6.5),col=colorRampPalette(c("Red","White","Blue"))(200)[101:200], labels=signif(seq(0,.98, length.out=11),1), at=sqrt(seq(0,.98, length.out=11)), cex=.7)

?colorlegend
?corrplot

#################################################################################
#### - HMSC - ####
#################################################################################
library(Hmsc)

allspecies_dat=subset(comb,select=c("MelCc","MelCe","CotPc","CotPe","PhoPc","PhoPe",
                                          "logMayPrec","logJunePrec","logJulyPrec","logAugPrec",
                                          "logArea","logPLM2","VS","logNH","logNH_Mil","logNH_Cot","Road","Year","Patch"))

drop=which(is.na(rowSums(allspecies_dat[,7:17])))
allspecies_dat=allspecies_dat[-drop,]
head(allspecies_dat)

Y=as.matrix(allspecies_dat[,1:6])
head(Y)

X=subset(allspecies_dat[,7:17])
head(X)

xList=list()
for(i in 1:6){xList[[i]]=X}
xList[[3]]$logNH=X$logNH_Cot
xList[[4]]$logNH=X$logNH_Cot
xList[[5]]$logNH=X$logNH_Mil
xList[[6]]$logNH=X$logNH_Mil


dfPi=subset(allspecies_dat,select=c("Year","Patch"))
dfPi$Year=factor(dfPi$Year)
dfPi$Patch=factor(dfPi$Patch)
head(dfPi,5)

rL1=HmscRandomLevel(units=unique(dfPi$Year))
rL2=HmscRandomLevel(units=unique(dfPi$Patch))
#rL2=HmscRandomLevel(data=xy,priors="default")

XFormula=~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec

TrData=data.frame(colext=factor(c("col","ext","col","ext","col","ext")))
TrFormula=~colext

m=Hmsc(Y=Y, XFormula = XFormula, XData = xList,
        studyDesign = dfPi, ranLevels = list(Year=rL1, Patch=rL2),
        TrFormula = TrFormula, TrData = TrData,
        distr = "probit")
head(m$Y)
head(m$X[[1]])
head(m$studyDesign)

dim(m$Y)
dim(m$X[[1]])
dim(m$studyDesign)

samples=100
thin=1
transient=.5*(thin*samples)
adaptNf=.4*(thin*samples)
nChains=1

m=sampleMcmc(hM = m, samples = samples, transient = transient, thin = thin, adaptNf=rep(adaptNf, m$nr), nChains=1, nParallel=1)


post=convertToCodaObject(out)
plot(post$Beta)

summary(post$Beta)
beta=getPostEstimate(hM=out,"Beta")
beta



#### - Predictions for the Plantago gradient - ####
outPL=list()
for(e in 1:length(alldat)){
  outPL[[e]]=matrix(NA,nrow=100,ncol=length(br)-1)
}

for(i in 1:length(alldat)){
  for(j in 1:100){
    outPL[[i]][j,]=tapply(alldat[[i]][,j+3],alldat[[i]]$PLcat,sum,na.rm=T)/tapply(alldat[[i]][,j+3]>-1,alldat[[i]]$PLcat,sum,na.rm=T)
  }}

outPL[[1]]

means=upper=lower=matrix(NA,nrow=18,ncol=length(br)-1)
for(e in 1:length(alldat)){
  means[e,]=apply(outPL[[e]],2,mean,na.rm=T)
  upper[e,]=apply(outPL[[e]],2,quantile,.025,na.rm=T)
  lower[e,]=apply(outPL[[e]],2,quantile,.975,na.rm=T)
}

Means=apply(means,2,mean,na.rm=T)
altCI=apply(means,2,quantile, c(.025, .975), na.rm=T)
Lower=apply(upper,2,mean,na.rm=T)
Upper=apply(lower,2,mean,na.rm=T)

br2=NULL #Midpoints
for(i in 1:(length(br)-1)){
  br2[i]=(br[i]+br[i+1])/2
}

xx=log(br2)

plot(xx,Means,pch=16,col="darkgrey",ylim=c(-.05,.7),type="p", xlab="Plantago cover", ylab="Butterfly occurrence")
arrows(xx, Lower, xx, Upper, length=0)
#polygon(c(xx,rev(xx)),c(Lower,rev(Upper)),col = colr, border = FALSE)
points(xx,Means,pch=16,col="darkblue",ylim=c(-.2,1))
points(xx,obs,pch=16)

cor(Means,obs)

