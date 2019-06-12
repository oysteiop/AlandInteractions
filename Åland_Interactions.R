##########################################################################
######################### - Åland Interactions - #########################
##########################################################################

rm(list=ls())
setwd("Z:/data/aland/")
list.files()

library(plyr)
library(reshape2)
library(corrplot)

memory.size()
memory.limit(50000)

##########################################
##### - Read and combine data files - ####
##########################################
indat=read.csv2("mergedpatchdata2016_Cot_Fix.csv",dec=".")
indat=indat[indat$Year!="2017",]
indat2017=read.csv2("mergedpatchdata2017.csv",dec=".")
indat2018=read.csv2("mergedpatchdata2018.csv",dec=".")
indat=rbind(indat,indat2017,indat2018)

oldindat=read.csv2("mergedpatchdata.csv",dec=".")
oldindat$SurveyType=substr(oldindat$Survey,1,2)
oldindat$Year=factor(oldindat$Year)
oldindat$Patch=factor(oldindat$Patch)

SiNdata=read.csv2("SiN.csv",dec=".")
SiNdata=subset(SiNdata,select=c("Name","SIN1","SIN2015"))
SiNdata$Patch=factor(SiNdata$Name)

indat=subset(indat, select= c("Locality","Dataset","Year","PLM2","VSM2","VS","PL","Mildew","Mildewcorrected","Absabund","Nest.Total","Total.Cotesia.Nests"))
indat$SurveyType=substr(indat$Dataset,1,2)
indat$Year=factor(indat$Year)
indat$Patch=factor(indat$Locality)
indat=merge(indat,SiNdata,by="Patch",all.x=T)

##############################
#### - Mildew formating - ####
##############################
tapply(indat$Mildew>-1,indat$Year,sum,na.rm=T)
tapply(indat$Absabund>-1,indat$Year,sum,na.rm=T)

#Replace Mildew with corrected value
md=indat$Mildew

w1=which(!is.na(md) & !is.na(indat$Mildewcorrected) & indat$Mildewcorrected!=md)
md[w1]=indat$Mildewcorrected[w1]

w2=which(is.na(md) & !is.na(indat$Mildewcorrected))
md[w2]=indat$Mildewcorrected[w2]
#check=cbind(indat$Mildew,indat$Mildewcorrected,md)
#View(check[which(check[,2]>-1),])
indat$Mildew=md

#Combine Mildew and Absabund
indat$Mildew2=apply(subset(indat, select=c("Mildew", "Absabund")), 1, sum, na.rm=T)
indat$Mildew2[which(is.na(indat$Mildew) & is.na(indat$Absabund))]=NA
head(subset(indat, select=c("Year","Mildew","Absabund","Mildew2")),20)

##########################################
#### - Compile yearly summary stats - ####
##########################################
AF=indat[indat$SurveyType=="AF",] #Butterfly autumn survey

#Insert 0 for NAs (blank cells in database) for Nest.Total
AF$Nest.Total[which(is.na(AF$Nest.Total))]=0

#Insert 0 for NAs for Mildew data.
AF$Mildew2[which(is.na(AF$Mildew2))]=0
AF$Absabund[which(is.na(AF$Absabund))]=0

years=2000:2018

AF=AF[AF$Year %in% years,]
AF$Year=factor(AF$Year)

yearly=ddply(AF,.(Year),summarize,
          np=sum(PLM2>-1,na.rm=T),
          nPatch=sum(Nest.Total>(-1),na.rm=T),
          nMildew=sum(Mildew>-1,na.rm=T),
          nMildew2=sum(Mildew2>-1,na.rm=T),
          nPopCinxia=sum(Nest.Total>0,na.rm=T),
          nPopMildew=sum(Mildew>0,na.rm=T),
          mPLM2=mean(PLM2,na.rm=T),
          mVSM2=mean(VSM2,na.rm=T),
          PlantSumM2=mean(PLM2,na.rm=T)+mean(VSM2,na.rm=T),
          mPL=mean(PL,na.rm=T),
          mVS=mean(VS,na.rm=T),
          PlantSumCat=mean(PL,na.rm=T)+mean(VS,na.rm=T),          
          MelC=sum(Nest.Total>0,na.rm=T)/sum(Nest.Total>(-1),na.rm=T),
          PhoP=sum(Mildew>0,na.rm=T)/sum(Mildew>-1,na.rm=T), #Based on Mildew variable
          PhoPA=sum(Absabund>0,na.rm=T)/sum(Absabund>-1,na.rm=T), #Based on Absabund variable
          PhoPcomb=sum(Mildew2>0,na.rm=T)/sum(Mildew2>-1,na.rm=T)) #Based on combination of both

#Mildew survey
MF=indat[indat$SurveyType=="MF",] #Mildew fall survey
MF=MF[MF$Year %in% years,]
MF$Year=factor(MF$Year)

#Insert 0 for NAs for Mildew data. NB TO BE CHECKED THAT THIS IS CORRECT
MF$Mildew2[which(is.na(MF$Mildew2))]=0
MF$Absabund[which(is.na(MF$Absabund))]=0

yearlyMF=ddply(MF,.(Year),summarize,
               nMildewMFA=sum(Absabund>-1,na.rm=T),
               PhoPAMF=sum(Absabund>0,na.rm=T)/sum(Absabund>-1,na.rm=T))

yearlyDat=merge(yearly,yearlyMF,by="Year",all=T)
head(subset(yearlyDat, select=c("Year","nMildew","nMildew2","PhoP","PhoPA","PhoPcomb","nMildewMFA","PhoPAMF")),20)

###Insert Mildew fall survey value in PhoPcomb (entered in different database that year)
yearlyDat$PhoPcomb[17]=yearlyDat$PhoPAMF[17]

#Cotesia survey
CS=indat[indat$SurveyType=="CS",] #Cotesia survey
CS=CS[CS$Year %in% years,]
CS$Year=factor(CS$Year)

yearlyC=ddply(CS,.(Year),summarize,
              nCot=sum(Total.Cotesia.Nests>-1,na.rm=T),
              CotSum=sum(Total.Cotesia.Nests>0,na.rm=T),
              CotP=sum(Total.Cotesia.Nests>0,na.rm=T)/sum(Total.Cotesia.Nests>-1,na.rm=T))

head(yearlyC,20)

yearlyDat=merge(yearlyDat,yearlyC,by="Year",all=T)
yearlyDat$CotSum[11]=37 #(value from SvN)


#Adjust Cotesia by one year
CotPAdj=NULL
for(i in 1:nrow(yearlyDat)){
  CotPAdj[i]=yearlyDat$CotSum[i+1]/yearlyDat$nPatch[i+1]
}
yearlyDat$CotPAdj=CotPAdj

head(yearlyDat)
#write.csv(yearlyDat, file="yearlyDat.csv",row.names=F)


#######################################
#### - Patchlevel data formating - ####
#######################################

Patchmeans=ddply(AF,.(Patch),summarize,
               PlaL=mean(PLM2,na.rm=T),
               Ver=mean(VS, na.rm=T),
               MelC=sum(Nest.Total>0,na.rm=T)/sum(Nest.Total>-1,na.rm=T),
               PhoP=sum(Mildew2>0, na.rm=T)/sum(Mildew2>-1, na.rm=T)) #Combined Mildew and Absabund

CotMeans=data.frame(tapply(CS$Total.Cotesia.Nests>0, CS$Patch, mean, na.rm=T)/tapply(CS$Total.Cotesia.Nests>-1, CS$Patch, mean, na.rm=T))
CotMeans$Patch=row.names(CotMeans)
names(CotMeans)=c("Cot","Patch")

Patchmeans=merge(patches, Patchmeans, by="Patch", all.x=T)
Patchmeans=merge(Patchmeans, CotMeans, by="Patch", all.x=T)
head(Patchmeans)

Patchmeans$MelC[which(is.na(Patchmeans$MelC))]=0
Patchmeans$PhoP[which(is.na(Patchmeans$PhoP))]=0
Patchmeans$Cot[which(is.na(Patchmeans$Cot))]=0
Patchmeans$PlaL[which(is.na(Patchmeans$PlaL))]=0
Patchmeans$Ver[which(is.na(Patchmeans$Ver))]=0

Patchmeans$PlaL=log(Patchmeans$PlaL+1)
Patchmeans$PlaL_scaled=Patchmeans$PlaL/max(Patchmeans$PlaL)
Patchmeans$Ver_scaled=Patchmeans$Ver/max(Patchmeans$Ver)
summary(Patchmeans)

#Shapefile
library(rgdal)
map=readOGR("aland_shapefile/åland bigland.shp")

pcex=.8

x11()
par(mfrow=c(2,3),mar=c(2,2,2,0))
plot(map, main=expression(paste(italic(Plantago)," ",italic(lanceolata))))
axis(2)
cols=rgb(red=0,green=100/255,blue=0,alpha=Patchmeans$PlaL_scaled*.65)
points(Patchmeans$X*1000, Patchmeans$Y*1000,pch=16,cex=pcex,col=cols)

plot(map, main=expression(paste(italic(Veronica)," ",italic(spicata))))
cols=rgb(red=0,green=100/255,blue=0,Patchmeans$Ver_scaled*.75)
points(Patchmeans$X*1000, Patchmeans$Y*1000,pch=16,cex=pcex,col=cols)

plot.new()

plot(map, main=expression("Butterfly"))
axis(1)
axis(2)
cols=rgb(red=0,green=0,blue=1,Patchmeans$MelC*.75)
points(Patchmeans$X*1000, Patchmeans$Y*1000,pch=16,cex=pcex,col=cols)

plot(map, main=expression("Parasitoid"))
axis(1)
cols=rgb(red=1,green=0,blue=0,Patchmeans$Cot*.75)
points(Patchmeans$X*1000, Patchmeans$Y*1000,pch=16,cex=pcex,col=cols)

plot(map, main=expression("Mildew"))
axis(1)
cols=rgb(red=0,green=1,blue=0,Patchmeans$PhoP*.75)
points(Patchmeans$X*1000, Patchmeans$Y*1000,pch=16,cex=pcex,col=cols)



library(ggplot2)
sp <- ggplot(data = Patchmeans, aes(x=X, y=Y, color=sqrt(PlaL)))+geom_point(size=3)
sp + ggtitle("Plantago lanceolata") + scale_color_gradient(low="white", high="blue")

sp <- ggplot(data = Patchmeans, aes(x=X, y=Y, color=Ver))+geom_point(size=3)
sp + ggtitle("Veronica spicata") + scale_color_gradient(low="white", high="blue")

sp <- ggplot(data = Patchmeans, aes(x=X, y=Y, color=MelC))+geom_point(size=3)
sp + ggtitle("Butterfly") + scale_color_gradient(low="white", high="blue")

sp <- ggplot(data = Patchmeans, aes(x=X, y=Y, color=PhoP))+geom_point(size=3)
sp + ggtitle("Mildew") + scale_color_gradient(low="white", high="blue")

sp <- ggplot(data = Patchmeans, aes(x=X, y=Y, color=Cot))+geom_point(size=3)
sp + ggtitle("Parasitoid") + scale_color_gradient(low="white", high="blue")



PatchDat=ddply(AF,.(Patch,Year),summarize,
          PlaL=PLM2,
          Ver=VSM2,
          PLcat=PL,
          VScat=VS,
          PLVSsum=PLM2+VSM2,
          PLVCcatsum=PL+VS,
          MelC=Nest.Total,
          PhoP=Mildew2) #Combined Mildew and Absabund
PatchDat$PatchYear=paste(PatchDat$Patch,PatchDat$Year,sep="_")
head(PatchDat,10)

CS$Total.Cotesia.Nests[which(is.na(CS$Total.Cotesia.Nests))]=999 #Prelim identifier for real NAs

CotesiaPatchDat=ddply(CS,.(Patch,Year),summarize,
                        CotP=Total.Cotesia.Nests)
CotesiaPatchDat$PatchYear=paste(CotesiaPatchDat$Patch,CotesiaPatchDat$Year,sep="_")
head(CotesiaPatchDat,10)

TotalPatchDat=merge(PatchDat,CotesiaPatchDat[,3:4],by="PatchYear",all.x=T)
head(TotalPatchDat,10)

#Put 0 for all patches not in data, no Cotesia ever recorded. Return NAs for real NAs
TotalPatchDat$CotP[which(is.na(TotalPatchDat$CotP))]=0
TotalPatchDat$CotP[which(TotalPatchDat$CotP>998)]=NA

####Mildew
MildewPatchDat=ddply(MF,.(Patch,Year),summarize,
                      PhoPAMF=Absabund)
MildewPatchDat$PatchYear=paste(MildewPatchDat$Patch,MildewPatchDat$Year,sep="_")
head(MildewPatchDat,10)

TotalPatchDat=merge(TotalPatchDat,MildewPatchDat[,3:4],by="PatchYear",all.x=T)
head(TotalPatchDat,10)

#Insert MF survey values for 2016 Mildew
#summary(TotalPatchDat$PhoP[TotalPatchDat$Year=="2016"])
TotalPatchDat$PhoP[which(TotalPatchDat$Year=="2016")]=TotalPatchDat$PhoPAMF[which(TotalPatchDat$Year=="2016")]

tapply(TotalPatchDat$PhoP>0,TotalPatchDat$Year,sum,na.rm=T)
tapply(TotalPatchDat$PhoP>-1,TotalPatchDat$Year,sum,na.rm=T)

#Adjust Cotesia by one year
CotPAdj=NULL
for(i in 1:nrow(TotalPatchDat)){
  CotPAdj[i]=TotalPatchDat$CotP[i+1]
}
TotalPatchDat$CotPAdj=CotPAdj

TotalPatchDat=subset(TotalPatchDat, select=c("Patch","Year","PatchYear","PlaL","Ver","PLcat","VScat","PLVCcatsum","MelC","PhoP","CotPAdj"))
names(TotalPatchDat)=c("Patch","Year","PatchYear","PLM2","VSM2","PL","VS","PLVS","MelC","PhoP","CotP")
head(TotalPatchDat)

TotalPatchDat$CotP[which(TotalPatchDat$MelC==0)]=NA  #Set Cotesia to NA when M. cinxia is absent
TotalPatchDat$CotP[which(is.na(TotalPatchDat$MelC))]=NA #Set Cotesia to NA when M. cinxia is NA

TotalPatchDat$State=paste(as.numeric(TotalPatchDat$MelC>0),as.numeric(TotalPatchDat$CotP>0),as.numeric(TotalPatchDat$PhoP>0),sep="_")
head(TotalPatchDat,20)
summary(TotalPatchDat)

#tapply(TotalPatchDat$CotP>0,TotalPatchDat$Year,sum,na.rm=T)
#tapply(TotalPatchDat$CotP>-1,TotalPatchDat$Year,sum,na.rm=T)

##################################################
#### - Wide-format summary for patch states - ####
##################################################
wide=dcast(TotalPatchDat,Year~State)

wide=wide[-c(19,20),colnames(wide)%in%c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1")]
wide[c(10),]=NA
rownames(wide)=2000:2017
wide

#write.csv(wide, file = "wide.csv", row.names=F)
#write.csv(TotalPatchDat,"TotalPatchDat.csv",row.names=F)


############################################
##### - Compute connectivity measures - ####
############################################

pdat=read.csv("patch_precipitation.csv")
xy=ddply(pdat,.(Patch),summarize,
      X=mean(Centroid_X,na.rm=T),
      Y=mean(Centroid_Y,na.rm=T))
xy[which(xy[,2]<100),2]=NA
xy[which(xy[,3]<100),3]=NA

patches=ddply(oldindat,.(Patch),summarize,
              Lat=mean(Latitude,na.rm=T),
              Lon=mean(Longitude,na.rm=T),
              AREA=mean(Area,na.rm=T),
              Road=max(BorderRoad,na.rm=T))
patches$Road[which(patches$Road<0)]=0
patches$Road=as.numeric(patches$Road>0)

patches=merge(patches,xy,by="Patch",all.x=T)
patches$X=patches$X/1000
patches$Y=patches$Y/1000
head(patches)

summary(patches)
#write.csv(patches,file="patch_data.csv",row.names=F)

#Butterfly
npatch=length(patches$Patch)
years=sort(unique(AF$Year))
years=2000:2018
outdat=matrix(NA,nrow=npatch,ncol=1)
outdat=as.data.frame(outdat)
outdat[,1]=patches$Patch[1:npatch]
colnames(outdat)="Patch"

a=Sys.time()
for(y in 1:length(years)){  #Years
  year=years[y]
  print(paste("Now computing year",year))
  
  sdat=AF[AF$Year==year,]
  ssdat=merge(patches,sdat,by="Patch",all.x=T)
  
  patches_y=ddply(ssdat,.(Patch),summarize,
                  MelC=Nest.Total>0,
                  MelCN=Nest.Total,
                  Lat=mean(Y,na.rm=T),
                  Lon=mean(X,na.rm=T),
                  AREA=mean(AREA,na.rm=T),
                  Host=mean(PLM2,na.rm=T))
 
  dmat=as.matrix(dist(patches_y[,4:5],diag=T,upper=T))
  Nvec=sqrt(patches_y$AREA)*as.numeric(patches_y$MelC)
  edist=exp(-dmat)
  p2=edist*Nvec
  
  aa=apply(p2,2,sum,na.rm=T)
  for(i in 1:length(aa)){
    aa[i]=aa[i]-p2[i,i]
  }
  
  out=data.frame(Patch=patches_y$Patch,w=aa)
  colnames(out)=c("Patch",paste(year))
  outdat=merge(outdat,out,by="Patch",all.x=T)
}
b=Sys.time()
b-a

head(outdat)
outdat2=outdat[order(as.numeric(outdat$Patch)),]
outdat2[which(outdat2[,2]==0),2:20]=NA

write.csv(outdat2,"NEWneighbourhooddata.csv",row.names=F)

#Total neighbourhood
hostsize=data.frame(tapply(TotalPatchDat$PLM2, TotalPatchDat$Patch, mean, na.rm=T))
hostsize$Patch=rownames(hostsize)
names(hostsize)=c("HostCover","Patch")

patches2=merge(patches, hostsize, by="Patch",all.x=T)
head(patches2)

wj=NULL
w=NULL
for(i in 1:5){
  print(paste("Patch",i))
  for(j in 1:nrow(patches2)){
    
    d=dist(patches2[c(i,j),6:7])
    wj[j]=exp(-d)*sqrt(patches2$HostCover[j])
  }   
  w[i]=sum(wj[-i],na.rm=T)
}

head(w)

dmat=as.matrix(dist(patches2[,6:7],diag=T,upper=T))
Nvec=as.vector(sqrt(patches2$HostCover))
edist=exp(-dmat)
p2=edist*Nvec

aa=apply(p2,2,sum,na.rm=T)
for(i in 1:length(aa)){
  aa[i]=aa[i]-p2[i,i]
}

out=data.frame(Patch=patches2$Patch,TotalNH=aa)
head(out)
#write.csv(out,"TotalNH.csv",row.names=F)

####Parasitoid
npatch=length(patches$Patch)
years=sort(unique(TotalPatchDat$Year))
years=2000:2018
outdat=matrix(NA,nrow=npatch,ncol=1)
outdat=as.data.frame(outdat)
outdat[,1]=patches$Patch[1:npatch]
colnames(outdat)="Patch"

a=Sys.time()
for(y in 1:length(years)){  #Years
  year=years[y]
  print(paste("Now computing year",year))
  
  sdat=TotalPatchDat[TotalPatchDat$Year==year,]
  ssdat=merge(patches,sdat,by="Patch",all.x=T)
  
  patches_y=ddply(ssdat,.(Patch),summarize,
                  CotP=CotP>0,
                  CotPM=CotP,
                  Lat=mean(Y,na.rm=T),
                  Lon=mean(X,na.rm=T),
                  AREA=mean(AREA,na.rm=T),
                  Host=mean(PLM2,na.rm=T))
  head(patches_y)
  
  dmat=as.matrix(dist(patches_y[,4:5],diag=T,upper=T))
  
  Nvec=sqrt(patches_y$AREA)*as.numeric(patches_y$CotP)
  Nvec[which(is.na(Nvec))]=0
  
  edist=exp(-dmat)
  
  p2=edist*Nvec
  
  aa=apply(p2,2,sum,na.rm=T)
  for(i in 1:length(aa)){
    aa[i]=aa[i]-p2[i,i]
  }
  
  out=data.frame(Patch=patches_y$Patch,w=aa)
  colnames(out)=c("Patch",paste(year))
  outdat=merge(outdat,out,by="Patch",all.x=T)
}
b=Sys.time()
b-a

head(outdat)
outdat2=outdat[order(as.numeric(outdat$Patch)),]
outdat2[which(outdat2[,2]==0),2:20]=NA
outdat2$`2009`=NA
head(outdat2)

#write.csv(outdat2,"NEWCotesia_neighbourhooddata.csv",row.names=F)

#Mildew

head(TotalPatchDat)

#Insert patch mean for missing PL data
pmeanPL=data.frame(tapply(TotalPatchDat$PLM2,TotalPatchDat$Patch,mean,na.rm=T))
pmeanPL$Patch=row.names(pmeanPL)
names(pmeanPL)=c("PL","Patch")
newPL=subset(TotalPatchDat,select=c("Patch","PLM2"))
head(newPL,20)
for(i in 1:nrow(newPL)){
  if(is.na(newPL$PLM2[i])){
    newPL$PLM2[i]=pmeanPL$PL[which(pmeanPL$Patch==newPL$Patch[i])]
  }
}
head(newPL,20)
TotalPatchDat$PLM2=newPL$PLM2

npatch=length(patches$Patch)
years=sort(unique(TotalPatchDat$Year))
outdat=matrix(NA,nrow=npatch,ncol=1)
outdat=as.data.frame(outdat)
outdat[,1]=patches$Patch[1:npatch]
colnames(outdat)="Patch"

a=Sys.time()
for(y in 1:length(years)){  #Years
  year=years[y]
  print(paste("Now computing year",year))
  
  sdat=TotalPatchDat[TotalPatchDat$Year==year,]
  ssdat=merge(patches,sdat,by="Patch",all.x=T)
  patches_y=ddply(ssdat,.(Patch),summarize,
                  PhoP=PhoP>0,
                  Lat=mean(Y,na.rm=T),
                  Lon=mean(X,na.rm=T),
                  AREA=mean(AREA,na.rm=T),
                  Host=mean(PLM2,na.rm=T))
  
  #patches2=patches[order(patches$Patch),]
  #patches_y2=patches_y[order(patches_y$Patch),]
  #which(patches2$Patch!=patches_y2$Patch)
  #which(patches2$Patch!=patches_y$Patch)
  
  dmat=as.matrix(dist(patches_y[,3:4],diag=T,upper=T))
  
  Nvec=sqrt(patches_y$AREA)*as.numeric(patches_y$PhoP)
  #Nvec=sqrt(patches_y$Host)*as.numeric(patches_y$PhoP)
  
  edist=exp(-dmat)
  p2=edist*Nvec
  
  aa=apply(p2,2,sum,na.rm=T)
  for(i in 1:length(aa)){
    aa[i]=aa[i]-p2[i,i]
  }
  
  out=data.frame(Patch=patches_y$Patch[1:npatch],w=aa)
  colnames(out)=c("Patch",paste(year))
  outdat=merge(outdat,out,by="Patch",all.x=T)
}
#colnames(outdat)=c("Patch",as.character(years))
b=Sys.time()
b-a

outdat2=outdat[order(as.numeric(outdat$Patch)),]
outdat2[which(outdat2[,2]==0),2:20]=NA

#write.csv(outdat2,"Mildew_neighbourhooddata.csv",row.names=F)

outdat3=outdat[order(as.numeric(outdat$Patch)),]
outdat3[which(outdat3[,2]==0),2:20]=NA

#write.csv(outdat3,"Mildew_neighbourhooddata_PatchArea.csv",row.names=F)


##################################################
#### Observed transition matrices in two ways ####
##################################################

head(TotalPatchDat,5)
#TotalPatchDat=TotalPatchDat[-which(is.na(TotalPatchDat$Patch)),]

TotalPatchDat$Year=factor(TotalPatchDat$Year)
TotalPatchDat$Patch=factor(TotalPatchDat$Patch)
table(TotalPatchDat$Year)
table(as.numeric(TotalPatchDat$Year))
TotalPatchDat$Year2=as.numeric(TotalPatchDat$Year)

out=list()
patches=levels(TotalPatchDat$Patch)
for(p in 1:length(patches)){
pa=TotalPatchDat[TotalPatchDat$Patch==patches[p],]  

MelClast=PhoPlast=CotPlast=Statelast=NULL
MelClast[1]=PhoPlast[1]=CotPlast[1]=Statelast[1]=NA
for(i in 2:nrow(pa)){
  MelClast[i]=ifelse(length(which(pa$Year2==pa$Year2[i]-1))==1,pa$MelC[which(pa$Year2==pa$Year2[i]-1)],NA)
  PhoPlast[i]=ifelse(length(which(pa$Year2==pa$Year2[i]-1))==1,pa$PhoP[which(pa$Year2==pa$Year2[i]-1)],NA)
  CotPlast[i]=ifelse(length(which(pa$Year2==pa$Year2[i]-1))==1,pa$CotP[which(pa$Year2==pa$Year2[i]-1)],NA)
  Statelast[i]=ifelse(length(which(pa$Year2==pa$Year2[i]-1))==1,pa$State[which(pa$Year2==pa$Year2[i]-1)],NA)
}
out[[p]]=data.frame(pa,MelClast=MelClast,PhoPlast=PhoPlast,CotPlast=CotPlast,Statelast=Statelast)
}

newdat=rbind.fill(out)
head(newdat,20)

d=newdat[which(newdat$CotPlast>0 & newdat$MelClast>0 & newdat$MelC==0),]

###Observed state transitions
s000=newdat[newdat$Statelast=="0_NA_0",]
sum=table(s000$State)["0_NA_0"]+table(s000$State)["0_NA_1"]+table(s000$State)["1_0_0"]+table(s000$State)["1_1_0"]+table(s000$State)["1_0_1"]+table(s000$State)["1_1_1"]
sum
o000_000=table(s000$State)["0_NA_0"]/sum
o000_001=table(s000$State)["0_NA_1"]/sum
o000_100=table(s000$State)["1_0_0"]/sum
o000_110=table(s000$State)["1_1_0"]/sum
o000_101=table(s000$State)["1_0_1"]/sum
o000_111=table(s000$State)["1_1_1"]/sum

s001=newdat[newdat$Statelast=="0_NA_1",] #NO 1-1-1, excluded from sum
sum=table(s001$State)["0_NA_0"]+table(s001$State)["0_NA_1"]+table(s001$State)["1_0_0"]+table(s001$State)["1_1_0"]+table(s001$State)["1_0_1"]
sum
o001_000=table(s001$State)["0_NA_0"]/sum
o001_001=table(s001$State)["0_NA_1"]/sum
o001_100=table(s001$State)["1_0_0"]/sum
o001_110=table(s001$State)["1_1_0"]/sum
o001_101=table(s001$State)["1_0_1"]/sum
o001_111=table(s001$State)["1_1_1"]/sum
o001_111=0

s100=newdat[newdat$Statelast=="1_0_0",]
sum=table(s100$State)["0_NA_0"]+table(s100$State)["0_NA_1"]+table(s100$State)["1_0_0"]+table(s100$State)["1_1_0"]+table(s100$State)["1_0_1"]+table(s100$State)["1_1_1"]
sum
o100_000=table(s100$State)["0_NA_0"]/sum
o100_001=table(s100$State)["0_NA_1"]/sum
o100_100=table(s100$State)["1_0_0"]/sum
o100_110=table(s100$State)["1_1_0"]/sum
o100_101=table(s100$State)["1_0_1"]/sum
o100_111=table(s100$State)["1_1_1"]/sum

s110=newdat[newdat$Statelast=="1_1_0",]
sum=table(s110$State)["0_NA_0"]+table(s110$State)["0_NA_1"]+table(s110$State)["1_0_0"]+table(s110$State)["1_1_0"]+table(s110$State)["1_0_1"]+table(s110$State)["1_1_1"]
sum
o110_000=table(s110$State)["0_NA_0"]/sum
o110_001=table(s110$State)["0_NA_1"]/sum
o110_100=table(s110$State)["1_0_0"]/sum
o110_110=table(s110$State)["1_1_0"]/sum
o110_101=table(s110$State)["1_0_1"]/sum
o110_111=table(s110$State)["1_1_1"]/sum

s101=newdat[newdat$Statelast=="1_0_1",]
sum=table(s101$State)["0_NA_0"]+table(s101$State)["0_NA_1"]+table(s101$State)["1_0_0"]+table(s101$State)["1_1_0"]+table(s101$State)["1_0_1"]+table(s101$State)["1_1_1"]
sum
o101_000=table(s101$State)["0_NA_0"]/sum
o101_001=table(s101$State)["0_NA_1"]/sum
o101_100=table(s101$State)["1_0_0"]/sum
o101_110=table(s101$State)["1_1_0"]/sum
o101_101=table(s101$State)["1_0_1"]/sum
o101_111=table(s101$State)["1_1_1"]/sum

s111=newdat[newdat$Statelast=="1_1_1",]#No 0_NA_1, excluded from sum
sum=table(s111$State)["0_NA_0"]+table(s111$State)["1_0_0"]+table(s111$State)["1_1_0"]+table(s111$State)["1_0_1"]+table(s111$State)["1_1_1"]
sum
o111_000=table(s111$State)["0_NA_0"]/sum
o111_001=table(s111$State)["0_NA_1"]/sum
o111_001=0
o111_100=table(s111$State)["1_0_0"]/sum
o111_110=table(s111$State)["1_1_0"]/sum
o111_101=table(s111$State)["1_0_1"]/sum
o111_111=table(s111$State)["1_1_1"]/sum


vals=c(o000_000,o000_001,o000_100,o000_110,o000_101,o000_111,
       o001_000,o001_001,o001_100,o001_110,o001_101,o001_111,
       o100_000,o100_001,o100_100,o100_110,o100_101,o100_111,
       o110_000,o110_001,o110_100,o110_110,o110_101,o110_111,
       o101_000,o101_001,o101_100,o101_110,o101_101,o101_111,
       o111_000,o111_001,o111_100,o111_110,o111_101,o111_111)

omat=matrix(vals,nrow=6,byrow=T)
rowSums(omat)
colnames(omat)=rownames(omat)=c("000","001","100","110","101","111")
omat

##############################################
####Defining colonisations and extinctions####
##############################################

###Cinxia col ext
diff=(as.numeric(newdat$MelC>0)-as.numeric(newdat$MelClast>0))
diff[1:10]

col=as.numeric(diff==1) #col
col[which(newdat$MelClast>0)]=NA

ext=as.numeric(diff==-1) #ext
ext[which(newdat$MelClast==0)]=NA

newdat$MelCc=col
newdat$MelCe=ext

###Mildew col ext
diff=(as.numeric(newdat$PhoP>0)-as.numeric(newdat$PhoPlast>0))

col=as.numeric(diff==1) #col
col[which(newdat$PhoPlast>0)]=NA

ext=as.numeric(diff==-1) #ext
ext[which(newdat$PhoPlast==0)]=NA

newdat$PhoPc=col
newdat$PhoPe=ext

###Cotesia col ext
diff=(as.numeric(newdat$CotP>0)-as.numeric(newdat$CotPlast>0))
diff[1:100]

col=as.numeric(diff==1) #col
col[which(newdat$CotPlast>0)]=NA
col[which(newdat$MelC<1)]=NA #Cinxia extinctions

ext=as.numeric(diff==-1) #ext
ext[which(newdat$CotPlast==0)]=NA

newdat$CotPc=col
newdat$CotPe=ext

#Special cases

#Joint colonisations
newdat[which(newdat$MelCc>0 & newdat$CotP>0),]
newdat$CotPc[which(newdat$MelCc>0 & newdat$CotP>0)]=1 #Joint colonisations

#Butterfly colonisation but no Cotesia colonisation
sss=newdat[which(newdat$MelCc>0 & newdat$CotP==0),]
head(sss,20)
newdat$CotPc[which(newdat$MelCc>0 & newdat$CotP==0)]=0 #Butterfly colonisation but no Cotesia colonisation

#Cotesia extinction because of butterfly extinction
newdat$CotPe[which(newdat$MelCe>0 & newdat$CotPlast>0)]=1


ss=subset(newdat,select=c("Patch","Year","CotP","CotPlast","CotPc","CotPe","MelC","MelClast"))
head(ss,20)
#ss[which(ss$CotPlast>0),]

#write.csv(newdat,"colextdat.csv",row.names=F)

######################
####State with lag####
######################
head(newdat)
newdat$stateLag=paste(as.numeric(newdat$MelC>0),as.numeric(newdat$CotP>0),as.numeric(newdat$PhoPlast>0),sep="_")
wideLag=dcast(newdat,Year~stateLag)
wideLag=wideLag[-c(19),colnames(wideLag)%in%c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1")]
wideLag[c(10),]=NA
rownames(wideLag)=2000:2017
wideLag

write.csv(wideLag, file = "wideLag.csv", row.names=F)


####
head(newdat,10)
ss=(subset(newdat,select=c("Patch","Year","MelC","CotP","CotPlast","CotPc","CotPe")))
ss2=ss[ss$MelC>0,]
head(ss2,20)

ss=(subset(newdat,select=c("Patch","Year","MelClast","MelC","MelCc","MelCe")))
ss=(subset(newdat,select=c("Patch","Year","PhoPlast","PhoP","PhoPc","PhoPe")))
head(ss)


c1=sum(newdat$MelCc>0,na.rm=T)/sum(newdat$MelCc>-1,na.rm=T)
e1=sum(newdat$MelCe>0,na.rm=T)/sum(newdat$MelCe>-1,na.rm=T)

c2=sum(newdat$CotPc>0,na.rm=T)/sum(newdat$CotPc>-1,na.rm=T)
e2=sum(newdat$CotPe>0,na.rm=T)/sum(newdat$CotPe>-1,na.rm=T)

c3=sum(newdat$PhoPc>0,na.rm=T)/sum(newdat$PhoPc>-1,na.rm=T)
e3=sum(newdat$PhoPe>0,na.rm=T)/sum(newdat$PhoPe>-1,na.rm=T)

c1
e1
c2
e2
c3
e3

p000_000=(1-c1)*(1-c3) #No Cot
p000_001=(1-c1)*(c3) #No Cot
p000_100=c1*(1-c2)*(1-c3)
p000_110=c1*c2*(1-c3)
p000_101=c1*(1-c2)*c3
p000_111=c1*c2*c3

p001_000=(1-c1)*e3
p001_001=(1-c1)*(1-e3)
p001_100=c1*(1-c2)*e3
p001_110=c1*c2*e3
p001_101=c1*(1-c2)*(1-e3)
p001_111=c1*c2*(1-e3)

p100_000=e1*(1-c3)
p100_001=e1*c3
p100_100=(1-e1)*(1-c2)*(1-c3)
p100_110=(1-e1)*c2*(1-c3)
p100_101=(1-e1)*(1-c2)*c3
p100_111=(1-e1)*c2*c3

p110_000=e1*(1-c3)
p110_001=e1*c3
p110_100=(1-e1)*e2*(1-c3)
p110_110=(1-e1)*(1-e2)*(1-c3)
p110_101=(1-e1)*e2*c3
p110_111=(1-e1)*(1-e2)*c3

p101_000=e1*e3
p101_001=e1*(1-e3)
p101_100=(1-e1)*(1-c2)*e3
p101_110=(1-e1)*c2*e3
p101_101=(1-e1)*(1-c2)*(1-e3)
p101_111=(1-e1)*c2*(1-e3)

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
par(mfrow=c(1,2))

corrplot(mat,is.corr=F,method="color",tl.col="black",title="Expected from c and e")
corrplot(omat,is.corr=F,method="color",tl.col="black",title="Observed")

plot(mat[-c(5,6),],omat[-c(5,6),])
lines(0:1,0:1)
cor(c(mat[-c(5,6),]),c(omat[-c(5,6),]))





##############################
#### - Yearly data SiNs - ####
##############################
str(AF)
AF$SIN1=as.factor(AF$SIN1)
yearlySIN=ddply(AF,.(Year),summarize,
                nPatch=sum(tapply(Nest.Total>(-1),SIN1,sum,na.rm=T),na.rm=T),
                nSIN=sum(tapply(Nest.Total>(-1),SIN1,sum,na.rm=T)>-1,na.rm=T),
                nPopCinxia=sum(tapply(Nest.Total>0,SIN1,sum,na.rm=T),na.rm=T),
                mPLM2=mean(tapply(PLM2,SIN1,sum,na.rm=T),na.rm=T),
                mVSM2=mean(tapply(VSM2,SIN1,mean,na.rm=T),na.rm=T),
                mPL=mean(tapply(PL,SIN1,mean,na.rm=T),na.rm=T),
                mVS=mean(tapply(VS,SIN1,mean,na.rm=T),na.rm=T),
                MelC=sum(tapply(Nest.Total>0,SIN1,sum,na.rm=T)>0,na.rm=T)/sum(tapply(Nest.Total>(-1),SIN1,sum,na.rm=T)>-1,na.rm=T),
                PhoP=sum(tapply(Mildew2>0,SIN1,sum,na.rm=T)>0,na.rm=T)/sum(tapply(Mildew2>(-1),SIN1,sum,na.rm=T)>-1,na.rm=T))

#Mildew
yearlyMFSIN=ddply(MF,.(Year),summarize,
                  nMildew=sum(tapply(Absabund>-1,SIN1,sum,na.rm=T)>-1,na.rm=T),
                  #CotSum=sum(tapply(Total.Cotesia.Nests>0,SIN1,sum,na.rm=T)>0,na.rm=T),
                  PhoPAMF=sum(tapply(Absabund>0,SIN1,sum,na.rm=T)>0,na.rm=T)/sum(tapply(Absabund>-1,SIN1,sum,na.rm=T)>-1,na.rm=T))
yearlyMFSIN

#Insert 2016 data from MF survey
yearlySIN$PhoP[17]=yearlyMFSIN$PhoPAMF[5]

#Cotesia
yearlyCSIN=ddply(CS,.(Year),summarize,
                 nCot=sum(tapply(Total.Cotesia.Nests>-1,SIN1,sum,na.rm=T)>-1,na.rm=T),
                 CotPsum=sum(tapply(Total.Cotesia.Nests>0,SIN1,sum,na.rm=T)>0,na.rm=T))
#CotP=sum(tapply(Total.Cotesia.Nests>0,SIN1,sum,na.rm=T)>0,na.rm=T)/sum(tapply(Total.Cotesia.Nests>-1,SIN1,sum,na.rm=T)>-1,na.rm=T))

yearlySINDat=merge(yearlySIN,yearlyCSIN,by="Year",all=T)

CotPAdj=NULL
for(i in 1:nrow(yearlySINDat)){
  CotPAdj[i]=yearlySINDat$CotPsum[i+1]
}
yearlySINDat$CotPAdj=CotPAdj

yearlySINDat$CotP=yearlySINDat$CotPAdj/yearlySINDat$nSIN


x11()
par(mfrow=c(1,1),mar=c(4,4,2,4))
plot(years,yearlyDat$MelC,type="l",lwd=2,pch=16,ylim=c(0,.3),xlab="",ylab="Proportion of patches occupied",las=1)
lines(years,yearlyDat$PhoPcomb,lwd=2,pch=16,col="red")
lines(years,yearlyDat$CotPAdj,lwd=2,pch=16,col="blue")
legend("topleft",lwd=2,bty="n",col=c("black","red","blue","darkgreen"),legend=c("M. cinxia", "Mildew", "Cotesia", "Plantago"))

par(new=T)
plot(years,yearlyDat$mPLM2,type="l",lwd=2,col="darkgreen",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4,c(5,10,15,20,25),c(5,10,15,20,25),las=1)
mtext(4,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)


#years=c(2000:2016,2018)
plot(years,yearlySINDat$MelC,type="l",lwd=2,pch=16,ylim=c(0,1),xlab="",ylab="Proportion of networks occupied",las=1)
lines(years,yearlySINDat$PhoP,lwd=2,pch=16,col="red")
lines(years,yearlySINDat$CotP,lwd=2,pch=16,col="darkgreen")
legend("topleft",lwd=2,bty="n",col=c("black","red","darkgreen"),legend=c("M. cinxia", "Mildew", "Cotesia"))




###More old crap
##################################
#### - Random distributions - ####
##################################
n_mel=apply(wide,1,function(x)sum(x[3:6]))
n_cot=apply(wide,1,function(x)sum(x[5:6]))
n_mil=apply(wide,1,function(x)sum(x[c(2,4,6)]))
nn=apply(wide,1,function(x)sum(x))
cbind(nn,n_mel,n_cot,n_mil)

dist=matrix(NA,ncol=6,nrow=1000)
medians=matrix(NA,ncol=6,nrow=18)
upper=matrix(NA,ncol=6,nrow=18)
lower=matrix(NA,ncol=6,nrow=18)

s=Sys.time() #1000 in 2.5min
for(i in c(1:9,11:18)){
  for(j in 1:1000){
    npa=nn[i]
    mel_pa=sample(1:npa,n_mel[i])
    cot_pa=sample(mel_pa,size=n_cot[i])
    mil_pa=sample(1:npa,n_mil[i])
    df=data.frame(patch=1:npa,mel=rep(0,npa),cot=rep(0,npa),mil=rep(0,npa))
    df$mel[mel_pa]=1
    df$cot[cot_pa]=1
    df$cot[which(df$mel==0)]=NA
    df$mil[mil_pa]=1
    df$state=paste(as.numeric(df$mel>0),as.numeric(df$cot>0),as.numeric(df$mil>0),sep="_")
    
    df$state=factor(df$state)
    levels(df$state)=colnames(wide)
    
    dist[j,]=table(df$state)
  }
  
  a=apply(dist,2,quantile,c(.025,.5,.975))
  medians[i,]=a[2,]
  upper[i,]=a[3,]
  lower[i,]=a[1,]
}
e=Sys.time()
e-s


xx=2000:2008
xx2=2010:2017
x11()
par(mfrow=c(2,3))
for(i in 1:6){
  #polygon(c(xx,rev(xx)),c(lower[1:9,i],rev(upper[1:9,i])),col = "grey75", border = FALSE)
  plot(xx, medians[1:9,i], type = "n",ylab="Count",xlab="",xlim=c(2000,2017),
       ylim=c(min(c(lower[,i],wide[,i]),na.rm=T),max(c(upper[,i],wide[,i]),na.rm=T)),main=colnames(wide)[i])
  polygon(c(xx,rev(xx)),c(lower[1:9,i],rev(upper[1:9,i])),col = "grey75", border = FALSE)
  polygon(c(xx2,rev(xx2)),c(lower[11:18,i],rev(upper[11:18,i])),col = "grey75", border = FALSE)
  
  #lines(xx,medians[1:9,i], lwd = 2)
  lines(xx, upper[1:9,i], col="red",lty=2)
  lines(xx, lower[1:9,i], col="red",lty=2)
  
  lines(xx2, upper[11:18,i], col="red",lty=2)
  lines(xx2, lower[11:18,i], col="red",lty=2)
  
  lines(c(2000:2017),wide[,i],lwd=2)
  
}

#################
wide2=t(apply(wide,1,function(x){x/sum(x)}))
medians2=t(apply(medians,1,function(x){x/sum(x)}))
upper2=t(apply(upper,1,function(x){x/sum(x)}))
lower2=t(apply(lower,1,function(x){x/sum(x)}))

xx=2000:2008
xx2=2010:2017
x11()
par(mfrow=c(2,3))
for(i in 1:6){
  #polygon(c(xx,rev(xx)),c(lower[1:9,i],rev(upper[1:9,i])),col = "grey75", border = FALSE)
  plot(xx, medians2[1:9,i], type = "n",ylab="Proportion of patches",xlab="",xlim=c(2000,2017),
       ylim=c(min(c(lower2[,i],wide2[,i]),na.rm=T),max(c(upper2[,i],wide2[,i]),na.rm=T)),main=colnames(wide2)[i])
  polygon(c(xx,rev(xx)),c(lower2[1:9,i],rev(upper2[1:9,i])),col = "grey75", border = FALSE)
  polygon(c(xx2,rev(xx2)),c(lower2[11:18,i],rev(upper2[11:18,i])),col = "grey75", border = FALSE)
  
  #lines(xx,medians[1:9,i], lwd = 2)
  lines(xx, upper2[1:9,i], col="red",lty=2)
  lines(xx, lower2[1:9,i], col="red",lty=2)
  
  lines(xx2, upper2[11:18,i], col="red",lty=2)
  lines(xx2, lower2[11:18,i], col="red",lty=2)
  
  lines(c(2000:2017),wide2[,i],lwd=2)
  
}

#####################################
#### - SINlevel data formating - ####
#####################################
head(AF)
SINDat=ddply(AF,.(SIN1,Year),summarize,
             PlaL=mean(PLM2,na.rm=T),
             Ver=mean(VSM2,na.rm=T),
             PLcat=mean(PL,na.rm=T),
             VScat=mean(VS,na.rm=T),
             PLVSsum=mean(PLM2+VSM2,na.rm=T),
             PLVCcatsum=mean(PL+VS,na.rm=T),
             MelC=sum(Nest.Total,na.rm=T),
             MelC_occ=sum(Nest.Total>0, na.rm=T)/sum(Nest.Total>-Inf, na.rm=T),
             PhoP=sum(Mildew2,na.rm=T),
             PhoP_occ=sum(Mildew2>0, na.rm=T)/sum(Mildew2>-Inf, na.rm=T))
SINDat$SINYear=paste(SINDat$SIN1,SINDat$Year,sep="_")
SINDat$SINYear=as.factor(SINDat$SINYear)
head(SINDat,10)

MF$SIN1=as.factor(MF$SIN1)
MildewSINDat=ddply(MF,.(SIN1,Year),summarize,
                   PhoPAMF=sum(Absabund,na.rm=T),
                   PhoPAMF_occ=sum(Absabund>0,na.rm=T)/sum(Absabund>-Inf,na.rm=T))
MildewSINDat$SINYear=paste(MildewSINDat$SIN,MildewSINDat$Year,sep="_")
MildewSINDat$SINYear=as.factor(MildewSINDat$SINYear)
head(MildewSINDat,10)

SINDat=merge(SINDat,MildewSINDat[,3:5],by="SINYear",all.x=T)

#Insert MF survey values for 2016 Mildew
summary(SINDat$PhoP[SINDat$Year=="2016"])
SINDat$PhoP[which(SINDat$Year=="2016")]=SINDat$PhoPAMF[which(SINDat$Year=="2016")]
SINDat$PhoP_occ[which(SINDat$Year=="2016")]=SINDat$PhoPAMF_occ[which(SINDat$Year=="2016")]

CS2=CS
CS2$Total.Cotesia.Nests[which(CS2$Total.Cotesia.Nests>998)]=NA

CotesiaSINDat=ddply(CS2,.(SIN1,Year),summarize,
                    CotP=sum(Total.Cotesia.Nests,na.rm=T))
CotesiaSINDat$SINYear=paste(CotesiaSINDat$SIN1,CotesiaSINDat$Year,sep="_")
head(CotesiaSINDat,10)

TotalSINDat=merge(SINDat,CotesiaSINDat,by="SINYear",all.x=T)
head(TotalSINDat,10)

TotalSINDat$CotP[which(is.na(TotalSINDat$CotP))]=0

CotPAdj=NULL
for(i in 1:nrow(TotalSINDat)){
  CotPAdj[i]=TotalSINDat$CotP[i+1]
}
TotalSINDat$CotPAdj=CotPAdj
names(TotalSINDat)
TotalSINDat=subset(TotalSINDat, select=c("SIN1.x","Year.x","SINYear","MelC_occ","PhoP_occ","MelC","PhoP","PhoPAMF","PlaL"))
names(TotalSINDat)=c("SIN","Year","SINYear","MelC_occ","PhoP_occ","MelC","PhoP","PhoPAMF","PlaL")
head(TotalSINDat)

TotalSINDat$PhoPmax=apply(subset(TotalSINDat, select=c("PhoP", "PhoPAMF")), 1, max, na.rm=T)

#### - Plots for A-L - ####

SINmeans=ddply(TotalSINDat, .(SIN), summarize,
               MelC_occ=mean(MelC_occ,na.rm=T),
               PhoP_occ=mean(PhoP_occ,na.rm=T),
               MelC=mean(MelC, na.rm=T),
               PhoP=mean(PhoPmax, na.rm=T),
               PlaL=mean(PlaL, na.rm=T))
head(SINmeans)
#SINmeans$SIN=as.numeric(SINmeans$SIN)

pdf("SIN_occurrence.pdf",height=4,width=7)
par(mfrow=c(1,2))
SINmeans2=SINmeans[order(SINmeans$MelC_occ),]
plot(SINmeans2$MelC_occ,ylab="Occurrence (proportion patches)", xlab="SIN",ylim=c(0,1))
points(SINmeans2$PhoP_occ, pch=16)
legend("topleft",pch=c(1,16),legend=c("Butterfly", "Mildew"),cex=.7)

SINmeans2=SINmeans[order(SINmeans$PhoP_occ),]
plot(SINmeans2$PhoP_occ,ylab="Occurrence (proportion patches)", xlab="SIN", pch=16, ylim=c(0,1))
points(SINmeans2$MelC_occ, pch=1)
#legend("topleft",pch=c(1,16),legend=c("Butterfly", "Mildew"))
dev.off()

years=unique(TotalSINDat$Year)

pdf("SIN_occurrence_yearly.pdf",height=4,width=7)
for(i in 1:length(years)){
year=years[i]
SIN_Y=TotalSINDat[TotalSINDat$Year==year,]
SINmeans2=SIN_Y[order(SIN_Y$MelC_occ),]

par(mfrow=c(1,2))
plot(SINmeans2$MelC_occ,ylab="Occurrence (proportion patches)", xlab="SIN",main=paste("Year = ", year),ylim=c(0,1))
points(SINmeans2$PhoP_occ, pch=16)
legend("topleft",pch=c(1,16),legend=c("Butterfly", "Mildew"),cex=.7,bg="white")

SINmeans2=SIN_Y[order(SIN_Y$PhoP_occ),]
plot(SINmeans2$PhoP_occ,ylab="Occurrence (proportion patches)", xlab="SIN", pch=16,ylim=c(0,1))
points(SINmeans2$MelC_occ, pch=1)
#legend("topleft",pch=c(1,16),legend=c("Butterfly", "Mildew"))
}
dev.off()

#Occurrence frequency per year
i=1
pdf("SIN_occurrence_freq_yearly.pdf",height=4,width=7)
for(i in 1:length(years)){
  year=years[i]
  SIN_Y=TotalSINDat[TotalSINDat$Year==year,]
  SINmeans2=SIN_Y[order(SIN_Y$MelC_occ),]
  par(mfrow=c(1,2))
  
  hist(SIN_Y$MelC_occ, xlab="Occurrence (proportion patches)", main=paste("Butterfly", "Year= ", paste(year)), las=1)
  hist(SIN_Y$PhoP_occ, xlab="Occurrence (proportion patches)", main="Mildew", las=1)
    }
dev.off()

#Abundance frequency per year
i=1
pdf("SIN_abundance_freq_yearly.pdf",height=4,width=7)
for(i in 1:length(years)){
  year=years[i]
  SIN_Y=TotalSINDat[TotalSINDat$Year==year,]
  SINmeans2=SIN_Y[order(SIN_Y$MelC_occ),]
  par(mfrow=c(1,2))
  
  hist(SIN_Y$MelC, xlab="Abundance (number of nests)", main=paste("Butterfly", "Year= ", paste(year)), las=1)
  hist(SIN_Y$PhoPmax, xlab="Abundance (number of ...)", main="Mildew", las=1)
}
dev.off()


#Abundance per network per species
head(SINmeans)

pdf("SIN_abundance.pdf",height=3,width=9)
par(mfrow=c(1,3))
plot(SINmeans$MelC[order(SINmeans$MelC)],xlab="SIN",ylab="Mean abundance over years (nests)", main ="Butterfly")
plot(SINmeans$PhoP[order(SINmeans$MelC)],xlab="SIN",ylab="Mean abundance over years (...)", main="Mildew")
plot(SINmeans$PlaL[order(SINmeans$MelC)],xlab="SIN",ylab="Mean abundance over years (m2 per patch)", main="Plantago")
dev.off()

#Per year
pdf("SIN_abundance_yearly.pdf",height=3,width=9)
for(i in c(1:9,12:19)){
  year=years[i]
  SIN_Y=TotalSINDat[TotalSINDat$Year==year,]
  
  par(mfrow=c(1,3))
  plot(SIN_Y$MelC[order(SIN_Y$MelC)],ylab="Abundance (number of nests)", xlab="SIN")
  legend("topleft","Butterfly")
  plot(SIN_Y$PhoP[order(SIN_Y$MelC)],ylab="Abundance (number of ...)", xlab="SIN",main=paste("Year = ", year))
  legend("topleft","Mildew")
  plot(SIN_Y$PlaL[order(SIN_Y$MelC)],ylab="Abundance (m2 per patch)", xlab="SIN")
  legend("topleft","Plantago")
}
dev.off()


####################################################
####################################################

TotalSINDat$CotP[which(TotalSINDat$MelC==0)]=NA  #Set Cotesia to NA when M. cinxia is absent
TotalSINDat$CotP[which(is.na(TotalSINDat$MelC))]=NA #Set Cotesia to NA when M. cinxia is NA

TotalSINDat$State=paste(as.numeric(TotalSINDat$MelC>0),as.numeric(TotalSINDat$CotP>0),as.numeric(TotalSINDat$PhoP>0),sep="_")
head(TotalSINDat,10)



wide=dcast(TotalSINDat,Year~State)

wide=wide[-c(19),colnames(wide)%in%c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1")]
wide[10,]=NA
wide

x11()
par(mfrow=c(1,1),mar=c(4,4,2,8),xpd=T)
matplot(2000:2017,apply(wide,2,function(x){log(x+1)}),type="l",lty=1,lwd=2,xlab="Year",ylab="Prevalence",yaxt="n")
axis(2,at=c(0,2,4,6,8),labels=round(exp(c(0,2,4,6,8))),las=1)
legend(x=2014,y=8,col=1:6,lty=1,lwd=2,legend=colnames(wide))

##################################
#### - Random distributions - ####
##################################
n_mel=apply(wide,1,function(x)sum(x[3:6]))
n_cot=apply(wide,1,function(x)sum(x[5:6]))
n_mil=apply(wide,1,function(x)sum(x[c(2,4,6)]))
nn=apply(wide,1,function(x)sum(x))
cbind(nn,n_mel,n_cot,n_mil)

dist=matrix(NA,ncol=6,nrow=1000)
medians=matrix(NA,ncol=6,nrow=18)
upper=matrix(NA,ncol=6,nrow=18)
lower=matrix(NA,ncol=6,nrow=18)

s=Sys.time() #1000 in 2.5min
for(i in c(1:9,11:18)){
  for(j in 1:1000){
    npa=nn[i]
    mel_pa=sample(1:npa,n_mel[i])
    cot_pa=sample(mel_pa,size=n_cot[i])
    mil_pa=sample(1:npa,n_mil[i])
    df=data.frame(patch=1:npa,mel=rep(0,npa),cot=rep(0,npa),mil=rep(0,npa))
    df$mel[mel_pa]=1
    df$cot[cot_pa]=1
    df$cot[which(df$mel==0)]=NA
    df$mil[mil_pa]=1
    df$state=paste(as.numeric(df$mel>0),as.numeric(df$cot>0),as.numeric(df$mil>0),sep="_")
    
    df$state=factor(df$state)
    levels(df$state)=colnames(wide)
    
    dist[j,]=table(df$state)
  }
  
  a=apply(dist,2,quantile,c(.025,.5,.975))
  medians[i,]=a[2,]
  upper[i,]=a[3,]
  lower[i,]=a[1,]
}
e=Sys.time()
e-s


xx=2000:2008
xx2=2010:2017
x11()
par(mfrow=c(2,3))
for(i in 1:6){
  #polygon(c(xx,rev(xx)),c(lower[1:9,i],rev(upper[1:9,i])),col = "grey75", border = FALSE)
  plot(xx, medians[1:9,i], type = "n",ylab="Count",xlab="",xlim=c(2000,2017),
       ylim=c(min(c(lower[,i],wide[,i]),na.rm=T),max(c(upper[,i],wide[,i]),na.rm=T)),main=colnames(wide)[i])
  polygon(c(xx,rev(xx)),c(lower[1:9,i],rev(upper[1:9,i])),col = "grey75", border = FALSE)
  polygon(c(xx2,rev(xx2)),c(lower[11:18,i],rev(upper[11:18,i])),col = "grey75", border = FALSE)
  
  #lines(xx,medians[1:9,i], lwd = 2)
  lines(xx, upper[1:9,i], col="red",lty=2)
  lines(xx, lower[1:9,i], col="red",lty=2)
  
  lines(xx2, upper[11:18,i], col="red",lty=2)
  lines(xx2, lower[11:18,i], col="red",lty=2)
  
  lines(c(2000:2017),wide[,i],lwd=2)
  
}


par(mfrow=c(1,1),mar=c(4,4,2,8),xpd=T)
matplot(2000:2013,apply(wide,2,function(x){log(x+1)}),type="l",lty=1,lwd=2,xlab="Year",ylab="Prevalence",yaxt="n")
axis(2,at=c(0,2,4,6,8),labels=round(exp(c(0,2,4,6,8))),las=1)
legend(x=2014,y=8,col=1:6,lty=1,lwd=2,legend=colnames(wide))
legend(x=2008,y=300,col=1:6,lty=1,lwd=2,legend=colnames(wide))



#Neighbourhood with for-loops
npatch=length(patches$Patch)
a=Sys.time()

years=sort(unique(AF$Year))
years=2000:2018
outdat=matrix(NA,nrow=npatch,ncol=1)
outdat=as.data.frame(outdat)
outdat[,1]=patches$Patch[1:npatch]
colnames(outdat)="Patch"

for(y in 1:1){  #Years
  year=years[y+1]
  print(paste("Now computing year",year))
  
  sdat=AF[AF$Year==year,]
  ssdat=merge(patches,sdat,by="Patch",all.x=T)
  
  patches_y=ddply(ssdat,.(Patch),summarize,
                  MelC=Nest.Total>0,
                  MelCN=Nest.Total,
                  Lat=mean(Y,na.rm=T),
                  Lon=mean(X,na.rm=T),
                  AREA=mean(AREA,na.rm=T))
  
  sdatlast=AF[AF$Year==year-1,]
  ssdatlast=merge(patches,sdatlast,by="Patch",all.x=T)
  
  patches_last=ddply(ssdatlast,.(Patch),summarize,
                     MelClast=Nest.Total>0,
                     MelCNlast=Nest.Total)
  patches_y=merge(patches_y,patches_last,by="Patch",all.x=T)
  
  wj=NULL
  w=NULL
  for(i in 1:npatch){
    #print(paste("Patch",i))
    for(j in 1:nrow(patches_y)){
      
      d=dist(patches_y[c(i,j),4:5])
      wj[j]=exp(-d)*patches_y$MelCNlast[j]
    }   
    w[i]=sum(wj[-i],na.rm=T)
  }
  out=data.frame(Patch=patches_y$Patch[1:npatch],w=w)
  colnames(out)=c("Patch",paste(year))
  outdat=merge(outdat,out,by="Patch",all.x=T)
}
#colnames(outdat)=c("Patch",as.character(years))
b=Sys.time()
b-a #1000 patches 14 years 1.29 hours
#All patches all years 5.74hours 

outdat2=outdat[order(as.numeric(outdat$Patch)),]

outdat2[which(outdat2[,2]==0),2:20]=NA
plot(outdat2[,2],outdat2[,3])
lines(0:1000,0:1000)
hist(sqrt(outdat2[,2]))
hist(sqrt(outdat2[,3]))

cor(outdat2[,2],outdat2[,3],"pairwise")


#write.csv(outdat2,"neighbourhooddata.csv",row.names=F)

#Parasitoid
npatch=length(patches$Patch)
a=Sys.time()

years=sort(unique(TotalPatchDat$Year))
outdat=matrix(NA,nrow=npatch,ncol=1)
outdat=as.data.frame(outdat)
outdat[,1]=patches$Patch[1:npatch]
colnames(outdat)="Patch"

for(y in 3:length(years)){  #Years
  year=years[y]
  print(paste("Now computing year",year))
  
  sdat=TotalPatchDat[TotalPatchDat$Year==year,]
  ssdat=merge(patches,sdat,by="Patch",all.x=T)
  patches_y=ddply(ssdat,.(Patch),summarize,
                  CotPp=CotP>0,
                  CotPN=CotP,
                  Lat=mean(Y,na.rm=T),
                  Lon=mean(X,na.rm=T),
                  AREA=mean(AREA,na.rm=T))
  
  wj=NULL
  w=NULL
  for(i in 1:npatch){
    #print(paste("Patch",i))
    for(j in 1:nrow(patches_y)){
      
      d=dist(patches_y[c(i,j),4:5])
      wj[j]=exp(-d)*patches_y$CotPN[j]
    }   
    w[i]=sum(wj[-i],na.rm=T)
  }
  out=data.frame(Patch=patches_y$Patch[1:npatch],w=w)
  colnames(out)=c("Patch",paste(year))
  outdat=merge(outdat,out,by="Patch",all.x=T)
}
#colnames(outdat)=c("Patch",as.character(years))
b=Sys.time()
b-a #1000 patches 14 years 1.29 hours
#All patches all years 5.74hours 

outdat2=outdat[order(as.numeric(outdat$Patch)),]

outdat2[which(outdat2[,2]==0),2:20]=NA
summary(outdat2)

hist(log(outdat[,2]))
hist(log(outdat[,13]))

plot(outdat2[,2],outdat2[,3])
lines(0:1000,0:1000)
hist(sqrt(outdat2[,2]))
hist(sqrt(outdat2[,3]))

cor(outdat2[,2],outdat2[,3],"pairwise")


write.csv(outdat2,"Cotesia_neighbourhooddata.csv",row.names=F)


npatch=length(patches$Patch)

years=sort(unique(TotalPatchDat$Year))
outdat=matrix(NA,nrow=npatch,ncol=1)
outdat=as.data.frame(outdat)
outdat[,1]=patches$Patch[1:npatch]
colnames(outdat)="Patch"

a=Sys.time()
for(y in 1:1){  #Years
  year=years[y]
  print(paste("Now computing year",year))
  
  sdat=TotalPatchDat[TotalPatchDat$Year==year,]
  ssdat=merge(patches,sdat,by="Patch",all.x=T)
  patches_y=ddply(ssdat,.(Patch),summarize,
                  PhoP=PhoP>0,
                  Lat=mean(Y,na.rm=T),
                  Lon=mean(X,na.rm=T),
                  Host=mean(PLM2,na.rm=T))
  
  wj=NULL
  w=NULL
  for(i in 1:npatch){
    #print(paste("Patch",i))
    for(j in 1:nrow(patches_y)){
      
      d=dist(patches_y[c(i,j),3:4])
      wj[j]=exp(-d)*sqrt(patches_y$Host[j])*as.numeric(patches_y$PhoP[j])
    }   
    w[i]=sum(wj[-i],na.rm=T)
  }
  out=data.frame(Patch=patches_y$Patch[1:npatch],w=w)
  colnames(out)=c("Patch",paste(year))
  outdat=merge(outdat,out,by="Patch",all.x=T)
}
#colnames(outdat)=c("Patch",as.character(years))
b=Sys.time()
b-a #1000 patches 14 years 1.29 hours
#All patches all years 5.74hours 

outdat2=outdat[order(as.numeric(outdat$Patch)),]

outdat2[which(outdat2[,2]==0),2:20]=NA


write.csv(outdat2,"Mildew_neighbourhooddata.csv",row.names=F)

####Matrix algebra
