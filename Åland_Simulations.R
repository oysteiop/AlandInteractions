##########################################################################
######################### - Åland Interactions - #########################
##########################################################################

###########################################################
#### - Model-based simulation of population dynamics - ####
###########################################################

setwd("Z:/data/aland")
rm(list=ls())
library(lme4)
library(reshape2)
library(plyr)

newdat=read.csv("colextdat.csv")
newdat$Patch=as.factor(newdat$Patch)

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
summary(newdat$VS)

years=unique(newdat$Year)
comb=read.csv("combdat.csv")
patches=read.csv("patch_data.csv")
yearlyDat=read.csv("yearlyDat.csv")

invlogit=function(x){1/(1+exp(-x))}
colr=rgb(0,0,1,alpha=.5)
colr2=rgb(1,0,0,alpha=.5)


#### - M. cinxia - ####
load("fitted_mods/MelCc_modList.RData")
mMelCc=MelCc_modList[["mMelCc"]]
load("fitted_mods/MelCe_modList.RData")
mMelCe14=MelCe_modList[["mMelCe14"]]

#Simulation setup
cmod=mMelCc
emod=mMelCe14
nsim=100
years=unique(newdat$Year)
predN=matrix(NA,nrow=length(years)-1,ncol=nsim)
predP=matrix(NA,nrow=length(years)-1,ncol=nsim)
alldat=list()
ns=tapply(newdat$Patch,newdat$Year,length)
for(y in 1:length(years)){
  alldat[[y]]=as.data.frame(matrix(NA,nrow=ns[1],ncol=3+nsim))
}

PLtrend=0  #All PL covers changed by this factor

a=Sys.time()
for(j in 1:nsim){
  print(paste("Run",j))
  
  rc=MASS::mvrnorm(1,mu=summary(cmod)$coef[,1],Sigma=vcov(cmod))
  re=MASS::mvrnorm(1,mu=summary(emod)$coef[,1],Sigma=vcov(emod))

  for(y in 1:(length(years)-1)){
    year=years[y+1]
    
    lastyear=subset(newdat[newdat$Year==year-1,],select=c("Patch","MelC","PLM2"))
    lastyear$MelC=as.numeric(lastyear$MelC>0)
    
    if(y>1){lastyear=subset(preddat,select=c("Patch","predBin","PLM2"))
    }
    
    colnames(lastyear)=c("Patch","MelC","PLM2last")
    
    thisyear=comb[comb$Year==year,]
    
    preddat=merge(lastyear,thisyear,by="Patch",all.x=T)
    #preddat$logPLM2=preddat$logPLM2 + PLtrend
    
    #Update neighbourhood
    ssdat=merge(patches,preddat,by="Patch",all.x=T)
    ssdat$Patch=as.factor(ssdat$Patch)
    ssdat=ssdat[-which(is.na(ssdat$X)),]
    
    patches_y=subset(ssdat,select=c("Patch","MelC","Y","X","AREA"))
    names(patches_y)=c("Patch","MelC","Lat","Lon","Area")
    
    dmat=as.matrix(dist(patches_y[,3:4],diag=T,upper=T))
    Nvec=sqrt(patches_y$Area)*as.numeric(patches_y$MelC>0)

    edist=exp(-dmat)
    p2=edist*Nvec
    
    aa=apply(p2,2,sum,na.rm=T)
    for(i in 1:length(aa)){
      aa[i]=aa[i]-p2[i,i]
    }
    
    out=data.frame(Patch=patches_y$Patch,w=aa)
    preddat=merge(preddat,out,by="Patch",all.x=T)
    ####
    
    alldat[[y]][,1:3]=data.frame(Year=rep(paste(year),nrow(preddat)),Patch=preddat$Patch,Pred=rep(NA,nrow(preddat)))
    
    rYear=ranef(cmod)$Year[which(row.names(ranef(cmod)$Year)==year),1]
    if(length(rYear)<1){rYear=rnorm(1,mean = 0, sd = attr(summary(cmod)$varcor$Year,"stddev"))}

    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(rc[-1])))
    pd$logNH=log(preddat$w)
    rcy=rc
    rcy[1]=rc[1]+rYear
    predc=invlogit(as.matrix(pd)%*%as.matrix(rcy))
    
    rYear=ranef(emod)$Year[which(row.names(ranef(emod)$Year)==year),1]
    if(length(rYear)<1){rYear=rnorm(1,mean = 0, sd = attr(summary(emod)$varcor$Year,"stddev"))}
    
    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(re[-1])))
    pd$logNH=log(preddat$w)
    rey=re
    rey[1]=re[1]+rYear
    prede=invlogit(as.matrix(pd)%*%as.matrix(rey))
    
    preddat$pred=rep(NA,nrow(preddat))
    preddat$pred[which(preddat$MelC==1)]=(1-prede[which(preddat$MelC==1)])
    preddat$pred[which(preddat$MelC==0)]=predc[which(preddat$MelC==0)]
    
    preddat$predBin=rbinom(nrow(preddat),1,p=preddat$pred)
    predN[y,j]=sum(preddat$predBin,na.rm=T)
    predP[y,j]=sum(preddat$predBin,na.rm=T)/(sum(preddat$predBin>-1,na.rm=T))
    
    alldat[[y]][,j+3]=preddat$predBin
    }
}
Sys.time()-a

allCinxiadat=alldat
save(allCinxiadat, file = "allCinxiadat_100.RData")

#### - Mildew - ####
load("fitted_mods/PhoPc_modList.RData")
mPhoPc9=PhoPc_modList[["mPhoPc9"]]
load("fitted_mods/PhoPe_modList.RData")
mPhoPe3=PhoPe_modList[["mPhoPe3"]]

#Simulation setup
cmod=mPhoPc9
emod=mPhoPe3
nsim=100
years=unique(newdat$Year)

predN=matrix(NA,nrow=length(years)-1,ncol=nsim)
predP=matrix(NA,nrow=length(years)-1,ncol=nsim)
alldat=list()
ns=tapply(newdat$Patch,newdat$Year,length)
for(y in 1:length(years)){
  alldat[[y]]=as.data.frame(matrix(NA,nrow=ns[1],ncol=3+nsim))
}

PLtrend=0  #All PL covers changed by this factor

#### - Run simulations - ####
a=Sys.time()
for(j in 1:nsim){
  print(paste("Run",j))
  
  #Regression params drawn from multivariate sampling distributions
  rc=MASS::mvrnorm(1,mu=summary(cmod)$coef[,1],Sigma=vcov(cmod))
  re=MASS::mvrnorm(1,mu=summary(emod)$coef[,1],Sigma=vcov(emod))
  
for(y in 1:(length(years)-1)){
  year=years[y+1]
  
  lastyear=subset(newdat[newdat$Year==year-1,],select=c("Patch","PhoP","PLM2"))
  lastyear$PhoP=as.numeric(lastyear$PhoP>0)
  
  if(y>1){lastyear=subset(preddat,select=c("Patch","predBin","PLM2"))
  }
  
  colnames(lastyear)=c("Patch","PhoP","PLM2last")
  thisyear=comb[comb$Year==year,]
  
  preddat=merge(lastyear,thisyear,by="Patch",all.x=T)
  preddat$logPLM2=preddat$logPLM2+PLtrend
  
  #Update neighbourhood
    ssdat=merge(patches,preddat,by="Patch",all.x=T)
    ssdat$Patch=as.factor(ssdat$Patch)
    ssdat=ssdat[-which(is.na(ssdat$X)),]
    
    patches_y=subset(ssdat,select=c("Patch","PhoP","Y","X","PLM2last"))
    patches_y$PhoP=patches_y$PhoP>0
    names(patches_y)=c("Patch","PhoP","Lat","Lon","Host")
    
    dmat=as.matrix(dist(patches_y[,3:4],diag=T,upper=T))
    Nvec=sqrt(patches_y$Host)*as.numeric(patches_y$PhoP)
    edist=exp(-dmat)
    p2=edist*Nvec
    
    aa=apply(p2,2,sum,na.rm=T)
    for(i in 1:length(aa)){
      aa[i]=aa[i]-p2[i,i]
    }
    
    out=data.frame(Patch=patches_y$Patch,w=aa)
    preddat=merge(preddat,out,by="Patch",all.x=T)
  ####
  
    alldat[[y]][,1:3]=data.frame(Year=rep(paste(year),nrow(preddat)),Patch=preddat$Patch,Pred=rep(NA,nrow(preddat)))
  
    rYear=ranef(cmod)$Year[which(row.names(ranef(cmod)$Year)==year),1]
    if(length(rYear)<1){rYear=rnorm(1,mean = 0, sd = attr(summary(cmod)$varcor$Year,"stddev"))}
  
    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(rc[-1])))
    pd$logNH_Mil=log(preddat$w)
    rcy=rc
    rcy[1]=rc[1]+rYear
    predc=invlogit(as.matrix(pd)%*%as.matrix(rcy))
    
    rYear=ranef(emod)$Year[which(row.names(ranef(emod)$Year)==year),1]
    if(length(rYear)<1){rYear=rnorm(1,mean = 0, sd = attr(summary(emod)$varcor$Year,"stddev"))}
    
    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(re[-1])))
    pd$logNH_Mil=log(preddat$w)
    rey=re
    rey[1]=re[1]+rYear
    prede=invlogit(as.matrix(pd)%*%as.matrix(rey))
    
    preddat$pred=rep(NA,nrow(preddat))
    preddat$pred[which(preddat$PhoP==1)]=(1-prede[which(preddat$PhoP==1)])
    preddat$pred[which(preddat$PhoP==0)]=predc[which(preddat$PhoP==0)]
    
    preddat$predBin=rbinom(nrow(preddat),1,p=preddat$pred)
    predN[y,j]=sum(preddat$predBin,na.rm=T)
    predP[y,j]=sum(preddat$predBin,na.rm=T)/(sum(preddat$predBin>-1,na.rm=T))
    
    alldat[[y]][,j+3]=preddat$predBin
    }
}

Sys.time()-a

allMildewdat=alldat
save(allMildewdat, file = "allMildewdat_100.RData")


#### - Cotesia - ####
load("fitted_mods/CotPc_modList.RData")
mCotPc2=CotPc_modList[["mCotPc2"]]
load("fitted_mods/CotPe_modList.RData")
mCotPe2=CotPe_modList[["mCotPe2"]]

#Simulation setup
cmod=mCotPc2
emod=mCotPe2
nsim=100
years=unique(newdat$Year)
predN=matrix(NA,nrow=length(years)-1,ncol=nsim)
predP=matrix(NA,nrow=length(years)-1,ncol=nsim)

ns=tapply(newdat$Patch,newdat$Year,length)
alldat=list()
for(y in 1:length(years)){
  alldat[[y]]=as.data.frame(matrix(NA,nrow=ns[1],ncol=3+nsim))
}

extdat=list()
for(y in 1:length(years)){
  extdat[[y]]=as.data.frame(matrix(NA,nrow=ns[1],ncol=3+nsim))
}

coldat=list()
for(y in 1:length(years)){
  coldat[[y]]=as.data.frame(matrix(NA,nrow=ns[1],ncol=3+nsim))
}

PLtrend=0  #All PL covers changed by this factor

a=Sys.time()
for(j in 1:nsim){
  print(paste("Run",j))
  
  rc=MASS::mvrnorm(1,mu=summary(cmod)$coef[,1],Sigma=vcov(cmod))
  re=MASS::mvrnorm(1,mu=summary(emod)$coef[,1],Sigma=vcov(emod))
  
  for(y in 1:(length(years)-1)){
    year=years[y+1]
    
    lastyear=subset(newdat[newdat$Year==year-1,],select=c("Patch","CotP","PLM2"))
    lastyear$CotP=as.numeric(lastyear$CotP>0)
    
    if(y>1){lastyear=subset(preddat,select=c("Patch","predBin","PLM2"))
    }
    
    colnames(lastyear)=c("Patch","CotP","PLM2last")
    
    thisyear=comb[comb$Year==year,]
    
    preddat=merge(lastyear,thisyear,by="Patch",all.x=T)
    preddat$logPLM2=preddat$logPLM2+PLtrend
    
    #Update neighbourhood
    ssdat=merge(patches,preddat,by="Patch",all.x=T)
    ssdat$Patch=as.factor(ssdat$Patch)
    ssdat=ssdat[-which(is.na(ssdat$X)),]
    
    patches_y=subset(ssdat,select=c("Patch","CotP","Y","X","AREA"))
    names(patches_y)=c("Patch","CotP","Lat","Lon","Area")
    
    dmat=as.matrix(dist(patches_y[,3:4],diag=T,upper=T))
    Nvec=sqrt(patches_y$Area)*as.numeric(patches_y$CotP>0)
    Nvec[which(is.na(Nvec))]=0
    
    edist=exp(-dmat)
    p2=edist*Nvec
    
    aa=apply(p2,2,sum,na.rm=T)
    for(i in 1:length(aa)){
      aa[i]=aa[i]-p2[i,i]
    }
    
    out=data.frame(Patch=patches_y$Patch,w=aa)
    preddat=merge(preddat,out,by="Patch",all.x=T)
    ####
    
    alldat[[y]][,1:3]=data.frame(Year=rep(paste(year),nrow(preddat)),Patch=preddat$Patch,Pred=rep(NA,nrow(preddat)))
    coldat[[y]][,1:3]=data.frame(Year=rep(paste(year),nrow(preddat)),Patch=preddat$Patch,Pred=rep(NA,nrow(preddat)))
    extdat[[y]][,1:3]=data.frame(Year=rep(paste(year),nrow(preddat)),Patch=preddat$Patch,Pred=rep(NA,nrow(preddat)))
    
    rYear=ranef(cmod)$Year[which(row.names(ranef(cmod)$Year)==year),1]
    if(length(rYear)<1){rYear=rnorm(1,mean = 0, sd = attr(summary(cmod)$varcor$Year,"stddev"))}

    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(rc[-1])))
    pd$logNH_Cot=sqrt(preddat$w)
    rcy=rc
    rcy[1]=rc[1]+rYear
    predc=invlogit(as.matrix(pd)%*%as.matrix(rcy))
    
    rYear=ranef(emod)$Year[which(row.names(ranef(emod)$Year)==year),1]
    if(length(rYear)<1){rYear=rnorm(1,mean = 0, sd = attr(summary(emod)$varcor$Year,"stddev"))}
    
    pd=data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select = names(re[-1])))
    pd$logNH_Cot=sqrt(preddat$w)
    rey=re
    rey[1]=re[1]+rYear
    prede=invlogit(as.matrix(pd)%*%as.matrix(rey))
    
    preddat$pred=rep(NA,nrow(preddat))
    preddat$pred[which(preddat$CotP==1)]=(1-prede[which(preddat$CotP==1)])
    preddat$pred[which(preddat$CotP==0)]=predc[which(preddat$CotP==0)]
    
    preddat$predBin=rbinom(nrow(preddat),1,p=preddat$pred)
    predN[y,j]=sum(preddat$predBin,na.rm=T)
    #predP[y,j]=sum(preddat$predBin,na.rm=T)/(sum(preddat$predBin>-1,na.rm=T))
    predP[y,j]=sum(preddat$predBin,na.rm=T)/length(preddat$predBin) #Scale by total patches surveyed
    
    alldat[[y]][,j+3]=preddat$predBin
    coldat[[y]][,j+3]=predc
    extdat[[y]][,j+3]=prede
    }
}
Sys.time()-a

allCotesiadat=alldat

allCotesiadat=alldat
save(allCotesiadat, file = "allCotesiadat_100.RData")
save(coldat, file = "Cotesia_coldat_100.RData")


######################################
#### - Plot population dynamics - ####
######################################
yearlyDat=read.csv("yearlyDat.csv")
load("allCinxiadat_100.RData")
load("allMildewdat_100.RData")
load("allCotesiadat_100.RData")

#load("allCinxiadat_100_PL_m20.RData")
#load("allMildewdat_100_PL_m20.RData")
#load("allCotesiadat_100_PL_m20.RData")

cols=c(rgb(0,0,1,.5), rgb(0,1,0,.5), rgb(1,0,0,.5))


#Cinxia
predP=matrix(NA,nrow=length(allCinxiadat)-1,ncol=100)
for(i in 1:(length(allCinxiadat)-1)){
  sums=apply(allCinxiadat[[i]][,4:103], 2, sum, na.rm=T)
  ns=apply(allCinxiadat[[i]][,4:103]>-1, 2, sum, na.rm=T)
  predP[i,]=sums/ns
}

df=data.frame(Year=years[1:18]+1,
              predP=apply(predP,1,mean,na.rm=T),
              predPupper=apply(predP,1,quantile,c(.025,.975),na.rm=T)[2,],
              predPlower=apply(predP,1,quantile,c(.025,.975),na.rm=T)[1,])

x11()
par(mfrow=c(1,3))
plot(years,yearlyDat$MelC,type="b",pch=16,ylim=c(0,.6),xlab="",las=1,ylab="",main=expression("Butterfly"))
xx=2001:2018
polygon(c(xx,rev(xx)),c(df$predPlower[1:18],rev(df$predPupper[1:18])),col = cols[1], border = FALSE)
points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years,yearlyDat$MelC,type="b",pch=16)

mtext("Proportion of patches occupied", 2, line=2.5)

#Mildew
nc=ncol(allMildewdat[[1]])
predP=matrix(NA,nrow=length(allMildewdat)-1,ncol=100)
for(i in 1:(length(allMildewdat)-1)){
  sums=apply(allMildewdat[[i]][,4:nc], 2, sum, na.rm=T)
  ns=apply(allMildewdat[[i]][,4:nc]>-1, 2, sum, na.rm=T)
  predP[i,]=sums/ns
}

df=data.frame(Year=years[1:18]+1,
              predP=apply(predP,1,mean,na.rm=T),
              predPupper=apply(predP,1,quantile,c(.025,.975),na.rm=T)[2,],
              predPlower=apply(predP,1,quantile,c(.025,.975),na.rm=T)[1,])

plot(years,yearlyDat$PhoPcomb,type="b",pch=16,ylim=c(0,.5),xlab="",las=1,ylab="",main=expression("Mildew"))
xx=2001:2018
polygon(c(xx,rev(xx)),c(df$predPlower[1:18],rev(df$predPupper[1:18])),col = cols[2], border = FALSE)
points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years,yearlyDat$PhoPcomb,type="b",pch=16)

#Cotesia
predP=matrix(NA,nrow=length(allCotesiadat)-1,ncol=100)
for(i in 1:(length(allCotesiadat)-1)){
  sums=apply(allCotesiadat[[i]][,4:103], 2, sum, na.rm=T)
  ns=apply(allCotesiadat[[i]][,4:103]>-1, 2, length)
  predP[i,]=sums/ns
}

df=data.frame(Year=years[1:18]+1,
              predP=apply(predP,1,mean,na.rm=T),
              predPupper=apply(predP,1,quantile,c(.025,.975),na.rm=T)[2,],
              predPlower=apply(predP,1,quantile,c(.025,.975),na.rm=T)[1,])

plot(years,yearlyDat$CotPAdj,type="b",pch=16,ylim=c(0,.03),xlab="",ylab="",las=1,main=expression("Parasitoid"))
xx=2001:2018
polygon(c(xx,rev(xx)),c(df$predPlower[1:18],rev(df$predPupper[1:18])),col = cols[3], border = FALSE)
points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years,yearlyDat$CotPAdj,type="b",pch=16)


###################################
#### - Community predictions - ####
###################################
load("allCinxiadat_100.RData")
load("allCotesiadat_100.RData")
load("allMildewdat_100_PLx1.RData")
load("Cotesia_coldat_100.RData")

years=unique(newdat$Year)
dist=matrix(NA,ncol=6,nrow=100)
medians=matrix(NA,ncol=6,nrow=18)
upper=matrix(NA,ncol=6,nrow=18)
lower=matrix(NA,ncol=6,nrow=18)
w=list()

for(y in 1:(length(years)-1)){
  for(i in 1:100){
    mc=allCinxiadat[[y]][,c(1:2,i+3)]
    names(mc)=c("Year","Patch","MelC")
    cm=allCotesiadat[[y]][,c(1:2,i+3)]
    names(cm)=c("Year","Patch","CotP")
    pp=allMildewdat[[y]][,c(1:2,i+3)]
    names(pp)=c("Year","Patch","PhoP")
    
    rep=cbind(mc,cm,pp)
    rep=rep[,c(1:3,6,9)]
    head(rep,20)
    
    rep$CotP[which(rep$MelC==0)]=NA #'Butterfly extinction'
    rep$CotP[which(is.na(rep$MelC))]=NA #Missing data
    
    Mel_col_patches=which(rep$MelC==1 & is.na(rep$CotP))
    rep$CotP[Mel_col_patches] = rbinom(length(Mel_col_patches),1,coldat[[y]][Mel_col_patches,c(i+3)]) #Joint colonisations
    
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

head(wide)
head(medians)

rs=rowSums(medians)
wide=t(apply(wide,1,function(x){x/sum(x)}))
medians=t(apply(medians,1,function(x){x/sum(x)}))

upper2=upper
for(i in 1:nrow(upper2)){upper2[i,]=upper[i,]/rs[i]}
upper=upper2

lower2=lower
for(i in 1:nrow(lower2)){lower2[i,]=lower[i,]/rs[i]}
lower=lower2


titles=c("No species", "Mildew only", "Butterfly only", "Butterfly & Mildew", "Butterfly & Parasitoid", "All species")

x11()
par(mfrow=c(2,3),mar=c(2,2,2,2), oma=c(3,3,0,0))
xx=2000:2017
for(i in 1:6){
  plot(xx+1, medians[,i], type = "n",ylab="",las=1,xlab="",xlim=c(2000,2018),
       ylim=c(min(c(lower[,i],wide[,i]),na.rm=T),max(c(upper[,i],wide[,i]),na.rm=T)),main=as.expression(titles[i]))
  polygon(c(xx+1,rev(xx+1)),c(lower[,i],rev(upper[,i])),col = colr, border = FALSE)

  lines(rownames(wide),wide[,i],lwd=2)
  lines(xx+1,medians[,i],lwd=2,col="darkblue")
  
}
#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=1, outer=T)
mtext(2,text="Proportion of patches",line=1, outer=T)



x11()
par(mfrow=c(3,3),mar=c(2,2,2,2), oma=c(3,3,0,0))
xx=2000:2017
for(i in 1:6){
  plot(xx+1, medians[,i], type = "n",ylab="",las=1,xlab="",xlim=c(2000,2018),
       ylim=c(min(c(lower[,i],wide[,i]),na.rm=T),max(c(upper[,i],wide[,i]),na.rm=T)),main=as.expression(titles[i]))
  polygon(c(xx+1,rev(xx+1)),c(lower[,i],rev(upper[,i])),col = colr, border = FALSE)
  
  lines(rownames(wide),wide[,i],lwd=2)
  lines(xx+1,medians[,i],lwd=2,col="darkblue")
  
}
#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=1, outer=T)
mtext(2,text="Proportion of patches",line=1, outer=T)



####Communities with Mildew lag####
dist=matrix(NA,ncol=6,nrow=100)
medians=matrix(NA,ncol=6,nrow=18)
upper=matrix(NA,ncol=6,nrow=18)
lower=matrix(NA,ncol=6,nrow=18)
w=list()

y=2
j=i=1
a=Sys.time()
for(y in 2:(length(years)-1)){
  print(paste("Year=",y))
  for(i in 1:100){
    mc=allCinxiadat[[y]][,c(1:2,i+3)]
    names(mc)=c("Year","Patch","MelC")
    cm=allCotesiadat[[y]][,c(1:2,i+3)]
    names(cm)=c("Year","Patch","CotP")
    pp=allMildewdat[[y-1]][,c(1:2,i+3)]
    names(pp)=c("Year","Patch","PhoP")
    
    rep=cbind(mc,cm)
    rep=rep[,c(1:3,6)]
    
    rep=merge(rep, pp, by="Patch", sort=F)
    rep=rep[,c(1:4,6)]
    colnames(rep)[2]="Year"
    
    rep=ddply(rep,.(Year,Patch),summarize,
               MelC=mean(MelC),
               CotP=mean(CotP),
               PhoP=mean(PhoP))
    rep$Patch=factor(rep$Patch)
    
    coldat[[y]]$V2=factor(coldat[[y]]$V2)
    coldat[[y]]$col=coldat[[y]][,i+3]
    colp=ddply(coldat[[y]],.(V2), summarize,
               colp=mean(col,na.rm=T))
    
    rep$CotP[which(rep$MelC==0)]=NA #'Butterfly extinction'
    rep$CotP[which(is.na(rep$MelC))]=NA #Missing data
    
    Mel_col_patches=which(rep$MelC==1 & is.na(rep$CotP))
    rep$CotP[Mel_col_patches] = rbinom(length(Mel_col_patches),1,colp$colp[Mel_col_patches]) #Joint colonisations
    
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
colnames(wideLag)=sub("X","",colnames(wideLag))

mm=match(colnames(wideLag),colnames(medians))
medians=medians[,mm]
upper=upper[,mm]
lower=lower[,mm]
dim(wideLag)
dim(medians)

head(wideLag,20) #Observed data, 2009=NA, 2000=NULL
head(medians,20) #Predicted data, 2001=NA
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
       ylim=c(min(c(lower[,i],wideLag[,i]),na.rm=T),max(c(upper[,i],wideLag[,i]),na.rm=T)),main=colnames(wideLag)[i])
  polygon(c(xx+1,rev(xx+1)),c(lower[,i],rev(upper[,i])),col = colr, border = FALSE)
  
  lines(rownames(wideLag),wideLag[,i],lwd=2)
  lines(xx+1,medians[,i],lwd=2,col="darkblue")
  
  
}

write.csv(medians, "mediansLag.csv")
write.csv(upper, "upperLag.csv")
write.csv(lower, "lowerLag.csv")

medians=read.csv("mediansLag.csv")
upper=read.csv("upperLag.csv")
lower=read.csv("lowerLag.csv")

medians=medians[,-1]
colnames(medians)=sub("X","",colnames(medians))
upper=upper[,-1]
colnames(upper)=sub("X","",colnames(upper))
lower=lower[,-1]
colnames(lower)=sub("X","",colnames(lower))
