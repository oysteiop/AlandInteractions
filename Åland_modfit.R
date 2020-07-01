##########################################################################
######################### - Aland Interactions - #########################
##########################################################################

###########################################
#### - Single-species col/ext models - ####
###########################################
rm(list=ls())
gc()
setwd("C:/data/aland/")

library(reshape2)
library(plyr)
library(mvtnorm)  
library(lme4)
library(MuMIn)

comb = read.csv("data/combdat.csv")

# Summary stats
names(comb)
apply(comb[,9:14], 2, sum, na.rm=T) #Number of events
apply(comb[,9:14]>-Inf, 2, sum, na.rm=T) #Number of possible events
signif(apply(comb[,9:14], 2, sum, na.rm=T)/apply(comb[,9:14]>-Inf, 2, sum, na.rm=T), 2)

d = comb[which(comb$CotPe>0),]
table(d$MelCe) #Joint extinctions

d = comb[which(comb$CotPc>0),]
table(d$MelCc) #Joint colonisations

# Weather
names(comb)
totalprec = apply(comb[,20:23], 1, sum)

par(mfrow=c(2,3))
yprec = tapply(totalprec, comb$Year, mean, na.rm=T)
yprec[18]/mean(yprec)*100
plot(2001:2018, yprec, type="b", main="Total")

yprec = tapply(comb$MayPrec, comb$Year, mean, na.rm=T)
yprec[18]/mean(yprec)*100
plot(2001:2018, yprec, type="b", main="May")

yprec = tapply(comb$junePrec, comb$Year, mean, na.rm=T)
yprec[18]/mean(yprec)*100
plot(2001:2018, yprec, type="b", main="June")

yprec = tapply(comb$julyPrec, comb$Year, mean, na.rm=T)
yprec[18]/mean(yprec)*100
plot(2001:2018, yprec, type="b", main="July")

yprec = tapply(comb$augustPrec, comb$Year, mean, na.rm=T)
yprec[18]/mean(yprec)*100
plot(2001:2018, yprec, type="b", main="August")

##########################################
#### - Fitting the candidate models - ####
##########################################

#### Plantago abundance models ####
# Remember to define this subset of data without inserting means for missing values
PlaL_dat = na.omit(subset(comb, select=c("logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","logArea","Year","Patch")))
sort(unique(PlaL_dat$Year))

mPlaL = lmer(logPLM2~logArea+logMayPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL2 = lmer(logPLM2~logArea+logJunePrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL3 = lmer(logPLM2~logArea+logJulyPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL4 = lmer(logPLM2~logArea+logAugPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL5 = lmer(logPLM2~logArea+logMayPrec+logJunePrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL6 = lmer(logPLM2~logArea+logMayPrec+logJulyPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL7 = lmer(logPLM2~logArea+logMayPrec+logAugPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL8 = lmer(logPLM2~logArea+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL9 = lmer(logPLM2~logArea+logJunePrec+logAugPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL10 = lmer(logPLM2~logArea+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL11 = lmer(logPLM2~logArea+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL12 = lmer(logPLM2~logArea+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL13 = lmer(logPLM2~logArea+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL14 = lmer(logPLM2~logArea+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), data=PlaL_dat)
mPlaL15 = lmer(logPLM2~logArea+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), data=PlaL_dat)

PlaL_modList = list(mPlaL,mPlaL2,mPlaL3,mPlaL4,mPlaL5,mPlaL6,mPlaL7,mPlaL8,mPlaL9,mPlaL10,mPlaL11,mPlaL12,mPlaL13,mPlaL14,mPlaL15)
names(PlaL_modList) = c("mPlaL", paste0("mPlaL", 2:15))
save(PlaL_modList, file="fitted_mods/PlaL_modList.RData")

#### Butterfly colonisation models ####
MelCc_dat = na.omit(subset(comb, select=c("MelCc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","VS","logNH","Road","Year","Patch")))

mMelCc = glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc2 = glmer(MelCc~logPLM2+VS+logNH+Road+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc3 = glmer(MelCc~logPLM2+VS+logNH+Road+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc4 = glmer(MelCc~logPLM2+VS+logNH+Road+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc5 = glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc6 = glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc7 = glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc8 = glmer(MelCc~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc9 = glmer(MelCc~logPLM2+VS+logNH+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc10 = glmer(MelCc~logPLM2+VS+logNH+Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc11 = glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc12 = glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc13 = glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc14 = glmer(MelCc~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)
mMelCc15 = glmer(MelCc~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCc_dat)

MelCc_modList = list(mMelCc,mMelCc2,mMelCc3,mMelCc4,mMelCc5,mMelCc6,mMelCc7,mMelCc8,mMelCc9,mMelCc10,mMelCc11,mMelCc12,mMelCc13,mMelCc14,mMelCc15)
names(MelCc_modList) = c("mMelCc", paste0("mMelCc", 2:15))
save(MelCc_modList, file="fitted_mods/MelCc_modList.RData")

#### Butterfly extinction models ####
MelCe_dat = na.omit(subset(comb, select=c("MelCe","logMayPrec","logJunePrec","logJulyPrec","julyPrec","logAugPrec","logPLM2","VS","logNH","Road","Year","Patch")))

mMelCe = glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe2 = glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe3 = glmer(MelCe~logPLM2+VS+logNH+Road+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe4 = glmer(MelCe~logPLM2+VS+logNH+Road+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe5 = glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe6 = glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe7 = glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe8 = glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe9 = glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe10 = glmer(MelCe~logPLM2+VS+logNH+Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe11 = glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe12 = glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe13 = glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe14 = glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
mMelCe15 = glmer(MelCe~logPLM2+VS+logNH+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=MelCe_dat)

MelCe_modList = list(mMelCe,mMelCe2,mMelCe3,mMelCe4,mMelCe5,mMelCe6,mMelCe7,mMelCe8,mMelCe9,mMelCe10,mMelCe11,mMelCe12,mMelCe13,mMelCe14,mMelCe15)
names(MelCe_modList) = c("mMelCe", paste0("mMelCe", 2:15))
save(MelCe_modList, file="fitted_mods/MelCe_modList.RData")

#### Mildew colonisation models ####
PhoPc_dat = na.omit(subset(comb, select=c("PhoPc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","logNH_Mil","logTotalNH","Road","Year","Patch")))

mPhoPc = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc2 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc3 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc4 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc5 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc6 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc7 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc8 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc9 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc10 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc11 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc12 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc13 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc14 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
mPhoPc15 = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)

PhoPc_modList = list(mPhoPc,mPhoPc2,mPhoPc3,mPhoPc4,mPhoPc5,mPhoPc6,mPhoPc7,mPhoPc8,mPhoPc9,mPhoPc10,mPhoPc11,mPhoPc12,mPhoPc13,mPhoPc14,mPhoPc15)
names(PhoPc_modList) = c("mPhoPc", paste0("mPhoPc", 2:15))
save(PhoPc_modList, file="fitted_mods/PhoPc_modList.RData")

#### Mildew extinction models ####
PhoPe_dat = na.omit(subset(comb, select=c("PhoPe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","logNH_Mil","logTotalNH","Road","Year","Patch")))

mPhoPe = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe2 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe3 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH +Road+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe4 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe5 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe6 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe7 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe8 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe9 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe10 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe11 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe12 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe13 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe14 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe15 = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)

PhoPe_modList = list(mPhoPe,mPhoPe2,mPhoPe3,mPhoPe4,mPhoPe5,mPhoPe6,mPhoPe7,mPhoPe8,mPhoPe9,mPhoPe10,mPhoPe11,mPhoPe12,mPhoPe13,mPhoPe14,mPhoPe15)
names(PhoPe_modList) = c("mPhoPe", paste0("mPhoPe", 2:15))
save(PhoPe_modList, file="fitted_mods/PhoPe_modList.RData")

#### Parasitoid colonisation models ####
CotPc_dat = na.omit(subset(comb, select=c("CotPc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPc_dat = CotPc_dat[CotPc_dat$Year!="2018",]

mCotPc = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc2 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc3 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc4 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc5 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc6 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc7 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc8 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc9 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc10 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc11 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc12 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc13 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc14 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc15 = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPc_dat)

CotPc_modList = list(mCotPc,mCotPc2,mCotPc3,mCotPc4,mCotPc5,mCotPc6,mCotPc7,mCotPc8,mCotPc9,mCotPc10,mCotPc11,mCotPc12,mCotPc13,mCotPc14,mCotPc15)
names(CotPc_modList) = c("mCotPc", paste0("mCotPc", 2:15))
save(CotPc_modList, file="fitted_mods/CotPc_modList.RData")

#### Parasitoid extinction models ####
CotPe_dat = na.omit(subset(comb, select=c("CotPe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPe_dat = CotPe_dat[CotPe_dat$Year!="2018",]

mCotPe = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe2 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe3 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe4 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe5 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe6 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe7 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe8 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe9 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe10 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe11 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe12 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe13 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe14 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe15 = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logMayPrec+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)

CotPe_modList = list(mCotPe,mCotPe2,mCotPe3,mCotPe4,mCotPe5,mCotPe6,mCotPe7,mCotPe8,mCotPe9,mCotPe10,mCotPe11,mCotPe12,mCotPe13,mCotPe14,mCotPe15)
names(CotPe_modList) = c("mCotPe", paste0("mCotPe", 2:15))
save(CotPe_modList, file="fitted_mods/CotPe_modList.RData")

############################################################
#### - Parameter estimates from highest ranked models - ####
############################################################

# Plantago abundance
# Remember to define this subset of data without inserting means for missing values
PlaL_dat = na.omit(subset(comb, select=c("logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logPLM2","logArea","Year","Patch")))
load("fitted_mods/PlaL_modList.RData")

AICtab = data.frame(AIC=unlist(lapply(PlaL_modList, AIC)))
AICtab$delta = AICtab$AIC-min(AICtab$AIC)
AICtab = AICtab[order(AICtab$delta),]
AICtab #Mod1 best, 7 good

mPlaL = PlaL_modList[["mPlaL"]]
summary(mPlaL)
r2 = r.squaredGLMM(mPlaL)
r2

# Varcomps
mod = mPlaL
dat = PlaL_dat

xdat = as.matrix(data.frame(Intercept=rep(1,nrow(dat)), subset(dat, select = rownames(summary(mod)$coef)[-1])))
fixed = t(summary(mod)$coef[,1])%*%cov(xdat)%*%summary(mod)$coef[,1]
yr = as.numeric(VarCorr(mod)$Year)
pa = as.numeric(VarCorr(mod)$Patch)
total = fixed+pa+yr
df = cbind(fixed, yr, pa)
round(df/sum(df)*100, 2)

area = (summary(mod)$coef[2,1]^2)*var(dat$logArea)
mp = (summary(mod)$coef[3,1]^2)*var(dat$logMayPrec)
round(area/total*100, 2)
round(mp/total*100, 2)

plot(comb$logMayPrec, comb$logPLM2)
xx = seq(min(comb$logMayPrec, na.rm=T), max(comb$logMayPrec, na.rm=T), 0.01)
la = comb$logArea[comb$logArea>-Inf]
yy = -1.79 + 0.397*mean(la, na.rm=T) + 0.1204*xx
lines(xx, yy, lwd=3, col="red")

# Butterfly col
MelCc_dat = na.omit(subset(comb, select=c("MelCc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH","logTotalNH","Road","Year","Patch")))
load("fitted_mods/MelCc_modList.RData")

AICtab = data.frame(AIC=unlist(lapply(MelCc_modList,AIC)))
AICtab$delta = AICtab$AIC-min(AICtab$AIC)
AICtab = AICtab[order(AICtab$delta),]
AICtab #Mod1 best, 11 good

mMelCc = MelCc_modList[["mMelCc"]]
summary(mMelCc)
r2 = r.squaredGLMM(mMelCc)
r2

# Varcomps
mod = mMelCc
dat = MelCc_dat

xdat = as.matrix(data.frame(Intercept=rep(1,nrow(dat)), subset(dat, select = rownames(summary(mod)$coef)[-1])))
fixed = t(summary(mod)$coef[,1])%*%cov(xdat)%*%summary(mod)$coef[,1]
yr = as.numeric(VarCorr(mod)$Year)
pa = as.numeric(VarCorr(mod)$Patch)
total = fixed+pa+yr
df = cbind(fixed, yr, pa)
round(df/sum(df)*100, 2)

pl = (summary(mod)$coef[2,1]^2)*var(dat$logPLM2)
vs = (summary(mod)$coef[3,1]^2)*var(dat$VS)
round(pl/total*100, 2)
round(vs/total*100, 2)

# Butterfly ext
MelCe_dat = na.omit(subset(comb, select=c("MelCe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH","logTotalNH","Road","Year","Patch")))
load("fitted_mods/MelCe_modList.RData")

AICtab = data.frame(AIC=unlist(lapply(MelCe_modList, AIC)))
AICtab$delta = AICtab$AIC-min(AICtab$AIC)
AICtab = AICtab[order(AICtab$delta),]
AICtab #Mod10 best but poor convergence, use mod14

mMelCe14 = MelCe_modList[["mMelCe14"]]
summary(mMelCe14)
r2 = r.squaredGLMM(mMelCe14)
r2

# Mildew col
PhoPc_dat = na.omit(subset(comb, select=c("PhoPc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Mil","logTotalNH","Road","Year","Patch")))
load("fitted_mods/PhoPc_modList.RData")

AICtab = data.frame(AIC=unlist(lapply(PhoPc_modList, AIC)))
AICtab$delta = AICtab$AIC-min(AICtab$AIC)
AICtab = AICtab[order(AICtab$delta),]
AICtab #Mod9 best

mPhoPc9 = PhoPc_modList[["mPhoPc9"]]
summary(mPhoPc9)
r2 = r.squaredGLMM(mPhoPc9)
r2

# Varcomps
mod = mPhoPc9
dat = PhoPc_dat

xdat = as.matrix(data.frame(Intercept=rep(1,nrow(dat)), subset(dat, select = rownames(summary(mod)$coef)[-1])))
fixed = t(summary(mod)$coef[,1])%*%cov(xdat)%*%summary(mod)$coef[,1]
yr = as.numeric(VarCorr(mod)$Year)
pa = as.numeric(VarCorr(mod)$Patch)
total = fixed+pa+yr
df = cbind(fixed, yr, pa)
round(df/sum(df)*100,2)

pl = (summary(mod)$coef[2,1]^2)*var(dat$logPLM2)
nhMil = (summary(mod)$coef[3,1]^2)*var(dat$logNH_Mil)
nh = (summary(mod)$coef[4,1]^2)*var(dat$logTotalNH)
june = (summary(mod)$coef[6,1]^2)*var(dat$logJunePrec)
prec = t((summary(mod)$coef[6:7,1])) %*% cov(subset(dat, select=c("logJunePrec", "logAugPrec"))) %*% (summary(mod)$coef[6:7,1])

round(pl/total*100, 2)
round(nhMil/total*100, 2)
round(nh/total*100, 2)
round(prec/total*100, 2)

# Mildew extinction
PhoPe_dat = na.omit(subset(comb, select=c("PhoPe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Mil","Road","Year","Patch")))
load("fitted_mods/PhoPe_modList.RData")

AICtab = data.frame(AIC=unlist(lapply(PhoPe_modList, AIC)))
AICtab$delta = AICtab$AIC-min(AICtab$AIC)
AICtab = AICtab[order(AICtab$delta),]
AICtab #Mod3 best

mPhoPe3 = PhoPe_modList[["mPhoPe3"]]
summary(mPhoPe3)
r2 = r.squaredGLMM(mPhoPe3)
r2

# Parasitoid col
CotPc_dat = na.omit(subset(comb, select=c("CotPc","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPc_dat = CotPc_dat[CotPc_dat$Year!="2018",]
load("fitted_mods/CotPc_modList.RData")

AICtab = data.frame(AIC=unlist(lapply(CotPc_modList, AIC)))
AICtab$delta = AICtab$AIC-min(AICtab$AIC)
AICtab = AICtab[order(AICtab$delta),]
AICtab #Mod 2 best, Mod 1 and Mod 4 equal

mCotPc2 = CotPc_modList[["mCotPc2"]]
summary(mCotPc2)
r2 = r.squaredGLMM(mCotPc2)
r2

# Varcomps
mod = mCotPc2
dat = CotPc_dat

xdat = as.matrix(data.frame(Intercept=rep(1,nrow(dat)), subset(dat, select = rownames(summary(mod)$coef)[-1])))
fixed = t(summary(mod)$coef[,1])%*%cov(xdat)%*%summary(mod)$coef[,1]
yr = as.numeric(VarCorr(mod)$Year)
pa = as.numeric(VarCorr(mod)$Patch)
total = fixed+pa+yr
df = cbind(fixed, yr, pa)
round(df/sum(df)*100, 2)

pl = (summary(mod)$coef[2,1]^2)*var(dat$logPLM2)
vs = (summary(mod)$coef[3,1]^2)*var(dat$VS)
round(pl/total*100, 1)
round(vs/total*100, 1)

# Parasitoid ext
CotPe_dat = na.omit(subset(comb, select=c("CotPe","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPe_dat = CotPe_dat[CotPe_dat$Year!="2018",]
load("fitted_mods/CotPe_modList.RData")

AICtab = data.frame(AIC=unlist(lapply(CotPe_modList, AIC)))
AICtab$delta = AICtab$AIC-min(AICtab$AIC)
AICtab = AICtab[order(AICtab$delta),]
AICtab #Mod2 best

mCotPe2 = CotPe_modList[["mCotPe2"]]
summary(mCotPe2)
r2 = r.squaredGLMM(mCotPe2)
r2

#############################################################
#### - Adding patch state to the highest ranked models - ####
#############################################################

comb2 = read.csv("data/combdat2.csv")

# Butterfly col
MelCc_dat = na.omit(subset(comb2, select=c("MelCc","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH","logTotalNH","Road","Year","Patch")))

mMelCca = glmer(MelCc~logPLM2 + VS + logNH + Road + logMayPrec + (1|Patch) + (1|Year), family="binomial", data=MelCc_dat)
mMelCcb = glmer(MelCc~logPLM2 + VS + logNH + Road + logMayPrec + Statelast + (1|Patch) + (1|Year), family="binomial", data=MelCc_dat)
AIC(mMelCca, mMelCcb)
AIC(mMelCca, mMelCcb)$AIC[2]-AIC(mMelCca, mMelCcb)$AIC[1]
logLik(mMelCca)
logLik(mMelCcb)

mMelCcc = glmer(MelCc~scale(logPLM2, scale=F) + scale(VS, scale=F) + scale(logNH, scale=F) + scale(Road, scale=F) +
                scale(logMayPrec, scale=F) + Statelast - 1 + (1|Patch) + (1|Year), family="binomial", data = MelCc_dat)
summary(mMelCcc)

#000
s = summary(mMelCcc)$coef[6,1]
sse = summary(mMelCcc)$coef[6,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#001
s = summary(mMelCcc)$coef[7,1]
sse = summary(mMelCcc)$coef[7,2]
round(1/(1+exp(-s))*100,2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

# Butterfly ext
MelCe_dat = na.omit(subset(comb2, select=c("MelCe","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH","logTotalNH","Road","Year","Patch")))

mMelCe14a = glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
mMelCe14b = glmer(MelCe~logPLM2+VS+logNH+Road+logJunePrec+logJulyPrec+logAugPrec+Statelast + (1|Patch)+(1|Year),family="binomial",data=MelCe_dat)
AIC(mMelCe14a, mMelCe14b)
AIC(mMelCe14a, mMelCe14b)$AIC[2]-AIC(mMelCe14a, mMelCe14b)$AIC[1]
logLik(mMelCe14a)
logLik(mMelCe14b)

mMelCe14c = glmer(MelCe~scale(logPLM2, scale=F) + scale(VS, scale=F) + scale(logNH, scale=F)+
                 scale(Road, scale=F) + scale(logJunePrec, scale=F) +scale(logJulyPrec, scale=F) + 
                 scale(logAugPrec, scale=F) + Statelast -1 + (1|Patch)+(1|Year), family="binomial", data=MelCe_dat)
summary(mMelCe14c)

#100
s = summary(mMelCe14c)$coef[8,1]
sse = summary(mMelCe14c)$coef[8,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#101
s = summary(mMelCe14c)$coef[9,1]
sse = summary(mMelCe14c)$coef[9,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#110
s = summary(mMelCe14c)$coef[10,1]
sse = summary(mMelCe14c)$coef[10,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#111
s = summary(mMelCe14c)$coef[11,1]
sse = summary(mMelCe14c)$coef[11,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

# Mildew col
PhoPc_dat = na.omit(subset(comb2, select=c("PhoPc","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Mil","logTotalNH","Road","Year","Patch")))

mPhoPc9a = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+logAugPrec+(1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
mPhoPc9b = glmer(PhoPc~logPLM2+logNH_Mil+logTotalNH+Road+logJunePrec+logAugPrec+Statelast + (1|Patch)+(1|Year),family="binomial",data=PhoPc_dat)
AIC(mPhoPc9a, mPhoPc9b)
AIC(mPhoPc9a, mPhoPc9b)$AIC[2]-AIC(mPhoPc9a, mPhoPc9b)$AIC[1]
logLik(mPhoPc9a)
logLik(mPhoPc9b)

mPhoPc9c = glmer(PhoPc~scale(logPLM2, scale=F) + scale(logNH_Mil, scale=F) + scale(logTotalNH, scale=F) + scale(Road, scale=F)+
                  scale(logJunePrec, scale=F) + scale(logAugPrec, scale=F) + Statelast -1 + (1|Patch)+(1|Year), family="binomial", data=PhoPc_dat)
summary(mPhoPc9c)

#000
s = summary(mPhoPc9c)$coef[7,1]
sse = summary(mPhoPc9c)$coef[7,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#100
s = summary(mPhoPc9c)$coef[8,1]
sse = summary(mPhoPc9c)$coef[8,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#110
s = summary(mPhoPc9c)$coef[9,1]
sse = summary(mPhoPc9c)$coef[9,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

# Mildew ext
PhoPe_dat = na.omit(subset(comb2, select=c("PhoPe","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Mil","logTotalNH","Road","Year","Patch")))

mPhoPe3a = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logJulyPrec+(1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
mPhoPe3b = glmer(PhoPe~logPLM2+logNH_Mil+logTotalNH+Road+logJulyPrec+ Statelast + (1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
AIC(mPhoPe3a, mPhoPe3b)
AIC(mPhoPe3a, mPhoPe3b)$AIC[2]-AIC(mPhoPe3a, mPhoPe3b)$AIC[1]
logLik(mPhoPe3a)
logLik(mPhoPe3b)

mPhoPe3c = glmer(PhoPe~scale(logPLM2, scale=F) + scale(logNH_Mil, scale=F) + scale(logTotalNH, scale=F) + scale(Road, scale=F) + 
                 scale(logJulyPrec, scale=F) + Statelast -1 + (1|Patch)+(1|Year), family="binomial", data=PhoPe_dat)
summary(mPhoPe3c)

#001
s = summary(mPhoPe3c)$coef[6,1]
sse = summary(mPhoPe3c)$coef[6,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#101
s = summary(mPhoPe3c)$coef[7,1]
sse = summary(mPhoPe3c)$coef[7,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#111
s = summary(mPhoPe3c)$coef[8,1]
sse = summary(mPhoPe3c)$coef[8,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

# Parasitoid col
CotPc_dat = na.omit(subset(comb2, select=c("CotPc","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPc_dat = CotPc_dat[CotPc_dat$Year!="2018",]

mCotPc2a = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec + (1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
mCotPc2b = glmer(CotPc~logPLM2+VS+logNH_Cot+Road+logJunePrec+Statelast + (1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
AIC(mCotPc2a, mCotPc2b)
AIC(mCotPc2a, mCotPc2b)$AIC[2]-AIC(mCotPc2a, mCotPc2b)$AIC[1]
logLik(mCotPc2a)
logLik(mCotPc2b)

mCotPc2c = glmer(CotPc~scale(logPLM2, scale=F) + scale(VS, scale=F) + scale(logNH_Cot, scale=F) + scale(Road, scale=F) + 
                scale(logJunePrec, scale=F) + Statelast - 1 + (1|Patch)+(1|Year), family="binomial", data=CotPc_dat)
summary(mCotPc2c)

#000
s = summary(mCotPc2c)$coef[6,1]
sse = summary(mCotPc2c)$coef[6,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#001
s = summary(mCotPc2c)$coef[7,1]
sse = summary(mCotPc2c)$coef[7,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#100
s = summary(mCotPc2c)$coef[8,1]
sse = summary(mCotPc2c)$coef[8,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#101
s = summary(mCotPc2c)$coef[9,1]
sse = summary(mCotPc2c)$coef[9,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

# Parasitoid ext
CotPe_dat = na.omit(subset(comb2, select=c("CotPe","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPe_dat = CotPe_dat[CotPe_dat$Year!="2018",]

mCotPe2aX = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
mCotPe2b = glmer(CotPe~logPLM2+VS+logNH_Cot+Road+logJunePrec+Statelast+(1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
AIC(mCotPe2a, mCotPe2b)
AIC(mCotPe2a, mCotPe2b)$AIC[2]-AIC(mCotPe2a, mCotPe2b)$AIC[1]
logLik(mCotPe2a)
logLik(mCotPe2b)

mCotPe2c = glmer(CotPe~scale(logPLM2, scale=F) + scale(VS, scale=F) + scale(logNH_Cot, scale=F)+
                 scale(Road, scale=F) + scale(logJunePrec, scale=F) +Statelast -1 + (1|Patch)+(1|Year), family="binomial", data=CotPe_dat)
summary(mCotPe2c)

#110
s = summary(mCotPe2c)$coef[6,1]
sse = summary(mCotPe2c)$coef[6,2]
round(1/(1+exp(-s))*100, 2) 
hround(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

#111
s = summary(mCotPe2c)$coef[7,1]
sse = summary(mCotPe2c)$coef[7,2]
round(1/(1+exp(-s))*100, 2) 
round(1/(1+exp(-s+1.96*sse))*100, 2)
round(1/(1+exp(-s-1.96*sse))*100, 2)

# Saving the 'b' models including patch state as a predictor
State_modList = list(mMelCcb, mMelCe14b, mPhoPc9b, mPhoPe3b, mCotPc2b, mCotPe2b)
names(State_modList) = c("mMelCcb", "mMelCe14b", "mPhoPc9b", "mPhoPe3b", "mCotPc2b", "mCotPe2b")
save(State_modList, file="fitted_mods/State_modList.RData")

######################################################
#### - Generating predicted population dynamics - ####
######################################################

TotalPatchDat = read.csv("data/TotalPatchDat.csv")
newdat = read.csv("data/colextdat.csv")
newdat$Patch = as.factor(newdat$Patch)

# Insert patch mean for missing PL data
pmeanPL = data.frame(tapply(newdat$PLM2, newdat$Patch, mean, na.rm=T))
pmeanPL$Patch = row.names(pmeanPL)
names(pmeanPL) = c("PL","Patch")
newPL = subset(newdat, select=c("Patch","PLM2"))
for(i in 1:nrow(newPL)){
  if(is.na(newPL$PLM2[i])){
    newPL$PLM2[i]=pmeanPL$PL[which(pmeanPL$Patch==newPL$Patch[i])]
  }
}
newdat$PLM2 = newPL$PLM2

# Insert patch mean for missing VS data
newdat$VS[which(newdat$VS>3)]=3

pmeanVS = data.frame(tapply(newdat$VS, newdat$Patch, mean, na.rm=T))
pmeanVS$Patch = row.names(pmeanVS)
names(pmeanVS) = c("VS","Patch")
newVS = subset(newdat, select=c("Patch","VS"))
for(i in 1:nrow(newVS)){
  if(is.na(newVS$VS[i])){
    newVS$VS[i]=pmeanVS$VS[which(pmeanVS$Patch==newVS$Patch[i])]
  }
}
newdat$VS = newVS$VS

comb = read.csv("data/combdat.csv")
comb = comb[,-which(colnames(comb)=="MelC")]

invlogit = function(x){1/(1+exp(-x))}

# Butterfly
load("fitted_mods/MelCc_modList.RData")
mMelCc = MelCc_modList[["mMelCc"]]
load("fitted_mods/MelCe_modList.RData")
mMelCe14 = MelCe_modList[["mMelCe14"]]

br = quantile(newdat$PLM2, c(seq(0,1,length.out=15)), na.rm=T)
test = cut(newdat$PLM2, breaks=br)
obs = tapply(newdat$MelC>0, test, sum, na.rm=T)/tapply(newdat$MelC>-1, test, sum, na.rm=T)

# Setup
years = unique(newdat$Year)
predN = matrix(NA, nrow=length(years)-1, ncol=100)
predP = matrix(NA, nrow=length(years)-1, ncol=100)
alldat = list()
cmod = mMelCc
emod = mMelCe14

for(y in 1:(length(years)-1)){
  year = years[y+1]
  print(paste("Year", year))
  
  lastyear = subset(newdat[newdat$Year==year-1,], select=c("Patch", "MelC"))
  lastyear$MelC = as.numeric(lastyear$MelC>0)
  thisyear = comb[comb$Year==year,]
  
  preddat = merge(lastyear, thisyear, by="Patch", all.x=T)
  preddat$PLcat = cut(preddat$PLM2, breaks=br)
  
  outdat = data.frame(Year=rep(paste(year), nrow(preddat)), Patch=preddat$Patch, PLcat=preddat$PLcat, Pred=rep(NA, nrow(preddat)))
  
  for(j in 1:100){
    
    #Col
    rr = MASS::mvrnorm(1, mu=summary(cmod)$coef[,1], Sigma=vcov(cmod))
    rYear = ranef(cmod)$Year[which(row.names(ranef(cmod)$Year)==year), 1]
    #rYear=0
    pd = data.frame(Intercept=rep(1, nrow(preddat)), subset(preddat, select = names(rr[-1])))
    rr[1] = rr[1] + rYear
    predc = invlogit(as.matrix(pd) %*% as.matrix(rr))
  
    #Ext
    rr = MASS::mvrnorm(1, mu=summary(emod)$coef[,1], Sigma=vcov(emod))
    rYear = ranef(emod)$Year[which(row.names(ranef(emod)$Year)==year), 1]
    #rYear=0
    pd = data.frame(Intercept=rep(1, nrow(preddat)), subset(preddat, select = names(rr[-1])))
    rr[1] = rr[1] + rYear
    prede = invlogit(as.matrix(pd) %*% as.matrix(rr))
  
  preddat$pred = rep(NA, nrow(preddat))
  preddat$pred[which(preddat$MelC==1)] = (1-prede[which(preddat$MelC==1)])
  preddat$pred[which(preddat$MelC==0)] = predc[which(preddat$MelC==0)]
  
  preddat$predBin = rbinom(nrow(preddat), 1, p=preddat$pred)
  predN[y,j] = sum(preddat$predBin, na.rm=T)
  predP[y,j] = sum(preddat$predBin, na.rm=T)/(sum(preddat$predBin>-1, na.rm=T))
  
  outdat[,j+3] = preddat$predBin
  }
  alldat[[y]] = outdat
  }

allCinxiaDat = alldat

save(allCinxiaDat, file = "allCinxiadat_pred.RData")

df = data.frame(Year = years[1:18]+1, 
                predN = apply(predN, 1, mean, na.rm=T),
                predNupper = apply(predN, 1, quantile,c(.025,.975), na.rm=T)[2,],
                predNlower = apply(predN, 1, quantile,c(.025,.975), na.rm=T)[1,],
                predP = apply(predP, 1, mean, na.rm=T),
                predPupper = apply(predP, 1, quantile,c(.025,.975), na.rm=T)[2,],
                predPlower = apply(predP, 1, quantile, c(.025,.975), na.rm=T)[1,])

colr = rgb(0,0,1, alpha=.5)

yearlyDat = read.csv("data/yearlyDat.csv")

x11()
par(mfrow=c(1,3))
plot(years, yearlyDat$MelC, type="b", pch=16, ylim=c(0,.5), xlab="", las=1, 
     ylab="Proportion of patches occupied", main="Melitaea cinxia")
xx = 2001:2018
polygon(c(xx, rev(xx)), c(df$predPlower[1:18], rev(df$predPupper[1:18])), col = colr, border = FALSE)
points(df$Year, df$predP, pch=16, type="b", col="darkblue")
points(years, yearlyDat$MelC, type="b", pch=16)

# Mildew
test = cut(newdat$PLM2, breaks=br)
obs = tapply(newdat$PhoP>0, test, sum, na.rm=T)/tapply(newdat$PhoP>-1, test, sum, na.rm=T)

load("fitted_mods/PhoPc_modList.RData")
mPhoPc9 = PhoPc_modList[["mPhoPc9"]]
load("fitted_mods/PhoPe_modList.RData")
mPhoPe3 = PhoPe_modList[["mPhoPe3"]]

# Setup
years = unique(newdat$Year)
predN = matrix(NA,nrow=length(years)-1, ncol=100)
predP = matrix(NA,nrow=length(years)-1, ncol=100)
alldat = list()
cmod = mPhoPc9
emod = mPhoPe3

a=Sys.time()
for(y in 1:(length(years)-1)){
  year = years[y+1]
  print(paste("Year", year))
  lastyear = subset(newdat[newdat$Year==year-1,], select=c("Patch","PhoP"))
  lastyear$PhoP = as.numeric(lastyear$PhoP>0)
  thisyear = comb[comb$Year==year,]
  
  preddat = merge(lastyear, thisyear, by="Patch", all.x=T)
  preddat$PLcat = cut(preddat$PLM2, breaks=br) 
  #preddat$PLcat=cut(preddat$logJulyPrec,breaks=br)
  outdat = data.frame(Year=rep(paste(year), nrow(preddat)), Patch=preddat$Patch, PLcat=preddat$PLcat, Pred=rep(NA, nrow(preddat)))
  
  
  for(j in 1:100){
    
    #Col
    rr = MASS::mvrnorm(1, mu=summary(cmod)$coef[,1], Sigma=vcov(cmod))
    rYear = ranef(cmod)$Year[which(row.names(ranef(cmod)$Year)==year),1]
    pd = data.frame(Intercept=rep(1, nrow(preddat)), subset(preddat, select=names(rr[-1])))
    rr[1] = rr[1]+rYear
    predc = invlogit(as.matrix(pd)%*%as.matrix(rr))
    
    #Ext
    rr = MASS::mvrnorm(1, mu=summary(emod)$coef[,1], Sigma=vcov(emod))
    rYear = ranef(emod)$Year[which(row.names(ranef(emod)$Year)==year), 1]
    pd = data.frame(Intercept=rep(1, nrow(preddat)), subset(preddat, select=names(rr[-1])))
    rr[1] = rr[1]+rYear
    prede = invlogit(as.matrix(pd)%*%as.matrix(rr))
    
    preddat$pred = rep(NA, nrow(preddat))
    preddat$pred[which(preddat$PhoP==1)] = (1-prede[which(preddat$PhoP==1)])
    preddat$pred[which(preddat$PhoP==0)] = predc[which(preddat$PhoP==0)]
    
    preddat$predBin = rbinom(nrow(preddat), 1, p=preddat$pred)
    predN[y,j] = sum(preddat$predBin, na.rm=T)
    predP[y,j]=sum(preddat$predBin, na.rm=T)/(sum(preddat$predBin>-1, na.rm=T))
    
    outdat[,j+3] = preddat$predBin
    }
  alldat[[y]] = outdat
}

allMildewDat = alldat
save(allMildewDat, file = "allMildewdat_pred.RData")

df=data.frame(Year = years[1:18]+1,
              predN = apply(predN, 1, mean, na.rm=T),
              predNupper = apply(predN, 1, quantile, c(.025,.975), na.rm=T)[2,],
              predNlower = apply(predN, 1, quantile, c(.025,.975), na.rm=T)[1,],
              predP = apply(predP, 1, mean, na.rm=T),
              predPupper = apply(predP, 1, quantile, c(.025,.975), na.rm=T)[2,],
              predPlower = apply(predP, 1, quantile, c(.025,.975), na.rm=T)[1,])

plot(years, yearlyDat$PhoPcomb, type="b", pch=16, ylim=c(0,.4), xlab="", 
     las=1, ylab="", main="Podosphaera plantaginis")
xx = 2001:2018
polygon(c(xx,rev(xx)), c(df$predPlower[1:18],rev(df$predPupper[1:18])), col = colr, border = FALSE)
points(df$Year, df$predP, pch=16, type="b", col="darkblue")
points(years, yearlyDat$PhoPcomb, type="b", pch=16)

# Cotesia
test = cut(newdat$PLM2, breaks=br)
obs = tapply(newdat$CotP>0, test, sum, na.rm=T)/tapply(newdat$CotP>-1, test, sum, na.rm=T)

load("fitted_mods/CotPc_modList.RData")
mCotPc2 = CotPc_modList[["mCotPc2"]]
load("fitted_mods/CotPe_modList.RData")
mCotPe2 = CotPe_modList[["mCotPe2"]]

# Setup
years = unique(newdat$Year)
predN = matrix(NA, nrow=length(years)-1, ncol=100)
predP = matrix(NA, nrow=length(years)-1, ncol=100)
alldat = list()
coldat = list()
extdat = list()
cmod = mCotPc2
emod = mCotPe2

for(y in 1:(length(years)-1)){
  year = years[y+1]
  print(paste("Year", year))
  lastyear = subset(newdat[newdat$Year==year-1,], select=c("Patch","CotP"))
  lastyear$CotP = as.numeric(lastyear$CotP>0)
  thisyear = comb[comb$Year==year,]
  
  preddat = merge(lastyear, thisyear, by="Patch", all.x=T)
  preddat$PLcat = cut(preddat$PLM2, breaks=br)
  
  outdat = data.frame(Year=rep(paste(year),nrow(preddat)), Patch=preddat$Patch, PLcat=preddat$PLcat,Pred=rep(NA,nrow(preddat)))
  cdat = data.frame(Year=rep(paste(year),nrow(preddat)), Patch=preddat$Patch, PLcat=preddat$PLcat,Pred=rep(NA,nrow(preddat)))
  edat = data.frame(Year=rep(paste(year),nrow(preddat)), Patch=preddat$Patch, PLcat=preddat$PLcat,Pred=rep(NA,nrow(preddat)))
  
  for(j in 1:100){

    #Col
    rr = MASS::mvrnorm(1, mu=summary(cmod)$coef[,1], Sigma=vcov(cmod))
    rYear = ranef(cmod)$Year[which(row.names(ranef(cmod)$Year)==year), 1]
    if(length(rYear)<1){rYear=0}
    pd = data.frame(Intercept=rep(1,nrow(preddat)), subset(preddat, select=names(rr[-1])))
    rr[1] = rr[1]+rYear
    predc = invlogit(as.matrix(pd) %*% as.matrix(rr))
    
    #Ext
    rr = MASS::mvrnorm(1, mu=summary(emod)$coef[,1], Sigma=vcov(emod))
    rYear = ranef(emod)$Year[which(row.names(ranef(emod)$Year)==year),1]
    if(length(rYear)<1){rYear=0}
    pd = data.frame(Intercept=rep(1, nrow(preddat)), subset(preddat, select=names(rr[-1])))
    rr[1] = rr[1] + rYear
    prede = invlogit(as.matrix(pd)%*%as.matrix(rr))
    
    preddat$pred = rep(NA, nrow(preddat))
    preddat$pred[which(preddat$CotP==1)] = (1-prede[which(preddat$CotP==1)])
    preddat$pred[which(preddat$CotP==0)] = predc[which(preddat$CotP==0)]
    
    preddat$predBin = rbinom(nrow(preddat), 1, p=preddat$pred)
    predN[y,j] = sum(preddat$predBin, na.rm=T)
    predP[y,j] = sum(preddat$predBin, na.rm=T)/length(preddat$predBin) #Scale by total patches surveyed
    
    outdat[,j+3] = preddat$predBin
    cdat[,j+3] = predc
    edat[,j+3] = prede
    
    #predN[y,j]=sum(preddat$pred,na.rm=T)
    #predP[y,j]=sum(preddat$pred,na.rm=T)/(sum(preddat$pred>-1,na.rm=T))
    #predP[y,j]=sum(preddat$pred,na.rm=T)/nrow(preddat) #Total Patches surveyed
      }
  alldat[[y]] = outdat
  coldat[[y]] = cdat
  extdat[[y]] = edat
  
}

allCotesiaDat = alldat
save(allCotesiaDat, file="allCotesiadat_pred.RData")
save(coldat, file="coldat_pred.RData")
save(extdat, file="extdat_pred.RData")

df=data.frame(Year = years[1:18]+1,
              predN = apply(predN, 1, mean, na.rm=T),
              predNupper = apply(predN, 1, quantile, c(.025,.975), na.rm=T)[2,],
              predNlower = apply(predN, 1, quantile, c(.025,.975), na.rm=T)[1,],
              predP = apply(predP, 1, mean, na.rm=T),
              predPupper = apply(predP, 1, quantile, c(.025,.975), na.rm=T)[2,],
              predPlower = apply(predP, 1, quantile, c(.025,.975), na.rm=T)[1,])

plot(years, yearlyDat$CotPAdj, type="b", pch=16, ylim=c(0,.06), xlab="", ylab="",
     las=1, main="Cotesia melitaearum")
xx = 2001:2018
polygon(c(xx,rev(xx)), c(df$predPlower[1:18],rev(df$predPupper[1:18])), col = colr, border = FALSE)
points(df$Year, df$predP, pch=16, type="b", col="darkblue")
points(years, yearlyDat$CotPAdj, type="b", pch=16)

######################################
#### - Plot population dynamics - ####
######################################
yearlyDat = read.csv("data/yearlyDat.csv")
newdat = read.csv("data/colextdat.csv")
newdat$Patch = as.factor(newdat$Patch)

load(file="allCinxiadat_pred.RData")
load(file="allMildewdat_pred.RData")
load(file="allCotesiadat_pred.RData")
load(file="coldat_pred.RData")
load(file="extdat_pred.RData")

cols = c(rgb(0,0,1,.5), rgb(0,1,0,.5), rgb(1,0,0,.5))
years = unique(newdat$Year)

# Cinxia
predP = matrix(NA, nrow=length(allCinxiaDat), ncol=100)
for(i in 1:(length(allCinxiaDat))){
  sums=apply(allCinxiaDat[[i]][,4:103], 2, sum, na.rm=T)
  ns=apply(allCinxiaDat[[i]][,4:103]>-1, 2, sum, na.rm=T)
  predP[i,]=sums/ns
}

df = data.frame(Year=years[1:18]+1,
                predP=apply(predP, 1, mean, na.rm=T),
                predPupper=apply(predP, 1, quantile, c(.025,.975), na.rm=T)[2,],
                predPlower=apply(predP, 1, quantile, c(.025,.975), na.rm=T)[1,])

x11(height=5, width=10)
par(mfrow=c(1,2), mar=c(4,4,2,4))

pdf("timeseriesfig_ab.pdf", height=5, width=10, family="Times")
par(mfrow=c(1,2), mar=c(4,4,2,4))

plot(years, yearlyDat$MelC, type="b", pch=16, ylim=c(0,.4), 
     xlab="", las=1, ylab="", main="")
xx=2001:2018
polygon(c(xx,rev(xx)), c(df$predPlower[1:18],rev(df$predPupper[1:18])), col = cols[1], border = FALSE)
#points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years, yearlyDat$MelC, type="b", pch=16)

mtext("Proportion of patches occupied", 2, line=2.5)

# Mildew
nc = ncol(allMildewDat[[1]])
predP = matrix(NA, nrow=length(allMildewDat), ncol=100)
for(i in 1:(length(allMildewDat))){
  sums=apply(allMildewDat[[i]][,4:nc], 2, sum, na.rm=T)
  ns=apply(allMildewDat[[i]][,4:nc]>-1, 2, sum, na.rm=T)
  predP[i,]=sums/ns
}

df = data.frame(Year=years[1:18]+1,
              predP=apply(predP, 1, mean, na.rm=T),
              predPupper=apply(predP, 1, quantile, c(.025,.975), na.rm=T)[2,],
              predPlower=apply(predP, 1, quantile, c(.025,.975), na.rm=T)[1,])

xx = 2001:2018
polygon(c(xx,rev(xx)), c(df$predPlower[1:18],rev(df$predPupper[1:18])), col = cols[2], border = FALSE)
#points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years, yearlyDat$PhoPcomb, type="b", pch=16)

# Cotesia
predP = matrix(NA, nrow=length(allCotesiaDat), ncol=100)
for(i in 1:(length(allCotesiaDat))){
  sums=apply(allCotesiaDat[[i]][,4:103], 2, sum, na.rm=T)
  ns=apply(allCotesiaDat[[i]][,4:103]>-1, 2, length)
  #nsCinxia = apply(allCinxiaDat[[i]][,4:103], 2, sum, na.rm=T)
  predP[i,]=sums/ns
  #predP[i,]=sums/nsCinxia
  }

df = data.frame(Year=years[1:18]+1,
                predP=apply(predP, 1, mean, na.rm=T),
                predPupper=apply(predP, 1, quantile, c(.025,.975), na.rm=T)[2,],
                predPlower=apply(predP, 1, quantile, c(.025,.975), na.rm=T)[1,])

xx1 = 2001:2008
xx2 = 2011:2018
polygon(c(xx1,rev(xx1)), c(df$predPlower[1:8],rev(df$predPupper[1:8])), col = cols[3], border = FALSE)
polygon(c(xx2,rev(xx2)), c(df$predPlower[11:18],rev(df$predPupper[11:18])), col = cols[3], border = FALSE)
#points(df$Year,df$predP,pch=16,type="b",col="darkblue")
points(years, yearlyDat$CotPAdj, type="b", pch=16)

text(1995.5, 0.4, labels="(a)", xpd=T)

# Plantago
par(new=T)
plot(yearlyDat$Year, yearlyDat$mPLM2, type="l", lwd=3, col="darkgreen", xaxt="n",
     yaxt="n", xlab="", ylab="")
axis(4, c(5,10,15,20,25), c(5,10,15,20,25), las=1)
mtext(4, text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")), line=2.5)

legend("topleft", col=c("darkgreen", cols), lwd=3, 
       legend=c("Plantago lanceolata", "Butterfly", "Mildew", "Parasitoid"), bty="n")

#plot(yearlyDat$Year,yearlyDat$mVS,type="l",lwd=3,col="darkgreen")
#dev.off()

# Cotesia relative to M. cinxia
predP = matrix(NA, nrow=length(allCotesiaDat), ncol=100)
for(i in 1:(length(allCotesiaDat))){
  sums=apply(allCotesiaDat[[i]][,4:103], 2, sum, na.rm=T)
  ns=apply(allCotesiaDat[[i]][,4:103]>-1, 2, length)
  nsCinxia=apply(allCinxiaDat[[i]][,4:103], 2, sum, na.rm=T)
  predP[i,]=sums/ns
  predP[i,]=sums/nsCinxia
}

df = data.frame(Year=years[1:18]+1,
                predP=apply(predP, 1, mean, na.rm=T),
                predPupper=apply(predP, 1, quantile, c(.025,.975), na.rm=T)[2,],
                predPlower=apply(predP, 1, quantile, c(.025,.975), na.rm=T)[1,])


obs = tapply(newdat$CotP>0, newdat$Year, sum, na.rm=T)/tapply(newdat$MelC>0, newdat$Year, sum, na.rm=T)
xx1 = 2001:2008
xx2 = 2011:2018

plot(years, yearlyDat$MelC, type="b", pch=16, ylim=c(0,.4), xlab="", las=1,
     ylab="", main="", lwd=2, col=cols[1])
xx = 2001:2018
mtext("Proportion of patches occupied", 2, line=2.5)

polygon(c(xx1,rev(xx1)), c(df$predPlower[1:8], rev(df$predPupper[1:8])),
        col = cols[3], border = FALSE)
polygon(c(xx2,rev(xx2)), c(df$predPlower[11:18], rev(df$predPupper[11:18])),
        col = cols[3], border = FALSE)
points(2000:2008, obs[1:9], pch=16, type="b")
points(2011:2017, obs[11:17], type="b", pch=16)

legend("topleft", col=c(cols[c(1,3)]), lwd=3, legend=c("Butterfly", "Parasitoid"), bty="n")
text(1995.5, 0.4, labels="(b)", xpd=T)
dev.off()

###################################
#### - Community predictions - ####
###################################

allCinxiaDat[[1]][1:10, 1:10]
allCotesiaDat[[1]][1:10, 1:10]
allMildewDat[[1]][1:10, 1:10]

dist = matrix(NA, ncol=6, nrow=100)
medians = matrix(NA, ncol=6, nrow=18)
upper = matrix(NA, ncol=6, nrow=18)
lower = matrix(NA, ncol=6, nrow=18)
w = list()

for(y in 1:(length(years)-1)){
  for(i in 1:100){
mc = allCinxiaDat[[y]][,c(1:2,i+3)]
names(mc) = c("Year","Patch","MelC")
cm = allCotesiaDat[[y]][,c(1:2,i+3)]
names(cm) = c("Year","Patch","CotP")
pp = allMildewDat[[y]][,c(1:2,i+3)]
names(pp) = c("Year","Patch","PhoP")

rep = cbind(mc,cm,pp)
rep = rep[,c(1:3,6,9)]
head(rep, 10)

rep$CotP[which(rep$MelC==0)]=NA #'Butterfly extinction'
rep$CotP[which(is.na(rep$MelC))]=NA #Missing data

Mel_col_patches=which(rep$MelC==1 & is.na(rep$CotP))
rep$CotP[Mel_col_patches] = rbinom(length(Mel_col_patches),1,coldat[[y]][Mel_col_patches,c(i+3)]) #Joint colonisations
head(rep, 20)

rep$State=paste(as.numeric(rep$MelC>0), as.numeric(rep$CotP>0), as.numeric(rep$PhoP>0), sep="_")

rep$State=factor(rep$State, levels=c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1"))

dist[i,]=table(rep$State)
}
  a = apply(dist, 2, quantile, c(.025,.5,.975))
  medians[y,] = a[2,]
  upper[y,] = a[3,]
  lower[y,] = a[1,]
}

colnames(medians)=colnames(upper)=colnames(lower)=levels(rep$State)
rownames(medians)=rownames(upper)=rownames(lower)=years[-1]
wide = read.csv("data/wide.csv")
rownames(wide) = 2000:2017
colnames(wide) = sub("X", "", colnames(wide))

mm = match(colnames(wide), colnames(medians))
medians = medians[,mm]
upper = upper[,mm]
lower = lower[,mm]
dim(wide)
dim(medians)

head(wide, 20) #Observed data, 2009=NA
head(medians, 20)

rs = rowSums(medians)
wide = t(apply(wide, 1, function(x){x/sum(x)}))
medians = t(apply(medians, 1, function(x){x/sum(x)}))

upper2 = upper
for(i in 1:nrow(upper2)){upper2[i,]=upper[i,]/rs[i]}
upper = upper2

lower2 = lower
for(i in 1:nrow(lower2)){lower2[i,]=lower[i,]/rs[i]}
lower = lower2

xx = 2000:2017
colr = rgb(0, 0, 1, alpha=.5)

x11()
par(mfrow=c(2,3))
for(i in 1:6){
  plot(xx+1, medians[,i], type = "n", ylab="Proportion of patches", las=1, xlab="", xlim=c(2000,2018),
       ylim=c(min(c(lower[,i], wide[,i]), na.rm=T), max(c(upper[,i], wide[,i]), na.rm=T)), main=colnames(wide)[i])
  polygon(c(xx+1,rev(xx+1)), c(lower[,i],rev(upper[,i])), col = colr, border = FALSE)
  
  lines(rownames(wide), wide[,i], lwd=2)
  lines(xx+1, medians[,i], lwd=2, col="darkblue")
  }

medians
ord = order(medians[,1])
barplot(t(medians[,]), legend=F, col=topo.colors(6))


############################################################
#### - Community predictions for the Plantago gradient- ####
############################################################

allCinxiaDat[[1]][1:10, 1:10]
allCotesiaDat[[1]][1:10, 1:10]
allMildewDat[[1]][1:10, 1:10]
dim(allCinxiaDat[[1]])
dim(allMildewDat[[1]])
dim(allCotesiaDat[[1]])

dist = matrix(NA, ncol=6, nrow=100)
medians = matrix(NA, ncol=6, nrow=18)
upper = matrix(NA, ncol=6, nrow=18)
lower = matrix(NA, ncol=6, nrow=18)
w = list()

medianlist = list()
lowerlist = list()
upperlist = list()

y=j=i=1
y=3
for(y in 1:(length(years)-1)){
  for(i in 1:100){
    mc = allCinxiaDat[[y]][,c(1:3,i+3)]
    names(mc) = c("Year","Patch","PLcat","MelC")
    cm = allCotesiaDat[[y]][,c(1:3,i+3)]
    names(cm) = c("Year","Patch","PLcat","CotP")
    pp = allMildewDat[[y]][,c(1:3,i+3)]
    names(pp) = c("Year","Patch","PLcat","PhoP")
    
    rep = cbind(mc,cm,pp)
    rep = rep[,c(1:4,8,12)]
    head(rep, 10)
    
    rep$CotP[which(rep$MelC==0)]=NA #'Butterfly extinction'
    rep$CotP[which(is.na(rep$MelC))]=NA #Missing data
    
    Mel_col_patches = which(rep$MelC==1 & is.na(rep$CotP))
    rep$CotP[Mel_col_patches] = rbinom(length(Mel_col_patches),1,coldat[[y]][Mel_col_patches,c(i+3)]) #Joint colonisations
    head(rep, 20)
    
    rep$State = paste(as.numeric(rep$MelC>0), as.numeric(rep$CotP>0), as.numeric(rep$PhoP>0), sep="_")
    
    rep$State = factor(rep$State, levels=c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1"))
    
    l = lapply(tapply(rep$State, rep$PLcat, table), function(x) as.data.frame(x))
    
    for(ii in 1:length(l)){
      if(dim(l[[ii]])[1]<1){
          l[[ii]]=data.frame(Var1=factor(c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1")), Freq=rep(NA,6))}}
    
    w[[i]] = sapply(l, function(x) x[,2])
    rownames(w[[i]]) = c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1")
    
    dist[i,] = table(rep$State)
  }
  a = apply(dist, 2, quantile, c(.025,.5,.975))
  medians[y,] = a[2,]
  upper[y,] = a[3,]
  lower[y,] = a[1,]
  medianlist[[y]] = apply(simplify2array(w), 1:2, mean, na.rm=T)
  upperlist[[y]] = apply(simplify2array(w), 1:2, quantile, .975, na.rm=T)
  lowerlist[[y]] = apply(simplify2array(w), 1:2, quantile, .025, na.rm=T)
  }


temp = medianlist
temp = lapply(medianlist, function(x) apply(x, 2, function(xx){xx/sum(xx)}))

tempmean = apply(simplify2array(temp), 1:2, mean, na.rm=T)
upper = apply(simplify2array(temp), 1:2, quantile, .975, na.rm=T)
lower = apply(simplify2array(temp), 1:2, quantile, .025, na.rm=T)

# Observed
newdat$PLcat = cut(newdat$PLM2, breaks=br) #br from TotalPatchDat
newdat$State = factor(newdat$State, levels=c("0_NA_0","1_0_0","1_1_0","1_1_1","0_NA_1","1_0_1"))
tapply(newdat$State, newdat$PLcat, table)
l = lapply(tapply(newdat$State, newdat$PLcat, table), function(x) as.data.frame(x))
obs = sapply(l, function(x) x[,2])
obs = apply(obs, 2, function(x){x/sum(x)})

obs = obs[c(1,5,2,6,3,4),]
tempmean = tempmean[c(1,5,2,6,3,4),]
upper = upper[c(1,5,2,6,3,4),]
lower = lower[c(1,5,2,6,3,4),]

titles = rownames(tempmean)
titles = c("No species", "Mildew only", "Butterfly only", "Butterfly & Mildew", "Butterfly & Parasitoid", "All species")

br2=NULL #Midpoints
for(i in 1:(length(br)-1)){
  br2[i]=(br[i]+br[i+1])/2
}

xx = br2
xx = log(br2)

x11()
par(mfrow=c(2,3), mar=c(2,2,2,2), oma=c(3,3,0,0))
for(i in 1:6){
  plot(xx, tempmean[i,], type = "n",ylab="",las=1, xlab="",xaxt="n",
       ylim=c(min(lower[i,]), max(upper[i,])), main=as.expression(titles[i]))
  axis(1, seq(-4,8,2), signif(exp(seq(-4,8,2)), 1))
  arrows(xx, upper[i,], xx, lower[i,], length=0, col="darkgrey")
  
  #polygon(c(xx+1,rev(xx+1)),c(lower[,i],rev(upper[,i])),col = colr, border = FALSE)
  points(xx, tempmean[i,], pch=16, col="darkblue")
  #points(xx, obs[i,],pch=16, col="black")
 }
mtext(1, text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")), line=1, outer=T)
mtext(2, text="Proportion of patches", line=1, outer=T)

pdf("barplotfig.pdf", height=4, width=7.25, family="Times")
x11()
par(mar=c(4,4,2,9))
barplot(tempmean, col=topo.colors(6), legend.text = titles, axisnames = F, las=1, ylab="",
        args.legend = list(x=24.25, bty="n", xpd=T), ylim=c(.3,1), xpd=F)
axis(1, seq(.75,16.25, length=14), labels=signif(exp(xx), 2))
mtext(1, text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")), line=2.5, cex=1.1)
mtext(2, text="Proporton of patches", line=2.5, cex=1.1)

dev.off()

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

# Datasets
#newdat = read.csv("colextdat.csv")
comb2 = read.csv("data/combdat2.csv")

MelCc_dat = na.omit(subset(comb2, select=c("MelCc","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH","logTotalNH","Road","Year","Patch")))
MelCe_dat = na.omit(subset(comb2, select=c("MelCe","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH","logTotalNH","Road","Year","Patch")))
PhoPc_dat = na.omit(subset(comb2, select=c("PhoPc","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Mil","logTotalNH","Road","Year","Patch")))
PhoPe_dat = na.omit(subset(comb2, select=c("PhoPe","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Mil","logTotalNH","Road","Year","Patch")))
CotPc_dat = na.omit(subset(comb2, select=c("CotPc","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPc_dat = CotPc_dat[CotPc_dat$Year!="2018",]
CotPe_dat = na.omit(subset(comb2, select=c("CotPe","Statelast","logMayPrec","logJunePrec","logJulyPrec","logAugPrec","logArea","logPLM2","VS","logNH_Cot","Road","Year","Patch")))
CotPe_dat = CotPe_dat[CotPe_dat$Year!="2018",]

# Models
load("fitted_mods/State_modList.RData")

#x11()
pdf("hostplanteffectfig.pdf", height=5, width=8, family="Times")

# Colonisation
par(mfrow=c(2,4),mar=c(2,2,2,2),oma=c(2,2,0,0))

# Butterly colonisation
head(MelCc_dat)
cmod = State_modList$mMelCcb
summary(cmod)

newdata = data.frame(logPLM2 = seq(min(MelCc_dat$logPLM2, na.rm=T), max(MelCc_dat$logPLM2, na.rm=T),length.out=1000),
                     VS = rep(mean(MelCc_dat$VS, na.rm=T), 1000),
                     logNH = rep(mean(MelCc_dat$logNH, na.rm=T), 1000),
                     Road = rep(mean(MelCc_dat$Road, na.rm=T), 1000),
                     logMayPrec = rep(mean(MelCc_dat$logMayPrec, na.rm=T), 1000),
                     Statelast = rep("0_NA_0", 1000))
levels(newdata$Statelast)=levels(factor(MelCc_dat$Statelast))

preds = predict(cmod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(cmod, newdata=newdata)

xx = newdata$logPLM2
plot(xx, preds, type="l", lwd=2, ylim=c(0,1), las=1, col="white", xaxt="n",
     xlab="", ylab="", main="Butterfly")
polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,0,1,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1, c(0,2,4,6,8), signif(exp(c(0,2,4,6,8)), 1))

newdata$Statelast = rep("0_NA_1", 1000)
newdata$Statelast = factor(newdata$Statelast)
levels(newdata$Statelast) = rev(levels(factor(MelCc_dat$Statelast)))

preds = predict(cmod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(cmod, newdata=newdata)

polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,1,0,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

legend("topleft", c("Mildew present", "Mildew absent"), col=c(rgb(0,1,0), 
                                        rgb(0,0,1)), lty=1, lwd=2, bty="n")

mtext(2, text="Colonisation probability", line=2.5, outer=F, xpd=T)

# Cotesia colonisation
head(CotPc_dat)
cmod = State_modList$mCotPc2b
summary(cmod)


newdata=data.frame(logPLM2 = seq(min(CotPc_dat$logPLM2, na.rm=T), max(CotPc_dat$logPLM2, na.rm=T), length.out=1000),
                   VS = rep(mean(CotPc_dat$VS, na.rm=T), 1000),
                   logNH_Cot = rep(mean(CotPc_dat$logNH_Cot, na.rm=T), 1000),
                   Road = rep(mean(CotPc_dat$Road, na.rm=T), 1000),
                   logJunePrec = rep(mean(CotPc_dat$logJunePrec, na.rm=T), 1000),
                   Statelast = rep("0_NA_0", 1000))
levels(newdata$Statelast) = levels(factor(CotPc_dat$Statelast))
preds = predict(cmod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(cmod, newdata=newdata)

xx = newdata$logPLM2
plot(xx, preds, type="l",lwd=2, ylim=c(0,.1), las=1, col="white", xaxt="n",
     xlab="", ylab="", main="Parasitoid")
polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,0,1,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1, c(0,2,4,6,8), signif(exp(c(0,2,4,6,8)), 1))

newdata$Statelast = rep("1_0_0", 1000)
newdata$Statelast = factor(newdata$Statelast)
levels(newdata$Statelast) = levels(factor(CotPc_dat$Statelast))[c(3,1,2,4)]

preds = predict(cmod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(cmod, newdata=newdata)

polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,1,0,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

legend("topleft", c("Butterfly present", "Butterfly absent"), col=c(rgb(0,1,0), rgb(0,0,1)), lty=1, lwd=2, bty="n")

# Cotesia colonisation, mildew effect
cmod = State_modList$mCotPc2b
summary(cmod)

newdata=data.frame(logPLM2 = seq(min(CotPc_dat$logPLM2, na.rm=T), max(CotPc_dat$logPLM2, na.rm=T), length.out=1000),
                   VS = rep(mean(CotPc_dat$VS, na.rm=T), 1000),
                   logNH_Cot = rep(mean(CotPc_dat$logNH_Cot, na.rm=T), 1000),
                   Road = rep(mean(CotPc_dat$Road, na.rm=T), 1000),
                   logJunePrec = rep(mean(CotPc_dat$logJunePrec, na.rm=T), 1000),
                   Statelast = rep("0_NA_0", 1000))
levels(newdata$Statelast) = levels(factor(CotPc_dat$Statelast))
preds = predict(cmod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(cmod, newdata=newdata)

xx = newdata$logPLM2
plot(xx, preds, type="l", lwd=2, ylim=c(0,.1), las=1, col="white", xaxt="n",
     xlab="", ylab="", main="Parasitoid")
polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,0,1,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1,c(0,2,4,6,8), signif(exp(c(0,2,4,6,8)), 1))

newdata$Statelast = rep("0_NA_1", 1000)
newdata$Statelast = factor(newdata$Statelast)
levels(newdata$Statelast) = levels(factor(CotPc_dat$Statelast))[c(2,1,3,4)]

preds = predict(cmod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(cmod, newdata=newdata)

polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,1,0,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

legend("topleft", c("Mildew present", "Mildew absent"), col=c(rgb(0,1,0), rgb(0,0,1)),
       lty=1, lwd=2, bty="n")

# Mildew colonisation
head(PhoPc_dat)
cmod = State_modList$mPhoPc9b
summary(cmod)

newdata=data.frame(logPLM2 = seq(min(PhoPc_dat$logPLM2, na.rm=T), max(PhoPc_dat$logPLM2, na.rm=T), length.out=1000),
                   logNH_Mil = rep(mean(PhoPc_dat$logNH_Mil, na.rm=T), 1000),
                   logTotalNH = rep(mean(PhoPc_dat$logTotalNH, na.rm=T), 1000),
                   Road = rep(mean(PhoPc_dat$Road, na.rm=T), 1000),
                   logJunePrec = rep(mean(PhoPc_dat$logJunePrec, na.rm=T), 1000),
                   logAugPrec = rep(mean(PhoPc_dat$logAugPrec, na.rm=T), 1000),
                   Statelast = rep("0_NA_0", 1000))
levels(newdata$Statelast) = levels(factor(PhoPc_dat$Statelast))
head(newdata)

preds = predict(cmod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(cmod, newdata=newdata)

xx = newdata$logPLM2
plot(xx, preds, type="l", lwd=2, ylim=c(0,1), las=1, col="white", xaxt="n",
     xlab="", ylab="",main="Mildew")
polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,0,1,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1, c(0,2,4,6,8), signif(exp(c(0,2,4,6,8)), 1))

newdata$Statelast = rep("1_0_0", 1000)
newdata$Statelast = factor(newdata$Statelast)
levels(newdata$Statelast) = levels(factor(PhoPc_dat$Statelast))[c(2,1,3)]

preds = predict(cmod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(cmod, newdata=newdata)

polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,1,0,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

legend("topleft", c("Butterfly present", "Butterfly absent"), col=c(rgb(0,1,0), rgb(0,0,1)), 
       lty=1, lwd=2, bty="n")

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=0.5, outer=T)
#mtext(2,text="Colonisation probability",line=0.5, outer=T)



# Extinction
#x11()
#par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,0,0))

# Butterly extinction: Mildew effects
head(MelCe_dat)
emod = State_modList$mMelCe14b
summary(emod)

newdata=data.frame(logPLM2 = seq(min(MelCe_dat$logPLM2, na.rm=T), max(MelCe_dat$logPLM2, na.rm=T), length.out=1000),
                   VS = rep(mean(MelCe_dat$VS, na.rm=T), 1000),
                   logNH = rep(mean(MelCe_dat$logNH, na.rm=T), 1000),
                   Road = rep(mean(MelCe_dat$Road, na.rm=T), 1000),
                   logJunePrec = rep(mean(MelCe_dat$logJunePrec, na.rm=T), 1000),
                   logJulyPrec = rep(mean(MelCe_dat$logJulyPrec, na.rm=T), 1000),
                   logAugPrec = rep(mean(MelCe_dat$logAugPrec, na.rm=T), 1000),
                   Statelast = rep("1_0_0", 1000))
levels(newdata$Statelast) = levels(factor(MelCe_dat$Statelast))
head(newdata)

preds = predict(emod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(emod, newdata=newdata)

xx = newdata$logPLM2
plot(xx, preds, type="l", lwd=2, ylim=c(0,1), las=1, col="white", xaxt="n",
     xlab="", ylab="", main="Butterfly")
polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,0,1,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1, c(0,2,4,6,8), signif(exp(c(0,2,4,6,8)), 1))

newdata$Statelast = rep("1_0_1", 1000)
newdata$Statelast = factor(newdata$Statelast)
levels(newdata$Statelast) = levels(factor(MelCe_dat$Statelast))[c(2,1,3,4)]

preds = predict(emod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(emod, newdata=newdata)

polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,1,0,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

legend("topleft", c("Mildew present", "Mildew absent"), col=c(rgb(0,1,0), rgb(0,0,1)),
       lty=1, lwd=2, bty="n")

mtext(2, text="Extinction probability", line=2.5, outer=F, xpd=T)

# Butterly extinction: Cotesia effects
emod = State_modList$mMelCe14b
summary(emod)

newdata=data.frame(logPLM2 = seq(min(MelCe_dat$logPLM2, na.rm=T), max(MelCe_dat$logPLM2, na.rm=T), length.out=1000),
                   VS = rep(mean(MelCe_dat$VS, na.rm=T), 1000),
                   logNH = rep(mean(MelCe_dat$logNH, na.rm=T), 1000),
                   Road = rep(mean(MelCe_dat$Road, na.rm=T), 1000),
                   logJunePrec = rep(mean(MelCe_dat$logJunePrec, na.rm=T), 1000),
                   logJulyPrec = rep(mean(MelCe_dat$logJulyPrec, na.rm=T), 1000),
                   logAugPrec = rep(mean(MelCe_dat$logAugPrec, na.rm=T), 1000),
                   Statelast = rep("1_0_0", 1000))
levels(newdata$Statelast) = levels(factor(MelCe_dat$Statelast))
head(newdata)

preds = predict(emod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(emod, newdata=newdata)

xx = newdata$logPLM2
plot(xx, preds, type="l", lwd=2, ylim=c(0,1), las=1, col="white", xaxt="n",
     xlab="", ylab="", main="Butterfly")
polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,0,1,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1, c(0,2,4,6,8), signif(exp(c(0,2,4,6,8)), 1))

newdata$Statelast = rep("1_1_0", 1000)
newdata$Statelast = factor(newdata$Statelast)
levels(newdata$Statelast) = levels(factor(MelCe_dat$Statelast))[c(3,1,2,4)]

preds = predict(emod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(emod, newdata=newdata)

polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,1,0,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

legend("topleft", c("Parasitoid present", "Parasitoid absent"), col=c(rgb(0,1,0), rgb(0,0,1)),
       lty=1, lwd=2, bty="n")

# Cotesia extinction: Mildew effect
head(CotPe_dat)
emod = State_modList$mCotPe2b
summary(emod)


newdata=data.frame(logPLM2 = seq(min(CotPe_dat$logPLM2, na.rm=T), max(CotPe_dat$logPLM2, na.rm=T), length.out=1000),
                   VS = rep(mean(CotPe_dat$VS, na.rm=T), 1000),
                   logNH_Cot = rep(mean(CotPe_dat$logNH_Cot, na.rm=T), 1000),
                   Road = rep(mean(CotPe_dat$Road, na.rm=T), 1000),
                   logJunePrec = rep(mean(CotPe_dat$logJunePrec, na.rm=T), 1000),
                   Statelast = rep("1_1_0", 1000))
levels(newdata$Statelast) = levels(factor(CotPe_dat$Statelast))
preds = predict(emod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(emod, newdata=newdata)

xx = newdata$logPLM2
plot(xx, preds, type="l", lwd=2, ylim=c(0,1), las=1, col="white", xaxt="n", main="Parasitoid",
     xlab="", ylab="")
polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,0,1,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1, c(0,2,4,6,8), signif(exp(c(0,2,4,6,8)), 1))

newdata$Statelast = rep("1_1_1", 1000)
newdata$Statelast = factor(newdata$Statelast)
levels(newdata$Statelast) = levels(factor(CotPe_dat$Statelast))[c(2,1)]

preds = predict(emod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(emod, newdata=newdata)

polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,1,0,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

legend("topleft", c("Mildew present", "Mildew absent"), col=c(rgb(0,1,0), rgb(0,0,1)),
       lty=1, lwd=2, bty="n")

# Mildew extinction: Butterfly effect
head(PhoPe_dat)
emod = State_modList$mPhoPe3b
summary(emod)

newdata=data.frame(logPLM2 = seq(min(PhoPe_dat$logPLM2, na.rm=T), max(PhoPe_dat$logPLM2, na.rm=T), length.out=1000),
                   logNH_Mil = rep(mean(PhoPe_dat$logNH_Mil, na.rm=T), 1000),
                   logTotalNH = rep(mean(PhoPe_dat$logTotalNH, na.rm=T), 1000),
                   Road = rep(mean(PhoPe_dat$Road, na.rm=T), 1000),
                   logJulyPrec = rep(mean(PhoPe_dat$logJulyPrec, na.rm=T), 1000),
                   Statelast = rep("0_NA_1", 1000))
levels(newdata$Statelast) = levels(factor(PhoPe_dat$Statelast))
head(newdata)

preds = predict(emod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(emod, newdata=newdata)

xx = newdata$logPLM2
plot(xx, preds, type="l", lwd=2, ylim=c(0,1), las=1, col="white", xaxt="n",
     xlab="", ylab="", main="Mildew")
polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,0,1,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

#mtext(1,text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")),line=2.5)
axis(1, c(0,2,4,6,8), signif(exp(c(0,2,4,6,8)), 1))

newdata$Statelast = rep("1_0_1", 1000)
newdata$Statelast = factor(newdata$Statelast)
levels(newdata$Statelast) = levels(factor(PhoPe_dat$Statelast))[c(2,1,3)]

preds = predict(emod, newdata=newdata, re.form=NA, type="response")
predCI = easyPredCI(emod, newdata=newdata)

polygon(c(xx,rev(xx)), c(predCI[,1],rev(predCI[,2])), col = rgb(0,1,0,.5), border = FALSE)
lines(xx, preds, type="l", lwd=2)

legend("topleft", c("Butterfly present", "Butterfly absent"), col=c(rgb(0,1,0), rgb(0,0,1)),
       lty=1, lwd=2, bty="n")

mtext(1, text=expression(paste(italic(Plantago)," ",italic(lanceolata)," cover (", m^2,")")), line=1, outer=T)
#mtext(2,text="Extinction probability",line=0.5, outer=T)
dev.off()
