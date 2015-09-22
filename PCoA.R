# Principal Coordinate Analysis
library(vegan)
library(rgl)
library(ape)
dat=read.csv("VSDandPVALS_lunar.csv")
names(dat)
head(dat)
data=dat[,2:23] # Change to fit your data
row.names(data)=data$X
head(data)

#-------------Set experimental conditions

names(data)
indiv= c("F1", "F1", "F1", "F1", "F1", "F1", "F1", "F2", "F2", "F2", "F2", "F2", "F2", "F2", "F3", "F3", "F3","F3", "F3", "F3")  
lunar=c("3Q", "FM", "NM", "1Q","3Q", "FM", "NM", "3Q", "FM", "NM", "1Q", "3Q", "FM", "NM", "FM", "NM", "1Q", "3Q", "FM", "NM")
dayTime=c("day","day","day", "night", "night", "night", "night", "day","day","day", "night", "night", "night", "night", "day","day", "night", "night", "night", "night")
conditions=data.frame(cbind(lunar,indiv,dayTime))
conditions

#-------------Calulate principal coordinates
dd.veg=vegdist(t(data), "manhattan")
div.dd.veg=dd.veg/1000
head(div.dd.veg)

dd.pcoa=pcoa(div.dd.veg) 
head(dd.pcoa)
scores=dd.pcoa$vectors

#-------------First and second axes
plot(scores[,1], scores[,2],col=as.numeric(conditions$lunar), pch=16)
ordihull(scores,lunar,label=T)

plot(scores[,1], scores[,2],col=as.numeric(conditions$dayTime), pch=16)
ordihull(scores,dayTime,label=T)

plot(scores[,1], scores[,2],col=as.numeric(conditions$indiv), pch=16)
ordihull(scores,indiv,label=T)

#-------------Second and third axes
plot(scores[,2], scores[,3],col=as.numeric(conditions$lunar),pch=16)
ordihull(scores[,2:3],lunar,label=T)

plot(scores[,2], scores[,3],col=as.numeric(conditions$dayTime),pch=16)
ordihull(scores[,2:3],dayTime,label=T)

plot(scores[,2], scores[,3],col=as.numeric(conditions$indiv),pch=16)
ordihull(scores[,2:3],indiv,label=T)

#-------------PERMANOVA
adonis(t(data)~lunar+indiv+dayTime,data=conditions,method="manhattan")  
#adonis(t(data)~lunar,data=conditions,method="manhattan")  

pco1=scores[,1]
TukeyHSD(aov(pco1~lunar))
TukeyHSD(aov(pco1~indiv))
TukeyHSD(aov(pco1~dayTime))

##########-----Unadjustd pvalue < 0.05 Only------##############

#-------------Subset for significant
head(dat)
AHdat=row.names(dat[dat$pvalAH<0.05 & !is.na(dat$pvalAH),])
DHdat=row.names(dat[dat$pvalDH<0.05 & !is.na(dat$pvalDH),])
DAdat=row.names(dat[dat$pvalDA<0.05 & !is.na(dat$pvalDA),])
NMdat=row.names(dat[dat$pvalNM<0.05 & !is.na(dat$pvalNM),])
DDdat=row.names(dat[dat$pvalDD<0.05 & !is.na(dat$pvalDD),])
DVdat=row.names(dat[dat$pvalDV<0.05 & !is.na(dat$pvalDV),])
DNdat=row.names(dat[dat$pvalDN<0.05 & !is.na(dat$pvalDN),])

sdat=union(DHdat,AHdat)
sdat=union(sdat,DAdat)
sdat=union(sdat,NMdat)
sdat=union(sdat,DDdat)
sdat=union(sdat,DVdat)
sdat=union(sdat,DNdat)

length(sdat)

sdata=dat[(row.names(dat) %in% sdat),]
head(sdata)
names(sdata)

data=sdata[,2:23]
row.names(data)=sdata$X
head(data)

#-------------Calulate principal coordinates
dd.veg=vegdist(t(data), "manhattan")
div.dd.veg=dd.veg/1000
head(div.dd.veg)

dd.pcoa=pcoa(div.dd.veg) 
head(dd.pcoa)
scores=dd.pcoa$vectors

#-------------First and second axes
plot(scores[,1], scores[,2],col=as.numeric(conditions$lunar), pch=16)
ordihull(scores,lunar,label=T)

plot(scores[,1], scores[,2],col=as.numeric(conditions$dayTime), pch=16)
ordihull(scores,dayTime,label=T)

plot(scores[,1], scores[,2],col=as.numeric(conditions$indiv), pch=16)
ordihull(scores,indiv,label=T)
```
#-------------Second and third axes 
plot(scores[,2], scores[,3],col=as.numeric(conditions$lunar),pch=16)
ordihull(scores[,2:3],lunar,label=T)

plot(scores[,2], scores[,3],col=as.numeric(conditions$dayTime),pch=16)
ordihull(scores[,2:3],dayTime,label=T)

plot(scores[,2], scores[,3],col=as.numeric(conditions$indiv),pch=16)
ordihull(scores[,2:3],indiv,label=T)

#-------------PERMANOVA
adonis(t(data)~lunar+indiv+dayTime,data=conditions,method="manhattan")
#adonis(t(data)~lunar,data=conditions,method="manhattan")  

pco1=scores[,1]
TukeyHSD(aov(pco1~lunar))
TukeyHSD(aov(pco1~dayTime))
TukeyHSD(aov(pco1~indiv))