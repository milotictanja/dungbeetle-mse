rm(list=ls())
library(vegan)

setwd("~/PhD/8 MSE Europe/Github/dungbeetle-mse/data")

####-NMDS analysis####
db_fg<-read.delim("dbFG_envSites.txt",header=T)
summary(db_fg)
db_fg$total<-db_fg$BT+db_fg$D+db_fg$SR+db_fg$ST

#remove empty samples
db_fg1<-subset(db_fg,db_fg$total>0)

#proportions of functional groups
db_fg1$pBT<-db_fg1$BT/db_fg1$total
db_fg1$pD<-db_fg1$D/db_fg1$total
db_fg1$pSR<-db_fg1$SR/db_fg1$total
db_fg1$pST<-db_fg1$ST/db_fg1$total

#different datasets for environmental factors and diversity
env<-db_fg1[,1:18]
div<-db_fg1[,19:22]
pdiv<-db_fg1[,24:27]

region<-env$biogeo
summary(region)
climate<-env$climate
summary(climate)
site<-env$SiteCode
summary(site)

#traptype T4= +/- T5 en T1 ongeveer T1b
env$trap[env$trap=="T4"]<-"T5"
env$trap[env$trap=="T1b"]<-"T1"

#correlation environmental variables
library(dplyr)

env<-group_by(env,SiteCode,start,year)
env_sum<-summarise(env,count=n(),lat=mean(lat),lon=mean(lon),alt=mean(alt),MMXT=mean(MMXT), MMNT=mean(MMNT),MNTM=mean(MNTM),EMXT=mean(EMXT),EMNT=mean(EMNT),TotalP=mean(TotalP),EMXP=mean(EMXP))

env_num<-env_sum[,5:14]
env_num<-env_num[complete.cases(env_num),]
cor(env_num)
cor.test(env_num$MMXT,env_num$MMNT)
cor.test(env_num$MMXT,env_num$TotalP)
cor.test(env_num$EMXP,env_num$TotalP)
cor.test(env_num$MNTM,env_num$TotalP)
cor.test(env_num$MNTM,env_num$MMXT)
cor.test(env_num$MNTM,env_num$MMNT)

env2<-group_by(env,SiteCode)
env2_sum<-summarise(env2,count=n(),lat=mean(lat),lon=mean(lon),alt=mean(alt))
env2_num<-env2_sum[,3:5]
cor(env2_num)

cor.test(env2_num$lat,env2_num$lon)
cor.test(env2_num$lat,env2_num$alt)
cor.test(env2_num$lon,env2_num$alt)

#mean temp highly correlated with temp vars, extreme min and max temp correlated with min and max temp
#x, y, z correlated, but yet needed in analysis
#keep lat, lon, alt, MNTM and TotalP


#NMDS with absolute values
NMDS1<-metaMDS(div,k=2,trymax=2000)#Bray-Curtis distance by default

#shepard plot
stressplot(NMDS1)
plot(NMDS1)

envfit1<-envfit(NMDS1~dung+trap+lat+lon+alt+MNTM+TotalP, env,na.rm=T)
envfit1

plot(NMDS1, type = "n")
ordiellipse(NMDS1, region, col="skyblue")
orditorp(NMDS1,display="species",col="red",air=0.01)
orditorp(NMDS1,display="sites",cex=1.25,air=0.01)
plot(envfit1,p.max=0.05)

#NMDS with proportions of functional groups
NMDS2<-metaMDS(pdiv,k=2,trymax=2000)#Bray-Curtis distance by default

#shepard plot
stressplot(NMDS2)
plot(NMDS2)

envfit2<-envfit(NMDS2~dung+trap+ lat+lon+alt+MMXT+MMNT+TotalP+EMXP, env,na.rm=T)
envfit2

plot(NMDS2, type = "n")
ordiellipse(NMDS2, region, col="skyblue")
orditorp(NMDS2,display="species",col="red",air=0.01)
#orditorp(NMDS2,display="sites",cex=1.25,air=0.01)
plot(envfit2,p.max=0.05)
#absolute values much better

#make a nice graph

library(ggplot2)
library(grid)

#data for plotting 
##NMDS points
NMDS.data<-env
NMDS.data$NMDS1<-NMDS1$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe
NMDS.data$NMDS2<-NMDS1$points[ ,2] 

##species data
spps <- data.frame(scores(NMDS1, display = "species")) #dataframe of species scores for plotting
spps$species <- row.names(spps) # making a column with species names
spps$species<-c("lT","D","sR","sT") #more intuitive abbreviations
spps<-spps[!is.na(spps$NMDS1) & !is.na(spps$NMDS2),] #removes NAs
spps$labelX<-spps$NMDS1
spps$labelY<-spps$NMDS2+0.05

# data for the envfit arrows
env.scores <- as.data.frame(scores(envfit1, display = "vectors")) #extracts relevant scores from envfit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) #and then gives them their names
#enkel significante variabelen plotten
envfit1

env.scores<-subset(env.scores,((env.variables!="alt")))
env.scores
env.scores$var<-c("y","x","mT","tP")

mult <- 1.5#multiplier for the arrows and text for envfit below
#plaats voor labels
env.scores$labelX<-ifelse(env.scores$NMDS1>0,(env.scores$NMDS1)*mult,(env.scores$NMDS1-0.075)*mult)
env.scores$labelY<-ifelse(env.scores$NMDS2>0,(env.scores$NMDS2)*mult,(env.scores$NMDS2-0.075)*mult)
#manueel aanpassen
env.scores$labelY<-ifelse(env.scores$var=="y",env.scores$labelY+0.1,(ifelse(env.scores$var=="x",env.scores$labelY+0.025,(ifelse(env.scores$var=="mT",env.scores$labelY-0.05,env.scores$labelY)))))

# function for ellipsess - just run this, is used later
#taken from the excellent stackoverflow Q+A: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#N_LY site has only 3 dung beetles -> same value for 3 samples -> impossible to plot ellipse
NMDS.data<-subset(NMDS.data,SiteCode!="N_LY")
#remove level Taiga
NMDS.data$biogeo<-factor(NMDS.data$biogeo)

#data for ellipse, in this case using region
df_ell.region <- data.frame() #sets up a data frame before running the function.
for(g in levels(NMDS.data$biogeo)){
  df_ell.region <- rbind(df_ell.region, cbind(as.data.frame(with(NMDS.data [NMDS.data$biogeo==g,],
                                                                 veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))),
                                              biogeo=g))
}

# data for labelling the ellipse
NMDS.mean<-aggregate(NMDS.data[,c("NMDS1", "NMDS2")], list(group = NMDS.data$biogeo), mean)

## finally plotting. 
NMDSplot <- ggplot()+
  theme_classic()+
  coord_cartesian(xlim = c(-1.5,2), ylim=c(-1.5,1))+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=24),axis.text.y=element_text(size=24))+
  theme(axis.title.x=element_text(size=28),axis.title.y=element_text(size=28))+
  ylab(expression(bold("NMDS2")))+
  xlab(expression(bold("NMDS1")))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))+
  geom_hline(aes(yintercept=0),lty=2)+
  geom_vline(aes(xintercept=0),lty=2)+
  geom_point(data=spps,aes(x=NMDS1,y=NMDS2),color="grey30",shape=19,size=14)+
  geom_path(data = df_ell.region, aes(x = NMDS1, y = NMDS2, group = biogeo, col=biogeo,lty=biogeo),size=1)+ #this is the ellipse, seperate ones by Site. If you didn't change the "alpha" (the shade) then you need to keep the "group 
  geom_segment(data = env.scores,
               aes(x = 0, xend = mult*NMDS1, y = 0, yend = mult*NMDS2),
               arrow = arrow(length = unit(0.25, "cm"),type="closed"),size=1, colour = "black") + #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
  geom_text(data = env.scores, #labels the environmental variable arrows * "mult" as for the arrows
            aes(x = labelX, y = labelY, label=var),size = 6,fontface=2,colour="black",vjust = 0,hjust=0)+
  geom_text(data=spps,aes(x=NMDS1,y=NMDS2,label=species),size = 6,colour="white",fontface=4,hjust = 0.5,vjust=0.5)
NMDSplot

NMDSplot<-NMDSplot+
  scale_color_manual(name=expression(bold("Biogeographic region")),
                     values=c("lightgoldenrod","thistle4","coral3","forestgreen","green3","goldenrod2","olivedrab3","aquamarine2","chocolate4","deeppink"))+
  scale_linetype_manual(name=expression(bold("Biogeographic region")),
                        values=c("solid","longdash","solid","longdash","solid","longdash","dotdash","solid","solid","longdash"))+
  theme(legend.key.width = unit(4, 'lines'),legend.key.size = unit(2.5, 'lines'),legend.title=element_text(size=28,face="bold"),legend.text=element_text(size=24))
NMDSplot


tiff(filename="NMDSplotBiogeoFG_v1.tiff",width=1200,height=800,unit="px")
NMDSplot
dev.off()



####ANOVA analysis of functional group composition and species richness####
rm(list=ls())
library(car)
library(agricolae)
source("HelperFunctions_summary.txt")
library(ggplot2)

#abundance of each functional group
db_fg<-read.delim("dbFG_envSites.txt",header=T)
summary(db_fg)
db_fg$total<-db_fg$BT+db_fg$D+db_fg$SR+db_fg$ST

#remove empty samples
db_fg<-subset(db_fg,db_fg$total>0)

#proportions of functional groups
db_fg$pBT<-db_fg$BT/db_fg$total
db_fg$pD<-db_fg$D/db_fg$total
db_fg$pSR<-db_fg$SR/db_fg$total
db_fg$pST<-db_fg$ST/db_fg$total

#traptype T4= +/- T5 en T1 ongeveer T1b
db_fg$trap[db_fg$trap=="T4"]<-"T5"
db_fg$trap[db_fg$trap=="T1b"]<-"T1"

logitTransform <- function(p) { log(p/(1-p)) }

#D
qqPlot(db_fg$D)
qqPlot(log(db_fg$D+1))

qqPlot(db_fg$pD)
qqPlot(log(db_fg$pD+1))

qqPlot(logitTransform(db_fg$pD+0.01))
hist(logitTransform(db_fg$pD+0.01))

it<-db_fg$biogeo

tapply(logitTransform(db_fg$pD+0.01),it,shapiro.test)

leveneTest(logitTransform(db_fg$pD+0.01),it)

mean<-tapply(logitTransform(db_fg$pD+0.01),it,mean,na.rm=T)
var<-tapply(logitTransform(db_fg$pD+0.01),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

aovD<-aov(logitTransform(db_fg$pD+0.01)~dung+trap+(biogeo/Code),data=db_fg)
anova(aovD)

D<-HSD.test(aovD,"biogeo",group=T)
D

D.sum<-summarySE(data=db_fg,measurevar = "D",groupvars = c("biogeo"),na.rm=T)
D.sum$FG<-rep("D",11)
D.groups<-D$groups[order(D$groups$trt),]
D.groups
D.groups<-rbind(D.groups,c("West Eurasian taiga",0,"NA"))
D.sum<-cbind(D.sum,D.groups)
D.sum

#BT
qqPlot(db_fg$BT)
qqPlot(log(db_fg$BT+1))

qqPlot(db_fg$pBT)
qqPlot(log(db_fg$pBT+1))

qqPlot(logitTransform(db_fg$pBT+0.01))
hist(logitTransform(db_fg$pBT+0.01))

it<-db_fg$biogeo

tapply(logitTransform(db_fg$pBT+0.01),it,shapiro.test)

leveneTest(logitTransform(db_fg$pBT+0.01),it)

mean<-tapply(logitTransform(db_fg$pBT+0.01),it,mean,na.rm=T)
var<-tapply(logitTransform(db_fg$pBT+0.01),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

aovBT<-aov(logitTransform(db_fg$pBT+0.01)~dung+trap+(biogeo/Code),data=db_fg)
anova(aovBT)

BT<-HSD.test(aovBT,"biogeo",group=T)
BT

BT.sum<-summarySE(data=db_fg,measurevar = "BT",groupvars = c("biogeo"),na.rm=T)
BT.sum$FG<-rep("BT",11)
BT.groups<-BT$groups[order(BT$groups$trt),]
BT.groups
BT.sum<-cbind(BT.sum,BT.groups)
BT.sum


#ST
qqPlot(db_fg$ST)
qqPlot(log(db_fg$ST+1))

qqPlot(db_fg$pST)
qqPlot(log(db_fg$pST+1))

qqPlot(logitTransform(db_fg$pST+0.01))
hist(logitTransform(db_fg$pST+0.01))

it<-db_fg$biogeo

tapply(logitTransform(db_fg$pST+0.01),it,shapiro.test)

leveneTest(logitTransform(db_fg$pST+0.01),it)

mean<-tapply(logitTransform(db_fg$pST+0.01),it,mean,na.rm=T)
var<-tapply(logitTransform(db_fg$pST+0.01),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

aovST<-aov(logitTransform(db_fg$pST+0.01)~dung+trap+(biogeo/Code),data=db_fg)
anova(aovST)

ST<-HSD.test(aovST,"biogeo",group=T)
ST

ST.sum<-summarySE(data=db_fg,measurevar = "ST",groupvars = c("biogeo"),na.rm=T)
ST.sum$FG<-rep("ST",11)
ST.groups<-ST$groups[order(ST$groups$trt),]
ST.groups
ST.sum<-cbind(ST.sum,ST.groups)
ST.sum


#SR
qqPlot(db_fg$SR)
qqPlot(log(db_fg$SR+1))

qqPlot(db_fg$pSR)
qqPlot(log(db_fg$pSR+1))

qqPlot(logitTransform(db_fg$pSR+0.01))
hist(logitTransform(db_fg$pSR+0.01))

it<-db_fg$biogeo

tapply(logitTransform(db_fg$pSR+0.01),it,shapiro.test)

leveneTest(logitTransform(db_fg$pSR+0.01),it)

mean<-tapply(logitTransform(db_fg$pSR+0.01),it,mean,na.rm=T)
var<-tapply(logitTransform(db_fg$pSR+0.01),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

aovSR<-aov(logitTransform(db_fg$pSR+0.01)~dung+trap+(biogeo/Code),data=db_fg)
anova(aovSR)

SR<-HSD.test(aovSR,"biogeo",group=T)
SR

SR.sum<-summarySE(data=db_fg,measurevar = "SR",groupvars = c("biogeo"),na.rm=T)
SR.sum$FG<-rep("SR",11)
SR.groups<-SR$groups[order(SR$groups$trt),]
SR.groups
SR.sum<-cbind(SR.sum,SR.groups)
SR.sum


#total
qqPlot(db_fg$total)
qqPlot(log(db_fg$total+1))

it<-db_fg$biogeo

tapply(log(db_fg$total+1),it,shapiro.test)

leveneTest(log(db_fg$total+1),it)

mean<-tapply(log(db_fg$total+1),it,mean,na.rm=T)
var<-tapply(log(db_fg$total+1),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

aovtotal<-aov(log(db_fg$total+1)~dung+trap+(biogeo/Code),data=db_fg)
anova(aovtotal)

total<-HSD.test(aovtotal,"biogeo",group=T)
total


total.sum<-summarySE(data=db_fg,measurevar = "total",groupvars = c("biogeo"),na.rm=T)
total.sum$FG<-rep("total",11)
total.groups<-total$groups[order(total$groups$trt),]
total.groups
total.sum<-cbind(total.sum,total.groups)
total.sum

#export of combined results D, BT, ST, SR, total
colnames(D.sum)[3]<-"mean"
colnames(BT.sum)[3]<-"mean"
colnames(ST.sum)[3]<-"mean"
colnames(SR.sum)[3]<-"mean"
colnames(total.sum)[3]<-"mean"

fg<-rbind(D.sum,BT.sum, ST.sum,SR.sum,total.sum)
fg$result<-paste(round(fg$mean,digits=1),"±",round(fg$se,digits=1))
write.csv(fg,file="fg_exportANOVA.csv")

####
#species richness of each functional group
db_fgR<-read.delim("dbFGrichness_envSites.txt",header=T)
summary(db_fgR)
db_fgR$total<-db_fgR$BT+db_fgR$D+db_fgR$SR+db_fgR$ST

#remove empty samples
db_fgR<-subset(db_fgR,db_fgR$total>0)

#traptype T4= +/- T5 en T1 ongeveer T1b
db_fgR$trap[db_fgR$trap=="T4"]<-"T5"
db_fgR$trap[db_fgR$trap=="T1b"]<-"T1"


#D
qqPlot(db_fgR$D)
qqPlot(log(db_fgR$D+1))

it<-db_fgR$biogeo

tapply(log(db_fgR$D+1),it,shapiro.test)

leveneTest(log(db_fgR$D+1),it)

mean<-tapply(log(db_fgR$D+1),it,mean,na.rm=T)
var<-tapply(log(db_fgR$D+1),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

aovDr<-aov(log(db_fgR$D+1)~dung+trap+(biogeo/Code),data=db_fgR)
anova(aovDr)

Dr<-HSD.test(aovDr,"biogeo",group=T)
Dr

Dr.sum<-summarySE(data=db_fgR,measurevar = "D",groupvars = c("biogeo"),na.rm=T)
Dr.sum$FG<-rep("D",11)
Dr.groups<-Dr$groups[order(Dr$groups$trt),]
Dr.groups
Dr.sum<-cbind(Dr.sum,Dr.groups)
Dr.sum

#BT
qqPlot(db_fgR$BT)
qqPlot(log(db_fgR$BT+1))

it<-db_fgR$biogeo

tapply(log(db_fgR$BT+1),it,shapiro.test)

leveneTest(log(db_fgR$BT+1),it)

mean<-tapply(log(db_fgR$BT+1),it,mean,na.rm=T)
var<-tapply(log(db_fgR$BT+1),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

aovBTr<-aov(log(db_fgR$BT+1)~dung+trap+(biogeo/Code),data=db_fgR)
anova(aovBTr)

BTr<-HSD.test(aovBTr,"biogeo",group=T)
BTr

BTr.sum<-summarySE(data=db_fgR,measurevar = "BT",groupvars = c("biogeo"),na.rm=T)
BTr.sum$FG<-rep("BT",11)
BTr.groups<-BTr$groups[order(BTr$groups$trt),]
BTr.groups
BTr.sum<-cbind(BTr.sum,BTr.groups)
BTr.sum

#ST
qqPlot(db_fgR$ST)
qqPlot(log(db_fgR$ST+1))

it<-db_fgR$biogeo

tapply(log(db_fgR$ST+1),it,shapiro.test)

leveneTest(log(db_fgR$ST+1),it)

mean<-tapply(log(db_fgR$ST+1),it,mean,na.rm=T)
var<-tapply(log(db_fgR$ST+1),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

aovSTr<-aov(log(db_fgR$ST+1)~dung+trap+(biogeo/Code),data=db_fgR)
anova(aovSTr)

STr<-HSD.test(aovSTr,"biogeo",group=T)
STr

STr.sum<-summarySE(data=db_fgR,measurevar = "ST",groupvars = c("biogeo"),na.rm=T)
STr.sum$FG<-rep("ST",11)
STr.groups<-STr$groups[order(STr$groups$trt),]
STr.groups
STr.sum<-cbind(STr.sum,STr.groups)
STr.sum

#SR
qqPlot(db_fgR$SR)
qqPlot(log(db_fgR$SR+1))

it<-db_fgR$biogeo

tapply(log(db_fgR$SR+1),it,shapiro.test)

leveneTest(log(db_fgR$SR+1),it)

mean<-tapply(log(db_fgR$SR+1),it,mean,na.rm=T)
var<-tapply(log(db_fgR$SR+1),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

aovSRr<-aov(log(db_fgR$SR+1)~dung+trap+(biogeo/Code),data=db_fgR)
anova(aovSRr)

SRr<-HSD.test(aovSRr,"biogeo",group=T)
SRr

SRr.sum<-summarySE(data=db_fgR,measurevar = "SR",groupvars = c("biogeo"),na.rm=T)
SRr.sum$FG<-rep("SR",11)
SRr.groups<-SRr$groups[order(SRr$groups$trt),]
SRr.groups
SRr.sum<-cbind(SRr.sum,SRr.groups)
SRr.sum

#SR
qqPlot(db_fgR$SR)
qqPlot(log(db_fgR$SR+1))

it<-db_fgR$biogeo

tapply(log(db_fgR$SR+1),it,shapiro.test)

leveneTest(log(db_fgR$SR+1),it)

mean<-tapply(log(db_fgR$SR+1),it,mean,na.rm=T)
var<-tapply(log(db_fgR$SR+1),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

aovSRr<-aov(log(db_fgR$SR+1)~dung+trap+(biogeo/Code),data=db_fgR)
anova(aovSRr)

SRr<-HSD.test(aovSRr,"biogeo",group=T)
SRr

SRr.sum<-summarySE(data=db_fgR,measurevar = "SR",groupvars = c("biogeo"),na.rm=T)
SRr.sum$FG<-rep("SR",11)
SRr.groups<-SRr$groups[order(SRr$groups$trt),]
SRr.groups
SRr.sum<-cbind(SRr.sum,SRr.groups)
SRr.sum


#total
qqPlot(db_fgR$total)
qqPlot(log(db_fgR$total+1))

it<-db_fgR$biogeo

tapply(log(db_fgR$total+1),it,shapiro.test)

leveneTest(log(db_fgR$total+1),it)

mean<-tapply(log(db_fgR$total+1),it,mean,na.rm=T)
var<-tapply(log(db_fgR$total+1),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

aovtotalr<-aov(log(db_fgR$total+1)~dung+trap+(biogeo/Code),data=db_fgR)
anova(aovtotalr)

totalr<-HSD.test(aovtotalr,"biogeo",group=T)
totalr

totalr.sum<-summarySE(data=db_fgR,measurevar = "total",groupvars = c("biogeo"),na.rm=T)
totalr.sum$FG<-rep("total",11)
totalr.groups<-totalr$groups[order(totalr$groups$trt),]
totalr.groups
totalr.sum<-cbind(totalr.sum,totalr.groups)
totalr.sum

#export of combined results D, BT, ST, SR, total
colnames(Dr.sum)[3]<-"mean"
colnames(BTr.sum)[3]<-"mean"
colnames(STr.sum)[3]<-"mean"
colnames(SRr.sum)[3]<-"mean"
colnames(totalr.sum)[3]<-"mean"

fgr<-rbind(Dr.sum,BTr.sum, STr.sum,SRr.sum,totalr.sum)
fgr$result<-paste(round(fgr$mean,digits=1),"±",round(fgr$se,digits=1))
write.csv(fgr,file="fgrichness_exportANOVA.csv")

#anova results
response<-c("dung","trap","biogeo","biogeo:Code","Residuals")
colnames<-c("response","Df","Fvalue","p","FG")
anovaD1<-anova(aovD)
FG<-rep("D",5)
anovaD<-data.frame(response,anovaD1$Df,anovaD1$`F value`,anovaD1$`Pr(>F)`,FG)
colnames(anovaD)<-colnames

anovaBT1<-anova(aovBT)
FG<-rep("BT",5)
anovaBT<-data.frame(response,anovaBT1$Df,anovaBT1$`F value`,anovaBT1$`Pr(>F)`,FG)
colnames(anovaBT)<-colnames

anovaST1<-anova(aovST)
FG<-rep("ST",5)
anovaST<-data.frame(response,anovaST1$Df,anovaST1$`F value`,anovaST1$`Pr(>F)`,FG)
colnames(anovaST)<-colnames

anovaSR1<-anova(aovSR)
FG<-rep("SR",5)
anovaSR<-data.frame(response,anovaSR1$Df,anovaSR1$`F value`,anovaSR1$`Pr(>F)`,FG)
colnames(anovaSR)<-colnames

anovatotal1<-anova(aovtotal)
FG<-rep("total",5)
anovatotal<-data.frame(response,anovatotal1$Df,anovatotal1$`F value`,anovatotal1$`Pr(>F)`,FG)
colnames(anovatotal)<-colnames

anovaDr1<-anova(aovDr)
FG<-rep("Dr",5)
anovaDr<-data.frame(response,anovaDr1$Df,anovaDr1$`F value`,anovaDr1$`Pr(>F)`,FG)
colnames(anovaDr)<-colnames

anovaBTr1<-anova(aovBTr)
FG<-rep("BTr",5)
anovaBTr<-data.frame(response,anovaBTr1$Df,anovaBTr1$`F value`,anovaBTr1$`Pr(>F)`,FG)
colnames(anovaBTr)<-colnames

anovaSTr1<-anova(aovSTr)
FG<-rep("STr",5)
anovaSTr<-data.frame(response,anovaSTr1$Df,anovaSTr1$`F value`,anovaSTr1$`Pr(>F)`,FG)
colnames(anovaSTr)<-colnames

anovaSRr1<-anova(aovSRr)
FG<-rep("SRr",5)
anovaSRr<-data.frame(response,anovaSRr1$Df,anovaSRr1$`F value`,anovaSRr1$`Pr(>F)`,FG)
colnames(anovaSRr)<-colnames

anovatotalr1<-anova(aovtotalr)
FG<-rep("totalr",5)
anovatotalr<-data.frame(response,anovatotalr1$Df,anovatotalr1$`F value`,anovatotalr1$`Pr(>F)`,FG)
colnames(anovatotalr)<-colnames

anova.results<-rbind(anovaD,anovaBT,anovaST,anovaSR, anovatotal,anovaDr,anovaBTr,anovaSTr,anovaSRr,anovatotalr)

anova.results

anova.results<-subset(anova.results,response!="Residuals")
anova.results$Fvalue<-round(anova.results$Fvalue,digits=3)
anova.results$p<-round(anova.results$p,digits=3)

write.csv(anova.results,file="anova_results.csv")
