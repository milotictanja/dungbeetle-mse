rm(list=ls())
library(vegan)


####NMDS analysis####

setwd("~/PhD/8 MSE Europe/Github/dungbeetle-mse/data")
dbspec<-read.delim("DBassemblageSUMAbs.txt",header=T)
head(dbspec)

dbenv<-read.delim("DBenvSUM.txt",header=T)
head(dbenv)
summary(dbenv)
#variabelen die niet gebruikt worden verwijderen?

#traptype T4= +/- T5 en T1 ongeveer T1b
dbenv$trap[dbenv$trap=="T4"]<-"T5"
dbenv$trap[dbenv$trap=="T1b"]<-"T1"

NMDS1<-metaMDS(dbspec,k=2,trymax=2000)#Bray-Curtis distance by default

#shepard plot
stressplot(NMDS1)
plot(NMDS1)


ordiplot(NMDS1,type="n")
orditorp(NMDS1,display="species",col="red",air=0.01)
orditorp(NMDS1,display="sites",cex=1.25,air=0.01)

region<-dbenv$biogeo.prov
summary(region)
climate<-dbenv$climate
summary(climate)

ordiplot(NMDS1,type="n")
ordihull(NMDS1,groups=region,draw="polygon",col="grey90",label=F)
orditorp(NMDS1,display="species",col="red",air=0.01)
orditorp(NMDS1,display="sites",air=0.01,cex=1.25)

ordiplot(NMDS1,type="n")
ordihull(NMDS1,groups=climate,draw="polygon",col="grey90",label=F)
orditorp(NMDS1,display="species",col="red",air=0.01)
orditorp(NMDS1,display="sites",air=0.01,cex=1.25)


#linear relationship between scores and variables

envfit<-envfit(NMDS1, dbenv,nperm=200,na.rm=T)
envfit

#formule gebruiken gelijkaardig aan CCA probeersel
envfit2<-envfit(NMDS1~dung+trap+ y+x+z+mmaxT+mminT+mT+emaxT+eminT+tP+
                  eP, dbenv,na.rm=T)
envfit2

plot(NMDS1, type = "n")
ordispider(NMDS1, region, col="skyblue")
points(NMDS1, display = "sites", col = as.numeric(region), pch=16)
orditorp(NMDS1,display="species",col="red",air=0.01)
plot(envfit2,p.max=0.05)

plot(NMDS1, type = "n")
ordispider(NMDS1, climate, col="skyblue")
points(NMDS1, display = "sites", col = as.numeric(climate), pch=16)
orditorp(NMDS1,display="species",col="red",air=0.01)
plot(envfit2,p.max=0.05)


#♣dung type G en D maar 1 keer gebruikt
levels(dbenv$dung)<-c("C","D","G","H","S","O")
dbenv$dung[dbenv$dung=="G"]<-"O"
dbenv$dung[dbenv$dung=="D"]<-"O"
envfit3<-envfit(NMDS1~dung+trap+ y+x+z+mmaxT+mminT+mT+emaxT+eminT+tP+
                  eP, dbenv,na.rm=T)
envfit3

plot(NMDS1, type = "n")
ordispider(NMDS1, region, col="skyblue")
points(NMDS1, display = "sites", col = as.numeric(region), pch=16)
orditorp(NMDS1,display="species",col="red",air=0.5)
plot(envfit3,p.max=0.05)

plot(NMDS1, type = "n")
ordiellipse(NMDS1, region, col="skyblue")
points(NMDS1, display = "sites", col = as.numeric(region), pch=16)
orditorp(NMDS1,display="species",col="red",air=0.5)
plot(envfit3,p.max=0.05)

plot(NMDS1, type = "n")
ordiellipse(NMDS1, climate, col="skyblue")
points(NMDS1, display = "sites", col = as.numeric(climate), pch=16)
orditorp(NMDS1,display="species",col="red",air=0.5)
plot(envfit3,p.max=0.05)

#deftige grafiek, volgens https://oliviarata.wordpress.com/2014/04/17/ordinations-in-ggplot2/

#plot1: by biogeographic zone

library(ggplot2)
library(grid)

#data for plotting 
##NMDS points
NMDS.data<-dbenv #there are other ways of doing this. But this is the way I do it for ease of plotting
NMDS.data$NMDS1<-NMDS1$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
NMDS.data$NMDS2<-NMDS1$points[ ,2] 

##species data
stems<-colSums(dbspec) #total abundances for each species
spps <- data.frame(scores(NMDS1, display = "species")) #dataframe of species scoes for plotting
spps$species <- row.names(spps) # making a column with species names
spps$colsums <- stems #adding the colSums from above
spps<-spps[!is.na(spps$NMDS1) & !is.na(spps$NMDS2),] #removes NAs
spps.colmedian <- median(spps$colsums) #create an object that is the median of the abundance of the measured species
spps.colmean <- mean(spps$colsums) #creates a mean instead if you wish to use
spps.NMDS1mean <- mean(spps$NMDS1)
spps.NMDS2mean <- mean(spps$NMDS2)
spps.NMDS1q95<-quantile(spps$NMDS1,0.95)
spps.NMDS1q05<-quantile(spps$NMDS1,0.05)
spps.NMDS2q95<-quantile(spps$NMDS2,0.95)
spps.NMDS2q05<-quantile(spps$NMDS2,0.05)
spps.colq95<-quantile(spps$colsums,0.95)
spps.colq05<-quantile(spps$colsums,0.05)

spps2 <- subset(spps,spps$colsums > spps.colmean) #select the most abundant species. Could discard fewer by going something like - spps$colsums>(spps.colmedian/2) instead
spps2$species <- factor(spps2$species) #otherwise factor doesn't drop unused levels and it will throw an error
#nog te veel op grafiek

#species tekenen indien NMDS1/2 in 0.05 of 0.95 quantiel zitten of abundantie in 0.95 quantiel
spps3<-subset(spps,((spps$NMDS1>spps.NMDS1q95)|(spps$NMDS2>spps.NMDS2q95)|
                      (spps$NMDS1<spps.NMDS1q05)|(spps$NMDS2<spps.NMDS2q05)|
                      (spps$colsums>spps.colq95)))


#functionele groepen toevoegen
abbr<-read.delim("abbr_list.txt",header=T)
spps<-cbind(abbr,spps)

# data for the envfit arrows
env.scores <- as.data.frame(scores(envfit3, display = "vectors")) #extracts relevant scores from envifit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) #and then gives them their names
#enkel significante variabelen plotten
envfit3

env.scores<-subset(env.scores,((env.variables!="eP")&(env.variables!="mmaxT")&(env.variables!="emaxT")))
env.scores

mult <- 2 #multiplier for the arrows and text for envfit below. You can change this and then rerun the plot command.
#plaats voor labels
rangeNMDS1<-max(NMDS.data$NMDS1)-min(NMDS.data$NMDS1)
rangeNMDS2<-max(NMDS.data$NMDS2)-min(NMDS.data$NMDS2)
env.scores$labelX<-ifelse(env.scores$NMDS1>0,(env.scores$NMDS1*mult)+0.1,(env.scores$NMDS1*mult)-0.25)
env.scores$labelY<-ifelse(env.scores$NMDS2>0,(env.scores$NMDS2*mult)+(rangeNMDS2/60),(env.scores$NMDS2*mult))
env.scores$labelX<-ifelse(env.scores$env.variables=="mminT",env.scores$labelX-0.75,env.scores$labelX)
env.scores$labelY<-ifelse(env.scores$env.variables=="mminT",env.scores$labelY-0.1,env.scores$labelY)
env.scores$labelX<-ifelse(env.scores$env.variables=="mT",env.scores$labelX-0.5,env.scores$labelX)
env.scores$labelY<-ifelse(env.scores$env.variables=="mT",env.scores$labelY+0.05,env.scores$labelY)
env.scores$labelX<-ifelse(env.scores$env.variables=="tP",env.scores$labelX-0.05,env.scores$labelX)
env.scores$labelY<-ifelse(env.scores$env.variables=="tP",env.scores$labelY+0.05,env.scores$labelY)

# function for ellipsess - just run this, is used later
#taken from the excellent stackoverflow Q+A: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using region
df_ell.region <- data.frame() #sets up a data frame before running the function.
for(g in levels(NMDS.data$biogeo.prov)){
  df_ell.region <- rbind(df_ell.region, cbind(as.data.frame(with(NMDS.data [NMDS.data$biogeo.prov==g,],
                                                                 veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))),
                                              biogeo.prov=g))
}

# data for labelling the ellipse
NMDS.mean=aggregate(NMDS.data[,c("NMDS1", "NMDS2")], 
                    list(group = NMDS.data$biogeo.prov), mean)

## finally plotting. 
NMDSplot <- ggplot()+
  geom_polygon(data = df_ell.region, aes(x = NMDS1, y = NMDS2, fill = biogeo.prov),alpha=.33,show.legend=F)+
  geom_path(data = df_ell.region, aes(x = NMDS1, y = NMDS2, group = biogeo.prov, col=biogeo.prov),size=1)+ #this is the ellipse, seperate ones by Site. If you didn't change the "alpha" (the shade) then you need to keep the "group 
  geom_point(data=NMDS.data,aes(x=NMDS1,y=NMDS2,col=biogeo.prov),show.legend=F)+
  geom_point(data=spps,aes(x=NMDS1,y=NMDS2,pch=FunctionalGroup,size=FunctionalGroup))+
  #annotate("text",x = NMDS.mean$NMDS1,y = NMDS.mean$NMDS2,label=NMDS.mean$group) + #labels for the centroids - I haven't used this since we have a legend. but you could also dithc the legend, but plot will get v messy
  #meest abundante soorten labellen
  #annotate("text",x=spps3$NMDS1,y=spps3$NMDS2,label=spps3$species)+ #te warrig
  geom_segment(data = env.scores,
               aes(x = 0, xend = mult*NMDS1, y = 0, yend = mult*NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),size=1, colour = "red4") + #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
  geom_text(data = env.scores, #labels the environmental variable arrows * "mult" as for the arrows
            aes(x = labelX, y = labelY, label=env.variables),
            size = 8,colour="red4",
            fontface=4,
            hjust = 0)
NMDSplot

NMDSplot<-NMDSplot+theme_classic()+
  #Vcoord_cartesian(xlim = c(-5,5), ylim=c(-2,2))+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=24),axis.text.y=element_text(size=24))+
  theme(axis.title.x=element_text(size=28),axis.title.y=element_text(size=28))+
  ylab(expression(bold("NMDS2")))+
  xlab(expression(bold("NMDS1")))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))+
  geom_hline(aes(yintercept=0),lty=2)+
  geom_vline(aes(xintercept=0),lty=2)
NMDSplot

NMDSplot<-NMDSplot+
  scale_color_manual(name=expression(bold("Biogeographic province")),
                     values=c("lightgoldenrod","thistle4","coral3","forestgreen","green3","goldenrod2","olivedrab3","aquamarine2","chocolate4","deeppink","steelblue4"))+
  scale_shape_manual(name=expression(bold("Functional group")),
                     breaks=c("D","ST","BT","SR"),
                     labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                     values=c(17,3,19,24))+
  scale_size_manual(name=expression(bold("Functional group")),
                     breaks=c("D","ST","BT","SR"),
                     labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                     values=c(6,2,8,4))+
  scale_fill_manual(name=expression(bold("Biogeographic province")),
                    values=c("lightgoldenrod","thistle4","coral3","forestgreen","green3","goldenrod2","olivedrab3","aquamarine2","chocolate4","deeppink","steelblue4"))+
  theme(legend.key.width = unit(2, 'lines'),legend.key.size = unit(2.5, 'lines'),legend.title=element_text(size=28,face="bold"),legend.text=element_text(size=24))
NMDSplot

#waarde stress toevoegen
NMDS1

NMDSplot<-NMDSplot+
  annotate(geom = "text",x=-4.5,y=-1.8,label="Stress: 0.108",size=10)
NMDSplot
NMDSplotbiogeo<-NMDSplot

tiff(filename="NMDSplotBiogeo_v2.tiff",width=1200,height=800,unit="px")
NMDSplotbiogeo
dev.off()

NMDSplotbiogeo_NL<-NMDSplotbiogeo+ theme(legend.position="none")

NMDSplot1 <- ggplot()+
  geom_polygon(data = df_ell.region, aes(x = NMDS1, y = NMDS2, fill = biogeo.prov),alpha=.33,show.legend=F)+
  geom_path(data = df_ell.region, aes(x = NMDS1, y = NMDS2, group = biogeo.prov, col=biogeo.prov),size=1)+ #this is the ellipse, seperate ones by Site. If you didn't change the "alpha" (the shade) then you need to keep the "group 
  geom_point(data=NMDS.data,aes(x=NMDS1,y=NMDS2,col=biogeo.prov),show.legend=F)+
  geom_point(data=spps,aes(x=NMDS1,y=NMDS2,pch=FunctionalGroup,size=FunctionalGroup),show.legend=F)+
  #annotate("text",x = NMDS.mean$NMDS1,y = NMDS.mean$NMDS2,label=NMDS.mean$group) + #labels for the centroids - I haven't used this since we have a legend. but you could also dithc the legend, but plot will get v messy
  #meest abundante soorten labellen
  #annotate("text",x=spps3$NMDS1,y=spps3$NMDS2,label=spps3$species)+ #te warrig
  geom_segment(data = env.scores,
               aes(x = 0, xend = mult*NMDS1, y = 0, yend = mult*NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),size=1, colour = "red4") + #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
  geom_text(data = env.scores, #labels the environmental variable arrows * "mult" as for the arrows
            aes(x = labelX, y = labelY, label=env.variables),
            size = 8,colour="red4",
            fontface=4,
            hjust = 0)
NMDSplot1

NMDSplot1<-NMDSplot1+theme_classic()+
  #Vcoord_cartesian(xlim = c(-5,5), ylim=c(-2,2))+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=24),axis.text.y=element_text(size=24))+
  theme(axis.title.x=element_text(size=28),axis.title.y=element_text(size=28))+
  ylab(expression(bold("NMDS2")))+
  xlab(expression(bold("NMDS1")))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))+
  geom_hline(aes(yintercept=0),lty=2)+
  geom_vline(aes(xintercept=0),lty=2)
NMDSplot1

NMDSplot1<-NMDSplot1+
  scale_color_manual(name=expression(bold("Biogeographic province")),
                     values=c("lightgoldenrod","thistle4","coral3","forestgreen","green3","goldenrod2","olivedrab3","aquamarine2","chocolate4","deeppink","steelblue4"))+
  scale_shape_manual(name=expression(bold("Functional group")),
                     breaks=c("D","ST","BT","SR"),
                     labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                     values=c(17,18,19,24))+
  scale_size_manual(name=expression(bold("Functional group")),
                    breaks=c("D","ST","BT","SR"),
                    labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                    values=c(6,6,8,4))+
  scale_fill_manual(name=expression(bold("Biogeographic province")),
                    values=c("lightgoldenrod","thistle4","coral3","forestgreen","green3","goldenrod2","olivedrab3","aquamarine2","chocolate4","deeppink","steelblue4"))+
  theme(legend.key.width = unit(2, 'lines'),legend.key.size = unit(2.5, 'lines'),legend.title=element_text(size=28,face="bold"),legend.text=element_text(size=24))
NMDSplot1

#waarde stress toevoegen
NMDS1

#NMDSplot1<-NMDSplot1+
#  annotate(geom = "text",x=-4.5,y=-1.8,label="Stress: 0.108",size=10)
NMDSplot1
NMDSplotbiogeo_NL<-NMDSplot1+ theme(legend.position="none")

NMDSplot1b <- ggplot()+
  geom_polygon(data = df_ell.region, aes(x = NMDS1, y = NMDS2, fill = biogeo.prov))+
  geom_path(data = df_ell.region, aes(x = NMDS1, y = NMDS2, group = biogeo.prov, col=biogeo.prov),size=1,show.legend=F)+ #this is the ellipse, seperate ones by Site. If you didn't change the "alpha" (the shade) then you need to keep the "group 
  geom_point(data=NMDS.data,aes(x=NMDS1,y=NMDS2,col=biogeo.prov),show.legend=F)+
  geom_point(data=spps,aes(x=NMDS1,y=NMDS2,pch=FunctionalGroup,size=FunctionalGroup),show.legend=F)+
  #annotate("text",x = NMDS.mean$NMDS1,y = NMDS.mean$NMDS2,label=NMDS.mean$group) + #labels for the centroids - I haven't used this since we have a legend. but you could also dithc the legend, but plot will get v messy
  #meest abundante soorten labellen
  #annotate("text",x=spps3$NMDS1,y=spps3$NMDS2,label=spps3$species)+ #te warrig
  geom_segment(data = env.scores,
               aes(x = 0, xend = mult*NMDS1, y = 0, yend = mult*NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),size=1, colour = "red4") + #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
  geom_text(data = env.scores, #labels the environmental variable arrows * "mult" as for the arrows
            aes(x = labelX, y = labelY, label=env.variables),
            size = 8,colour="red4",
            fontface=4,
            hjust = 0)
NMDSplot1b
NMDSplot1b<-NMDSplot1b+theme_classic()+
  #Vcoord_cartesian(xlim = c(-5,5), ylim=c(-2,2))+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=24),axis.text.y=element_text(size=24))+
  theme(axis.title.x=element_text(size=28),axis.title.y=element_text(size=28))+
  ylab(expression(bold("NMDS2")))+
  xlab(expression(bold("NMDS1")))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))+
  geom_hline(aes(yintercept=0),lty=2)+
  geom_vline(aes(xintercept=0),lty=2)
NMDSplot1b

NMDSplot1b<-NMDSplot1b+
  scale_color_manual(name=expression(bold("Biogeographic province")),
                     values=c("lightgoldenrod","thistle4","coral3","forestgreen","green3","goldenrod2","olivedrab3","aquamarine2","chocolate4","deeppink","steelblue4"))+
  scale_shape_manual(name=expression(bold("Functional group")),
                     breaks=c("D","ST","BT","SR"),
                     labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                     values=c(17,18,19,24))+
  scale_size_manual(name=expression(bold("Functional group")),
                    breaks=c("D","ST","BT","SR"),
                    labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                    values=c(6,6,8,4))+
  scale_fill_manual(name=expression(bold("Biogeographic province")),
                    values=c("lightgoldenrod","thistle4","coral3","forestgreen","green3","goldenrod2","olivedrab3","aquamarine2","chocolate4","deeppink","steelblue4"))+
  theme(legend.key.width = unit(2, 'lines'),legend.key.size = unit(2.5, 'lines'),legend.title=element_text(size=28,face="bold"),legend.text=element_text(size=24))
NMDSplot1b

library(gtable)
# extract the different legends
leg <- NMDSplot1b + guides(group=FALSE)
legend <- gtable_filter(ggplot_gtable(ggplot_build(leg)), "guide-box") 
legGrob <- grobTree(legend)

# final plot were you can place the legends were you want them
legend.plot <- qplot(1:10, 1:10, geom = "blank") + theme_classic()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())+
  annotation_custom(grob = legGrob, xmin = 0, xmax = 5, ymin = 1, ymax = 10)
legend.plot



#plot2: by climate zone

#data for ellipse, in this case using climate
df_ell.climate <- data.frame() #sets up a data frame before running the function.
for(g in levels(NMDS.data$climate)){
  df_ell.climate <- rbind(df_ell.climate, cbind(as.data.frame(with(NMDS.data [NMDS.data$climate==g,],
                                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                                ,climate=g))
}

# data for labelling the ellipse
NMDS.mean=aggregate(NMDS.data[ ,c("NMDS1", "NMDS2")], 
                    list(group = NMDS.data$climate), mean)

## finally plotting. 
NMDSplot2 <- ggplot()+
  geom_polygon(data = df_ell.climate, aes(x = NMDS1, y = NMDS2, fill = climate),alpha=.33,show.legend=F)+
  geom_path(data = df_ell.climate, aes(x = NMDS1, y = NMDS2, group = climate, col=climate),size=1)+ #this is the ellipse, seperate ones by Site. If you didn't change the "alpha" (the shade) then you need to keep the "group 
  geom_point(data=NMDS.data,aes(x=NMDS1,y=NMDS2,col=climate),show.legend=F)+
  geom_point(data=spps,aes(x=NMDS1,y=NMDS2,pch=FunctionalGroup,size=FunctionalGroup),show.legend=F)+
  geom_segment(data = env.scores,
               aes(x = 0, xend = mult*NMDS1, y = 0, yend = mult*NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),size=1, colour = "red4") + #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
  geom_text(data = env.scores, #labels the environmental variable arrows * "mult" as for the arrows
            aes(x = labelX, y = labelY, label=env.variables),
            size = 8,colour="red4",
            fontface=4,
            hjust = 0)
NMDSplot2

NMDSplot2<-NMDSplot2+theme_classic()+
  #coord_cartesian(xlim = c(-7.5,7.5), ylim=c(-9,9))+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=24),axis.text.y=element_text(size=24))+
  theme(axis.title.x=element_text(size=28),axis.title.y=element_text(size=28))+
  ylab(expression(bold("NMDS2")))+
  xlab(expression(bold("NMDS1")))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))+
  geom_hline(aes(yintercept=0),lty=2)+
  geom_vline(aes(xintercept=0),lty=2)
NMDSplot2

NMDSplot2<-NMDSplot2+
  scale_color_manual(name=expression(bold("Climate zone")),
                     values=c("#FF9999","#66FF00","forestgreen","goldenrod3","#33FFFF","#33CCCC"))+
  scale_shape_manual(name=expression(bold("Functional group")),
                     breaks=c("D","ST","BT","SR"),
                     labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                     values=c(17,18,19,24))+
  scale_size_manual(name=expression(bold("Functional group")),
                    breaks=c("D","ST","BT","SR"),
                    labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                    values=c(6,6,8,4))+
  scale_fill_manual(name=expression(bold("Climate zone")),
                    values=c("#FF9999","#66FF00","forestgreen","goldenrod3","#33FFFF","#33CCCC"))+
  theme(legend.key.width = unit(2, 'lines'),legend.key.size = unit(2.5, 'lines'),legend.title=element_text(size=28,face="bold"),legend.text=element_text(size=24))
NMDSplot2

#waarde stress toevoegen
NMDS1

#NMDSplot2<-NMDSplot2+
#  annotate(geom = "text",x=-4.5,y=-1.8,label="Stress: 0.108",size=10)
#NMDSplot2

tiff(filename="NMDSplotClimate_v2.tiff",width=1200,height=800,unit="px")
NMDSplot2
dev.off()

NMDSplot2_NL<-NMDSplot2+ theme(legend.position="none")

NMDSplot2b <- ggplot()+
  geom_polygon(data = df_ell.climate, aes(x = NMDS1, y = NMDS2, fill = climate))+
  geom_path(data = df_ell.climate, aes(x = NMDS1, y = NMDS2, group = climate, col=climate),size=1,show.legend=F)+ #this is the ellipse, seperate ones by Site. If you didn't change the "alpha" (the shade) then you need to keep the "group 
  geom_point(data=NMDS.data,aes(x=NMDS1,y=NMDS2,col=climate),show.legend=F)+
  geom_point(data=spps,aes(x=NMDS1,y=NMDS2,pch=FunctionalGroup,size=FunctionalGroup),show.legend=F)+
  geom_segment(data = env.scores,
               aes(x = 0, xend = mult*NMDS1, y = 0, yend = mult*NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),size=1, colour = "red4") + #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
  geom_text(data = env.scores, #labels the environmental variable arrows * "mult" as for the arrows
            aes(x = labelX, y = labelY, label=env.variables),
            size = 8,colour="red4",
            fontface=4,
            hjust = 0)
NMDSplot2b

NMDSplot2b<-NMDSplot2b+theme_classic()+
  #coord_cartesian(xlim = c(-7.5,7.5), ylim=c(-9,9))+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=24),axis.text.y=element_text(size=24))+
  theme(axis.title.x=element_text(size=28),axis.title.y=element_text(size=28))+
  ylab(expression(bold("NMDS2")))+
  xlab(expression(bold("NMDS1")))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))+
  geom_hline(aes(yintercept=0),lty=2)+
  geom_vline(aes(xintercept=0),lty=2)
NMDSplot2b

NMDSplot2b<-NMDSplot2b+
  scale_color_manual(name=expression(bold("Climate zone")),
                     values=c("#FF9999","#66FF00","forestgreen","goldenrod3","#33FFFF","#33CCCC"))+
  scale_shape_manual(name=expression(bold("Functional group")),
                     breaks=c("D","ST","BT","SR"),
                     labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                     values=c(17,18,19,24))+
  scale_size_manual(name=expression(bold("Functional group")),
                    breaks=c("D","ST","BT","SR"),
                    labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                    values=c(6,6,8,4))+
  scale_fill_manual(name=expression(bold("Climate zone")),
                    values=c("#FF9999","#66FF00","forestgreen","goldenrod3","#33FFFF","#33CCCC"))+
  theme(legend.key.width = unit(2, 'lines'),legend.key.size = unit(2.5, 'lines'),legend.title=element_text(size=28,face="bold"),legend.text=element_text(size=24))
NMDSplot2b

NMDSplot3 <- ggplot()+
  geom_polygon(data = df_ell.climate, aes(x = NMDS1, y = NMDS2, fill = climate),alpha=.33,show.legend=F)+
  geom_path(data = df_ell.climate, aes(x = NMDS1, y = NMDS2, group = climate, col=climate),size=1,show.legend=F)+ #this is the ellipse, seperate ones by Site. If you didn't change the "alpha" (the shade) then you need to keep the "group 
  geom_point(data=NMDS.data,aes(x=NMDS1,y=NMDS2,col=climate),show.legend=F)+
  geom_point(data=spps,aes(x=NMDS1,y=NMDS2,pch=FunctionalGroup,size=FunctionalGroup))
NMDSplot3

NMDSplot3<-NMDSplot3+theme_classic()+
    scale_shape_manual(name=expression(bold("Functional group")),
                     breaks=c("D","ST","BT","SR"),
                     labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                     values=c(17,18,19,24))+
  scale_size_manual(name=expression(bold("Functional group")),
                    breaks=c("D","ST","BT","SR"),
                    labels=c("dwellers","small tunnelers","large tunnelers","small rollers"),
                    values=c(6,6,8,4))+
  theme(legend.key.width = unit(2, 'lines'),legend.key.size = unit(2.5, 'lines'),legend.title=element_text(size=28,face="bold"),legend.text=element_text(size=24))
NMDSplot3


library(gtable)
# extract the different legends
leg2 <- NMDSplot2b + guides(group=FALSE)
legend2 <- gtable_filter(ggplot_gtable(ggplot_build(leg2)), "guide-box") 
legGrob2 <- grobTree(legend2)

# final plot were you can place the legends were you want them
legend.plot2 <- qplot(1:10, 1:10, geom = "blank") + theme_classic()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())+
  annotation_custom(grob = legGrob2, xmin = 0, xmax = 1, ymin = 7, ymax = 9)
legend.plot2

leg3 <- NMDSplot3 + guides(group=FALSE)
legend3 <- gtable_filter(ggplot_gtable(ggplot_build(leg3)), "guide-box") 
legGrob3 <- grobTree(legend3)

# final plot were you can place the legends were you want them
legend.plot3 <- qplot(1:10, 1:10, geom = "blank") + theme_classic()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())+
  annotation_custom(grob = legGrob3, xmin = 0, xmax = 3, ymin = 7, ymax = 9)
legend.plot3


#samenvoegen#
library(cowplot)
tiff(filename="NMDS_v3.tiff",width=1200,height=1400,unit="px")
ggdraw() +
  draw_plot(NMDSplotbiogeo_NL, 0, 0.5, 0.6, 0.5) +
  draw_plot(NMDSplot2_NL, 0, 0, 0.6, 0.5) +
  draw_plot(legend.plot, 0.7, 0.5, 0.3, 0.5) +
  draw_plot(legend.plot2, 0.7, 0.2, 0.3, 0.25) +
  draw_plot(legend.plot3, 0.7, 0, 0.3, 0.25) +
  draw_plot_label(c("a","b"),c(0,0),c(1,0.5),size=36)
dev.off()



####Functional group composition per biogeographic region and climate zone####
rm(list=ls())
library(car)
library(agricolae)
source("HelperFunctions_summary.txt")
library(ggplot2)

#aparte analyse per functionele groep
FGratio2<-read.delim("FGratio2.txt",header=T)
head(FGratio2)
summary(FGratio2)

qqPlot(FGratio2$D)
qqPlot(asin(sqrt(FGratio2$D)))
qqPlot((sqrt(FGratio2$D)))
qqPlot(log(FGratio2$D+1))

it<-FGratio2$Biogeographic.provinceUdvardy1975:FGratio2$Traptype2:FGratio2$Dungtype

tapply(FGratio2$D,it,shapiro.test)
tapply(asin(sqrt(FGratio2$D)),it,shapiro.test)
tapply(asin(FGratio2$D),it,shapiro.test)
tapply(sqrt(FGratio2$D),it,shapiro.test)
tapply(log(FGratio2$D+1),it,shapiro.test)

#

leveneTest(FGratio2$D,it)
leveneTest(asin(sqrt(FGratio2$D)),it)
leveneTest(asin(FGratio2$D),it)
leveneTest(sqrt(FGratio2$D),it)
leveneTest(log(FGratio2$D+1),it)


mean<-tapply(asin(sqrt(FGratio2$D)),it,mean,na.rm=T)
var<-tapply(asin(sqrt(FGratio2$D)),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)
#ok

qqPlot(FGratio2$BT)
qqPlot(asin(sqrt(FGratio2$BT)))
qqPlot((sqrt(FGratio2$BT)))
qqPlot(log(FGratio2$BT+1))

it<-FGratio2$Biogeographic.provinceUdvardy1975:FGratio2$Traptype2:FGratio2$Dungtype

tapply(FGratio2$BT,it,shapiro.test)
tapply(asin(sqrt(FGratio2$BT)),it,shapiro.test)
tapply(asin(FGratio2$BT),it,shapiro.test)
tapply(sqrt(FGratio2$BT),it,shapiro.test)
tapply(log(FGratio2$BT+1),it,shapiro.test)

#

leveneTest(FGratio2$BT,it)
leveneTest(asin(sqrt(FGratio2$BT)),it)
leveneTest(asin(FGratio2$BT),it)
leveneTest(sqrt(FGratio2$BT),it)
leveneTest(log(FGratio2$BT+1),it)


mean<-tapply(asin(sqrt(FGratio2$BT)),it,mean,na.rm=T)
var<-tapply(asin(sqrt(FGratio2$BT)),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)
#niet ok

mean<-tapply(sqrt(FGratio2$BT),it,mean,na.rm=T)
var<-tapply(sqrt(FGratio2$BT),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

mean<-tapply(asin(FGratio2$BT),it,mean,na.rm=T)
var<-tapply(asin(FGratio2$BT),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)

#geen verbetering


qqPlot(FGratio2$ST)
qqPlot(asin(sqrt(FGratio2$ST)))
qqPlot((sqrt(FGratio2$ST)))
qqPlot(log(FGratio2$ST+1))

it<-FGratio2$Biogeographic.provinceUdvardy1975:FGratio2$Traptype2:FGratio2$Dungtype

tapply(FGratio2$ST,it,shapiro.test)
tapply(asin(sqrt(FGratio2$ST)),it,shapiro.test)
tapply(asin(FGratio2$ST),it,shapiro.test)
tapply(sqrt(FGratio2$ST),it,shapiro.test)
tapply(log(FGratio2$ST+1),it,shapiro.test)

#

leveneTest(FGratio2$ST,it)
leveneTest(asin(sqrt(FGratio2$ST)),it)
leveneTest(asin(FGratio2$ST),it)
leveneTest(sqrt(FGratio2$ST),it)
leveneTest(log(FGratio2$ST+1),it)


mean<-tapply(asin(sqrt(FGratio2$ST)),it,mean,na.rm=T)
var<-tapply(asin(sqrt(FGratio2$ST)),it,var,na.rm=T)
plot(mean~var)
cor.test(mean,var)
#ok

FGratio2$climate<-FGratio2$Koppen.climate.zonePeel2007
FGratio2$bioprov<-FGratio2$Biogeographic.provinceUdvardy1975

#no need to transform abundance data (RM)

#dwellers
aovD1<-aov(asin(sqrt(D))~Dungtype+Traptype2+(climate/SiteCode)+
             (bioprov/SiteCode),data=FGratio2)
anova(aovD1)

D1<-HSD.test(aovD1,"Dungtype",group=T)
D1

D2<-HSD.test(aovD1,"Traptype2",group=T)
D2

D3<-HSD.test(aovD1,"bioprov",group=T)
D3

D4<-HSD.test(aovD1,"climate",group=T)
D4

#gemiddelden per variabele, SE en Tukey-groep
dung.D.sum<-summarySE(FGratio2,measurevar="D",groupvars=c("Dungtype"),na.rm=T)
dung.D.sum
dung.D.sum$FG<-rep("D",5)
dung.D.sum

D1.res<-D1$groups
D1.res<-D1.res[order(D1.res$trt),]
dung.D.sum<-cbind(dung.D.sum,D1.res)
dung.D.sum

trap.D.sum<-summarySE(FGratio2,measurevar="D",groupvars=c("Traptype2"),na.rm=T)
trap.D.sum
trap.D.sum$FG<-rep("D",2)
trap.D.sum

D2.res<-D2$groups
D2.res<-D2.res[order(D2.res$trt),]
trap.D.sum<-cbind(trap.D.sum,D2.res)
trap.D.sum
trap.D.sum$M<-c("","") #enkel letters toevoegen indien sign verschil
trap.D.sum

biogeo.D.sum<-summarySE(FGratio2,measurevar="D",groupvars=c("bioprov"),na.rm=T)
biogeo.D.sum
biogeo.D.sum$FG<-rep("D",11)
biogeo.D.sum

D3.res<-D3$groups
D3.res<-D3.res[order(D3.res$trt),]
biogeo.D.sum<-cbind(biogeo.D.sum,D3.res)
biogeo.D.sum$M<-c("ab","bc","c","bc","ad","d","ad","d","abd","d","c")
biogeo.D.sum

climate.D.sum<-summarySE(FGratio2,measurevar="D",groupvars=c("climate"),na.rm=T)
climate.D.sum
climate.D.sum$FG<-rep("D",6)
climate.D.sum

D4.res<-D4$groups
D4.res<-D4.res[order(D4.res$trt),]
climate.D.sum<-cbind(climate.D.sum,D4.res)
climate.D.sum


#♠large tunnelers
aovBT1<-aov(asin(sqrt(BT))~Dungtype+Traptype2+(climate/SiteCode)+
              (bioprov/SiteCode),data=FGratio2)
anova(aovBT1)

BT1<-HSD.test(aovBT1,"Dungtype",group=T)
BT1

BT2<-HSD.test(aovBT1,"Traptype2",group=T)
BT2

BT3<-HSD.test(aovBT1,"bioprov",group=T)
BT3

BT4<-HSD.test(aovBT1,"climate",group=T)
BT4

#gemiddelden per variabele, SE en Tukey-groep
dung.BT.sum<-summarySE(FGratio2,measurevar="BT",groupvars=c("Dungtype"),na.rm=T)
dung.BT.sum
dung.BT.sum$FG<-rep("BT",5)
dung.BT.sum

BT1.res<-BT1$groups
BT1.res<-BT1.res[order(BT1.res$trt),]
dung.BT.sum<-cbind(dung.BT.sum,BT1.res)
dung.BT.sum

trap.BT.sum<-summarySE(FGratio2,measurevar="BT",groupvars=c("Traptype2"),na.rm=T)
trap.BT.sum
trap.BT.sum$FG<-rep("BT",2)
trap.BT.sum

BT2.res<-BT2$groups
BT2.res<-BT2.res[order(BT2.res$trt),]
trap.BT.sum<-cbind(trap.BT.sum,BT2.res)
trap.BT.sum
trap.BT.sum$M<-c("","") 
trap.BT.sum

biogeo.BT.sum<-summarySE(FGratio2,measurevar="BT",groupvars=c("bioprov"),na.rm=T)
biogeo.BT.sum
biogeo.BT.sum$FG<-rep("BT",11)
biogeo.BT.sum

BT3.res<-BT3$groups
BT3.res<-BT3.res[order(BT3.res$trt),]
biogeo.BT.sum<-cbind(biogeo.BT.sum,BT3.res)
biogeo.BT.sum
biogeo.BT.sum$M<-c("a","ab","b","b","c","ab","ab","ab","ab","b","b")
biogeo.BT.sum

climate.BT.sum<-summarySE(FGratio2,measurevar="BT",groupvars=c("climate"),na.rm=T)
climate.BT.sum
climate.BT.sum$FG<-rep("BT",6)
climate.BT.sum

BT4.res<-BT4$groups
BT4.res<-BT4.res[order(BT4.res$trt),]
climate.BT.sum<-cbind(climate.BT.sum,BT4.res)
climate.BT.sum
climate.BT.sum$M<-c("a","b","a","ab","ab","c")
climate.BT.sum


#small tunnelers
aovST1<-aov(asin(sqrt(ST))~Dungtype+Traptype2+(climate/SiteCode)+
              (bioprov/SiteCode),data=FGratio2)
anova(aovST1)
ST1<-HSD.test(aovST1,"Dungtype",group=T)
ST1

ST2<-HSD.test(aovST1,"Traptype2",group=T)
ST2

ST3<-HSD.test(aovST1,"bioprov",group=T)
ST3

ST4<-HSD.test(aovST1,"climate",group=T)
ST4

#gemiddelden per variabele, SE en Tukey-groep
dung.ST.sum<-summarySE(FGratio2,measurevar="ST",groupvars=c("Dungtype"),na.rm=T)
dung.ST.sum
dung.ST.sum$FG<-rep("ST",5)
dung.ST.sum

ST1.res<-ST1$groups
ST1.res<-ST1.res[order(ST1.res$trt),]
dung.ST.sum<-cbind(dung.ST.sum,ST1.res)
dung.ST.sum
dung.ST.sum$M<-c("","","","","")
dung.ST.sum

trap.ST.sum<-summarySE(FGratio2,measurevar="ST",groupvars=c("Traptype2"),na.rm=T)
trap.ST.sum
trap.ST.sum$FG<-rep("ST",2)
trap.ST.sum

ST2.res<-ST2$groups
ST2.res<-ST2.res[order(ST2.res$trt),]
trap.ST.sum<-cbind(trap.ST.sum,ST2.res)
trap.ST.sum

biogeo.ST.sum<-summarySE(FGratio2,measurevar="ST",groupvars=c("bioprov"),na.rm=T)
biogeo.ST.sum
biogeo.ST.sum$FG<-rep("ST",11)
biogeo.ST.sum

ST3.res<-ST3$groups
ST3.res<-ST3.res[order(ST3.res$trt),]
biogeo.ST.sum<-cbind(biogeo.ST.sum,ST3.res)
biogeo.ST.sum
biogeo.ST.sum$M<-c("ab","ab","a","abc","a","c","c","c","bc","c","a")
biogeo.ST.sum

climate.ST.sum<-summarySE(FGratio2,measurevar="ST",groupvars=c("climate"),na.rm=T)
climate.ST.sum
climate.ST.sum$FG<-rep("ST",6)
climate.ST.sum

ST4.res<-ST4$groups
ST4.res<-ST4.res[order(ST4.res$trt),]
climate.ST.sum<-cbind(climate.ST.sum,ST4.res)
climate.ST.sum
climate.ST.sum$M<-c("a","a","a","b","b","a")
climate.ST.sum

#small rollers

aovSR1<-aov(asin(sqrt(SR))~Dungtype+Traptype2+(climate/SiteCode)+
              (bioprov/SiteCode),data=FGratio2)
anova(aovSR1)

SR1<-HSD.test(aovSR1,"Dungtype",group=T)
SR1

SR2<-HSD.test(aovSR1,"Traptype2",group=T)
SR2

SR3<-HSD.test(aovSR1,"bioprov",group=T)
SR3

SR4<-HSD.test(aovSR1,"climate",group=T,alpha=0.1)
SR4

#gemiddelden per variabele, SE en Tukey-groep
dung.SR.sum<-summarySE(FGratio2,measurevar="SR",groupvars=c("Dungtype"),na.rm=T)
dung.SR.sum
dung.SR.sum$FG<-rep("SR",5)
dung.SR.sum

SR1.res<-SR1$groups
SR1.res<-SR1.res[order(SR1.res$trt),]
dung.SR.sum<-cbind(dung.SR.sum,SR1.res)
dung.SR.sum
dung.SR.sum$M<-c("","","","","")
dung.SR.sum


trap.SR.sum<-summarySE(FGratio2,measurevar="SR",groupvars=c("Traptype2"),na.rm=T)
trap.SR.sum
trap.SR.sum$FG<-rep("SR",2)
trap.SR.sum

SR2.res<-SR2$groups
SR2.res<-SR2.res[order(SR2.res$trt),]
trap.SR.sum<-cbind(trap.SR.sum,SR2.res)
trap.SR.sum
trap.SR.sum$M<-c("","")
trap.SR.sum

biogeo.SR.sum<-summarySE(FGratio2,measurevar="SR",groupvars=c("bioprov"),na.rm=T)
biogeo.SR.sum
biogeo.SR.sum$FG<-rep("SR",11)
biogeo.SR.sum

SR3.res<-SR3$groups
SR3.res<-SR3.res[order(SR3.res$trt),]
biogeo.SR.sum<-cbind(biogeo.SR.sum,SR3.res)
biogeo.SR.sum
biogeo.SR.sum$M<-c("a","a","a","a","a","a","b","a","a","a","a")
biogeo.SR.sum

climate.SR.sum<-summarySE(FGratio2,measurevar="SR",groupvars=c("climate"),na.rm=T)
climate.SR.sum
climate.SR.sum$FG<-rep("SR",6)
climate.SR.sum

SR4.res<-SR4$groups
SR4.res<-SR4.res[order(SR4.res$trt),]
climate.SR.sum<-cbind(climate.SR.sum,SR4.res)
climate.SR.sum
climate.SR.sum$M<-c("a","b","b","c","b","b")
climate.SR.sum

#gecombineerde dataset maken + ruimte voor grote rollers (lege dataset)

#samenvoegen
colnames(dung.D.sum)[3]<-"avg"
colnames(dung.BT.sum)[3]<-"avg"
colnames(dung.ST.sum)[3]<-"avg"
colnames(dung.SR.sum)[3]<-"avg"

dung.sum<-rbind(dung.D.sum,dung.ST.sum,dung.BT.sum,dung.SR.sum)
dung.sum

colnames(trap.D.sum)[3]<-"avg"
colnames(trap.BT.sum)[3]<-"avg"
colnames(trap.ST.sum)[3]<-"avg"
colnames(trap.SR.sum)[3]<-"avg"

trap.sum<-rbind(trap.D.sum,trap.ST.sum,trap.BT.sum,trap.SR.sum)
trap.sum

colnames(biogeo.D.sum)[3]<-"avg"
colnames(biogeo.BT.sum)[3]<-"avg"
colnames(biogeo.ST.sum)[3]<-"avg"
colnames(biogeo.SR.sum)[3]<-"avg"

biogeo.sum<-rbind(biogeo.D.sum,biogeo.ST.sum,biogeo.BT.sum,biogeo.SR.sum)
biogeo.sum

colnames(climate.D.sum)[3]<-"avg"
colnames(climate.BT.sum)[3]<-"avg"
colnames(climate.ST.sum)[3]<-"avg"
colnames(climate.SR.sum)[3]<-"avg"

climate.sum<-rbind(climate.D.sum,climate.ST.sum,climate.BT.sum,climate.SR.sum)
climate.sum

#grafieken
dung.plot<-ggplot(data=dung.sum,aes(x=FG,y=avg*100,fill=Dungtype,label=M))+
  geom_bar(position="dodge",colour="black",stat="identity")+
  geom_errorbar(aes(ymin=avg*100,ymax=(avg+se)*100),colour="black",position=position_dodge(.9),width=.5)
dung.plot
dung.plot<-dung.plot+theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Relative abundance (%)")))+
  xlab(expression(bold(" ")))+
  scale_x_discrete(limits=c("D","ST","BT","SR"),
                   breaks=c("D","ST","BT","SR"),
                   labels=c("dwellers","small tunnelers","large tunnelers","small rollers"))+
  theme(axis.text.x=element_blank())+
  theme(axis.ticks=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(dung.sum)/4), y = -4, label = c("dwellers","small tunnelers","large tunnelers","small rollers","large rollers"), size = 5)+ 
  geom_text(aes(y=(100*(avg+se+(0.05*max(avg+se))))),position=position_dodge(width=1),size=4)
dung.plot
dung.plot<-dung.plot+
  scale_fill_manual(name=expression(bold("Dung type")),
                    breaks=c("C","D","G","H","S"),
                    labels=c("cattle","deer","goat","horse","sheep"),
                    values=c("black","grey30","grey60","grey90","white"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))+
  theme(legend.position=c(0.9,0.8))
dung.plot

tiff(filename="DungFGplot.tiff",width=600,height=400,unit="px")
dung.plot
dev.off()


trap.plot<-ggplot(data=trap.sum,aes(x=FG,y=avg*100,fill=Traptype2,label=M))+
  geom_bar(position="dodge",colour="black",stat="identity")+
  geom_errorbar(aes(ymin=avg*100,ymax=(avg+se)*100),colour="black",position=position_dodge(.9),width=.5)
trap.plot
trap.plot<-trap.plot+theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Relative abundance (%)")))+
  xlab(expression(bold(" ")))+
  scale_x_discrete(limits=c("D","ST","BT","SR"),
                   breaks=c("D","ST","BT","SR"),
                   labels=c("dwellers","small tunnelers","large tunnelers","small rollers"))+
  theme(axis.text.x=element_blank())+
  theme(axis.ticks=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(trap.sum)/2), y = -4, label = c("dwellers","small tunnelers","large tunnelers","small rollers"), size = 5)+ 
  geom_text(aes(y=(100*(avg+se+(0.05*max(avg+se))))),position=position_dodge(width=1),size=4)
trap.plot
trap.plot<-trap.plot+
  scale_fill_manual(name=expression(bold("Trap type")),
                    breaks=c("T1","T5"),
                    labels=c("T1","T5"),
                    values=c("black","white"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))+
  theme(legend.position=c(0.9,0.8))
trap.plot

tiff(filename="trapFGplot.tiff",width=600,height=400,unit="px")
trap.plot
dev.off()

summary(biogeo.sum)
biogeo.sum[is.na(biogeo.sum)] <- 0

biogeo.plot<-ggplot(data=biogeo.sum,aes(x=FG,y=avg*100,fill=trt,label=M))+
  geom_bar(position="dodge",colour="black",stat="identity")+
  geom_errorbar(aes(ymin=avg*100,ymax=(avg+se)*100),colour="black",position=position_dodge(.9),width=.5)
biogeo.plot
biogeo.plot<-biogeo.plot+theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Relative abundance (%)")))+
  xlab(expression(bold(" ")))+
  scale_x_discrete(limits=c("D","ST","BT","SR"),
                   breaks=c("D","ST","BT","SR"),
                   labels=c("dwellers","small tunnelers","large tunnelers","small rollers"))+
  theme(axis.text.x=element_blank())+
  theme(axis.ticks.x=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(biogeo.sum)/10), y = -4, label = c("dwellers","small tunnelers","large tunnelers","small rollers"), size = 5)+ 
  geom_text(aes(y=(100*(avg+se+(0.05*max(avg+se)))),angle=90),position=position_dodge(width=.9),vjust=.5,hjust=.5,size=4)
biogeo.plot
biogeo.plot<-biogeo.plot+
  scale_fill_manual(name=expression(bold("Biogeographical province")),
                    values=c("lightgoldenrod","thistle4","coral3","forestgreen","green3","goldenrod2","olivedrab3","aquamarine2","chocolate4","deeppink","steelblue4"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))+
  theme(legend.justification=c(0,1), legend.position=c(0.6,.9))
biogeo.plot

tiff(filename="biogeoFGplot.tiff",width=600,height=400,unit="px")
biogeo.plot
dev.off()

summary(climate.sum)
climate.sum[is.na(climate.sum)] <- 0

climate.plot<-ggplot(data=climate.sum,aes(x=FG,y=avg*100,fill=trt,label=M))+
  geom_bar(position="dodge",colour="black",stat="identity")+
  geom_errorbar(aes(ymin=avg*100,ymax=(avg+se)*100),colour="black",position=position_dodge(.9),width=.5)
climate.plot
climate.plot<-climate.plot+theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Relative abundance (%)")))+
  xlab(expression(bold(" ")))+
  scale_x_discrete(limits=c("D","ST","BT","SR"),
                   breaks=c("D","ST","BT","SR"),
                   labels=c("dwellers","small tunnelers","large tunnelers","small rollers"))+
  theme(axis.text.x=element_blank())+
  theme(axis.ticks.x=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(climate.sum)/6), y = -4, label = c("dwellers","small tunnelers","large tunnelers","small rollers"), size = 5)+ 
  geom_text(aes(y=(100*(avg+se+(0.05*max(avg+se)))),angle=90),position=position_dodge(width=.9),vjust=.5,hjust=.5,size=4)
climate.plot
climate.plot<-climate.plot+
  scale_fill_manual(name=expression(bold("Köppen-Geiger \nclimate zone")),
                    values=c("#FF9999","#66FF00","forestgreen","goldenrod3","#33FFFF","#33CCCC"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))+
  theme(legend.justification=c(0,1), legend.position=c(0.6,.9))
climate.plot

tiff(filename="climateFGplot.tiff",width=600,height=400,unit="px")
climate.plot
dev.off()


#samenvoegen#
library(cowplot)
tiff(filename="RelAbundFG.tiff",width=600,height=800,unit="px")
ggdraw() +
  draw_plot(biogeo.plot, 0, 0.5, 1, 0.5) +
  draw_plot(climate.plot, 0, 0, 1, 0.5) +
  draw_plot_label(c("a","b"),c(0,0),c(1,0.5),size=18)
dev.off()



#stacked graphs functional groups
summary(biogeo.sum)
levels(biogeo.sum$bioprov)
biogeo.sum$bioprov<-revalue(biogeo.sum$bioprov,
                            c("Atlantic"="atl","Boreonemoral"="bor",
                              "British islands"="brit","Caucaso-Iranian highlands" ="cau",
                              "Central European highlands"="ceh","Iberian highlands"="ibh",
                              "Mediterranean sclerophyl" ="med","Middle European forests"="mef",
                              "Pannonian"="pan","Pontian Steppe"="pon","West Eurasian taiga"="taig"))
biogeo.sum$FG<-as.factor(biogeo.sum$FG)

summary(biogeo.sum)
biogeo.sum$FG<-as.factor(biogeo.sum$FG)
levels(biogeo.sum$FG)
biogeo.sum$avg<-as.numeric(biogeo.sum$avg)
  
biogeo.sum$FG2<-ifelse(biogeo.sum$FG=="SR",1,(ifelse(biogeo.sum$FG=="ST",2,(ifelse(biogeo.sum$FG=="BT",3,4)))))
biogeo.sum$FG3<-ifelse(biogeo.sum$FG=="SR",4,(ifelse(biogeo.sum$FG=="ST",3,(ifelse(biogeo.sum$FG=="BT",2,1)))))
biogeo.sum$FG<-reorder(biogeo.sum$FG,biogeo.sum$FG3)

levels(biogeo.sum$bioprov)
biogeo.sum$bioprov<- ordered(biogeo.sum$bioprov, levels = c("atl" , "bor"  ,"brit" ,"cau" , "ceh" ,"ibh"  , "med" , "mef" , "pan","pon" , "taig"))

biogeo.sum.plot2<-ggplot(biogeo.sum[order(biogeo.sum$FG2),],aes(x=bioprov,y=avg*100,fill=FG))+
  geom_bar(position="stack",stat="identity")
biogeo.sum.plot2
biogeo.sum.plot2<-biogeo.sum.plot2+theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Relative abundance (%)")))+
  xlab(expression(bold("Biogeographic province")))+
  theme(axis.ticks.x=element_blank())+
  theme(axis.text.x=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(biogeo.sum)/4), y = -4, label = c("atl","bor","brit","cau","ceh","ibh","med","mef","pan","pont","taig"), size = 5)
biogeo.sum.plot2
biogeo.sum.plot2<-biogeo.sum.plot2+
  scale_fill_manual(name=expression(bold("Functional group")),
                    values=c("royalblue4","red4","orange","yellowgreen"),
                    labels=c("dwellers","tunnelers (large)","tunnelers (small)","rollers (small)"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))
biogeo.sum.plot2

#semi-manueel lettertjes zetten
biogeo.sum #via excel nieuwe dataset maken
biogeo.relFG.label<-read.delim("biogeo.relFG.labels.txt",header=T)
head(biogeo.relFG.label)
biogeo.sum.plot2<-biogeo.sum.plot2+
  geom_text(data=biogeo.relFG.label,aes(x=bioprov,y=labelY*100,color=color,label=M),size=5,show.legend = F)+
  scale_colour_manual(values=c("white","black"))
biogeo.sum.plot2

tiff(filename="relFGBiogeo_v2.tiff",width=600,height=400,unit="px")
biogeo.sum.plot2
dev.off()


#climate graph
summary(climate.sum)
climate.sum$FG<-as.factor(climate.sum$FG)

levels(climate.sum$FG)
climate.sum$FG2<-ifelse(climate.sum$FG=="SR",1,(ifelse(climate.sum$FG=="BR",2,(ifelse(climate.sum$FG=="ST",3,(ifelse(climate.sum$FG=="BT",4,5)))))))
climate.sum$FG3<-ifelse(climate.sum$FG=="SR",5,(ifelse(climate.sum$FG=="BR",4,(ifelse(climate.sum$FG=="ST",3,(ifelse(climate.sum$FG=="BT",2,1)))))))
climate.sum$FG<-reorder(climate.sum$FG,climate.sum$FG3)

climate.sum.plot2<-ggplot(climate.sum[order(climate.sum$FG2),],aes(x=climate,y=avg*100,fill=FG))+
  geom_bar(position="stack",stat="identity")
climate.sum.plot2
climate.sum.plot2<-climate.sum.plot2+theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Relative abundance (%)")))+
  xlab(expression(bold("Climate zone")))+
  theme(axis.ticks.x=element_blank())+
  theme(axis.text.x=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(climate.sum)/4), y = -4, label = c("BWk","Cfb","Cfc","Csb","Dfa","Dfb"), size = 5)
climate.sum.plot2
climate.sum.plot2<-climate.sum.plot2+
  scale_fill_manual(name=expression(bold("Functional group")),
                    values=c("royalblue4","red4","orange","yellowgreen","grey"),
                    labels=c("dwellers","tunnelers (large)","tunnelers (small)","rollers (small)","soil fauna"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))
climate.sum.plot2

#semi-manueel sterretjes zetten
climate.sum #via excel nieuwe dataset maken
climate.relFG.labels<-read.delim("climate.relFG.labels.txt",header=T)
head(climate.relFG.labels)
climate.sum.plot2<-climate.sum.plot2+
  geom_text(data=climate.relFG.labels,aes(x=climate,y=labelY*100,color=color,label=M),size=5,show.legend = F)+
  scale_colour_manual(values=c("white","black"))
climate.sum.plot2

tiff(filename="relFGClimate_v2.tiff",width=600,height=400,unit="px")
climate.sum.plot2
dev.off()

#wegschrijven, samengevoegd zonder legende bij biogeo?
biogeo.sum.plot2<-biogeo.sum.plot2+theme(legend.position="none")
biogeo.sum.plot2

library(cowplot)
tiff(filename="RelFG_v2.tiff",width=600,height=700,unit="px")
ggdraw() +
  draw_plot(biogeo.sum.plot2, 0, 0.5, 1, 0.5) +
  draw_plot(climate.sum.plot2, 0, 0, 1, 0.5) +
  draw_plot_label(c("a","b"),c(0,0),c(1,0.5),size=18)
dev.off()

#biogeo and climate in 1 graph

#dataset
biogeoclimate<-read.delim("biogeo.relFG.labels.txt",header=T)
head(biogeoclimate)
biogeoclimate$FG4<-paste(biogeoclimate$FG3,biogeoclimate$FG)

plot3<-ggplot(biogeoclimate[order(biogeoclimate$FG2),],aes(x=bioprov2,y=avg*100,fill=FG4))+
  geom_bar(position="stack",stat="identity")
plot3
plot3<-plot3+theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Relative abundance (%)")))+
  xlab(expression(bold("Biogeographic province and climate zone")))+
  theme(axis.ticks.x=element_blank())+
  theme(axis.text.x=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(biogeoclimate)/4), y = -3.5, label = c("taig","bor","ceh","pan","brit","atl","mef","pont","med","ibh","cau"), size = 5)+
  annotate(geom="text",x=c(1,2,5,7,9,11),y=-10,label=c("Cfc","Dfb","Cfb","Dfa","Csb","Bwk"),color="gray50",fontface="italic", size=5.5)+
  geom_segment(x=0.6,xend=1.4,y=-7,yend=-7,color="gray75",size=1)+
  geom_segment(x=1.6,xend=4.4,y=-7,yend=-7,color="gray75",size=1)+
  geom_segment(x=4.6,xend=6.4,y=-7,yend=-7,color="gray75",size=1)+
  geom_segment(x=6.6,xend=8.4,y=-7,yend=-7,color="gray75",size=1)+
  geom_segment(x=8.6,xend=10.4,y=-7,yend=-7,color="gray75",size=1)+
  geom_segment(x=10.6,xend=11.4,y=-7,yend=-7,color="gray75",size=1)
plot3
plot3<-plot3+
  scale_fill_manual(name=expression(bold("Functional group")),
                    values=c("royalblue4","red4","orange","yellowgreen"),
                    labels=c("dwellers","tunnelers (large)","tunnelers (small)","rollers (small)"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))
plot3

plot3<-plot3+
  geom_text(data=biogeoclimate,aes(x=bioprov2,y=labelY*100,color=color,label=M),size=5,show.legend = F)+
  scale_colour_manual(values=c("white","black"))
plot3

tiff(filename="relFGBiogeoclimate_v2.tiff",width=600,height=400,unit="px")
plot3
dev.off()




####Biodiversity indices####

rm(list=ls())
library(vegan)

setwd('D:/Tmilotic/Documents/Doctoraat/8 MSE Europe/R scripts')
dbspec_env<-read.delim("dbspec_envFull.txt",header=T)
summary(dbspec_env)

dbspec<-read.delim("dbspecFull2.txt",header=T)
summary(dbspec)

dbenv<-read.delim("dbenvFull2.txt",header=T)
summary(dbenv)
dbenv$climate<-dbenv$Koppen.climate.zonePeel2007
dbenv$biogeo.prov<-dbenv$Biogeographic.provinceUdvardy1975

#biodiversity indices
H <- diversity(dbspec)
simp <- diversity(dbspec, "simpson")
invsimp <- diversity(dbspec, "inv")
## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
unbias.simp <- rarefy(dbspec, 2) - 1
## Fisher alpha
alpha <- fisher.alpha(dbspec)
## Plot all
pairs(cbind(H, simp, invsimp, unbias.simp, alpha), pch="+", col="blue")
## Species richness (S) and Pielou's evenness (J):

S <- specnumber(dbspec) ## rowSums(BCI > 0) does the same...
J <- H/log(S)
## beta diversity defined as gamma/alpha - 1:

#shannon wiener diversity (H), species richness (S) and eveness (J) per climate zone, biogeo and site
specClimate<-rowsum(dbspec,dbenv$climate)
Hclimate <- diversity(specClimate)
Sclimate<-specnumber(specClimate)
Jclimate<-Hclimate/log(Sclimate)

specbiogeo<-rowsum(dbspec,dbenv$biogeo.prov)
Hbiogeo <- diversity(specbiogeo)
Sbiogeo<-specnumber(specbiogeo)
Jbiogeo<-Hbiogeo/log(Sbiogeo)

specSite<-rowsum(dbspec,dbenv$SiteCode)
HSite <- diversity(specSite)
SSite<-specnumber(specSite)
JSite<-HSite/log(SSite)



#per site, period and dung type
head(dbenv)
dbenv$group<-paste(dbenv$SiteCode,dbenv$period,dbenv$Dungtype)
dbenv$group<-factor(dbenv$group)
levels(dbenv$group)

alpha <- with(dbenv, tapply(specnumber(dbspec), group, mean))
#test SE berekenen
alphaSD <- with(dbenv, tapply(specnumber(dbspec), group, sd))
alphaN <- with(dbenv, tapply(specnumber(dbspec), group, length))
alphaSE<-alphaSD/sqrt(alphaN)

gamma <- with(dbenv, specnumber(dbspec, group))
gamma/alpha - 1


#specpool Extrapolated Species Richness in a Species Pool
pool <- specpool(dbspec, dbenv$group)
pool

pool.site <- specpool(dbspec, dbenv$SiteCode)
pool.site

pool.region <- specpool(dbspec, dbenv$biogeo.prov)
pool.region

pool.climate <- specpool(dbspec, dbenv$climate)
pool.climate

#sampling effort, specaccum, by site
summary(dbenv)
dbspecSite<-cbind(dbspec,dbenv$SiteCode)
summary(dbspecSite)
colnames(dbspecSite)[95]<-"site"
levels(dbspecSite$site)

dbB_KH<-subset(dbspecSite,site=="B_KH")
dbB_ZW<-subset(dbspecSite,site=="B_ZW")
dbD_BAV<-subset(dbspecSite,site=="D_BAV")
dbD_BR<-subset(dbspecSite,site=="D_BR")
dbD_SB<-subset(dbspecSite,site=="D_SB")
dbEST_TP<-subset(dbspecSite,site=="EST_TP")
dbF_CH<-subset(dbspecSite,site=="F_CH")
dbF_MP<-subset(dbspecSite,site=="F_MP")
dbH_BUG<-subset(dbspecSite,site=="H_BUG")
dbH_VAC<-subset(dbspecSite,site=="H_VAC")
dbIR_SR<-subset(dbspecSite,site=="IR_SR")
dbIR_TS<-subset(dbspecSite,site=="IR_TS")
dbN_LY<-subset(dbspecSite,site=="N_LY")
dbSP_CV<-subset(dbspecSite,site=="SP_CV")
dbUK_MH<-subset(dbspecSite,site=="UK_MH")
dbUK_WW<-subset(dbspecSite,site=="UK_WW")
dbRO_BI<-subset(dbspecSite,site=="RO_BI")

dbB_KH<-dbB_KH[,1:94]
dbB_ZW<-dbB_ZW[,1:94]
dbD_BAV<-dbD_BAV[,1:94]
dbD_BR<-dbD_BR[,1:94]
dbD_SB<-dbD_SB[,1:94]
dbEST_TP<-dbEST_TP[,1:94]
dbF_CH<-dbF_CH[,1:94]
dbF_MP<-dbF_MP[,1:94]
dbH_VAC<-dbH_VAC[,1:94]
dbIR_SR<-dbIR_SR[,1:94]
dbIR_TS<-dbIR_TS[,1:94]
dbN_LY<-dbN_LY[,1:94]
dbSP_CV<-dbSP_CV[,1:94]
dbH_BUG<-dbH_BUG[,1:94]
dbUK_MH<-dbUK_MH[,1:94]
dbUK_WW<-dbUK_WW[,1:94]
dbRO_BI<-dbRO_BI[,1:94]

#by region
dbAtl<-rbind(dbB_KH,dbB_ZW,dbF_CH,dbD_SB)
dbcont<-rbind(dbD_BR,dbD_BAV)
dbBrit<-rbind(dbUK_WW,dbUK_MH)
dbIr<-rbind(dbIR_TS,dbIR_SR)

#atlantic
atl1.sp1 <- specaccum(dbB_KH)
atl1.sp2 <- specaccum(dbB_KH, "random")
atl1.sp2
summary(atl1.sp2)

atl1.sp2$sites
atl1.sp2$method
atl1.sp2$call
atl1.sp2$richness
atl1.sp2$sd
atl1.sp2$perm

plot(atl1.sp2$sites,atl1.sp2$richness)


atl2.sp1 <- specaccum(dbB_ZW)
atl2.sp2 <- specaccum(dbB_ZW, "random")
atl2.sp2
summary(atl2.sp2)

plot(atl2.sp2$sites,atl2.sp2$richness)

atl3.sp1 <- specaccum(dbF_CH)
atl3.sp2 <- specaccum(dbF_CH, "random")
atl3.sp2
summary(atl3.sp2)

plot(atl3.sp2$sites,atl3.sp2$richness)

atl4.sp1 <- specaccum(dbD_SB)
atl4.sp2 <- specaccum(dbD_SB, "random")
atl4.sp2
summary(atl4.sp2)

plot(atl4.sp2$sites,atl4.sp2$richness)

atl.sp1 <- specaccum(dbAtl)
atl.sp2 <- specaccum(dbAtl, "random")
atl.sp2
summary(atl.sp2)

plot(atl.sp2$sites,atl.sp2$richness)



#central european highlands
cont1.sp1 <- specaccum(dbD_BAV)
cont1.sp2 <- specaccum(dbD_BAV, "random")
cont1.sp2
summary(cont1.sp2)

plot(cont1.sp2$sites,cont1.sp2$richness)


cont2.sp1 <- specaccum(dbD_BR)
cont2.sp2 <- specaccum(dbD_BR, "random")
cont2.sp2
summary(cont2.sp2)

plot(cont2.sp2$sites,cont2.sp2$richness)

cont.sp1 <- specaccum(dbcont)
cont.sp2 <- specaccum(dbcont, "random")
cont.sp2
summary(cont.sp2)

plot(cont.sp2$sites,cont.sp2$richness)


#pannonian
pan1.sp1 <- specaccum(dbH_BUG)
pan1.sp2 <- specaccum(dbH_BUG, "random")
pan1.sp2
summary(pan1.sp2)

plot(pan1.sp2$sites,pan1.sp2$richness)


#British islands
brit1.sp1 <- specaccum(dbUK_WW)
brit1.sp2 <- specaccum(dbUK_WW, "random")
brit1.sp2
summary(brit1.sp2)

plot(brit1.sp2$sites,brit1.sp2$richness)

brit2.sp1 <- specaccum(dbUK_MH)
brit2.sp2 <- specaccum(dbUK_MH, "random")
brit2.sp2
summary(brit2.sp2)

plot(brit2.sp2$sites,brit2.sp2$richness)

brit.sp1 <- specaccum(dbBrit)
brit.sp2 <- specaccum(dbBrit, "random")
brit.sp2
summary(brit.sp2)

plot(brit.sp2$sites,brit.sp2$richness)


#Caucaso-Iranian highlands
ir1.sp1 <- specaccum(dbIR_TS)
ir1.sp2 <- specaccum(dbIR_TS, "random")
ir1.sp2
summary(ir1.sp2)

plot(ir1.sp2$sites,ir1.sp2$richness)

ir2.sp1 <- specaccum(dbIR_SR)
ir2.sp2 <- specaccum(dbIR_SR, "random")
ir2.sp2
summary(ir2.sp2)

plot(ir2.sp2$sites,ir2.sp2$richness)

ir.sp1 <- specaccum(dbIr)
ir.sp2 <- specaccum(dbIr, "random")
ir.sp2
summary(ir.sp2)

plot(ir.sp2$sites,ir.sp2$richness)


#Boreonemoral
est1.sp1 <- specaccum(dbEST_TP)
est1.sp2 <- specaccum(dbEST_TP, "random")
est1.sp2
summary(est1.sp2)

plot(est1.sp2$sites,est1.sp2$richness)


#iberian highlands
ib1.sp1 <- specaccum(dbSP_CV)
ib1.sp2 <- specaccum(dbSP_CV, "random")
ib1.sp2
summary(ib1.sp2)

plot(ib1.sp2$sites,ib1.sp2$richness)


#mediterranean sclerophyl
med1.sp1 <- specaccum(dbF_MP)
med1.sp2 <- specaccum(dbF_MP, "random")
med1.sp2
summary(med1.sp2)

plot(med1.sp2$sites,med1.sp2$richness)


#middle european forest
mef1.sp1 <- specaccum(dbH_VAC)
mef1.sp2 <- specaccum(dbH_VAC, "random")
mef1.sp2
summary(mef1.sp2)

plot(mef1.sp2$sites,mef1.sp2$richness)

#pontian steppe
pon1.sp1 <- specaccum(dbRO_BI)
pon1.sp2 <- specaccum(dbRO_BI, "random")
pon1.sp2
summary(pon1.sp2)

plot(pon1.sp2$sites,pon1.sp2$richness)


#west eurasian taiga
wet1.sp1 <- specaccum(dbN_LY)
wet1.sp2 <- specaccum(dbN_LY, "random")
wet1.sp2
summary(wet1.sp2)

plot(wet1.sp2$sites,wet1.sp2$richness)


#grafieken maken per regio
#atlantic
atl1<-cbind(atl1.sp2$sites,atl1.sp2$richness,atl1.sp2$sd)
colnames(atl1)<-c("samples","richness","sd")
atl1<-data.frame(atl1)
dim(atl1)
atl1$site<-c(rep("Kalmthout",83))
atl2<-cbind(atl2.sp2$sites,atl2.sp2$richness,atl2.sp2$sd)
colnames(atl2)<-c("samples","richness","sd")
atl2<-data.frame(atl2)
dim(atl2)
atl2$site<-c(rep("Zwin",66))
atl3<-cbind(atl3.sp2$sites,atl3.sp2$richness,atl3.sp2$sd)
colnames(atl3)<-c("samples","richness","sd")
atl3<-data.frame(atl3)
dim(atl3)
atl3$site<-c(rep("Le Chesnoy",49))

atl4<-cbind(atl4.sp2$sites,atl4.sp2$richness,atl4.sp2$sd)
colnames(atl4)<-c("samples","richness","sd")
atl4<-data.frame(atl4)
dim(atl4)
atl4$site<-c(rep("Steinbuhl",22))

atl<-cbind(atl.sp2$sites,atl.sp2$richness,atl.sp2$sd)
colnames(atl)<-c("samples","richness","sd")
atl<-data.frame(atl)
dim(atl)
atl$site<-c(rep("all Atlantic sites",220))

atl<-rbind(atl1,atl2,atl3,atl4,atl)
atl$lower<-atl$richness-atl$sd
atl$upper<-atl$richness+atl$sd

library(ggplot2)
levels(atl$site)
atl$site<-ordered(atl$site,levels=c("all Atlantic sites","Kalmthout","Le Chesnoy","Steinbuhl","Zwin"))

atl.plot<-ggplot()+
  geom_ribbon(data=atl,aes(x=samples,ymin=lower,ymax=upper,fill=site),alpha=.3)+
  geom_line(data=atl,aes(x=samples,y=richness,colour=site,lty=site),size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Species richness")))+
  xlab(expression(bold(" ")))+
  ggtitle(expression(bold("Atlantic")))+
  theme(plot.title=element_text(size=14,hjust=0))+
  scale_fill_discrete(name="Study site")+
  scale_linetype_discrete(name="Study site")+
  scale_color_discrete(name="Study site")+
  theme(legend.title=element_text(size=12,face="bold"))+
  theme(legend.text=element_text(size=10))+
  theme(legend.position=c(0.75,0.3))+
  coord_cartesian(ylim=c(0,45))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
atl.plot


#Central European highlands
cont1<-cbind(cont1.sp2$sites,cont1.sp2$richness,cont1.sp2$sd)
colnames(cont1)<-c("samples","richness","sd")
cont1<-data.frame(cont1)
dim(cont1)
cont1$site<-c(rep("Bavarian Forest NP",22))
cont2<-cbind(cont2.sp2$sites,cont2.sp2$richness,cont2.sp2$sd)
colnames(cont2)<-c("samples","richness","sd")
cont2<-data.frame(cont2)
dim(cont2)
cont2$site<-c(rep("Bayreuth",14))
cont<-cbind(cont.sp2$sites,cont.sp2$richness,cont.sp2$sd)
colnames(cont)<-c("samples","richness","sd")
cont<-data.frame(cont)
dim(cont)
cont$site<-c(rep("all c-Eur highland sites",36))

cont<-rbind(cont1,cont2,cont)
cont$lower<-cont$richness-cont$sd
cont$upper<-cont$richness+cont$sd

levels(cont$site)
cont$site<-ordered(cont$site,levels=c("all c-Eur highland sites","Bayreuth","Bavarian Forest NP"))

cont.plot<-ggplot()+
  geom_ribbon(data=cont,aes(x=samples,ymin=lower,ymax=upper,fill=site),alpha=.3)+
  geom_line(data=cont,aes(x=samples,y=richness,colour=site,lty=site),size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Species richness")))+
  xlab(expression(bold(" ")))+
  ggtitle(expression(bold("Central European highlands")))+
  theme(plot.title=element_text(size=14,hjust=0))+
  scale_fill_discrete(name="Study site")+
  scale_linetype_discrete(name="Study site")+
  scale_color_discrete(name="Study site")+
  theme(legend.title=element_text(size=12,face="bold"))+
  theme(legend.text=element_text(size=10))+
  theme(legend.position=c(0.75,0.75))+
  coord_cartesian(ylim=c(0,45))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
cont.plot

#Pannonian
pan1<-cbind(pan1.sp2$sites,pan1.sp2$richness,pan1.sp2$sd)
colnames(pan1)<-c("samples","richness","sd")
pan1<-data.frame(pan1)
dim(pan1)
pan1$site<-c(rep("Bugac",72))

pan<-pan1
pan$lower<-pan$richness-pan$sd
pan$upper<-pan$richness+pan$sd

pan.plot<-ggplot()+
  geom_ribbon(data=pan,aes(x=samples,ymin=lower,ymax=upper,fill=site),alpha=.3)+
  geom_line(data=pan,aes(x=samples,y=richness,colour=site,lty=site),size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold(" ")))+
  xlab(expression(bold("Sample number")))+
  ggtitle(expression(bold("Pannonian")))+
  theme(plot.title=element_text(size=14,hjust=0))+
  scale_fill_discrete(name="Study site")+
  scale_linetype_discrete(name="Study site")+
  scale_color_discrete(name="Study site")+
  theme(legend.title=element_text(size=12,face="bold"))+
  theme(legend.text=element_text(size=10))+
  theme(legend.position=c(0.75,0.2))+
  coord_cartesian(ylim=c(0,45))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
pan.plot


#British islands
brit1<-cbind(brit1.sp2$sites,brit1.sp2$richness,brit1.sp2$sd)
colnames(brit1)<-c("samples","richness","sd")
brit1<-data.frame(brit1)
dim(brit1)
brit1$site<-c(rep("Wytham woods",9))
brit2<-cbind(brit2.sp2$sites,brit2.sp2$richness,brit2.sp2$sd)
colnames(brit2)<-c("samples","richness","sd")
brit2<-data.frame(brit2)
dim(brit2)
brit2$site<-c(rep("Moor House NP",34))
brit<-cbind(brit.sp2$sites,brit.sp2$richness,brit.sp2$sd)
colnames(brit)<-c("samples","richness","sd")
brit<-data.frame(brit)
dim(brit)
brit$site<-c(rep("all British sites",43))

brit<-rbind(brit1,brit2,brit)
brit$lower<-brit$richness-brit$sd
brit$upper<-brit$richness+brit$sd

levels(brit$site)
brit$site<-ordered(brit$site,levels=c("all British sites","Moor House NP","Wytham woods"))

brit.plot<-ggplot()+
  geom_ribbon(data=brit,aes(x=samples,ymin=lower,ymax=upper,fill=site),alpha=.3)+
  geom_line(data=brit,aes(x=samples,y=richness,colour=site,lty=site),size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Species richness")))+
  xlab(expression(bold(" ")))+
  ggtitle(expression(bold("British islands")))+
  theme(plot.title=element_text(size=14,hjust=0))+
  scale_fill_discrete(name="Study site")+
  scale_linetype_discrete(name="Study site")+
  scale_color_discrete(name="Study site")+
  theme(legend.title=element_text(size=12,face="bold"))+
  theme(legend.text=element_text(size=10))+
  theme(legend.position=c(0.75,0.75))+
  coord_cartesian(ylim=c(0,45))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
brit.plot


#Caucaso-Iranian highlands
ir1<-cbind(ir1.sp2$sites,ir1.sp2$richness,ir1.sp2$sd)
colnames(ir1)<-c("samples","richness","sd")
ir1<-data.frame(ir1)
dim(ir1)
ir1$site<-c(rep("Tange Sayad",3))
ir2<-cbind(ir2.sp2$sites,ir2.sp2$richness,ir2.sp2$sd)
colnames(ir2)<-c("samples","richness","sd")
ir2<-data.frame(ir2)
dim(ir2)
ir2$site<-c(rep("Shahrekord",7))
ir<-cbind(ir.sp2$sites,ir.sp2$richness,ir.sp2$sd)
colnames(ir)<-c("samples","richness","sd")
ir<-data.frame(ir)
dim(ir)
ir$site<-c(rep("all Caucaso-Iranian sites",10))

ir<-rbind(ir1,ir2,ir)
ir$lower<-ir$richness-ir$sd
ir$upper<-ir$richness+ir$sd

levels(ir$site)
ir$site<-ordered(ir$site,levels=c("all Caucaso-Iranian sites","Shahrekord","Tange Sayad"))

ir.plot<-ggplot()+
  geom_ribbon(data=ir,aes(x=samples,ymin=lower,ymax=upper,fill=site),alpha=.3)+
  geom_line(data=ir,aes(x=samples,y=richness,colour=site,lty=site),size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold(" ")))+
  xlab(expression(bold(" ")))+
  ggtitle(expression(bold("Caucaso-Iranian highlands")))+
  theme(plot.title=element_text(size=14,hjust=0))+
  scale_fill_discrete(name="Study site")+
  scale_linetype_discrete(name="Study site")+
  scale_color_discrete(name="Study site")+
  theme(legend.title=element_text(size=12,face="bold"))+
  theme(legend.text=element_text(size=10))+
  theme(legend.position=c(0.75,0.75))+
  coord_cartesian(ylim=c(0,45))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
ir.plot


#Boreonemoral
est<-cbind(est1.sp2$sites,est1.sp2$richness,est1.sp2$sd)
colnames(est)<-c("samples","richness","sd")
est<-data.frame(est)
dim(est)
est$site<-c(rep("Tähtvere parish",32))

bor<-est
bor$lower<-bor$richness-bor$sd
bor$upper<-bor$richness+bor$sd

bor.plot<-ggplot()+
  geom_ribbon(data=bor,aes(x=samples,ymin=lower,ymax=upper,fill=site),alpha=.3)+
  geom_line(data=bor,aes(x=samples,y=richness,colour=site,lty=site),size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold(" ")))+
  xlab(expression(bold(" ")))+
  ggtitle(expression(bold("Boreonemoral")))+
  theme(plot.title=element_text(size=14,hjust=0))+
  scale_fill_discrete(name="Study site")+
  scale_linetype_discrete(name="Study site")+
  scale_color_discrete(name="Study site")+
  theme(legend.title=element_text(size=12,face="bold"))+
  theme(legend.text=element_text(size=10))+
  theme(legend.position=c(0.75,0.75))+
  coord_cartesian(ylim=c(0,45))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
bor.plot


#Iberian
ib<-cbind(ib1.sp2$sites,ib1.sp2$richness,ib1.sp2$sd)
colnames(ib)<-c("samples","richness","sd")
ib<-data.frame(ib)
dim(ib)
ib$site<-c(rep("Castillo de Viñuelas",6))

ib$lower<-ib$richness-ib$sd
ib$upper<-ib$richness+ib$sd

ib.plot<-ggplot()+
  geom_ribbon(data=ib,aes(x=samples,ymin=lower,ymax=upper,fill=site),alpha=.3)+
  geom_line(data=ib,aes(x=samples,y=richness,colour=site,lty=site),size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold(" ")))+
  xlab(expression(bold(" ")))+
  ggtitle(expression(bold("Iberian highlands")))+
  theme(plot.title=element_text(size=14,hjust=0))+
  scale_fill_discrete(name="Study site")+
  scale_linetype_discrete(name="Study site")+
  scale_color_discrete(name="Study site")+
  theme(legend.title=element_text(size=12,face="bold"))+
  theme(legend.text=element_text(size=10))+
  theme(legend.position=c(0.75,0.75))+
  coord_cartesian(ylim=c(0,45))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
ib.plot


#Mediterranean sclerofyl
med<-cbind(med1.sp2$sites,med1.sp2$richness,med1.sp2$sd)
colnames(med)<-c("samples","richness","sd")
med<-data.frame(med)
dim(med)
med$site<-c(rep("INRA, La Fage",10))

med$lower<-med$richness-med$sd
med$upper<-med$richness+med$sd

med.plot<-ggplot()+
  geom_ribbon(data=med,aes(x=samples,ymin=lower,ymax=upper,fill=site),alpha=.3)+
  geom_line(data=med,aes(x=samples,y=richness,colour=site,lty=site),size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Species richness")))+
  xlab(expression(bold(" ")))+
  ggtitle(expression(bold("Mediterranean sclerofyl")))+
  theme(plot.title=element_text(size=14,hjust=0))+
  scale_fill_discrete(name="Study site")+
  scale_linetype_discrete(name="Study site")+
  scale_color_discrete(name="Study site")+
  theme(legend.title=element_text(size=12,face="bold"))+
  theme(legend.text=element_text(size=10))+
  theme(legend.position=c(0.75,0.75))+
  coord_cartesian(ylim=c(0,45))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
med.plot

#Middle European forest
mef<-cbind(mef1.sp2$sites,mef1.sp2$richness,mef1.sp2$sd)
colnames(mef)<-c("samples","richness","sd")
mef<-data.frame(mef)
dim(mef)
mef$site<-c(rep("Vácrátót",36))

mef$lower<-mef$richness-mef$sd
mef$upper<-mef$richness+mef$sd

mef.plot<-ggplot()+
  geom_ribbon(data=mef,aes(x=samples,ymin=lower,ymax=upper,fill=site),alpha=.3)+
  geom_line(data=mef,aes(x=samples,y=richness,colour=site,lty=site),size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold(" ")))+
  xlab(expression(bold(" ")))+
  ggtitle(expression(bold("Middle European forest")))+
  theme(plot.title=element_text(size=14,hjust=0))+
  scale_fill_discrete(name="Study site")+
  scale_linetype_discrete(name="Study site")+
  scale_color_discrete(name="Study site")+
  theme(legend.title=element_text(size=12,face="bold"))+
  theme(legend.text=element_text(size=10))+
  theme(legend.position=c(0.75,0.2))+
  coord_cartesian(ylim=c(0,45))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
mef.plot

#Pontian steppe
pon<-cbind(pon1.sp2$sites,pon1.sp2$richness,pon1.sp2$sd)
colnames(pon)<-c("samples","richness","sd")
pon<-data.frame(pon)
dim(pon)
pon$site<-c(rep("Braila Islands",32))

pon$lower<-pon$richness-pon$sd
pon$upper<-pon$richness+pon$sd

pon.plot<-ggplot()+
  geom_ribbon(data=pon,aes(x=samples,ymin=lower,ymax=upper,fill=site),alpha=.3)+
  geom_line(data=pon,aes(x=samples,y=richness,colour=site,lty=site),size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Species richness")))+
  xlab(expression(bold("Sample number")))+
  ggtitle(expression(bold("Pontian steppe")))+
  theme(plot.title=element_text(size=14,hjust=0))+
  scale_fill_discrete(name="Study site")+
  scale_linetype_discrete(name="Study site")+
  scale_color_discrete(name="Study site")+
  theme(legend.title=element_text(size=12,face="bold"))+
  theme(legend.text=element_text(size=10))+
  theme(legend.position=c(0.75,0.75))+
  coord_cartesian(ylim=c(0,45))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
pon.plot


#combinatiegrafiek maken
library(gridExtra)

tiff(filename="SamplingEffort_v2.tiff",width=800,height=1500,units="px")
grid.arrange(atl.plot,bor.plot,brit.plot,
             ir.plot,cont.plot,ib.plot,
             med.plot, mef.plot,pon.plot,pan.plot,ncol=2)
dev.off()



####Ecological functions####
rm(list=ls())

setwd('D:/Tmilotic/Documents/Doctoraat/8 MSE Europe/R scripts')
db<-read.delim("dungSeedRemovalRelFG.txt",header=T)
head(db)


####hierarchical partitioning####
library(hier.part)

colnames(db)[names(db)== "Biogeographic.provinceUdvardy1975"] <- "biogeo"
colnames(db)[names(db)== "Koppen.climate.zonePeel2007"] <- "climate"
summary(db)

#removal as proportion in stead of %
db$DungRemoval<-db$DungRemovalCorr/100


####per biogeographic region, alle mestsoorten samen####
head(db)
levels(db$biogeo)
atl.db<-subset(db,biogeo=="Atlantic")
bor.db<-subset(db,biogeo=="Boreonemoral")
brit.db<-subset(db,biogeo=="British islands")
ir.db<-subset(db,biogeo=="Caucaso-Iranian highlands")
cont.db<-subset(db,biogeo=="Central European highlands")
ib.db<-subset(db,biogeo=="Iberian highlands")
med.db<-subset(db,biogeo=="Mediterranean sclerophyl")
mef.db<-subset(db,biogeo=="Middle European forests")
pan.db<-subset(db,biogeo=="Pannonian")
pont.db<-subset(db,biogeo=="Pontian Steppe")
taig.db<-subset(db,biogeo=="West Eurasian taiga")

#atl
atl.vars<-c("D","Tl","ts","rs","S")
atl.env<-atl.db[atl.vars]
atl.hier<-hier.part(atl.db$DungRemoval,atl.env,family = "quasibinomial")
atl.hier
atl.hier$I.perc #values for graph

atl.prob<-rand.hp(atl.db$DungRemoval,atl.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
atl.prob$Iprobs #Z scores for graph

#Rdev model berekenen
atl.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=atl.db,family="quasibinomial")
atl.glm
atl.sumglm<-summary.glm(atl.glm)
atl.sumglm$deviance
atl.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
atl.Rdev<-(atl.sumglm$null.deviance-atl.sumglm$deviance)/atl.sumglm$null.deviance
atl.Rdev

#no rolling species captured and not to be expected
atl2.vars<-c("Drel","Tlrel","tsrel","Srel")
atl2.env<-atl.db[atl2.vars]
atl2.hier<-hier.part(atl.db$DungRemoval,atl2.env,family = "quasibinomial")
atl2.hier
atl2.hier$I.perc #values for graph

atl2.prob<-rand.hp(atl.db$DungRemoval,atl2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
atl2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
atl2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Srel,data=atl.db,family="quasibinomial")
atl2.glm<-glm(DungRemoval~Drel,data=atl.db,family="quasibinomial")
atl2.glm<-glm(DungRemoval~Tlrel,data=atl.db,family="quasibinomial")
atl2.glm<-glm(DungRemoval~tsrel,data=atl.db,family="quasibinomial")
atl2.glm<-glm(DungRemoval~Srel,data=atl.db,family="quasibinomial")

atl2.sumglm<-summary.glm(atl2.glm)
atl2.sumglm$deviance
atl2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
atl2.Rdev<-(atl2.sumglm$null.deviance-atl2.sumglm$deviance)/atl2.sumglm$null.deviance
atl2.Rdev


#bor
bor.vars<-c("D","Tl","ts","S")
bor.env<-bor.db[bor.vars]
bor.hier<-hier.part(bor.db$DungRemoval,bor.env,family = "quasibinomial")
bor.hier
bor.hier$I.perc #values for graph

bor.prob<-rand.hp(bor.db$DungRemoval,bor.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
bor.prob$Iprobs #Z scores for graph

#Rdev model berekenen
bor.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=bor.db,family="quasibinomial")
bor.glm
bor.sumglm<-summary.glm(bor.glm)
bor.sumglm$deviance
bor.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
bor.Rdev<-(bor.sumglm$null.deviance-bor.sumglm$deviance)/bor.sumglm$null.deviance
bor.Rdev

#no rolling species captured and not to be expected
bor2.vars<-c("Drel","Tlrel","tsrel","Srel")
bor2.env<-bor.db[bor2.vars]
bor2.hier<-hier.part(bor.db$DungRemoval,bor2.env,family = "quasibinomial")
bor2.hier
bor2.hier$I.perc #values for graph

bor2.prob<-rand.hp(bor.db$DungRemoval,bor2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
bor2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
bor2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Srel,data=bor.db,family="quasibinomial")
bor2.glm<-glm(DungRemoval~Drel,data=bor.db,family="quasibinomial")
bor2.glm<-glm(DungRemoval~tsrel,data=bor.db,family="quasibinomial")
bor2.glm<-glm(DungRemoval~Srel,data=bor.db,family="quasibinomial")

bor2.glm
bor2.sumglm<-summary.glm(bor2.glm)
bor2.sumglm$deviance
bor2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
bor2.Rdev<-(bor2.sumglm$null.deviance-bor2.sumglm$deviance)/bor2.sumglm$null.deviance
bor2.Rdev


#brit
brit.vars<-c("D","Tl","ts","Rl","rs","S")
brit.env<-brit.db[brit.vars]
brit.hier<-hier.part(brit.db$DungRemoval,brit.env,family = "quasibinomial")
brit.hier
brit.hier$I.perc #values for graph

brit.prob<-rand.hp(brit.db$DungRemoval,brit.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
brit.prob$Iprobs #Z scores for graph

#Rdev model berekenen
brit.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=brit.db,family="quasibinomial")
brit.glm
brit.sumglm<-summary.glm(brit.glm)
brit.sumglm$deviance
brit.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
brit.Rdev<-(brit.sumglm$null.deviance-brit.sumglm$deviance)/brit.sumglm$null.deviance
brit.Rdev

#no rolling species captured and not to be expected
brit2.vars<-c("Drel","Tlrel","tsrel","Srel")
brit2.env<-brit.db[brit2.vars]
brit2.hier<-hier.part(brit.db$DungRemoval,brit2.env,family = "quasibinomial")
brit2.hier
brit2.hier$I.perc #values for graph

brit2.prob<-rand.hp(brit.db$DungRemoval,brit2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
brit2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
brit2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Srel,data=brit.db,family="quasibinomial")
brit2.glm<-glm(DungRemoval~Drel,data=brit.db,family="quasibinomial")
brit2.glm<-glm(DungRemoval~tsrel,data=brit.db,family="quasibinomial")
brit2.glm<-glm(DungRemoval~Srel,data=brit.db,family="quasibinomial")

brit2.sumglm<-summary.glm(brit2.glm)
brit2.sumglm$deviance
brit2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
brit2.Rdev<-(brit2.sumglm$null.deviance-brit2.sumglm$deviance)/brit2.sumglm$null.deviance
brit2.Rdev



#cont
cont.vars<-c("D","Tl","ts","Rl","rs","S")
cont.env<-cont.db[cont.vars]
cont.hier<-hier.part(cont.db$DungRemoval,cont.env,family = "quasibinomial")
cont.hier
cont.hier$I.perc #values for graph

cont.prob<-rand.hp(cont.db$DungRemoval,cont.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
cont.prob$Iprobs #Z scores for graph

#Rdev model berekenen
cont.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=cont.db,family="quasibinomial")
cont.glm
cont.sumglm<-summary.glm(cont.glm)
cont.sumglm$deviance
cont.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
cont.Rdev<-(cont.sumglm$null.deviance-cont.sumglm$deviance)/cont.sumglm$null.deviance
cont.Rdev

#no rolling species captured MAAR EFFECT???
summary(cont.db)
cont.db2a<- subset(cont.db,SiteCode=="D_BAV")
cont.db2b<- subset(cont.db,((SiteCode=="D_BR")&(Year==2015)))
cont.db2<-rbind(cont.db2a,cont.db2b)
cont2.vars<-c("Drel","Tlrel","tsrel","Srel")
cont2.env<-cont.db2[cont2.vars]
cont2.hier<-hier.part(cont.db2$DungRemoval,cont2.env,family = "quasibinomial")
cont2.hier
cont2.hier$I.perc #values for graph

cont2.prob<-rand.hp(cont.db2$DungRemoval,cont2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
cont2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
cont2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Srel,data=cont.db2,family="quasibinomial")
cont2.glm<-glm(DungRemoval~Drel,data=cont.db2,family="quasibinomial")
cont2.glm<-glm(DungRemoval~Tlrel,data=cont.db2,family="quasibinomial")
cont2.glm<-glm(DungRemoval~tsrel,data=cont.db2,family="quasibinomial")
cont2.glm<-glm(DungRemoval~Srel,data=cont.db2,family="quasibinomial")

cont2.sumglm<-summary.glm(cont2.glm)
cont2.sumglm$deviance
cont2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
cont2.Rdev<-(cont2.sumglm$null.deviance-cont2.sumglm$deviance)/cont2.sumglm$null.deviance
cont2.Rdev


#ir
ir.vars<-c("D","Tl","ts","Rl","rs","S")
ir.env<-ir.db[ir.vars]
ir.hier<-hier.part(ir.db$DungRemoval,ir.env,family = "quasibinomial")
ir.hier
ir.hier$I.perc #values for graph

ir.prob<-rand.hp(ir.db$DungRemoval,ir.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
ir.prob$Iprobs #Z scores for graph

#Rdev model berekenen
ir.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=ir.db,family="quasibinomial")
ir.glm
ir.sumglm<-summary.glm(ir.glm)
ir.sumglm$deviance
ir.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
ir.Rdev<-(ir.sumglm$null.deviance-ir.sumglm$deviance)/ir.sumglm$null.deviance
ir.Rdev

#rel abundance
summary(ir.db)
ir.db2<-subset(ir.db,SiteCode=="IR_SR")
ir2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
ir2.env<-ir.db2[ir2.vars]
ir2.hier<-hier.part(ir.db2$DungRemoval,ir2.env,family = "quasibinomial")
ir2.hier
ir2.hier$I.perc #values for graph

ir2.prob<-rand.hp(ir.db2$DungRemoval,ir2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
ir2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
ir2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+rsrel+Srel,data=ir.db2,family="quasibinomial")
ir2.glm<-glm(DungRemoval~Drel,data=ir.db2,family="quasibinomial")
ir2.glm<-glm(DungRemoval~Tlrel,data=ir.db2,family="quasibinomial")
ir2.glm<-glm(DungRemoval~tsrel,data=ir.db2,family="quasibinomial")
ir2.glm<-glm(DungRemoval~rsrel,data=ir.db2,family="quasibinomial")
ir2.glm<-glm(DungRemoval~Srel,data=ir.db2,family="quasibinomial")

ir2.sumglm<-summary.glm(ir2.glm)
ir2.sumglm$deviance
ir2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
ir2.Rdev<-(ir2.sumglm$null.deviance-ir2.sumglm$deviance)/ir2.sumglm$null.deviance
ir2.Rdev


#pan
pan.vars<-c("D","Tl","ts","Rl","rs","S")
pan.env<-pan.db[pan.vars]
pan.hier<-hier.part(pan.db$DungRemoval,pan.env,family = "quasibinomial")
pan.hier
pan.hier$I.perc #values for graph

pan.prob<-rand.hp(pan.db$DungRemoval,pan.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
pan.prob$Iprobs #Z scores for graph

#Rdev model berekenen
pan.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=pan.db,family="quasibinomial")
pan.glm
pan.sumglm<-summary.glm(pan.glm)
pan.sumglm$deviance
pan.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
pan.Rdev<-(pan.sumglm$null.deviance-pan.sumglm$deviance)/pan.sumglm$null.deviance
pan.Rdev

#no rolling species captured, but these results suggest something else
pan2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
pan2.env<-pan.db[pan2.vars]
pan2.hier<-hier.part(pan.db$DungRemoval,pan2.env,family = "quasibinomial")
pan2.hier
pan2.hier$I.perc #values for graph

pan2.prob<-rand.hp(pan.db$DungRemoval,pan2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
pan2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
pan2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Srel,data=pan.db,family="quasibinomial")
pan2.glm<-glm(DungRemoval~Drel,data=pan.db,family="quasibinomial")
pan2.glm<-glm(DungRemoval~Tlrel,data=pan.db,family="quasibinomial")
pan2.glm<-glm(DungRemoval~tsrel,data=pan.db,family="quasibinomial")
pan2.glm<-glm(DungRemoval~Srel,data=pan.db,family="quasibinomial")

pan2.sumglm<-summary.glm(pan2.glm)
pan2.sumglm$deviance
pan2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
pan2.Rdev<-(pan2.sumglm$null.deviance-pan2.sumglm$deviance)/pan2.sumglm$null.deviance
pan2.Rdev


#pont
pont.vars<-c("D","Tl","ts","Rl","rs","S")
pont.env<-pont.db[pont.vars]
summary(pont.db)
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

pont.db<-completeFun(pont.db, "D")
pont.db<-completeFun(pont.db, "DungRemovalCorr")

pont.vars<-c("D","Tl","ts","rs","S")
pont.env<-pont.db[pont.vars]
pont.hier<-hier.part(pont.db$DungRemoval,pont.env,family = "quasibinomial")
pont.hier
pont.hier$I.perc #values for graph

pont.prob<-rand.hp(pont.db$DungRemoval,pont.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
pont.prob$Iprobs #Z scores for graph

#Rdev model berekenen
pont.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=pont.db,family="quasibinomial")
pont.glm<-glm(DungRemoval~D,data=pont.db,family="quasibinomial")
pont.glm<-glm(DungRemoval~Tl,data=pont.db,family="quasibinomial")
pont.glm<-glm(DungRemoval~ts,data=pont.db,family="quasibinomial")
pont.glm<-glm(DungRemoval~Rl,data=pont.db,family="quasibinomial")
pont.glm<-glm(DungRemoval~rs,data=pont.db,family="quasibinomial")
pont.glm<-glm(DungRemoval~S,data=pont.db,family="quasibinomial")

pont.sumglm<-summary.glm(pont.glm)
pont.sumglm$deviance
pont.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
pont.Rdev<-(pont.sumglm$null.deviance-pont.sumglm$deviance)/pont.sumglm$null.deviance
pont.Rdev

#no rolling species captured, but these results suggest something else
pont2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
pont2.env<-pont.db[pont2.vars]
pont2.hier<-hier.part(pont.db$DungRemoval,pont2.env,family = "quasibinomial")
pont2.hier
pont2.hier$I.perc #values for graph

pont2.prob<-rand.hp(pont.db$DungRemoval,pont2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
pont2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
pont2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Srel,data=pont.db,family="quasibinomial")
pont2.glm<-glm(DungRemoval~Drel,data=pont.db,family="quasibinomial")
pont2.glm<-glm(DungRemoval~Tlrel,data=pont.db,family="quasibinomial")
pont2.glm<-glm(DungRemoval~tsrel,data=pont.db,family="quasibinomial")
pont2.glm<-glm(DungRemoval~Srel,data=pont.db,family="quasibinomial")

pont2.sumglm<-summary.glm(pont2.glm)
pont2.sumglm$deviance
pont2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
pont2.Rdev<-(pont2.sumglm$null.deviance-pont2.sumglm$deviance)/pont2.sumglm$null.deviance
pont2.Rdev



#taig
taig.vars<-c("D","Tl","ts","Rl","rs","S")
taig.env<-taig.db[taig.vars]
taig.hier<-hier.part(taig.db$DungRemoval,taig.env,family = "quasibinomial")
taig.hier
taig.hier$I.perc #values for graph

taig.prob<-rand.hp(taig.db$DungRemoval,taig.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
taig.prob$Iprobs #Z scores for graph

#Rdev model berekenen
taig.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=taig.db,family="quasibinomial")
taig.glm
taig.sumglm<-summary.glm(taig.glm)
taig.sumglm$deviance
taig.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
taig.Rdev<-(taig.sumglm$null.deviance-taig.sumglm$deviance)/taig.sumglm$null.deviance
taig.Rdev

#no rolling species captured and not to be expected
taig2.vars<-c("Drel","Tlrel","tsrel","Srel")
taig2.env<-taig.db[taig2.vars]
taig2.hier<-hier.part(taig.db$DungRemoval,taig2.env,family = "quasibinomial")
taig2.hier
taig2.hier$I.perc #values for graph

taig2.prob<-rand.hp(taig.db$DungRemoval,taig2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
taig2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
taig2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Srel,data=taig.db,family="quasibinomial")
taig2.glm<-glm(DungRemoval~Drel,data=taig.db,family="quasibinomial")
taig2.glm<-glm(DungRemoval~Srel,data=taig.db,family="quasibinomial")

taig2.sumglm<-summary.glm(taig2.glm)
taig2.sumglm$deviance
taig2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
taig2.Rdev<-(taig2.sumglm$null.deviance-taig2.sumglm$deviance)/taig2.sumglm$null.deviance
taig2.Rdev


#ib
ib.vars<-c("D","Tl","ts","Rl","rs","S")
ib.env<-ib.db[ib.vars]
ib.hier<-hier.part(ib.db$DungRemoval,ib.env,family = "quasibinomial")
ib.hier
ib.hier$I.perc #values for graph

ib.prob<-rand.hp(ib.db$DungRemoval,ib.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
ib.prob$Iprobs #Z scores for graph

#Rdev model berekenen
ib.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=ib.db,family="quasibinomial")
ib.glm
ib.sumglm<-summary.glm(ib.glm)
ib.sumglm$deviance
ib.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
ib.Rdev<-(ib.sumglm$null.deviance-ib.sumglm$deviance)/ib.sumglm$null.deviance
ib.Rdev

#no rolling species captured, but these results suggest something else
ib2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
ib2.env<-ib.db[ib2.vars]
ib2.hier<-hier.part(ib.db$DungRemoval,ib2.env,family = "quasibinomial")
ib2.hier
ib2.hier$I.perc #values for graph

ib2.prob<-rand.hp(ib.db$DungRemoval,ib2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
ib2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
ib2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Srel,data=ib.db,family="quasibinomial")
ib2.glm<-glm(DungRemoval~Drel,data=ib.db,family="quasibinomial")
ib2.glm<-glm(DungRemoval~Tlrel,data=ib.db,family="quasibinomial")
ib2.glm<-glm(DungRemoval~tsrel,data=ib.db,family="quasibinomial")
ib2.glm<-glm(DungRemoval~Srel,data=ib.db,family="quasibinomial")

ib2.sumglm<-summary.glm(ib2.glm)
ib2.sumglm$deviance
ib2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
ib2.Rdev<-(ib2.sumglm$null.deviance-ib2.sumglm$deviance)/ib2.sumglm$null.deviance
ib2.Rdev


#med
med.vars<-c("D","Tl","ts","Rl","rs","S")
med.env<-med.db[med.vars]
med.hier<-hier.part(med.db$DungRemoval,med.env,family = "quasibinomial")
med.hier
med.hier$I.perc #values for graph

med.prob<-rand.hp(med.db$DungRemoval,med.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
med.prob$Iprobs #Z scores for graph

#Rdev model berekenen
med.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=med.db,family="quasibinomial")
med.glm
med.sumglm<-summary.glm(med.glm)
med.sumglm$deviance
med.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
med.Rdev<-(med.sumglm$null.deviance-med.sumglm$deviance)/med.sumglm$null.deviance
med.Rdev

#relative values
med2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
med2.env<-med.db[med2.vars]
med2.hier<-hier.part(med.db$DungRemoval,med2.env,family = "quasibinomial")
med2.hier
med2.hier$I.perc #values for graph

med2.prob<-rand.hp(med.db$DungRemoval,med2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
med2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
med2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Rlrel+rsrel+Srel,data=med.db,family="quasibinomial")
med2.glm<-glm(DungRemoval~Drel,data=med.db,family="quasibinomial")
med2.glm<-glm(DungRemoval~Tlrel,data=med.db,family="quasibinomial")
med2.glm<-glm(DungRemoval~tsrel,data=med.db,family="quasibinomial")
med2.glm<-glm(DungRemoval~rsrel,data=med.db,family="quasibinomial")
med2.glm<-glm(DungRemoval~Srel,data=med.db,family="quasibinomial")

med2.sumglm<-summary.glm(med2.glm)
med2.sumglm$deviance
med2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
med2.Rdev<-(med2.sumglm$null.deviance-med2.sumglm$deviance)/med2.sumglm$null.deviance
med2.Rdev



#mef
mef.vars<-c("D","Tl","ts","Rl","rs","S")
mef.env<-mef.db[mef.vars]
mef.hier<-hier.part(mef.db$DungRemoval,mef.env,family = "quasibinomial")
mef.hier
mef.hier$I.perc #values for graph

mef.prob<-rand.hp(mef.db$DungRemoval,mef.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
mef.prob$Iprobs #Z scores for graph

#Rdev model berekenen
mef.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=mef.db,family="quasibinomial")
mef.glm
mef.sumglm<-summary.glm(mef.glm)
mef.sumglm$deviance
mef.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
mef.Rdev<-(mef.sumglm$null.deviance-mef.sumglm$deviance)/mef.sumglm$null.deviance
mef.Rdev

#no rolling species captured 
mef2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
mef2.env<-mef.db[mef2.vars]
mef2.hier<-hier.part(mef.db$DungRemoval,mef2.env,family = "quasibinomial")
mef2.hier
mef2.hier$I.perc #values for graph

mef2.prob<-rand.hp(mef.db$DungRemoval,mef2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
mef2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
mef2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Srel,data=mef.db,family="quasibinomial")
mef2.glm<-glm(DungRemoval~Drel,data=mef.db,family="quasibinomial")
mef2.glm<-glm(DungRemoval~Tlrel,data=mef.db,family="quasibinomial")
mef2.glm<-glm(DungRemoval~tsrel,data=mef.db,family="quasibinomial")
mef2.glm<-glm(DungRemoval~Srel,data=mef.db,family="quasibinomial")

mef2.sumglm<-summary.glm(mef2.glm)
mef2.sumglm$deviance
mef2.sumglm$null.deviance
mef2.Rdev<-(mef2.sumglm$null.deviance-mef2.sumglm$deviance)/mef2.sumglm$null.deviance
mef2.Rdev


#grafiek maken
atl<-atl2.hier$I.perc
atl<-rbind(atl,c(0,0))
atl<-rbind(atl,c(0,0))
atl$FG<-c("D","Tl","ts","S","Rl","rs")
atl$biogeo<-rep("atl",6)
atl$biogeo2<-rep("2atl",6)
atl
atl.probs<-atl2.prob$Iprobs
atl.probs<-rbind(atl.probs,c(0,0,""))
atl.probs<-rbind(atl.probs,c(0,0,""))
atl.probs
atl<-cbind(atl,atl.probs)
atl

bor<-bor2.hier$I.perc
bor<-rbind(bor,c(0,0))
bor<-rbind(bor,c(0,0))
bor$FG<-c("D","Tl","ts","S","Rl","rs")
bor$biogeo<-rep("bor",6)
bor$biogeo2<-rep("4bor",6)
bor
bor.probs<-bor2.prob$Iprobs
bor.probs<-rbind(bor.probs,c(0,0,""))
bor.probs<-rbind(bor.probs,c(0,0,""))
bor.probs
bor<-cbind(bor,bor.probs)
bor

brit<-brit2.hier$I.perc
brit<-rbind(brit,c(0,0))
brit<-rbind(brit,c(0,0))
brit$FG<-c("D","Tl","ts","S","Rl","rs")
brit$biogeo<-rep("brit",6)
brit$biogeo2<-rep("1brit",6)
brit
brit.probs<-brit2.prob$Iprobs
brit.probs<-rbind(brit.probs,c(0,0,""))
brit.probs<-rbind(brit.probs,c(0,0,""))
brit.probs
brit<-cbind(brit,brit.probs)
brit

cont<-cont2.hier$I.perc
cont<-rbind(cont,c(0,0))
cont<-rbind(cont,c(0,0))
cont$FG<-c("D","Tl","ts","S","Rl","rs")
cont$biogeo<-rep("cont",6)
cont$biogeo2<-rep("1cont",6)
cont
cont.probs<-cont2.prob$Iprobs
cont.probs<-rbind(cont.probs,c(0,0,""))
cont.probs<-rbind(cont.probs,c(0,0,""))
cont.probs
cont<-cbind(cont,cont.probs)
cont


ir<-ir2.hier$I.perc
ir$FG<-c("D","Tl","ts","Rl","rs","S")
ir$biogeo<-rep("cau",6)
ir$biogeo2<-rep("8cau",6)
ir
ir.probs<-ir2.prob$Iprobs
ir<-cbind(ir,ir.probs)
ir

pan<-pan2.hier$I.perc
pan$FG<-c("D","Tl","ts","Rl","rs","S")
pan$biogeo<-rep("pan",6)
pan$biogeo2<-rep("pan",6)
pan
pan.probs<-pan2.prob$Iprobs
pan<-cbind(pan,pan.probs)
pan

pont<-pont2.hier$I.perc
pont$FG<-c("D","Tl","ts","Rl","rs","S")
pont$biogeo<-rep("pont",6)
pont$biogeo2<-rep("7pont",6)
pont
pont.probs<-pont2.prob$Iprobs
pont<-cbind(pont,pont.probs)
pont

taig<-taig2.hier$I.perc
taig<-rbind(taig,c(0,0))
taig<-rbind(taig,c(0,0))
taig$FG<-c("D","Tl","ts","S","Rl","rs")
taig$biogeo<-rep("taig",6)
taig$biogeo2<-rep("3taig",6)
taig
taig.probs<-taig2.prob$Iprobs
taig.probs<-rbind(taig.probs,c(0,0,""))
taig.probs<-rbind(taig.probs,c(0,0,""))
taig.probs
taig<-cbind(taig,taig.probs)
taig

ib<-ib2.hier$I.perc
ib$FG<-c("D","Tl","ts","Rl","rs","S")
ib$biogeo<-rep("ib",6)
ib$biogeo2<-rep("ib",6)
ib
ib.probs<-ib2.prob$Iprobs
ib<-cbind(ib,ib.probs)
ib

med<-med2.hier$I.perc
med$FG<-c("D","Tl","ts","Rl","rs","S")
med$biogeo<-rep("med",6)
med$biogeo2<-rep("med",6)
med
med.probs<-med2.prob$Iprobs
med<-cbind(med,med.probs)
med

mef<-mef2.hier$I.perc
mef$FG<-c("D","Tl","ts","Rl","rs","S")
mef$biogeo<-rep("mef",6)
mef$biogeo2<-rep("mef",6)
mef
mef.probs<-mef2.prob$Iprobs
mef<-cbind(mef,mef.probs)
mef


hier<-rbind(atl,bor,brit,cont,ir,pan,pont,taig,ib,med,mef)
hier

hier$sig95b<-ifelse(hier$I==0,"",ifelse(hier$sig95=="*","*",""))


library(ggplot2)
hier.plot<-ggplot(data=hier,aes(x=FG,y=I,fill=biogeo,label=sig95b))+
  geom_bar(position="dodge",colour="black",stat="identity")
hier.plot
hier.plot<-hier.plot+theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Independent effect (%)")))+
  xlab(expression(bold(" ")))+
  scale_x_discrete(limits=c("D","ts","Tl","rs","Rl","S"),
                   breaks=c("D","ts","Tl","rs","Rl","S"),
                   labels=c("dwellers","small tunnelers","large tunnelers","small rollers","large rollers","soil fauna"))+
  theme(axis.text.x=element_blank())+
  theme(axis.ticks=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(hier)/10), y = -4, label = c("dwellers","small tunnelers","large tunnelers","small rollers","large rollers","soil fauna"), size = 5)+ 
  geom_text(aes(y=(I+0.02*(max(I)))),position=position_dodge(width=.9),size=5)
hier.plot
hier.plot<-hier.plot+
  scale_fill_manual(name=expression(bold("Biogeographical province")),
                    values=c("lightgoldenrod","thistle4","coral3","forestgreen","green3","goldenrod2","olivedrab3","aquamarine2","chocolate4","maroon1","steelblue4"),
                    labels=c("Atlantic","Boreonemoral","British islands",
                              "Caucaso-Iranian highlands","Central European highlands",
                              "Iberian highlands","Mediterranean sclerophyl","Middle European forests",
                              "Pannonian","Pontian steppe","West Eurasian taiga"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))
hier.plot

tiff(filename="biogeo.hier.plot.tiff",width=900,height=400,unit="px")
hier.plot
dev.off()

hier$FG<-as.factor(hier$FG)
levels(hier$FG)
hier$FG2<-ifelse(hier$FG=="S",1,(ifelse(hier$FG=="rs",2,(ifelse(hier$FG=="Rl",3,(ifelse(hier$FG=="ts",4,(ifelse(hier$FG=="Tl",5,6)))))))))
hier$FG3<-ifelse(hier$FG=="S",6,(ifelse(hier$FG=="rs",5,(ifelse(hier$FG=="Rl",4,(ifelse(hier$FG=="ts",3,(ifelse(hier$FG=="Tl",2,1)))))))))
hier$FG<-reorder(hier$FG,hier$FG3)

hier<-subset(hier,hier$FG!="Rl")

hier.plot2<-ggplot(hier[order(hier$FG2),],aes(x=biogeo,y=I,fill=FG))+
  geom_bar(position="stack",stat="identity")
hier.plot2
hier.plot2<-hier.plot2+theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Independent effect (%)")))+
  xlab(expression(bold("Biogeographic province")))+
  theme(axis.ticks.x=element_blank())+
  theme(axis.text.x=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(hier)/5), y = -4, label = c("atl","bor","brit","cau","ceh","ibh","med","mef","pan","pont","taig"), size = 5)
hier.plot2
hier.plot2<-hier.plot2+
  scale_fill_manual(name=expression(bold("Functional group")),
                    values=c("royalblue4","red4","orange","yellowgreen","grey"),
                    labels=c("dwellers","tunnelers (large)","tunnelers (small)","rollers (small)","soil macro-invertebrates"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))
hier.plot2

#semi-manueel sterretjes zetten
hier #via excel nieuwe dataset maken
biogeo.hier.label<-read.delim("biogeo.hier.labels.txt",header=T)
head(biogeo.hier.label)
hier.plot2<-hier.plot2+
  geom_text(data=biogeo.hier.label,aes(x=biogeo,y=labelY,color=color,label=sig95b),size=8,show.legend = F)+
  scale_colour_manual(values=c("white","black"))
hier.plot2


####hierarchical partitioning per climate zone####
head(db)
levels(db$climate)
BWk.db<-subset(db,climate=="BWk")
Cfb.db<-subset(db,climate=="Cfb")
Cfc.db<-subset(db,climate=="Cfc")
Csb.db<-subset(db,climate=="Csb")
Dfa.db<-subset(db,climate=="Dfa")
Dfb.db<-subset(db,climate=="Dfb")


#BWk
BWk.vars<-c("D","Tl","ts","Rl","rs","S")
BWk.env<-BWk.db[BWk.vars]
BWk.hier<-hier.part(BWk.db$DungRemoval,BWk.env,family = "quasibinomial")
BWk.hier
BWk.hier$I.perc #values for graph

BWk.prob<-rand.hp(BWk.db$DungRemoval,BWk.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
BWk.prob$Iprobs #Z scores for graph

#Rdev model berekenen
BWk.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=BWk.db,family="quasibinomial")
BWk.glm
BWk.sumglm<-summary.glm(BWk.glm)
BWk.sumglm$deviance
BWk.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
BWk.Rdev<-(BWk.sumglm$null.deviance-BWk.sumglm$deviance)/BWk.sumglm$null.deviance
BWk.Rdev

#data rel abund
summary(BWk.db)
BWk.db2<-subset(BWk.db,SiteCode=="IR_SR")
BWk2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
BWk2.env<-BWk.db2[BWk2.vars]
BWk2.hier<-hier.part(BWk.db2$DungRemoval,BWk2.env,family = "quasibinomial")
BWk2.hier
BWk2.hier$I.perc #values for graph

BWk2.prob<-rand.hp(BWk.db2$DungRemoval,BWk2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
BWk2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
BWk2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Srel,data=BWk.db2,family="quasibinomial")
BWk2.glm
BWk2.sumglm<-summary.glm(BWk2.glm)
BWk2.sumglm$deviance
BWk2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
BWk2.Rdev<-(BWk2.sumglm$null.deviance-BWk2.sumglm$deviance)/BWk2.sumglm$null.deviance
BWk2.Rdev


#Cfb
Cfb.vars<-c("D","Tl","ts","Rl","rs","S")
Cfb.env<-Cfb.db[Cfb.vars]
Cfb.hier<-hier.part(Cfb.db$DungRemoval,Cfb.env,family = "quasibinomial")
Cfb.hier
Cfb.hier$I.perc #values for graph

Cfb.prob<-rand.hp(Cfb.db$DungRemoval,Cfb.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
Cfb.prob$Iprobs #Z scores for graph

#Rdev model berekenen
Cfb.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=Cfb.db,family="quasibinomial")
Cfb.glm
Cfb.sumglm<-summary.glm(Cfb.glm)
Cfb.sumglm$deviance
Cfb.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
Cfb.Rdev<-(Cfb.sumglm$null.deviance-Cfb.sumglm$deviance)/Cfb.sumglm$null.deviance
Cfb.Rdev

#relative abundance
Cfb2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
Cfb2.env<-Cfb.db[Cfb2.vars]
Cfb2.hier<-hier.part(Cfb.db$DungRemoval,Cfb2.env,family = "quasibinomial")
Cfb2.hier
Cfb2.hier$I.perc #values for graph

Cfb2.prob<-rand.hp(Cfb.db$DungRemoval,Cfb2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
Cfb2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
Cfb2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Rlrel+rsrel+Srel,data=Cfb.db,family="quasibinomial")
Cfb2.glm<-glm(DungRemoval~Drel,data=Cfb.db,family="quasibinomial")
Cfb2.glm<-glm(DungRemoval~Tlrel,data=Cfb.db,family="quasibinomial")
Cfb2.glm<-glm(DungRemoval~tsrel,data=Cfb.db,family="quasibinomial")
Cfb2.glm<-glm(DungRemoval~Srel,data=Cfb.db,family="quasibinomial")

Cfb2.sumglm<-summary.glm(Cfb2.glm)
Cfb2.sumglm$deviance
Cfb2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
Cfb2.Rdev<-(Cfb2.sumglm$null.deviance-Cfb2.sumglm$deviance)/Cfb2.sumglm$null.deviance
Cfb2.Rdev


#Cfc
Cfc.vars<-c("D","Tl","ts","Rl","rs","S")
Cfc.env<-Cfc.db[Cfc.vars]
Cfc.hier<-hier.part(Cfc.db$DungRemoval,Cfc.env,family = "quasibinomial")
Cfc.hier
Cfc.hier$I.perc #values for graph

Cfc.prob<-rand.hp(Cfc.db$DungRemoval,Cfc.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
Cfc.prob$Iprobs #Z scores for graph

#Rdev model berekenen
Cfc.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=Cfc.db,family="quasibinomial")
Cfc.glm
Cfc.sumglm<-summary.glm(Cfc.glm)
Cfc.sumglm$deviance
Cfc.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
Cfc.Rdev<-(Cfc.sumglm$null.deviance-Cfc.sumglm$deviance)/Cfc.sumglm$null.deviance
Cfc.Rdev

#relative abundance
Cfc2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
Cfc2.env<-Cfc.db[Cfc2.vars]
Cfc2.hier<-hier.part(Cfc.db$DungRemoval,Cfc2.env,family = "quasibinomial")
Cfc2.hier
Cfc2.hier$I.perc #values for graph

Cfc2.prob<-rand.hp(Cfc.db$DungRemoval,Cfc2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
Cfc2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
Cfc2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Rlrel+rsrel+Srel,data=Cfc.db,family="quasibinomial")
Cfc2.glm
Cfc2.sumglm<-summary.glm(Cfc2.glm)
Cfc2.sumglm$deviance
Cfc2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
Cfc2.Rdev<-(Cfc2.sumglm$null.deviance-Cfc2.sumglm$deviance)/Cfc2.sumglm$null.deviance
Cfc2.Rdev


#Csb
Csb.vars<-c("D","Tl","ts","Rl","rs","S")
Csb.env<-Csb.db[Csb.vars]
Csb.hier<-hier.part(Csb.db$DungRemoval,Csb.env,family = "quasibinomial")
Csb.hier
Csb.hier$I.perc #values for graph

Csb.prob<-rand.hp(Csb.db$DungRemoval,Csb.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
Csb.prob$Iprobs #Z scores for graph

#Rdev model berekenen
Csb.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=Csb.db,family="quasibinomial")
Csb.glm
Csb.sumglm<-summary.glm(Csb.glm)
Csb.sumglm$deviance
Csb.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
Csb.Rdev<-(Csb.sumglm$null.deviance-Csb.sumglm$deviance)/Csb.sumglm$null.deviance
Csb.Rdev

#relative abundance
Csb2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
Csb2.env<-Csb.db[Csb2.vars]
Csb2.hier<-hier.part(Csb.db$DungRemoval,Csb2.env,family = "quasibinomial")
Csb2.hier
Csb2.hier$I.perc #values for graph

Csb2.prob<-rand.hp(Csb.db$DungRemoval,Csb2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
Csb2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
Csb2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Rlrel+rsrel+Srel,data=Csb.db,family="quasibinomial")
Csb2.glm<-glm(DungRemoval~Drel,data=Csb.db,family="quasibinomial")
Csb2.glm<-glm(DungRemoval~Tlrel,data=Csb.db,family="quasibinomial")
Csb2.glm<-glm(DungRemoval~tsrel,data=Csb.db,family="quasibinomial")
Csb2.glm<-glm(DungRemoval~rsrel,data=Csb.db,family="quasibinomial")
Csb2.glm<-glm(DungRemoval~Srel,data=Csb.db,family="quasibinomial")

Csb2.sumglm<-summary.glm(Csb2.glm)
Csb2.sumglm$deviance
Csb2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
Csb2.Rdev<-(Csb2.sumglm$null.deviance-Csb2.sumglm$deviance)/Csb2.sumglm$null.deviance
Csb2.Rdev


#Dfa
Dfa.db<-completeFun(Dfa.db, "D")
Dfa.db<-completeFun(Dfa.db, "DungRemovalCorr")
Dfa.vars<-c("D","Tl","ts","Rl","rs","S")
Dfa.env<-Dfa.db[Dfa.vars]
Dfa.hier<-hier.part(Dfa.db$DungRemoval,Dfa.env,family = "quasibinomial")
Dfa.hier
Dfa.hier$I.perc #values for graph

Dfa.prob<-rand.hp(Dfa.db$DungRemoval,Dfa.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
Dfa.prob$Iprobs #Z scores for graph

#Rdev model berekenen
Dfa.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=Dfa.db,family="quasibinomial")
Dfa.glm
Dfa.sumglm<-summary.glm(Dfa.glm)
Dfa.sumglm$deviance
Dfa.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
Dfa.Rdev<-(Dfa.sumglm$null.deviance-Dfa.sumglm$deviance)/Dfa.sumglm$null.deviance
Dfa.Rdev

#relative abundance
Dfa2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
Dfa2.env<-Dfa.db[Dfa2.vars]
Dfa2.hier<-hier.part(Dfa.db$DungRemoval,Dfa2.env,family = "quasibinomial")
Dfa2.hier
Dfa2.hier$I.perc #values for graph

Dfa2.prob<-rand.hp(Dfa.db$DungRemoval,Dfa2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
Dfa2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
Dfa2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Rlrel+rsrel+Srel,data=Dfa.db,family="quasibinomial")
Dfa2.glm<-glm(DungRemoval~Drel,data=Dfa.db,family="quasibinomial")
Dfa2.glm<-glm(DungRemoval~Tlrel,data=Dfa.db,family="quasibinomial")
Dfa2.glm<-glm(DungRemoval~tsrel,data=Dfa.db,family="quasibinomial")
Dfa2.glm<-glm(DungRemoval~Srel,data=Dfa.db,family="quasibinomial")
Dfa2.glm
Dfa2.sumglm<-summary.glm(Dfa2.glm)
Dfa2.sumglm$deviance
Dfa2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
Dfa2.Rdev<-(Dfa2.sumglm$null.deviance-Dfa2.sumglm$deviance)/Dfa2.sumglm$null.deviance
Dfa2.Rdev


#Dfb
Dfb.vars<-c("D","Tl","ts","Rl","rs","S")

summary(Dfb.db)
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

Dfb.db<-completeFun(Dfb.db, "D")
Dfb.db<-completeFun(Dfb.db, "DungRemovalCorr")

Dfb.env<-Dfb.db[Dfb.vars]

Dfb.hier<-hier.part(Dfb.db$DungRemoval,Dfb.env,family = "quasibinomial")
Dfb.hier
Dfb.hier$I.perc #values for graph

Dfb.prob<-rand.hp(Dfb.db$DungRemoval,Dfb.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
Dfb.prob$Iprobs #Z scores for graph

#Rdev model berekenen
Dfb.glm<-glm(DungRemoval~D+Tl+ts+Rl+rs+S,data=Dfb.db,family="quasibinomial")
Dfb.glm
Dfb.sumglm<-summary.glm(Dfb.glm)
Dfb.sumglm$deviance
Dfb.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
Dfb.Rdev<-(Dfb.sumglm$null.deviance-Dfb.sumglm$deviance)/Dfb.sumglm$null.deviance
Dfb.Rdev

#relative abundance
summary(Dfb.db)
Dfb2.db<-Dfb.db

Dfb2.db<-completeFun(Dfb2.db, "Drel")
Dfb2.db<-completeFun(Dfb2.db, "DungRemovalCorr")
summary(Dfb2.db)

Dfb2.vars<-c("Drel","Tlrel","tsrel","Rlrel","rsrel","Srel")
Dfb2.env<-Dfb2.db[Dfb2.vars]
Dfb2.hier<-hier.part(Dfb2.db$DungRemoval,Dfb2.env,family = "quasibinomial")
Dfb2.hier
Dfb2.hier$I.perc #values for graph

Dfb2.prob<-rand.hp(Dfb2.db$DungRemoval,Dfb2.env,family = "quasibinomial",gof = "Rsqu",num.reps=100)
Dfb2.prob$Iprobs #Z scores for graph

#Rdev model berekenen
Dfb2.glm<-glm(DungRemoval~Drel+Tlrel+tsrel+Rlrel+rsrel+Srel,data=Dfb2.db,family="quasibinomial")
Dfb2.glm<-glm(DungRemoval~Drel,data=Dfb2.db,family="quasibinomial")
Dfb2.glm<-glm(DungRemoval~Tlrel,data=Dfb2.db,family="quasibinomial")
Dfb2.glm<-glm(DungRemoval~tsrel,data=Dfb2.db,family="quasibinomial")
Dfb2.glm<-glm(DungRemoval~Srel,data=Dfb2.db,family="quasibinomial")

Dfb2.sumglm<-summary.glm(Dfb2.glm)
Dfb2.sumglm$deviance
Dfb2.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
Dfb2.Rdev<-(Dfb2.sumglm$null.deviance-Dfb2.sumglm$deviance)/Dfb2.sumglm$null.deviance
Dfb2.Rdev


#grafiek maken
BWk<-BWk2.hier$I.perc
BWk$FG<-c("D","Tl","ts","Rl","rs","S")
BWk$climate<-rep("BWk",6)
BWk
BWk.probs<-BWk2.prob$Iprobs
BWk.probs
BWk<-cbind(BWk,BWk.probs)
BWk

Cfb<-Cfb2.hier$I.perc
Cfb$FG<-c("D","Tl","ts","Rl","rs","S")
Cfb$climate<-rep("Cfb",6)
Cfb
Cfb.probs<-Cfb2.prob$Iprobs
Cfb.probs
Cfb<-cbind(Cfb,Cfb.probs)
Cfb

Cfc<-Cfc2.hier$I.perc
Cfc$FG<-c("D","Tl","ts","Rl","rs","S")
Cfc$climate<-rep("Cfc",6)
Cfc
Cfc.probs<-Cfc2.prob$Iprobs
Cfc.probs
Cfc<-cbind(Cfc,Cfc.probs)
Cfc

Csb<-Csb2.hier$I.perc
Csb$FG<-c("D","Tl","ts","Rl","rs","S")
Csb$climate<-rep("Csb",6)
Csb
Csb.probs<-Csb2.prob$Iprobs
Csb.probs
Csb<-cbind(Csb,Csb.probs)
Csb

Dfa<-Dfa2.hier$I.perc
Dfa$FG<-c("D","Tl","ts","Rl","rs","S")
Dfa$climate<-rep("Dfa",6)
Dfa
Dfa.probs<-Dfa2.prob$Iprobs
Dfa.probs
Dfa<-cbind(Dfa,Dfa.probs)
Dfa

Dfb<-Dfb2.hier$I.perc
Dfb$FG<-c("D","Tl","ts","Rl","rs","S")
Dfb$climate<-rep("Dfb",6)
Dfb
Dfb.probs<-Dfb2.prob$Iprobs
Dfb.probs
Dfb<-cbind(Dfb,Dfb.probs)
Dfb

hier.climate<-rbind(BWk,Cfb,Cfc,Csb,Dfa,Dfb)
hier.climate

hier.climate$sig95b<-ifelse(hier.climate$I==0,"",ifelse(hier.climate$sig95=="*","*",""))
hier.climate

library(ggplot2)

hier.climate$FG<-as.factor(hier.climate$FG)
levels(hier.climate$FG)
hier.climate$FG2<-ifelse(hier.climate$FG=="S",1,(ifelse(hier.climate$FG=="rs",2,(ifelse(hier.climate$FG=="Rl",3,(ifelse(hier.climate$FG=="ts",4,(ifelse(hier.climate$FG=="Tl",5,6)))))))))
hier.climate$FG3<-ifelse(hier.climate$FG=="S",6,(ifelse(hier.climate$FG=="rs",5,(ifelse(hier.climate$FG=="Rl",4,(ifelse(hier.climate$FG=="ts",3,(ifelse(hier.climate$FG=="Tl",2,1)))))))))
hier.climate$FG<-reorder(hier.climate$FG,hier.climate$FG3)
hier.climate<-subset(hier.climate,hier.climate$FG!="Rl")

hier.climate.plot2<-ggplot(hier.climate[order(hier.climate$FG2),],aes(x=climate,y=I,fill=FG))+
  geom_bar(position="stack",stat="identity")
hier.climate.plot2
hier.climate.plot2<-hier.climate.plot2+theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Independent effect (%)")))+
  xlab(expression(bold("Climate zone")))+
  theme(axis.ticks.x=element_blank())+
  theme(axis.text.x=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(hier.climate)/5), y = -4, label = c("BWk","Cfb","Cfc","Csb","Dfa","Dfb"), size = 5)
hier.climate.plot2
hier.climate.plot2<-hier.climate.plot2+
  scale_fill_manual(name=expression(bold("Functional group")),
                    values=c("royalblue4","red4","orange","yellowgreen","grey"),
                    labels=c("dwellers","tunnelers (large)","tunnelers (small)","rollers (small)","soil macro-invertebrates"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))
hier.climate.plot2

#semi-manueel sterretjes zetten
hier.climate #via excel nieuwe dataset maken
climate.hier.label<-read.delim("climate.hier.labels.txt",header=T)
head(climate.hier.label)
hier.climate.plot2<-hier.climate.plot2+
  geom_text(data=climate.hier.label,aes(x=climate,y=labelY,color=color,label=sig95b),size=8,show.legend = F)+
  scale_colour_manual(values=c("white","black"))
hier.climate.plot2

#wegschrijven, samengevoegd zonder legende bij biogeo?
hier.biogeo<-hier.plot2+theme(legend.position="none")
hier.biogeo

library(cowplot)
tiff(filename="hierpart.tiff",width=600,height=700,unit="px")
ggdraw() +
  draw_plot(hier.biogeo, 0, 0.5, 1, 0.5) +
  draw_plot(hier.climate.plot2, 0, 0, 1, 0.5) +
  draw_plot_label(c("a","b"),c(0,0),c(1,0.5),size=18)
dev.off()

#biogeo and climate in 1 graph

#dataset
biogeoclimate<-read.delim("biogeo.hier.labels.txt",header=T)
head(biogeoclimate)
biogeoclimate$FG4<-paste(biogeoclimate$FG3,biogeoclimate$FG)

plot3<-ggplot(biogeoclimate[order(biogeoclimate$FG2),],aes(x=biogeo2,y=I,fill=FG4))+
  geom_bar(position="stack",stat="identity")
plot3
plot3<-plot3+theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Independent effect (%)")))+
  xlab(expression(bold("Biogeographic province and climate zone")))+
  theme(axis.ticks.x=element_blank())+
  theme(axis.text.x=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(biogeoclimate)/5), y = -3.5, label = c("taig","bor","ceh","pan","brit","atl","mef","pont","med","ibh","cau"), size = 5)+
  annotate(geom="text",x=c(1,2,5,7,9,11),y=-10,label=c("Cfc","Dfb","Cfb","Dfa","Csb","Bwk"),color="gray50",fontface="italic", size=5.5)+
  annotate(geom='text',x=1:11,y=rep(104,11),label=c(0.167, 0.158, 0.062, "0.310", "0.100", 0.668, 0.587, 0.162, 0.159, 0.209, 0.387),size=4.5,fontface="italic")+
  geom_segment(x=0.6,xend=1.4,y=-7,yend=-7,color="gray75",size=1)+
  geom_segment(x=1.6,xend=4.4,y=-7,yend=-7,color="gray75",size=1)+
  geom_segment(x=4.6,xend=6.4,y=-7,yend=-7,color="gray75",size=1)+
  geom_segment(x=6.6,xend=8.4,y=-7,yend=-7,color="gray75",size=1)+
  geom_segment(x=8.6,xend=10.4,y=-7,yend=-7,color="gray75",size=1)+
  geom_segment(x=10.6,xend=11.4,y=-7,yend=-7,color="gray75",size=1)
plot3
plot3<-plot3+
  scale_fill_manual(name=expression(bold("Functional group")),
                    values=c("royalblue4","red4","orange","yellowgreen","grey"),
                    labels=c("dwellers","tunnelers (large)","tunnelers (small)","rollers (small)","soil macro- \ninvertebrates"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))
plot3

plot3<-plot3+
  geom_text(data=biogeoclimate,aes(x=biogeo2,y=labelY,color=color,label=sig95b),size=8,show.legend = F)+
  scale_colour_manual(values=c("white","black"))
plot3

tiff(filename="hierBiogeoclimate_v2.tiff",width=600,height=400,unit="px")
plot3
dev.off()


####relatie aantal functionele groepen en mestverwijdering per biogeografische regio####
#lme proberen
rm(list=ls())
library(lme4)
library(car)
library(multcomp)

setwd('D:/Tmilotic/Documents/Doctoraat/8 MSE Europe/R scripts')
db<-read.delim("DungRemFGno.txt",header=T)
head(db)

#function correlation between the fitted and the observed values
#http://stats.stackexchange.com/questions/95054/how-to-get-the-overall-effect-for-linear-mixed-model-in-lme4-in-r
r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

colnames(db)[names(db)== "Biogeographic.provinceUdvardy1975"] <- "biogeo"
colnames(db)[names(db)== "Koppen.climate.zonePeel2007"] <- "climate"

summary(db)
db$DungRemoval<-db$DungRemovalCorr/100
#db<-subset(db,Dung.type=="C")
#speciale treatments eruit
db<-subset(db,((TreatmentNo!="t12")&(TreatmentNo!="t13")&(TreatmentNo!="t14")&(TreatmentNo!="t15")&(TreatmentNo!="t16")))
db$TreatmentNo<-droplevels(db$TreatmentNo)
levels(db$TreatmentNo)

db$Study.site<-as.factor(db$Study.site)
levels(db$Study.site)

#linear models
#model met FG naargelang vangst en treatment
m1<-lmer(DungRemoval~FGcorr+biogeo+TreatmentNo+(1|Dung.type)+(1|Study.site),
         data=db)
summary(m1)
plot(m1)
Anova(m1)
r2.corr.mer(m1)

#of hier en daar lijkt het alsof er wel rollers voorkomen, maar niet gevangen zijn (waar theoretisch gezien rollers kunnen zitten rollers toegevoegd ook al zijn er geen gevangen)
m2<-lmer(DungRemoval~FGcorr2+biogeo+TreatmentNo+(1|Dung.type)+(1|Study.site),
         data=db)
summary(m2)
plot(m2)
Anova(m2)
r2.corr.mer(m2)

anova(m1,m2)

#toch maar met m1 verder werken, want niet zoveel verschil in fit én geen ondersteunende waarnemingen

summary(m1)
Anova(m1)
#model zonder onderscheid treatments
m1b<-lmer(DungRemoval~FGcorr+biogeo+(1|Dung.type)+(1|Study.site),
          data=db)
summary(m1b)
plot(m1b)
Anova(m1b)
r2.corr.mer(m1b)
anova(m1,m1b)

#m2 juister
m<-m2
summary(glht(m, mcp(biogeo="Tukey")))
summary(glht(m, mcp(TreatmentNo="Tukey")))


#grafiek maken
library(ggplot2)
source("HelperFunctions_summary.txt")

biogeolist<-levels(db$biogeo)
treatmentlist<-levels(db$TreatmentNo)
newdat<-data.frame(FGcorr2=rep((c(rep(0,11),rep(6,11))),11),biogeo=c(rep(biogeolist,22)),TreatmentNo=c(rep(treatmentlist,22)))
newdat$pred<-predict(m2,newdat,type="response",re.form=NA)

#means
avg<-summarySE(data=db,measurevar = "DungRemoval",groupvars=c("FGcorr2","biogeo","TreatmentNo"))
avg


plot<-ggplot()+
  geom_line(data=newdat,aes(x=FGcorr2,y=pred*100,colour=biogeo),size=1)+
  geom_pointrange(data=avg,aes(x=FGcorr2, y=DungRemoval*100,ymin=(DungRemoval-se)*100,ymax=(DungRemoval+se)*100,colour=biogeo),show.legend = F)+
  theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  theme(axis.line.x = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=12),axis.title.y=element_text(size=12))+
  ylab(expression(bold("Dung removal (%)")))+
  xlab(expression(bold("Functional group richness")))+
  scale_colour_manual(name=expression(bold("Biogeographic province")),
                      values=c("lightgoldenrod","thistle4","coral3","forestgreen","green3","goldenrod2","olivedrab3","aquamarine2","chocolate4","maroon1","steelblue4"))+
  theme(legend.title=element_text(size=12),
        legend.text=element_text(size=12))+
  coord_cartesian(ylim=c(0,100),xlim=c(0,5.2))+
  scale_x_continuous(expand = c(0.02, 0)) + scale_y_continuous(expand = c(0, 0))
plot

tiff(filename="FGnoRemovalGLM_v2.tiff",width=600,height=400,unit="px")
plot
dev.off()


#dung removal = ratio between 0 and 1 -> binomial distribution

#model met FG naargelang vangst en treatment
bm1<-glmer(DungRemoval~FGcorr2+biogeo+TreatmentNo+(1|Dung.type)+(1|Study.site),
           family=binomial,data=db)
summary(bm1)
plot(bm1)
Anova(bm1)
r2.corr.mer(bm1)

#better fit with linear model

library(lsmeans)
lsmeans(m1, pairwise ~ biogeo)
lsmeans(m1, pairwise ~ TreatmentNo)


####☺functional group number by climate zone####
mclimate1<-lmer(DungRemoval~FGcorr+climate+TreatmentNo+(1|Dung.type)+(FGcorr|Study.site),
         data=db)
summary(mclimate1)
plot(mclimate1)
Anova(mclimate1)
r2.corr.mer(mclimate1)

#of hier en daar lijkt het alsof er wel rollers voorkomen, maar niet gevangen zijn (waar theoretisch gezien rollers kunnen zitten rollers toegevoegd ook al zijn er geen gevangen)
mclimate2<-lmer(DungRemoval~FGcorr2+climate+TreatmentNo+(1|Dung.type)+(1|Study.site),
         data=db)
summary(mclimate2)
plot(mclimate2)
Anova(mclimate2)
r2.corr.mer(mclimate2)

anova(mclimate1,mclimate2)

#toch maar met mclimate1 verder werken, want niet zoveel verschil in fit én geen ondersteunende waarnemingen

summary(mclimate1)
Anova(mclimate1)
#model zonder onderscheid treatments
mclimate1b<-lmer(DungRemoval~FGcorr+climate+(1|Dung.type)+(1|Study.site),
          data=db)
summary(mclimate1b)
plot(mclimate1b)
Anova(mclimate1b)
r2.corr.mer(mclimate1b)
anova(mclimate1,mclimate1b)

climatelist<-levels(db$climate)
treatmentlist<-levels(db$TreatmentNo)
newdat<-data.frame(FGcorr=rep((c(rep(0,6),rep(6,6))),11),climate=c(rep(climatelist,22)),TreatmentNo=c(rep(treatmentlist,12)))
newdat$pred<-predict(mclimate1b,newdat,type="response",re.form=NA)


#means
avg<-summarySE(data=db,measurevar = "DungRemoval",groupvars=c("FGcorr","climate","TreatmentNo"))
avg


climate.plot<-ggplot()+
  geom_line(data=newdat,aes(x=FGcorr,y=pred*100,colour=climate),size=1)+
  geom_pointrange(data=avg,aes(x=FGcorr, y=DungRemoval*100,ymin=(DungRemoval-se)*100,ymax=(DungRemoval+se)*100,colour=climate),show.legend = F)+
  theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  theme(axis.line.x = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Dung removal (%)")))+
  xlab(expression(bold("Functional group number")))+
  scale_colour_manual(name=expression(bold("Köppen-Geiger climate zone")),
                      values=c("#FF9999","#66FF00","forestgreen","goldenrod3","#33FFFF","#33CCCC"))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=12))
climate.plot



####Dung removal vs relatieve abundantie van functionele groepen####
#mixed model niet relevant indien klimaatsvariabelen, relatieve abundantie en biogeo in model zitten (enkel relevant indien iedereen meerdere keren het experiment uitgevoerd zou hebben)
rm(list=ls())
library(car)
setwd('D:/Tmilotic/Documents/Doctoraat/8 MSE Europe/R scripts')
db<-read.delim("dungSeedRemovalRelFGenv2.txt",header=T)
head(db)

colnames(db)[names(db)== "Biogeographic.provinceUdvardy1975"] <- "biogeo"
colnames(db)[names(db)== "Koppen.climate.zonePeel2007"] <- "climate"

summary(db)
db$DungRemoval<-db$DungRemovalCorr/100
db$timing<-paste(db$Month_startExp,db$Year)
db$timing<-as.factor(db$timing)

#remove NA's
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

db<-completeFun(db, "Drel")
summary(db)
db$S<-db$Worms

dungrem.glm<-glm(DungRemoval~Drel+BTrel+STrel+SRrel+S+
                   biogeo+MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP+
                   Dung.type,data=db,family="quasibinomial")
dungrem.glm
Anova(dungrem.glm)
dungrem.sumglm<-summary.glm(dungrem.glm)
dungrem.sumglm$deviance
dungrem.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
dungrem.Rdev<-(dungrem.sumglm$null.deviance-dungrem.sumglm$deviance)/dungrem.sumglm$null.deviance
dungrem.Rdev

plot(dungrem.glm)

drop1(dungrem.glm)


#dung en biogeo genest want sommige soorten niet overal gebruikt
dungrem.glm2<-glm(DungRemoval~Drel+BTrel+STrel+SRrel+S+
                   Dung.type/biogeo+MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                  data=db,family="quasibinomial")
dungrem.glm2
Anova(dungrem.glm2)
dungrem.sumglm<-summary.glm(dungrem.glm2)
dungrem.sumglm$deviance
dungrem.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
dungrem.Rdev<-(dungrem.sumglm$null.deviance-dungrem.sumglm$deviance)/dungrem.sumglm$null.deviance
dungrem.Rdev

plot(dungrem.glm2)

drop1(dungrem.glm2)

#model3
dungrem.glm3<-glm(DungRemoval~Drel+BTrel+STrel+SRrel+S+biogeo+
                    Dung.type/biogeo+MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                  data=db,family="quasibinomial")
dungrem.glm3
summary(dungrem.glm3)
Anova(dungrem.glm3,test.statistic = "Wald")

dungrem.sumglm<-summary.glm(dungrem.glm3)
dungrem.sumglm$deviance
dungrem.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
dungrem.Rdev<-(dungrem.sumglm$null.deviance-dungrem.sumglm$deviance)/dungrem.sumglm$null.deviance
dungrem.Rdev

plot(dungrem.glm3)

drop1(dungrem.glm3)


sum(residuals(dungrem.glm3, type = "deviance")^2)

plot(residuals(dungrem.glm3))

library(multcomp)
summary(glht(dungrem.glm, mcp(Dung.type="Tukey")))
plot(glht(dungrem.glm, mcp(Dung.type="Tukey")))

summary(glht(dungrem.glm, mcp(biogeo="Tukey")))
plot(glht(dungrem.glm, mcp(biogeo="Tukey")))

#'zeldzame' mestsoorten enkel verschillend, anderen geen effect
#ofwel D en G verwijderen, ofwel dung type uit model verwijderen
db2<-subset(db,((Dung.type!="G")&(Dung.type!="D")))
#model 4 zonder D en G dung
dungrem.glm4<-glm(DungRemoval~Drel+BTrel+STrel+SRrel+S+biogeo+
                    Dung.type/biogeo+MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                  data=db2,family="quasibinomial")
dungrem.glm4
summary(dungrem.glm4)
Anova(dungrem.glm4,test.statistic = "Wald")

dungrem.sumglm<-summary.glm(dungrem.glm4)
dungrem.sumglm$deviance
dungrem.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
dungrem.Rdev<-(dungrem.sumglm$null.deviance-dungrem.sumglm$deviance)/dungrem.sumglm$null.deviance
dungrem.Rdev

plot(dungrem.glm4)

drop1(dungrem.glm4)

dungrem.glm4b<-glm(DungRemoval~Drel+BTrel+STrel+SRrel+S+biogeo+
                    Dung.type/SiteCode +MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                  data=db2,family="quasibinomial")
dungrem.glm4b
summary(dungrem.glm4b)
Anova(dungrem.glm4b,test.statistic = "Wald")

dungrem.sumglm<-summary.glm(dungrem.glm4b)
dungrem.sumglm$deviance
dungrem.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
dungrem.Rdev<-(dungrem.sumglm$null.deviance-dungrem.sumglm$deviance)/dungrem.sumglm$null.deviance
dungrem.Rdev

plot(dungrem.glm4b)

drop1(dungrem.glm4b)

dungrem.glm4c<-glm(DungRemoval~Drel+BTrel+STrel+SRrel+S+biogeo+biogeo/SiteCode+
                     Dung.type/SiteCode +MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                   data=db2,family="quasibinomial")
dungrem.glm4c
summary(dungrem.glm4c)
Anova(dungrem.glm4c,test.statistic = "Wald")

dungrem.sumglm<-summary.glm(dungrem.glm4c)
dungrem.sumglm$deviance
dungrem.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
dungrem.Rdev<-(dungrem.sumglm$null.deviance-dungrem.sumglm$deviance)/dungrem.sumglm$null.deviance
dungrem.Rdev

plot(dungrem.glm4c)

drop1(dungrem.glm4c)

anova(dungrem.glm4,dungrem.glm4b,dungrem.glm4c)

library(multcomp)
dungrem.glm1<-glm(DungRemoval~Drel+BTrel+STrel+SRrel+S+biogeo+
                    Dung.type+MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                  data=db2,family="quasibinomial")
summary(glht(dungrem.glm1, mcp(Dung.type="Tukey")))
plot(glht(dungrem.glm1, mcp(Dung.type="Tukey")))

summary(glht(dungrem.glm1, mcp(biogeo="Tukey")))
plot(glht(dungrem.glm1, mcp(biogeo="Tukey")))


#MODELKEUZE: dung genest in site, volledige dataset
dungrem.glm4d<-glm(DungRemoval~Drel+BTrel+STrel+SRrel+S+biogeo+biogeo/SiteCode+
                     Dung.type/SiteCode +MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                   data=db,family="quasibinomial")
dungrem.glm4d
summary(dungrem.glm4d)
Anova(dungrem.glm4d,test.statistic = "Wald")

dungrem.sumglm<-summary.glm(dungrem.glm4d)
dungrem.sumglm$deviance
dungrem.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
dungrem.Rdev<-(dungrem.sumglm$null.deviance-dungrem.sumglm$deviance)/dungrem.sumglm$null.deviance
dungrem.Rdev

plot(dungrem.glm4d)

drop1(dungrem.glm4d)


#model 5 zonder dung
dungrem.glm5<-glm(DungRemoval~Drel+BTrel+STrel+SRrel+S+biogeo+
                    MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                  data=db,family="quasibinomial")
dungrem.glm5
summary(dungrem.glm5)
Anova(dungrem.glm5,test.statistic = "Wald")

dungrem.sumglm<-summary.glm(dungrem.glm5)
dungrem.sumglm$deviance
dungrem.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
dungrem.Rdev<-(dungrem.sumglm$null.deviance-dungrem.sumglm$deviance)/dungrem.sumglm$null.deviance
dungrem.Rdev

plot(dungrem.glm5)

drop1(dungrem.glm5)

anova(dungrem.glm3,dungrem.glm5)
#slechter model


#keuze model 4d
dungrem.glm4d<-glm(DungRemoval~Drel+BTrel+STrel+SRrel+S+biogeo+biogeo/SiteCode+
                     Dung.type/SiteCode +MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                   data=db,family="quasibinomial")
dungrem.glm4d
summary(dungrem.glm4d)
Anova(dungrem.glm4d,test.statistic = "Wald")

dungrem.sumglm<-summary.glm(dungrem.glm4d)
dungrem.sumglm$deviance
dungrem.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
dungrem.Rdev<-(dungrem.sumglm$null.deviance-dungrem.sumglm$deviance)/dungrem.sumglm$null.deviance
dungrem.Rdev

drop1(dungrem.glm4d)

summary(glht(dungrem.glm, mcp(Dung.type="Tukey")))
plot(glht(dungrem.glm, mcp(Dung.type="Tukey")))
tapply(db$DungRemoval,db$Dung.type,mean)

summary(glht(dungrem.glm1, mcp(biogeo="Tukey")))
plot(glht(dungrem.glm1, mcp(biogeo="Tukey")))
tapply(db$DungRemoval,db$biogeo,mean)




#of zelfde model met abs abundantie?
db2<-completeFun(db2,"Dabs")
summary(db2)

dungrem.glm4e<-glm(DungRemoval~Dabs+BTabs+STabs+SRabs+S+biogeo+biogeo/SiteCode+
                     Dung.type/SiteCode +MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                   data=db2,family="quasibinomial")
dungrem.glm4e
summary(dungrem.glm4e)
Anova(dungrem.glm4e,test.statistic = "Wald")

plot(dungrem.glm4e)


dungrem.sumglm<-summary.glm(dungrem.glm4e)
dungrem.sumglm$deviance
dungrem.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
dungrem.Rdev<-(dungrem.sumglm$null.deviance-dungrem.sumglm$deviance)/dungrem.sumglm$null.deviance
dungrem.Rdev


#4f met relatieve abundanties MAAR totale mestkeverdichtheid
summary(db)
db3<-completeFun(db,"TotalAbund")
dungrem.glm4f<-glm(DungRemoval~TotalAbund+Drel+BTrel+STrel+SRrel+S+biogeo+biogeo/SiteCode+
                     Dung.type/SiteCode +MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                   data=db3,family="quasibinomial")
dungrem.glm4f
summary(dungrem.glm4f)
Anova(dungrem.glm4f,test.statistic = "Wald")

plot(dungrem.glm4f)

dungrem.sumglm<-summary.glm(dungrem.glm4f)
dungrem.sumglm$deviance
dungrem.sumglm$null.deviance

#Rdev² berekenen
#(null dev-res dev)/null dev
dungrem.Rdev<-(dungrem.sumglm$null.deviance-dungrem.sumglm$deviance)/dungrem.sumglm$null.deviance
dungrem.Rdev

#MAAR werkelijke abundantie helemaal andere grootte-orde dan soil faun (0/1) -> foute interpretatie
#+ niet alle aanwezige mestkevers komen naar elke plot -> relatieve abundantie veiliger!

m<-dungrem.glm4d
summary(m)

library(ggplot2)
source("HelperFunctions_summary.txt")
emaxT.pred<-plot_fit(m,focal_var = "EMXT")
summary(emaxT.pred)

emaxT.sum<-summarySE(data=db,measurevar = "DungRemoval",groupvars = c("EMXT","timing","Dung.type","TreatmentNo","biogeo"))
summary(emaxT.sum)

emaxT.plot<-ggplot()+geom_point(data=emaxT.sum,aes(x=EMXT,y=DungRemoval*100),colour="lightgrey")+
  geom_line(data=emaxT.pred,aes(x=EMXT,y=Pred*100),size=1)+
  geom_line(data=emaxT.pred,aes(x=EMXT,y=LC*100),lty=2,size=1)+
  geom_line(data=emaxT.pred,aes(x=EMXT,y=UC*100),lty=2,size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Dung removal (%)")))+
  xlab(expression(bold("emaxT (°C)")))+
  scale_y_continuous(expand = c(0, 0))
emaxT.plot

eminT.pred<-plot_fit(m,focal_var = "EMNT")
summary(eminT.pred)

eminT.sum<-summarySE(data=db,measurevar = "DungRemoval",groupvars = c("EMNT","timing","Dung.type","TreatmentNo","biogeo"))
summary(eminT.sum)

eminT.plot<-ggplot()+geom_point(data=eminT.sum,aes(x=EMNT,y=DungRemoval*100),colour="lightgrey")+
  geom_line(data=eminT.pred,aes(x=EMNT,y=Pred*100),size=1)+
  geom_line(data=eminT.pred,aes(x=EMNT,y=LC*100),lty=2,size=1)+
  geom_line(data=eminT.pred,aes(x=EMNT,y=UC*100),lty=2,size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Dung removal (%)")))+
  xlab(expression(bold("eminT (°C)")))+
  scale_y_continuous(expand = c(0, 0))
eminT.plot

mT.pred<-plot_fit(m,focal_var = "MNTM")
summary(mT.pred)
mT.sum<-summarySE(data=db,measurevar = "DungRemoval",groupvars = c("MNTM","timing","Dung.type","TreatmentNo","biogeo"))
summary(mT.sum)

mT.plot<-ggplot()+geom_point(data=mT.sum,aes(x=MNTM,y=DungRemoval*100),colour="lightgrey")+
  geom_line(data=mT.pred,aes(x=MNTM,y=Pred*100),size=1)+
  geom_line(data=mT.pred,aes(x=MNTM,y=LC*100),lty=2,size=1)+
  geom_line(data=mT.pred,aes(x=MNTM,y=UC*100),lty=2,size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Dung removal (%)")))+
  xlab(expression(bold("mT (°C)")))+
  scale_y_continuous(expand = c(0, 0))
mT.plot
 

mminT.pred<-plot_fit(m,focal_var = "MMNT")
summary(mminT.pred)
mminT.sum<-summarySE(data=db,measurevar = "DungRemoval",groupvars = c("MMNT","timing","Dung.type","TreatmentNo","biogeo"))
summary(mminT.sum)

mminT.plot<-ggplot()+geom_point(data=mminT.sum,aes(x=MMNT,y=DungRemoval*100),colour="lightgrey")+
  geom_line(data=mminT.pred,aes(x=MMNT,y=Pred*100),size=1)+
  geom_line(data=mminT.pred,aes(x=MMNT,y=LC*100),lty=2,size=1)+
  geom_line(data=mminT.pred,aes(x=MMNT,y=UC*100),lty=2,size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Dung removal (%)")))+
  xlab(expression(bold("mminT (°C)")))+
  scale_y_continuous(expand = c(0, 0))
mminT.plot

emaxT.pred<-plot_fit(m,focal_var = "EMXT")
summary(emaxT.pred)
emaxT.sum<-summarySE(data=db,measurevar = "DungRemoval",groupvars = c("EMXT","timing","Dung.type","TreatmentNo","biogeo"))
summary(emaxT.sum)

emaxT.plot<-ggplot()+geom_point(data=emaxT.sum,aes(x=EMXT,y=DungRemoval*100),colour="lightgrey")+
  geom_line(data=emaxT.pred,aes(x=EMXT,y=Pred*100),size=1)+
  geom_line(data=emaxT.pred,aes(x=EMXT,y=LC*100),lty=2,size=1)+
  geom_line(data=emaxT.pred,aes(x=EMXT,y=UC*100),lty=2,size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Dung removal (%)")))+
  xlab(expression(bold("emaxT (°C)")))+
  scale_y_continuous(expand = c(0, 0))
emaxT.plot

eP.pred<-plot_fit(m,focal_var = "EMXP")
summary(eP.pred)
eP.sum<-summarySE(data=db,measurevar = "DungRemoval",groupvars = c("EMXP","timing","Dung.type","TreatmentNo","biogeo"))
summary(eP.sum)

eP.plot<-ggplot()+geom_point(data=eP.sum,aes(x=EMXP,y=DungRemoval*100),colour="lightgrey")+
  geom_line(data=eP.pred,aes(x=EMXP,y=Pred*100),size=1)+
  geom_line(data=eP.pred,aes(x=EMXP,y=LC*100),lty=2,size=1)+
  geom_line(data=eP.pred,aes(x=EMXP,y=UC*100),lty=2,size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Dung removal (%)")))+
  xlab(expression(bold("eP (mm)")))+
  coord_cartesian(xlim=c(0,70))+
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))
eP.plot



D.pred<-plot_fit(m,focal_var = "Drel")
summary(D.pred)
dim(D.pred)
D.pred$FG<-rep("D",20)
colnames(D.pred)[1]<-"relAbund"

ST.pred<-plot_fit(m,focal_var = "STrel")
summary(ST.pred)
dim(ST.pred)
ST.pred$FG<-rep("ST",20)
colnames(ST.pred)[1]<-"relAbund"

SR.pred<-plot_fit(m,focal_var = "SRrel")
summary(SR.pred)
dim(SR.pred)
SR.pred$FG<-rep("SR",20)
colnames(SR.pred)[1]<-"relAbund"

FG.pred<-rbind(D.pred,ST.pred,SR.pred)

D.db<-data.frame(db$Drel,db$DungRemoval)
dim(db)
colnames(D.db)[1]<-"relAbund"
D.db$FG<-rep("D",2739)

ST.db<-data.frame(db$STrel,db$DungRemoval)
colnames(ST.db)[1]<-"relAbund"
ST.db$FG<-rep("ST",2739)

SR.db<-data.frame(db$SRrel,db$DungRemoval)
colnames(SR.db)[1]<-"relAbund"
SR.db$FG<-rep("SR",2739)

FG.db<-rbind(D.db,ST.db,ST.db)
FG.db$FG<-as.factor(FG.db$FG)

FG.plot<-ggplot()+geom_point(data=FG.db,aes(x=relAbund*100,y=db.DungRemoval*100,colour=FG))+
  geom_line(data=FG.pred,aes(x=relAbund*100,y=Pred*100,colour=FG),size=1)+
  geom_line(data=FG.pred,aes(x=relAbund*100,y=LC*100,colour=FG),lty=2,size=1)+
  geom_line(data=FG.pred,aes(x=relAbund*100,y=UC*100,colour=FG),lty=2,size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Dung removal (%)")))+
  xlab(expression(bold("Relative abundance (%)")))+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
FG.plot

library(multcomp)
dungrem.glm1<-glm(DungRemoval~Drel+BTrel+STrel+SRrel+S+biogeo+
                    Dung.type+MMXT+MMNT+MNTM+EMXT+EMNT+TotalPrecipitation+EMXP,
                  data=db,family="quasibinomial")
summary(glht(dungrem.glm1, mcp(Dung.type="Tukey")))
plot(glht(dungrem.glm1, mcp(Dung.type="Tukey")))

summary(glht(dungrem.glm1, mcp(biogeo="Tukey")))
plot(glht(dungrem.glm1, mcp(biogeo="Tukey")))

biogeo.sum<-summarySE(data=db,measurevar = "DungRemoval",groupvars = c("biogeo2"),na.rm = T)
biogeo.sum

biogeo.plot<-ggplot()+geom_bar(data=biogeo.sum,aes(x=biogeo2,y=DungRemoval*100),position="dodge",colour="black",stat="identity",fill="lightgrey")+
  geom_errorbar(data=biogeo.sum,aes(x=biogeo2,ymin=DungRemoval*100,ymax=(DungRemoval+se)*100),colour="black",position=position_dodge(.9),width=.5)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Dung removal (%)")))+
  xlab(expression(bold("Biogeographic province")))+
  scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim=c(0,90))+
  scale_x_discrete(labels=c("taig","bor","ceh","pan","brit","atl","mef","pont","med","ibh","cau"))+
  annotate(geom = "text",x=1:11,y=((biogeo.sum$DungRemoval+biogeo.sum$se)*100)+5,label=c("a","b","cd","ad","b","a","b","b","abcd","abcd","c"),size=4)
biogeo.plot



dung.sum<-summarySE(data=db,measurevar = "DungRemoval",groupvars = c("Dung.type"),na.rm = T)
dung.sum$dung<-c("cattle","red deer","goat","horse", "sheep")
dung.sum$dung<-as.factor(dung.sum$dung)
dung.sum
dung.sum$ymax<-(dung.sum$DungRemoval+dung.sum$se)*100

dung.plot<-ggplot()+geom_bar(data=dung.sum,aes(x=dung,y=DungRemoval*100),position="dodge",colour="black",stat="identity",fill="lightgrey")+
  geom_errorbar(data=dung.sum,aes(x=dung,ymin=DungRemoval*100,ymax=(DungRemoval+se)*100),colour="black",position=position_dodge(.9),width=.5)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Dung removal (%)")))+
  xlab(expression(bold("Dung type")))+
  scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim=c(0,45))+
  annotate(geom = "text",x=1:5,y=c(35,37,35,23,41),label=c("a","a","a","b","a"),size=4)
dung.plot

#combined graph
library(cowplot)
tiff(filename="DungRemGLM.tiff",width=800,height=800,unit="px")
ggdraw() +
  draw_plot(emaxT.plot, 0, 0.5, 0.5, 0.5) +
  draw_plot(eP.plot, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(dung.plot, 0, 0, 0.5, 0.5) +
  draw_plot(biogeo.plot, 0.5, 0, 0.5, 0.5) +
  draw_plot_label(c("a","b","c","d"),c(0,0.5,0,0.5),c(1,1,0.5,0.5),size=18)
dev.off()



####seed removal####

rm(list=ls())

library(ggplot2)

setwd('D:/Tmilotic/Documents/Doctoraat/8 MSE Europe/R scripts')
db<-read.delim("DungSeedSMA.txt",header=T)
head(db)


colnames(db)[names(db)== "Biogeographic.provinceUdvardy1975"] <- "biogeo"
colnames(db)[names(db)== "Koppen.climate.zonePeel2007"] <- "climate"
colnames(db)[names(db)== "T"] <- "Tl"
colnames(db)[names(db)== "R"] <- "Rl"
colnames(db)[names(db)== "Worms"] <- "S"

summary(db)

#in Iran other seed mimics were used
db<-subset(db, ((SiteCode!="IR_SR")&(SiteCode!="IR_TS")))
summary(db)
levels(db$Dung.type)
db$Dung.type<-factor(db$Dung.type)
db$SiteCode<-factor(db$SiteCode)

summary(db)
db$DungRemoval<-db$DungRemovalCorr/100
db$SeedRemoval<-db$SeedRem/100

#seed size as ordinal variable
db$Seed.size2<-db$Seed.size
db$Seed.size2<-revalue(db$Seed.size2, c("small"="1", "medium"="2","large"="3"))
db$Seed.size2<-as.numeric(db$Seed.size2)
summary(db)

#distributie bepalen
#1. Specifying fixed and random effects
head(db)
(reptab <- with(db, table(SiteCode, Dung.type)))
(reptab <- with(db, table(SiteCode, Seed.size)))


#2.Choose an error distribution; graphical checks#
## use within() to make multiple changes within a data frame
db <- within(db,
             {
               gna <- interaction(Dung.type,Seed.size,SiteCode)
               gna <- reorder(gna, SeedRemoval)
             })
#bwplot
print(bwplot(log(SeedRemoval + 1) ~ gna,
             data=db,
             scales=list(x=list(rot=90)) ## rotate x-axis labels
))

#We could also calculate the variance for each combination
#and provide a statistical summary of these variances.
grpVarL <- with(db, tapply(log(1+SeedRemoval),
                           list(gna), var) )
summary(grpVarL)

#We could calculate  for each genotype  treatment combination
#and provide a statistical summary of each group's .
grpMeans <- with(db, tapply(SeedRemoval, list(gna), mean))
summary(grpMeans)
#lambda<5 -> geen PQL

#poisson verdeling?
grpVars <- with(db,
                tapply(SeedRemoval,
                       list(gna), var) ) ## variance of UNTRANSFORMED data

lm1 <- lm(grpVars~grpMeans-1) ## `quasi-Poisson' fit
phi.fit <- coef(lm1)

lm2 <- lm(grpVars~I(grpMeans^2)+offset(grpMeans)-1) #negative binomial
k.fit <- 1/coef(lm2)

plot( grpVars ~ grpMeans, xlab="group means", ylab="group variances" )
abline(c(0,1), lty=2)
text(0.8,0.15,"Poisson")
curve(phi.fit*x, col=2,add=TRUE)
## bquote() is used to substitute numeric values
## in equations with symbols
text(0.8,0.1,
     bquote(paste("QP: ",sigma^2==.(round(phi.fit,1))*mu)),
     col=2)
curve(x*(1+x/k.fit),col=4,add=TRUE)
text(0.8,0.125,paste("NB: k=",round(k.fit,1),sep=""),col=4)
## loess fit
Lfit <- loess(grpVars~grpMeans)
mvec <- 0:2
lines(mvec,predict(Lfit,mvec),col=5)
text(0.8,0.075,"loess",col=5)

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  (rdf <- nrow(model@frame)-model.df)
  rp <- residuals(model)
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE,log.p=TRUE)
  c(chisq=Pearson.chisq,ratio=prat,p=exp(pval),logp=pval)
}

DungSeed.glm<-glm(log(SeedRemoval+1)~log(DungRemoval+1)*Dung.type+Seed.size,data=db,family=quasibinomial)
summary(DungSeed.glm)
plot(DungSeed.glm)
drop1(DungSeed.glm,test= "Chi")
#every term is highly significant

#different effects of study sites (searching effort for seeds varies for every researcher)
library(MASS)
DungSeedP.PQL <- glmmPQL(SeedRemoval~DungRemoval+Dung.type+Seed.size,
                         random = ~ 1 | SiteCode/Dung.type, family = poisson,
                         data=db)
summary(DungSeedP.PQL)

DungSeedQP.PQL <- glmmPQL(SeedRemoval~DungRemoval+Dung.type+Seed.size,
                          random = ~ 1 | SiteCode/Dung.type, family = quasipoisson,
                          data=db)
summary(DungSeedQP.PQL)

DungSeedB.PQL <- glmmPQL(SeedRemoval~DungRemoval+Dung.type+Seed.size,
                         random = ~ 1 | SiteCode/Dung.type, family = binomial,
                         data=db)
summary(DungSeedB.PQL)

DungSeedQB.PQL <- glmmPQL(SeedRemoval~DungRemoval+Dung.type+Seed.size,
                          random = ~ 1 | SiteCode/Dung.type, family = quasibinomial,
                          data=db)
summary(DungSeedQB.PQL)


library(lme4)
DungSeed.lmer <- glmer(SeedRemoval~DungRemoval+Dung.type+Seed.size+
                         (1 | Dung.type:SiteCode), family = binomial,
                       data=db)
summary(DungSeed.lmer)
drop1(DungSeed.lmer)

DungSeed.lmer1 <- glmer(SeedRemoval~DungRemoval+Seed.size+
                          (1 | Dung.type:SiteCode), family = binomial,
                        data=db)
summary(DungSeed.lmer1)

overdisp_fun(DungSeed.lmer) #p>0.05 -> OK, maar ratio laag?
overdisp_fun(DungSeed.lmer1) #eerder onderdispersie!

anova(DungSeed.lmer,DungSeed.lmer1)


DungSeed.lmer2 <- glmer(SeedRemoval~DungRemoval+Seed.size+
                          (1 | Dung.type/SiteCode), family = binomial,
                        data=db)
summary(DungSeed.lmer2)

anova(DungSeed.lmer1,DungSeed.lmer2)

DungSeed.lmer3 <- lmer(SeedRemoval~DungRemoval+Seed.size+
                         (1 | SiteCode), 
                       data=db)
summary(DungSeed.lmer3)
anova(DungSeed.lmer1,DungSeed.lmer3)
plot(DungSeed.lmer3)

overdisp_fun(DungSeed.lmer3) #underdispersed?

anova(DungSeed.lmer3,test="Chisq")


DungSeed.lmer4 <- glmer(SeedRemoval~DungRemoval+Dung.type+Seed.size+
                          biogeo+
                          (1 | Dung.type:SiteCode)+
                          (1|biogeo:SiteCode), family = binomial,
                        data=db)
summary(DungSeed.lmer4)

drop1(DungSeed.lmer4) #niets weglaten

DungSeed.lmer4b <- glmer(SeedRemoval~DungRemoval+Dung.type+Seed.size+
                           Dung.type:Seed.size+biogeo+
                           (1 | Dung.type:SiteCode)+
                           (1|biogeo:SiteCode), family = binomial,
                         data=db)
summary(DungSeed.lmer4b)

anova(DungSeed.lmer4,DungSeed.lmer4b)
overdisp_fun(DungSeed.lmer4b)

#4b lijkt beter


DungSeed.lmer5 <- glmer(SeedRemoval~DungRemoval+Dung.type+Seed.size+
                          biogeo+climate3+
                          (1 | Dung.type:SiteCode)+
                          (1|biogeo:SiteCode)+
                          (1|climate3:SiteCode), family = binomial,
                        data=db)
summary(DungSeed.lmer5)

anova(DungSeed.lmer,DungSeed.lmer1,DungSeed.lmer2,DungSeed.lmer3,DungSeed.lmer4,DungSeed.lmer4b,DungSeed.lmer5)

overdisp_fun(DungSeed.lmer4) #eerder onderdispersed

#model 4, 5 of 1 (voor de eenvoud) kiezen?

plot(fitted(DungSeed.lmer1), residuals(DungSeed.lmer1), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(DungSeed.lmer1), residuals(DungSeed.lmer1)))

plot(fitted(DungSeed.lmer4), residuals(DungSeed.lmer4), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(DungSeed.lmer4), residuals(DungSeed.lmer4)))

plot(fitted(DungSeed.lmer4b), residuals(DungSeed.lmer4b), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(DungSeed.lmer4b), residuals(DungSeed.lmer4b)))

plot(fitted(DungSeed.lmer5), residuals(DungSeed.lmer5), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(DungSeed.lmer5), residuals(DungSeed.lmer5)))

summary(DungSeed.lmer1)
summary(DungSeed.lmer4)
summary(DungSeed.lmer4b)
summary(DungSeed.lmer5)

anova(DungSeed.lmer4,DungSeed.lmer5,DungSeed.lmer1)

#AIC model4b best
Anova(DungSeed.lmer4b)
Anova(DungSeed.lmer4)
Anova(DungSeed.lmer5)
Anova(DungSeed.lmer1)

#model 4b lijkt best, wat verschilt?
m<-DungSeed.lmer4b

#function correlation between the fitted and the observed values
#http://stats.stackexchange.com/questions/95054/how-to-get-the-overall-effect-for-linear-mixed-model-in-lme4-in-r
r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

r2.corr.mer(m)

library(lsmeans)

lsm <- lsmeans(m, ~ Dung.type*Seed.size+biogeo)
summary(lsm, type = "response")
#To see the results graphically
plot(lsm, by = "Seed.size", intervals = TRUE, type = "response")
plot(lsm, by = "biogeo", intervals = TRUE, type = "response")

lsmip(lsm, Dung.type ~ Seed.size, type = "response")

summary(pairs(regrid(lsm)), type = "response")

summary(glht(m, mcp(Dung.type="Tukey")))
summary(glht(m, mcp(Seed.size="Tukey")))
summary(glht(m, mcp(biogeo="Tukey")))

#graphs
summary(m)
Anova(m)
drop1(m)

m2 <- glmer(SeedRemoval~DungRemoval+Dung.type+Seed.size+
                           Dung.type:Seed.size+
                           (1 | Dung.type:SiteCode), family = binomial,
                         data=db)
summary(m2)
r2.corr.mer(m2)
Anova(m2)

summary(glht(m2,linfct=mcp(Seed.size="Tukey")))

#cattle
dbC<-subset(db,Dung.type=="C")
m2C <- glmer(SeedRemoval~DungRemoval+Seed.size+
              (1 | SiteCode), family = binomial,
            data=dbC)
summary(glht(m2C,linfct=mcp(Seed.size="Tukey")))

#horse
dbH<-subset(db,Dung.type=="H")
m2H <- glmer(SeedRemoval~DungRemoval+Seed.size+
               (1 | SiteCode), family = binomial,
             data=dbH)
summary(glht(m2H,linfct=mcp(Seed.size="Tukey")))

#sheep
dbS<-subset(db,Dung.type=="S")
m2S <- glmer(SeedRemoval~DungRemoval+Seed.size+
               (1 | SiteCode), family = binomial,
             data=dbS)
summary(glht(m2S,linfct=mcp(Seed.size="Tukey")))

library(ggplot2)
source("HelperFunctions_summary.txt")



library("effects")
# obtain a fit at different estimates of the predictor
ef.1=effect(c("DungRemoval","Seed.size"),m2)
Seed.size.pred=data.frame(ef.1)
Seed.size.pred


Seed.size.sum<-summarySE(data=db,measurevar = "DungRemoval",groupvars = c("SeedRemoval","Seed.size","TreatmentNo"))
summary(Seed.size.sum)

Seed.size.plot<-ggplot()+
  geom_line(data=Seed.size.pred,aes(x=DungRemoval*100,y=fit*100,colour=Seed.size),size=1)+
  geom_line(data=Seed.size.pred,aes(x=DungRemoval*100,y=lower*100,colour=Seed.size),lty=2,size=1)+
  geom_line(data=Seed.size.pred,aes(x=DungRemoval*100,y=upper*100,colour=Seed.size),lty=2,size=1)+
  geom_point(data=Seed.size.sum,aes(x=DungRemoval*100,y=SeedRemoval*100,colour=Seed.size),pch=16,alpha=0.50,show.legend=F)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Seed removal (%)")))+
  xlab(expression(bold("Dung removal (%)")))+
  scale_colour_manual(name=expression(bold("Seed size")),
                    breaks=c("small","medium","large"),
                    values=c("#330033","#FF3300","#009900"))+
  coord_cartesian(xlim=c(0,105),ylim=c(0,105))+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position=c(0.8,0.25))
Seed.size.plot

#○dungtype:seedsize plot

dungseed.sum<-summarySE(data=db,measurevar="SeedRemoval",groupvars = c("Dung.type","Seed.size"))
dungseed.sum
dungseed.sum$label<-c("b","b","a","g","f","f","z","y","x")
dungseed.sum$y=dungseed.sum$SeedRemoval*100

dungseed.sum$order<-ifelse(dungseed.sum$Seed.size=="small",1,(ifelse(dungseed.sum$Seed.size=="medium",2,3)))
dungseed.sum$Seed.size<-reorder(dungseed.sum$Seed.size,dungseed.sum$order)

dungseed.plot<-ggplot(data=dungseed.sum,aes(x=Dung.type,y=y,label=label,fill=Seed.size,colour=Seed.size))+
  geom_bar(position="dodge",stat="identity")+
  geom_errorbar(aes(ymin=SeedRemoval*100,ymax=(SeedRemoval+se)*100,colour=Seed.size),position=position_dodge(.9),width=.5)+
  geom_text(aes(y = y * 1.1), position = position_dodge(0.9), vjust = 0,size=4)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Seed removal (%)")))+
  xlab(expression(bold("Dung type")))+
  coord_cartesian(ylim=c(0,100))+
  scale_colour_manual(name=expression(bold("Seed size")),
                      breaks=c("small","medium","large"),
                      values=c("#009900","#FF3300","#330033"))+
  scale_fill_manual(name=expression(bold("Seed size")),
                      breaks=c("small","medium","large"),
                    values=c("#009900","#FF3300","#330033"))+
  scale_x_discrete(breaks=c("C", "H", "S"),
                   labels=c("cattle", "horse", "sheep"))+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position="none")

dungseed.plot

tiff(filename="SeedRemGLM_v2.tiff",width=600,height=300,unit="px")
ggdraw() +
  draw_plot(Seed.size.plot, 0, 0, 0.5, 1) +
  draw_plot(dungseed.plot, 0.5, 0, 0.5, 1) +
  draw_plot_label(c("a","b"),c(0,0.5),c(1,1),size=18)
dev.off()


#simplified graph without dots
Seed.size.plot2<-ggplot()+
  geom_ribbon(data=Seed.size.pred, aes(x=DungRemoval*100,ymin=lower*100,ymax=upper*100,fill=Seed.size), alpha=0.5,show.legend=F)+
  geom_line(data=Seed.size.pred,aes(x=DungRemoval*100,y=fit*100,colour=Seed.size),size=1)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(legend.text=element_text(size=14),legend.title=element_text(size=14))+
  ylab(expression(bold("Seed removal (%)")))+
  xlab(expression(bold("Dung removal (%)")))+
  scale_colour_manual(name=expression(bold("Seed size")),
                      breaks=c("small","medium","large"),
                      values=c("#330033","#FF3300","#009900"))+
  scale_fill_manual(name=expression(bold("Seed size")),
                      breaks=c("small","medium","large"),
                      values=c("#330033","#FF3300","#009900"))+
  coord_cartesian(xlim=c(0,105),ylim=c(0,105))+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position=c(0.8,0.25))
Seed.size.plot2

tiff(filename="Seed.size.plot2.tiff",width=300,height=300,unit="px")
Seed.size.plot2
dev.off()


####Habitat specialization vs lat and lon####
setwd('D:/Tmilotic/Documents/Doctoraat/8 MSE Europe/R scripts')
dbHS<-read.delim("HabitatDiffSites.txt",header=T)
head(dbHS)
dbHS$specialist100<-dbHS$specialist*100

lat <- lm(latitude ~ specialist100, data = dbHS)
summary(lat)

lon <- lm(longitude ~ specialist100, data = dbHS)
summary(lon)


#graphs
library(ggplot2)
lat.plot<-ggplot(data=dbHS,aes(x=latitude,y=specialist100))+geom_point(colour="grey")+
  geom_smooth(method=lm, se=F,color="black")+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Specialist abundance (%)")))+
  xlab(expression(bold("Latitude")))+
  scale_y_continuous(expand = c(0, 0))+
  annotate(geom = "text",x=38,y=5,label="R²= 0.323, p= 0.017",size=4)
lat.plot

lon.plot<-ggplot(data=dbHS,aes(x=longitude,y=specialist100))+geom_point(colour="grey")+
  geom_smooth(method=lm, se=F,color="black")+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = .5),axis.line.y = element_line(color="black", size = .5))+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Specialist abundance (%)")))+
  xlab(expression(bold("Longitude")))+
  scale_y_continuous(expand = c(0, 0))+
  annotate(geom = "text",x=7.5,y=5,label="R²= 0.265, p= 0.034",size=4)
lon.plot

#♠graphs with abundance and richness per region
richn<-read.delim("biogeoRichn.txt",header=T)
head(richn)

abund.plot<-ggplot()+geom_bar(data=richn,aes(x=biogeo2,y=total),position="dodge",colour="black",stat="identity",fill="lightgrey")+
  theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Abundance (specimens/trap)")))+
  xlab(expression(bold("Biogeographic province and climate zone")))+
  theme(axis.ticks.x=element_blank())+
  theme(axis.text.x=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(richn)), y = -10.5, label = c("taig","bor","ceh","pan","brit","atl","mef","pont","med","ibh","cau"), size = 5)+
  annotate(geom="text",x=c(1,2,5,7,9,11),y=-30,label=c("Cfc","Dfb","Cfb","Dfa","Csb","Bwk"),color="gray50",fontface="italic", size=5.5)+
  geom_segment(aes(x=0.6,xend=1.4,y=-21,yend=-21),color="gray75",size=1)+
  geom_segment(aes(x=1.6,xend=4.4,y=-21,yend=-21),color="gray75",size=1)+
  geom_segment(aes(x=4.6,xend=6.4,y=-21,yend=-21),color="gray75",size=1)+
  geom_segment(aes(x=6.6,xend=8.4,y=-21,yend=-21),color="gray75",size=1)+
  geom_segment(aes(x=8.6,xend=10.4,y=-21,yend=-21),color="gray75",size=1)+
  geom_segment(aes(x=10.6,xend=11.4,y=-21,yend=-21),color="gray75",size=1)
abund.plot


richn.plot<-ggplot()+geom_bar(data=richn,aes(x=biogeo2,y=richness),position="dodge",colour="black",stat="identity",fill="lightgrey")+
  theme_classic()+
  theme(axis.line.y = element_line(color="black", size = .5))+
  geom_hline(yintercept=0,size=.5)+
  theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  ylab(expression(bold("Species richness")))+
  xlab(expression(bold("Biogeographic province and climate zone")))+
  theme(axis.ticks.x=element_blank())+
  theme(axis.text.x=element_blank())+
  annotate(geom = "text", x = seq_len(nrow(richn)), y = -1.575, label = c("taig","bor","ceh","pan","brit","atl","mef","pont","med","ibh","cau"), size = 5)+
  annotate(geom="text",x=c(1,2,5,7,9,11),y=-4.5,label=c("Cfc","Dfb","Cfb","Dfa","Csb","Bwk"),color="gray50",fontface="italic", size=5.5)+
  geom_segment(aes(x=0.6,xend=1.4,y=-3.15,yend=-3.15),color="gray75",size=1)+
  geom_segment(aes(x=1.6,xend=4.4,y=-3.15,yend=-3.15),color="gray75",size=1)+
  geom_segment(aes(x=4.6,xend=6.4,y=-3.15,yend=-3.15),color="gray75",size=1)+
  geom_segment(aes(x=6.6,xend=8.4,y=-3.15,yend=-3.15),color="gray75",size=1)+
  geom_segment(aes(x=8.6,xend=10.4,y=-3.15,yend=-3.15),color="gray75",size=1)+
  geom_segment(aes(x=10.6,xend=11.4,y=-3.15,yend=-3.15),color="gray75",size=1)
richn.plot

library(cowplot)
tiff(filename="FigSpRich.tiff",width=1000,height=800,unit="px")
ggdraw() +
  draw_plot(lat.plot, 0, 0.5, 0.4, 0.5) +
  draw_plot(abund.plot, 0.4, 0.5, 0.6, 0.5) +
  draw_plot(lon.plot, 0, 0, 0.4, 0.5) +
  draw_plot(richn.plot, 0.4, 0, 0.6, 0.5) +
  draw_plot_label(c("a","c","b","d"),c(0,0.4,0,0.4),c(1,1,0.5,0.5),size=18)
dev.off()

#separately
tiff(filename="LatLon.tiff",width=600,height=300,unit="px")
ggdraw() +
  draw_plot(lat.plot, 0, 0, 0.5, 1) +
  draw_plot(lon.plot, 0.5, 0, 0.5, 1) +
  draw_plot_label(c("a","b"),c(0,0.5),c(1,1),size=18)
dev.off()

tiff(filename="RichAbund.tiff",width=600,height=800,unit="px")
ggdraw() +
  draw_plot(abund.plot, 0, 0, 1, 0.5) +
  draw_plot(richn.plot, 0.0125, 0.5, 1, 0.5) +
  draw_plot_label(c("a","b"),c(0,0),c(1,0.5),size=18)
dev.off()


