####analysis of plant data, Quaternary Sept. 2021####
####load libraries####
library(tidyverse) # data handling and plotting
library(vegan) # numerical ecology
library(ca) # correspondence analysis
library(factoextra) # correspondence analysis
library(ggrepel) # plotting
library(ggpubr) # plotting
library(cowplot)
####references####
citation("tidyverse")
citation("vegan")
citation("ca")
####load data####
data.quat <- read.csv("data_v4.csv", header=TRUE)
data.quat<-as_tibble(data.quat)
data.quat<-mutate(data.quat,SampleID=row_number()) # add id number to each entry
data.quat$SampleID<-as.character(data.quat$SampleID)
data.quat[, 9:36][is.na(data.quat[, 9:36])] <- 0 # replace NA by zero in species columns
data.quat$culture<-factor(data.quat$culture,levels=c("SKC","vinča","butmir","sopot","LBK"),labels=c("SKC","Vinča","Butmir","Sopot","LBK"))
write.csv(data.quat,"/home/marc/Dropbox/m_pers/quaternary20/submission/data_map_01.csv") # for map all    
####species richness and sample numbers####
data.quat$richness<-as.numeric(rowSums(data.quat[,c(9:36)]))
data.quat<-filter(data.quat,richness!=1)
write.csv(data.quat,"/home/marc/Dropbox/m_pers/quaternary20/submission/data_map_02.csv") # for map all                  
data.quat$logrichness<-as.numeric(log(data.quat$richness))
data.quat$richnessagro<-as.numeric(rowSums(data.quat[,c(9:22)]))
data.quat$logrichnessagro<-as.numeric(log(data.quat$richnessagro))
data.quat$richnesswild<-as.numeric(rowSums(data.quat[,c(23:36)]))
data.quat$logrichnesswild<-as.numeric(log(data.quat$richnesswild))
data.quat$logrichnesswild<-data.quat$logrichnesswild%>%na_if(-Inf)
data.quat$lognumber<-as.numeric(log(data.quat$sample.number))
data.quat$logvolume<-as.numeric(log(data.quat$sample.volume))
data.quat$sample.cat<-as.factor(ifelse(data.quat$sample.number<=5,1,
                                       ifelse(data.quat$sample.number>5 & data.quat$sample.number<=20,2,
                                              ifelse(data.quat$sample.number>20 & data.quat$sample.number<=50,3,
                                                     ifelse(data.quat$sample.number>50 & data.quat$sample.number<=100,4,5)))))
####diversity analysis using Simpson's index of diversity (1-D)####
data.quat$simpson<-diversity(data.quat[,c(9:36)],index="simpson") # all assemblages
data.quat$simpsonagro<-diversity(data.quat[,c(9:22)],index="simpson") # crop/pulses only
data.quat$simpsonwild<-diversity(data.quat[,c(23:36)],index="simpson") #for wild only
####correspondence analysis####
#ca for entire assemblages but excluding species where species richness <3
data.quat.ca<-select(data.quat,c(9:15,17:20,22:24,27:36))
ca.quat <- ca(data.quat.ca)
#ca for crops only but excluding species where species richness <3
data.quat$richnesscropca<-as.numeric(rowSums(data.quat[,c(9:15,17:20)]))
data.quat.ca.crops<-filter(data.quat,richnesscropca!="0")
data.quat.ca.crops.2<-select(data.quat.ca.crops,c(9:15,17:20))
ca.quat.crops <- ca(data.quat.ca.crops.2)
#ca for wild only (species richness >3, removing sites where species richness>1 and acorn as appears as total outlier)
data.quat$richnesswildca<-as.numeric(rowSums(data.quat[,c(23:24,27:35)]))
data.quat.ca.wild<-filter(data.quat,richnesswildca!="0")
data.quat.ca.wild.2<-select(data.quat.ca.wild,c(23:24,27:35))
ca.quat.wild<-ca(data.quat.ca.wild.2)
#extract necessary ca data and create specific datasets for species scores + attach all site data to original dataset
#extract percentage of variance and round-up to two decimals to use for automatic legend of CA plot
eigenvalues <- get_eigenvalue(ca.quat)
var.percent<-data.frame(eigenvalues$variance.percent)
Dim1.raw<-var.percent[1,1]
Dim1<-round(Dim1.raw, digits=2)
Dim2.raw<-var.percent[2,1]
Dim2<-round(Dim2.raw, digits=2)
Dim3.raw<-var.percent[3,1]
Dim3<-round(Dim3.raw, digits=2)
#get column variables
col<-get_ca_col(ca.quat)
row<-get_ca_row(ca.quat)
# create data.frame combining CA results and necessary information for plotting (e.g. culture, bioregion)
x1<-row$coord[,1]
data.quat$dim1<-x1
x2<-row$coord[,2]
data.quat$dim2<-x2
x3<-row$coord[,3]
y1<-col$coord[,1]
y2<-col$coord[,2]
y3<-col$coord[,3]
colnames <- data.frame(names(data.quat.ca))
colscores <- data.frame(cbind(y1,y2,y3))
species.scores <- cbind(colnames,colscores)
names(species.scores) <- c("species", "y1", "y2","y3")
#crops/pulses assemblages
#extract percentage of variance and round-up to two decimals to use for automatic legend of CA plot
eigenvalues.crops <- get_eigenvalue(ca.quat.crops)
var.percent.crops<-data.frame(eigenvalues.crops$variance.percent)
Dim1.raw.crops<-var.percent.crops[1,1]
Dim1.crops<-round(Dim1.raw.crops, digits=2)
Dim2.raw.crops<-var.percent.crops[2,1]
Dim2.crops<-round(Dim2.raw.crops, digits=2)
Dim3.raw.crops<-var.percent.crops[3,1]
Dim3.crops<-round(Dim3.raw.crops, digits=2)
#get column variables
col.crops<-get_ca_col(ca.quat.crops)
row.crops<-get_ca_row(ca.quat.crops)
# create data.frame combining CA results and necessary information for plotting (e.g. culture, bioregion)
x1.crops<-row.crops$coord[,1]
x2.crops<-row.crops$coord[,2]
data.quat$dim1.crops<-x1.crops
data.quat$dim2.crops<-x2.crops
x3.crops<-row.crops$coord[,3]
y1.crops<-col.crops$coord[,1]
y2.crops<-col.crops$coord[,2]
y3.crops<-col.crops$coord[,3]
colnames.crops <- data.frame(names(data.quat.ca.crops.2))
colscores.crops <- data.frame(cbind(y1.crops,y2.crops,y3.crops))
species.scores.crops <- cbind(colnames.crops,colscores.crops)
names(species.scores.crops) <- c("species", "y1", "y2","y3")
#wild assemblages
#extract percentage of variance and round-up to two decimals to use for automatic legend of CA plot
eigenvalues.wild <- get_eigenvalue(ca.quat.wild)
var.percent.wild<-data.frame(eigenvalues.wild$variance.percent)
Dim1.raw.wild<-var.percent.wild[1,1]
Dim1.wild<-round(Dim1.raw.wild, digits=2)
Dim2.raw.wild<-var.percent.wild[2,1]
Dim2.wild<-round(Dim2.raw.wild, digits=2)
#get column variables
col.wild<-get_ca_col(ca.quat.wild)
row.wild<-get_ca_row(ca.quat.wild)
# create data.frame combining CA results and necessary information for plotting (e.g. culture, bioregion)
x1.wild<-row.wild$coord[,1]
x2.wild<-row.wild$coord[,2]
y1.wild<-col.wild$coord[,1]
y2.wild<-col.wild$coord[,2]
idwild<-data.frame(data.quat.ca.wild$SampleID)
dim1wild<-data.frame(x1.wild)
dim2wild<-data.frame(x2.wild)
wild.ca<-data.frame(cbind(idwild,dim1wild,dim2wild))
names(wild.ca) <- c("SampleID", "dim1.wild", "dim2.wild")
data.quat<-full_join(data.quat,wild.ca,by="SampleID")
#write.csv(data.quat.2,"/home/marc/Dropbox/m_pers/quaternary20/20210211/data_feb21.csv")
colnames.wild <- data.frame(names(data.quat.ca.wild.2))
colscores.wild <- data.frame(cbind(y1.wild,y2.wild))
species.scores.wild <- cbind(colnames.wild,colscores.wild)
names(species.scores.wild) <- c("species", "y1", "y2")


####data analysis and figures####
####data robustness####
#richness vs. number / volume of samples (crops)
rich.sample.agro<-ggplot(data.quat,aes(x=lognumber,y=logrichnessagro))+
  geom_point()+
  stat_smooth(method=lm)+
  labs(x = "Log number of samples", y="Log species richness (crops)")+
  theme_bw()
rich.volume.agro<-ggplot(data.quat,aes(x=logvolume,y=logrichnessagro))+
  geom_point()+
  stat_smooth(method=lm)+
  labs(x = "Log volume of samples", y="Log species richness (crops)")+
  theme_bw()
plot_grid(rich.sample.agro,rich.volume.agro,ncol=2) #figure 2
#correlation and regression analysis (number of samples vs. crops richness)
shapiro.test(data.quat$lognumber)
shapiro.test(data.quat$logrichnessagro)
cor.test(data.quat$logrichnessagro, data.quat$lognumber, method='spearman')
lm.rich.sample.agro<-lm(data.quat$logrichnessagro ~ data.quat$lognumber)
summary(lm.rich.sample.agro)
#correlation and regression analysis (volume of samples vs. crop richness)
shapiro.test(data.quat$logvolume)
cor.test(data.quat$logrichnessagro, data.quat$logvolume, method='spearman')
lm.rich.volume.agro<-lm(data.quat$logrichnessagro ~ data.quat$logvolume)
summary(lm.rich.volume.agro)
# richness vs. number/volume of samples (wild)
rich.sample.wild<-ggplot(data.quat,aes(x=lognumber,y=logrichnesswild))+
  geom_point()+
  stat_smooth(method=lm)+
  labs(x = "Log number of samples", y="Log species richness (wild)")+
  theme_bw()
rich.volume.wild<-ggplot(data.quat,aes(x=logvolume,y=logrichnesswild))+
  geom_point()+
  stat_smooth(method=lm)+
  labs(x = "Log volume of samples", y="Log species richness (wild)")+
  theme_bw()
plot_grid(rich.sample.wild,rich.volume.wild,ncol=2) # figure 3
#correlation and regression analysis (number of samples vs. wild richness)
shapiro.test(data.quat$lognumber)
shapiro.test(data.quat$logrichnesswild)
cor.test(data.quat$logrichnesswild, data.quat$lognumber, method='spearman')
lm.rich.sample.wild<-lm(data.quat$logrichnesswild ~ data.quat$lognumber)
summary(lm.rich.sample.wild)
#correlation and regression analysis (volume of samples vs. wild richness)
dq<-filter(data.quat,!is.na(logvolume))
dq<-filter(dq,!is.na(logrichnesswild))
cor.test(dq$logrichnesswild, dq$logvolume, method='spearman')
lm.rich.volume.wild<-lm(dq$logrichnesswild ~ dq$logvolume)
summary(lm.rich.volume.wild)

#histogram number of sample per culture (figure 4)
sample.hist<-ggplot(filter(data.quat,!is.na(sample.cat)),aes(x=sample.cat))+
  geom_bar()+
  scale_x_discrete(name="Number of samples", breaks=c(1,2,3,4,5),labels=c("1-5","6-20","21-50","51-100","100+"))+
  scale_y_continuous(limits=c(0,10),breaks=seq(0,10,2))+
  ylab("Number of sites")+
  facet_wrap(~culture,ncol=2)

# violin plot diversity per culture (figure 5)
plot.diversity.culture<-ggplot(filter(data.quat,simpsonagro!="0"), aes(x=culture, y=simpsonagro))+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  geom_jitter(aes(colour=sample.number, shape=biogroup),height=0,width=0.1)+
  scale_colour_gradient(low="grey87", high="black", guide="none")+
  scale_shape_manual(values=c(17,16), guide="none")+
  scale_y_continuous(limits=c(.5,1))+
  labs(x="Archaeological cultures",y="Simpson's index of diversity (crops)")+
  theme_bw()
plot.diversitywild.culture<-ggplot(filter(data.quat,simpsonwild!="0"), aes(x=culture, y=simpsonwild))+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  geom_jitter(aes(colour=sample.number, shape=biogroup),height=0,width=0.1)+
  scale_colour_gradient(low="grey87", high="black",name="Number of samples")+
  scale_shape_manual(values=c(17,16),name="Bioregions")+
  scale_y_continuous(limits=c(.5,1))+
  labs(x="Archaeological cultures",y="Simpson's index of diversity (wild)")+
  theme_bw()
plot_grid(plot.diversity.culture,plot.diversitywild.culture,rel_widths = c(1.45,2))

# violin plot diversity per bioregion (figure 6)
plot.diversity.bioregion<-ggplot(filter(data.quat,simpsonagro!="0"), aes(x=biogroup, y=simpsonagro))+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  geom_jitter(aes(colour=sample.number),height=0,width=0.1)+
  scale_colour_gradient(low="grey87", high="black", guide="none")+
  scale_y_continuous(limits=c(.5,1))+
  labs(x="Bioregions",y="Simpson's index of diversity (crops)")+
  theme_bw()
plot.diversitywild.bioregion<-ggplot(filter(data.quat,simpsonwild!="0"), aes(x=biogroup, y=simpsonwild))+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  geom_jitter(aes(colour=sample.number),height=0,width=0.1)+
  scale_colour_gradient(low="grey87", high="black",name="Number of samples")+
  scale_y_continuous(limits=c(.5,1))+
  labs(x="Bioregions",y="Simpson's index of diversity (wild)")+
  theme_bw()
plot_grid(plot.diversity.bioregion,plot.diversitywild.bioregion,rel_widths = c(1.4,2))


#statistical comparisons of richness and diversity, for both crops/pulses and wild species only#
#test for normality
shapiro.test(data.quat$richnessagro[data.quat$culture=="SKC"])
shapiro.test(data.quat$richnesswild[data.quat$culture=="SKC"])
shapiro.test(data.quat$simpsonagro[data.quat$culture=="SKC"])
shapiro.test(data.quat$simpsonwild[data.quat$simpsonwild!="0" & data.quat$culture=="SKC"])
shapiro.test(data.quat$richnessagro[data.quat$culture=="vinča"])
shapiro.test(data.quat$richnesswild[data.quat$culture=="vinča"])
shapiro.test(data.quat$simpsonagro[data.quat$culture=="vinča"])
shapiro.test(data.quat$simpsonwild[data.quat$simpsonwild!="0" & data.quat$culture=="vinča"])
shapiro.test(data.quat$richnessagro[data.quat$culture=="butmir"])
shapiro.test(data.quat$richnesswild[data.quat$culture=="butmir"])
shapiro.test(data.quat$simpsonagro[data.quat$culture=="butmir"])
shapiro.test(data.quat$simpsonwild[data.quat$simpsonwild!="0" & data.quat$culture=="butmir"])
shapiro.test(data.quat$richnessagro[data.quat$culture=="sopot"])
shapiro.test(data.quat$richnesswild[data.quat$culture=="sopot"])
shapiro.test(data.quat$simpsonagro[data.quat$culture=="sopot"])
shapiro.test(data.quat$simpsonwild[data.quat$simpsonwild!="0" & data.quat$culture=="sopot"])
shapiro.test(data.quat$richnessagro[data.quat$culture=="LBK"])
shapiro.test(data.quat$richnesswild[data.quat$culture=="LBK"])
shapiro.test(data.quat$simpsonagro[data.quat$culture=="LBK"])
shapiro.test(data.quat$simpsonwild[data.quat$simpsonwild!="0" & data.quat$culture=="LBK"])

#statistical comparisons
skc.vinca.richnessagro<-wilcox.test(data.quat$richnessagro[data.quat$culture=="SKC"],
                                    data.quat$richnessagro[data.quat$culture=="vinča"])
skc.butmir.richnessagro<-wilcox.test(data.quat$richnessagro[data.quat$culture=="SKC"],
                                     data.quat$richnessagro[data.quat$culture=="butmir"])
skc.sopot.richnessagro<-wilcox.test(data.quat$richnessagro[data.quat$culture=="SKC"],
                                    data.quat$richnessagro[data.quat$culture=="sopot"])
skc.lbk.richnessagro<-wilcox.test(data.quat$richnessagro[data.quat$culture=="SKC"],
                                  data.quat$richnessagro[data.quat$culture=="LBK"])
vinca.butmir.richnessagro<-t.test(data.quat$richnessagro[data.quat$culture=="vinča"],
                                  data.quat$richnessagro[data.quat$culture=="butmir"])
vinca.sopot.richnessagro<-wilcox.test(data.quat$richnessagro[data.quat$culture=="vinča"],
                                      data.quat$richnessagro[data.quat$culture=="sopot"])
vinca.lbk.richnessagro<-wilcox.test(data.quat$richnessagro[data.quat$culture=="vinča"],
                                    data.quat$richnessagro[data.quat$culture=="LBK"])
sopot.lbk.richnessagro<-wilcox.test(data.quat$richnessagro[data.quat$culture=="sopot"],
                                    data.quat$richnessagro[data.quat$culture=="LBK"])
butmir.sopot.richnessagro<-wilcox.test(data.quat$richnessagro[data.quat$culture=="butmir"],
                                       data.quat$richnessagro[data.quat$culture=="sopot"])
butmir.lbk.richnessagro<-wilcox.test(data.quat$richnessagro[data.quat$culture=="butmir"],
                                     data.quat$richnessagro[data.quat$culture=="LBK"])
skc.vinca.richnesswild<-wilcox.test(data.quat$richnesswild[data.quat$culture=="SKC"],
                                    data.quat$richnesswild[data.quat$culture=="vinča"])
skc.butmir.richnesswild<-wilcox.test(data.quat$richnesswild[data.quat$culture=="SKC"],
                                     data.quat$richnesswild[data.quat$culture=="butmir"])
skc.sopot.richnesswild<-wilcox.test(data.quat$richnesswild[data.quat$culture=="SKC"],
                                    data.quat$richnesswild[data.quat$culture=="sopot"])
skc.lbk.richnesswild<-wilcox.test(data.quat$richnesswild[data.quat$culture=="SKC"],
                                  data.quat$richnesswild[data.quat$culture=="LBK"])
vinca.butmir.richnesswild<-wilcox.test(data.quat$richnesswild[data.quat$culture=="vinča"],
                                       data.quat$richnesswild[data.quat$culture=="butmir"])
vinca.sopot.richnesswild<-wilcox.test(data.quat$richnesswild[data.quat$culture=="vinča"],
                                      data.quat$richnesswild[data.quat$culture=="sopot"])
vinca.lbk.richnesswild<-wilcox.test(data.quat$richnesswild[data.quat$culture=="vinča"],
                                    data.quat$richnesswild[data.quat$culture=="LBK"])
sopot.lbk.richnesswild<-wilcox.test(data.quat$richnesswild[data.quat$culture=="sopot"],
                                    data.quat$richnesswild[data.quat$culture=="LBK"])
butmir.sopot.richnesswild<-wilcox.test(data.quat$richnesswild[data.quat$culture=="butmir"],
                                       data.quat$richnesswild[data.quat$culture=="sopot"])
butmir.lbk.richnesswild<-wilcox.test(data.quat$richnesswild[data.quat$culture=="butmir"],
                                     data.quat$richnesswild[data.quat$culture=="LBK"])
skc.vinca.diversityagro<-wilcox.test(data.quat$simpsonagro[data.quat$culture=="SKC"],
                                     data.quat$simpsonagro[data.quat$culture=="vinča"])
skc.butmir.diversityagro<-wilcox.test(data.quat$simpsonagro[data.quat$culture=="SKC"],
                                      data.quat$simpsonagro[data.quat$culture=="butmir"])
skc.sopot.diversityagro<-wilcox.test(data.quat$simpsonagro[data.quat$culture=="SKC"],
                                     data.quat$simpsonagro[data.quat$culture=="sopot"])
skc.lbk.diversityagro<-wilcox.test(data.quat$simpsonagro[data.quat$culture=="SKC"],
                                   data.quat$simpsonagro[data.quat$culture=="LBK"])
vinca.butmir.diversityagro<-wilcox.test(data.quat$simpsonagro[data.quat$culture=="vinča"],
                                        data.quat$simpsonagro[data.quat$culture=="butmir"])
vinca.sopot.diversityagro<-wilcox.test(data.quat$simpsonagro[data.quat$culture=="vinča"],
                                       data.quat$simpsonagro[data.quat$culture=="sopot"])
vinca.lbk.diversityagro<-wilcox.test(data.quat$simpsonagro[data.quat$culture=="vinča"],
                                     data.quat$simpsonagro[data.quat$culture=="LBK"])
sopot.lbk.diversityagro<-wilcox.test(data.quat$simpsonagro[data.quat$culture=="sopot"],
                                     data.quat$simpsonagro[data.quat$culture=="LBK"])
butmir.sopot.diversityagro<-wilcox.test(data.quat$simpsonagro[data.quat$culture=="butmir"],
                                        data.quat$simpsonagro[data.quat$culture=="sopot"])
butmir.lbk.diversityagro<-wilcox.test(data.quat$simpsonagro[data.quat$culture=="butmir"],
                                      data.quat$simpsonagro[data.quat$culture=="LBK"])
skc.vinca.diversitywild<-wilcox.test(data.quat$simpsonwild[data.quat$culture=="SKC" & data.quat$simpsonwild!="0"],
                                     data.quat$simpsonwild[data.quat$culture=="vinča" & data.quat$simpsonwild!="0"])
skc.butmir.diversitywild<-wilcox.test(data.quat$simpsonwild[data.quat$culture=="SKC"& data.quat$simpsonwild!="0"],
                                      data.quat$simpsonwild[data.quat$culture=="butmir"& data.quat$simpsonwild!="0"])
skc.sopot.diversitywild<-wilcox.test(data.quat$simpsonwild[data.quat$culture=="SKC"& data.quat$simpsonwild!="0"],
                                     data.quat$simpsonwild[data.quat$culture=="sopot"& data.quat$simpsonwild!="0"])
skc.lbk.diversitywild<-wilcox.test(data.quat$simpsonwild[data.quat$culture=="SKC"& data.quat$simpsonwild!="0"],
                                   data.quat$simpsonwild[data.quat$culture=="LBK"& data.quat$simpsonwild!="0"])
vinca.butmir.diversitywild<-wilcox.test(data.quat$simpsonwild[data.quat$culture=="vinča"& data.quat$simpsonwild!="0"],
                                        data.quat$simpsonwild[data.quat$culture=="butmir"& data.quat$simpsonwild!="0"])
vinca.sopot.diversitywild<-wilcox.test(data.quat$simpsonwild[data.quat$culture=="vinča"& data.quat$simpsonwild!="0"],
                                       data.quat$simpsonwild[data.quat$culture=="sopot"& data.quat$simpsonwild!="0"])
vinca.lbk.diversitywild<-wilcox.test(data.quat$simpsonwild[data.quat$culture=="vinča"& data.quat$simpsonwild!="0"],
                                     data.quat$simpsonwild[data.quat$culture=="LBK"& data.quat$simpsonwild!="0"])
sopot.lbk.diversitywild<-wilcox.test(data.quat$simpsonwild[data.quat$culture=="sopot" & data.quat$simpsonwild!="0"],
                                     data.quat$simpsonwild[data.quat$culture=="LBK" & data.quat$simpsonwild!="0"])
butmir.sopot.diversitywild<-wilcox.test(data.quat$simpsonwild[data.quat$culture=="butmir"& data.quat$simpsonwild!="0"],
                                        data.quat$simpsonwild[data.quat$culture=="sopot"& data.quat$simpsonwild!="0"])
butmir.lbk.diversitywild<-wilcox.test(data.quat$simpsonwild[data.quat$culture=="butmir"& data.quat$simpsonwild!="0"],
                                      data.quat$simpsonwild[data.quat$culture=="LBK"& data.quat$simpsonwild!="0"])




####correspondence analysis: sensitivity analysis and plotting####
#plot CA crops: species (left) and CA (right) # figure 7
plot.ca.species.crops <- ggplot(data.quat, aes(x=dim1.crops, y=dim2.crops))+
  geom_blank()+
  xlab(label=paste("Dimension 1 (",Dim1.crops,"%)"))+
  ylab(label=paste("Dimension 2 (",Dim2.crops,"%)"))+
  geom_point(data=species.scores.crops, aes(x=y1, y=y2), shape=3, size=2)+
  geom_label_repel(data=species.scores.crops, aes(x=y1, y=y2), label=species.scores.crops$species)+
  theme_bw()+
  coord_fixed()+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=0,linetype="dotted")
plot.ca.crops <- ggplot(data.quat, aes(x=dim1.crops, y=dim2.crops))+
  geom_blank()+
  xlab(label=paste("Dimension 1 (",Dim1.crops,"%)"))+
  ylab(label=paste("Dimension 2 (",Dim2.crops,"%)"))+
  geom_point(data=data.quat,aes(x=dim1.crops, y=dim2.crops, colour=culture,shape=biogroup), size=2.6)+
  scale_colour_manual(values=c(SKC="#e41a1c", Vinča="#377eb8", Butmir="#4daf4a", Sopot="#984ea3", LBK="#ff7f00"),
                      guide=guide_legend(title="Archaeological\ncultures"))+
  scale_shape_manual(values=c(17,16), guide="none")+
  geom_point(data=species.scores.crops, aes(x=y1, y=y2), shape=3, size=2)+
  theme_bw()+
  coord_fixed(ratio=1)+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=0,linetype="dotted")+
  theme(
    legend.position=c(0.69,0.99),
    legend.justification=c(0,1),
    legend.background=element_blank(),
    legend.box.background=element_rect(colour="black"))
plot_grid(plot.ca.species.crops,plot.ca.crops)

# species contribution to the CA (crops) # figure 8
ca.crops.dim1<-fviz_contrib(ca.quat.crops, choice = "col", axes = 1,fill="lightgrey",color="lightgrey")+
  theme(axis.text.x=element_text(angle=30,size=rel(2)),axis.text.y=element_text(size=rel(2)),
        plot.title=element_text(size=rel(2)),axis.title=element_text(size=rel(2)))
ca.crops.dim2<-fviz_contrib(ca.quat.crops, choice = "col", axes = 2,fill="lightgrey",color="lightgrey")+
  theme(axis.text.x=element_text(angle=30,size=rel(2)),axis.text.y=element_text(size=rel(2)),
        plot.title=element_text(size=rel(2)),axis.title=element_text(size=rel(2)))
plot_grid(ca.crops.dim1,ca.crops.dim2)

#plot CA wild: species (left) and CA (right) # figure 9
plot.ca.species.wild <- ggplot(data.quat, aes(x=dim1.wild, y=dim2.wild))+
  geom_blank()+
  xlab(label=paste("Dimension 1 (",Dim1.wild,"%)"))+
  ylab(label=paste("Dimension 2 (",Dim2.wild,"%)"))+
  geom_point(data=species.scores.wild, aes(x=y1, y=y2), shape=3, size=2)+
  geom_label_repel(data=species.scores.wild, aes(x=y1, y=y2), label=species.scores.wild$species)+
  theme_bw()+
  coord_fixed()+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=0,linetype="dotted")
plot.ca.wild <- ggplot(data.quat, aes(x=dim1.wild, y=dim2.wild))+
  geom_blank()+
  xlab(label=paste("Dimension 1 (",Dim1.wild,"%)"))+
  ylab(label=paste("Dimension 2 (",Dim2.wild,"%)"))+
  geom_point(data=data.quat,aes(x=dim1.wild, y=dim2.wild, colour=culture,shape=biogroup), size=2.6)+
  scale_colour_manual(values=c(SKC="#e41a1c", Vinča="#377eb8", Butmir="#4daf4a", Sopot="#984ea3", LBK="#ff7f00"),
                      guide=guide_legend(title="Archaeological\ncultures"))+
  scale_shape_manual(values=c(17,16), guide="none")+
  geom_point(data=species.scores.wild, aes(x=y1, y=y2), shape=3, size=2)+
  theme_bw()+
  coord_fixed(ratio=1)+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=0,linetype="dotted")+
  theme(
    legend.position=c(0.02,0.98),
    legend.justification=c(0,1),
    legend.background=element_blank(),
    legend.box.background=element_rect(colour="black"))
plot_grid(plot.ca.species.wild,plot.ca.wild)

# species contribution to the CA (wild) # figure 10
ca.wild.dim1<-fviz_contrib(ca.quat.wild, choice = "col", axes = 1,fill="lightgrey",color="lightgrey")+
  theme(axis.text.x=element_text(angle=30,size=rel(2)),axis.text.y=element_text(size=rel(2)),
        plot.title=element_text(size=rel(2)),axis.title=element_text(size=rel(2)))
ca.wild.dim2<-fviz_contrib(ca.quat.wild, choice = "col", axes = 2,fill="lightgrey",color="lightgrey")+
  theme(axis.text.x=element_text(angle=30,size=rel(2)),axis.text.y=element_text(size=rel(2)),
        plot.title=element_text(size=rel(2)),axis.title=element_text(size=rel(2)))
plot_grid(ca.wild.dim1,ca.wild.dim2)


#plot CA results according to arbitrary bins of sample number # figure 11
plot.ca.crops.sample <- ggplot(data.quat, aes(x=dim1.crops, y=dim2.crops))+
  geom_blank()+
  xlab(label=paste("Dimension 1 (",Dim1.crops,"%)"))+
  ylab(label=paste("Dimension 2 (",Dim2.crops,"%)"))+
  geom_point(data=filter(data.quat,!is.na(sample.cat)),aes(x=dim1.crops, y=dim2.crops, colour=sample.cat,shape=biogroup), size=2.6)+
  scale_colour_manual(values=c("1"="#e41a1c", "2"="#377eb8", "3"="#4daf4a", "4"="#984ea3", "5"="#ff7f00"),
                      labels=c("1-5","6-20","21-50","51-100","100+"),
                      guide=guide_legend(title="Number of samples"))+
  scale_shape_manual(values=c(17,16), guide="none")+
  geom_point(data=species.scores.crops, aes(x=y1, y=y2), shape=3, size=2)+
  theme_bw()+
  coord_fixed(ratio=1)+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=0,linetype="dotted")+
  theme(
    legend.position=c(0.6,0.99),
    legend.justification=c(0,1),
    legend.background=element_blank(),
    legend.box.background=element_rect(colour="black"))
plot.ca.wild.sample <- ggplot(data.quat, aes(x=dim1.wild, y=dim2.wild))+
  geom_blank()+
  xlab(label=paste("Dimension 1 (",Dim1.wild,"%)"))+
  ylab(label=paste("Dimension 2 (",Dim2.wild,"%)"))+
  geom_point(data=filter(data.quat,!is.na(sample.cat)),aes(x=dim1.wild, y=dim2.wild, colour=sample.cat,shape=biogroup), size=2.6)+
  scale_colour_manual(values=c("1"="#e41a1c", "2"="#377eb8", "3"="#4daf4a", "4"="#984ea3", "5"="#ff7f00"),
                      labels=c("1-5","6-20","21-50","51-100","100+"),
                      guide=guide_legend(title="Number of samples"))+
  scale_shape_manual(values=c(17,16), guide="none")+
  geom_point(data=species.scores.wild, aes(x=y1, y=y2), shape=3, size=2)+
  theme_bw()+
  coord_fixed(ratio=1)+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=0,linetype="dotted")+
  theme(
    legend.position=c(0.02,0.98),
    legend.justification=c(0,1),
    legend.background=element_blank(),
    legend.box.background=element_rect(colour="black"))

####various statistical tests#### 
#test for normality
shapiro.test(data.quat$dim1.crops[data.quat$culture=="SKC"])
shapiro.test(data.quat$dim2.crops[data.quat$culture=="SKC"])
shapiro.test(data.quat$dim1.wild[data.quat$culture=="SKC"])
shapiro.test(data.quat$dim2.wild[data.quat$culture=="SKC"])
shapiro.test(data.quat$dim1.crops[data.quat$culture=="vinča"])
shapiro.test(data.quat$dim2.crops[data.quat$culture=="vinča"])
shapiro.test(data.quat$dim1.wild[data.quat$culture=="vinča"])
shapiro.test(data.quat$dim2.wild[data.quat$culture=="vinča"])
shapiro.test(data.quat$dim1.crops[data.quat$culture=="butmir"])
shapiro.test(data.quat$dim2.crops[data.quat$culture=="butmir"])
shapiro.test(data.quat$dim1.wild[data.quat$culture=="butmir"])
shapiro.test(data.quat$dim2.wild[data.quat$culture=="butmir"])
shapiro.test(data.quat$dim1.crops[data.quat$culture=="sopot"])
shapiro.test(data.quat$dim2.crops[data.quat$culture=="sopot"])
shapiro.test(data.quat$dim1.wild[data.quat$culture=="sopot"])
shapiro.test(data.quat$dim2.wild[data.quat$culture=="sopot"])
shapiro.test(data.quat$dim1.crops[data.quat$culture=="LBK"])
shapiro.test(data.quat$dim2.crops[data.quat$culture=="LBK"])
shapiro.test(data.quat$dim1.wild[data.quat$culture=="LBK"])
shapiro.test(data.quat$dim2.wild[data.quat$culture=="LBK"])
shapiro.test(data.quat$dim1.crops[data.quat$biogroup=="pannonian"])
shapiro.test(data.quat$dim2.crops[data.quat$biogroup=="pannonian"])
shapiro.test(data.quat$dim1.wild[data.quat$biogroup=="pannonian"])
shapiro.test(data.quat$dim2.wild[data.quat$biogroup=="pannonian"])
shapiro.test(data.quat$dim1.crops[data.quat$biogroup=="continental"])
shapiro.test(data.quat$dim2.crops[data.quat$biogroup=="continental"])
shapiro.test(data.quat$dim1.wild[data.quat$biogroup=="continental"])
shapiro.test(data.quat$dim2.wild[data.quat$biogroup=="continental"])

#testing for difference
ca.dim1.crops.skc.vinca<-t.test(data.quat$dim1.crops[data.quat$culture=="SKC"],
                                data.quat$dim1.crops[data.quat$culture=="vinča"])
ca.dim1.crops.skc.butmir<-t.test(data.quat$dim1.crops[data.quat$culture=="SKC"],
                                 data.quat$dim1.crops[data.quat$culture=="butmir"])
ca.dim1.crops.skc.sopot<-t.test(data.quat$dim1.crops[data.quat$culture=="SKC"],
                                data.quat$dim1.crops[data.quat$culture=="sopot"])
ca.dim1.crops.skc.lbk<-wilcox.test(data.quat$dim1.crops[data.quat$culture=="SKC"],
                                   data.quat$dim1.crops[data.quat$culture=="LBK"])
ca.dim1.crops.vinca.butmir<-t.test(data.quat$dim1.crops[data.quat$culture=="vinča"],
                                   data.quat$dim1.crops[data.quat$culture=="butmir"])
ca.dim1.crops.vinca.sopot<-t.test(data.quat$dim1.crops[data.quat$culture=="vinča"],
                                  data.quat$dim1.crops[data.quat$culture=="sopot"])
ca.dim1.crops.vinca.lbk<-wilcox.test(data.quat$dim1.crops[data.quat$culture=="vinča"],
                                     data.quat$dim1.crops[data.quat$culture=="LBK"])
ca.dim1.crops.butmir.sopot<-t.test(data.quat$dim1.crops[data.quat$culture=="butmir"],
                                   data.quat$dim1.crops[data.quat$culture=="sopot"])
ca.dim1.crops.butmir.lbk<-wilcox.test(data.quat$dim1.crops[data.quat$culture=="butmir"],
                                      data.quat$dim1.crops[data.quat$culture=="LBK"])
ca.dim1.crops.sopot.lbk<-wilcox.test(data.quat$dim1.crops[data.quat$culture=="sopot"],
                                     data.quat$dim1.crops[data.quat$culture=="LBK"])
ca.dim2.crops.skc.vinca<-wilcox.test(data.quat$dim2.crops[data.quat$culture=="SKC"],
                                     data.quat$dim2.crops[data.quat$culture=="vinča"])
ca.dim2.crops.skc.butmir<-wilcox.test(data.quat$dim2.crops[data.quat$culture=="SKC"],
                                      data.quat$dim2.crops[data.quat$culture=="butmir"])
ca.dim2.crops.skc.sopot<-t.test(data.quat$dim2.crops[data.quat$culture=="SKC"],
                                data.quat$dim2.crops[data.quat$culture=="sopot"])
ca.dim2.crops.skc.lbk<-wilcox.test(data.quat$dim2.crops[data.quat$culture=="SKC"],
                                   data.quat$dim2.crops[data.quat$culture=="LBK"])
ca.dim2.crops.vinca.butmir<-wilcox.test(data.quat$dim2.crops[data.quat$culture=="vinča"],
                                        data.quat$dim2.crops[data.quat$culture=="butmir"])
ca.dim2.crops.vinca.sopot<-wilcox.test(data.quat$dim2.crops[data.quat$culture=="vinča"],
                                       data.quat$dim2.crops[data.quat$culture=="sopot"])
ca.dim2.crops.vinca.lbk<-wilcox.test(data.quat$dim2.crops[data.quat$culture=="vinča"],
                                     data.quat$dim2.crops[data.quat$culture=="LBK"])
ca.dim2.crops.butmir.sopot<-wilcox.test(data.quat$dim2.crops[data.quat$culture=="butmir"],
                                        data.quat$dim2.crops[data.quat$culture=="sopot"])
ca.dim2.crops.butmir.lbk<-wilcox.test(data.quat$dim2.crops[data.quat$culture=="butmir"],
                                      data.quat$dim2.crops[data.quat$culture=="LBK"])
ca.dim2.crops.sopot.lbk<-wilcox.test(data.quat$dim2.crops[data.quat$culture=="sopot"],
                                     data.quat$dim2.crops[data.quat$culture=="LBK"])
ca.dim1.wild.skc.vinca<-wilcox.test(data.quat$dim1.wild[data.quat$culture=="SKC"],
                                    data.quat$dim1.wild[data.quat$culture=="vinča"])
ca.dim1.wild.skc.butmir<-t.test(data.quat$dim1.wild[data.quat$culture=="SKC"],
                                data.quat$dim1.wild[data.quat$culture=="butmir"])
ca.dim1.wild.skc.sopot<-t.test(data.quat$dim1.wild[data.quat$culture=="SKC"],
                               data.quat$dim1.wild[data.quat$culture=="sopot"])
ca.dim1.wild.skc.lbk<-t.test(data.quat$dim1.wild[data.quat$culture=="SKC"],
                             data.quat$dim1.wild[data.quat$culture=="LBK"])
ca.dim1.wild.vinca.butmir<-wilcox.test(data.quat$dim1.wild[data.quat$culture=="vinča"],
                                       data.quat$dim1.wild[data.quat$culture=="butmir"])
ca.dim1.wild.vinca.sopot<-wilcox.test(data.quat$dim1.wild[data.quat$culture=="vinča"],
                                      data.quat$dim1.wild[data.quat$culture=="sopot"])
ca.dim1.wild.vinca.lbk<-wilcox.test(data.quat$dim1.wild[data.quat$culture=="vinča"],
                                    data.quat$dim1.wild[data.quat$culture=="LBK"])
ca.dim1.wild.butmir.sopot<-t.test(data.quat$dim1.wild[data.quat$culture=="butmir"],
                                  data.quat$dim1.wild[data.quat$culture=="sopot"])
ca.dim1.wild.butmir.lbk<-t.test(data.quat$dim1.wild[data.quat$culture=="butmir"],
                                data.quat$dim1.wild[data.quat$culture=="LBK"])
ca.dim1.wild.sopot.lbk<-t.test(data.quat$dim1.wild[data.quat$culture=="sopot"],
                               data.quat$dim1.wild[data.quat$culture=="LBK"])
ca.dim2.wild.skc.vinca<-t.test(data.quat$dim2.wild[data.quat$culture=="SKC"],
                               data.quat$dim2.wild[data.quat$culture=="vinča"])
ca.dim2.wild.skc.butmir<-t.test(data.quat$dim2.wild[data.quat$culture=="SKC"],
                                data.quat$dim2.wild[data.quat$culture=="butmir"])
ca.dim2.wild.skc.sopot<-t.test(data.quat$dim2.wild[data.quat$culture=="SKC"],
                               data.quat$dim2.wild[data.quat$culture=="sopot"])
ca.dim2.wild.skc.lbk<-t.test(data.quat$dim2.wild[data.quat$culture=="SKC"],
                             data.quat$dim2.wild[data.quat$culture=="LBK"])
ca.dim2.wild.vinca.butmir<-t.test(data.quat$dim2.wild[data.quat$culture=="vinča"],
                                  data.quat$dim2.wild[data.quat$culture=="butmir"])
ca.dim2.wild.vinca.sopot<-t.test(data.quat$dim2.wild[data.quat$culture=="vinča"],
                                 data.quat$dim2.wild[data.quat$culture=="sopot"])
ca.dim2.wild.vinca.lbk<-t.test(data.quat$dim2.wild[data.quat$culture=="vinča"],
                               data.quat$dim2.wild[data.quat$culture=="LBK"])
ca.dim2.wild.butmir.sopot<-t.test(data.quat$dim2.wild[data.quat$culture=="butmir"],
                                  data.quat$dim2.wild[data.quat$culture=="sopot"])
ca.dim2.wild.butmir.lbk<-t.test(data.quat$dim2.wild[data.quat$culture=="butmir"],
                                data.quat$dim2.wild[data.quat$culture=="LBK"])
ca.dim2.wild.sopot.lbk<-t.test(data.quat$dim2.wild[data.quat$culture=="sopot"],
                               data.quat$dim2.wild[data.quat$culture=="LBK"])
ca.dim1.crops.biogroup<-wilcox.test(data.quat$dim1.crops[data.quat$biogroup=="pannonian"],
                                    data.quat$dim1.crops[data.quat$culture!="pannonian"])
ca.dim2.crops.biogroup<-wilcox.test(data.quat$dim2.crops[data.quat$biogroup=="pannonian"],
                                    data.quat$dim2.crops[data.quat$culture!="pannonian"])
ca.dim1.wild.biogroup<-wilcox.test(data.quat$dim1.wild[data.quat$biogroup=="pannonian"],
                                   data.quat$dim1.wild[data.quat$culture!="pannonian"])
ca.dim2.wild.biogroup<-wilcox.test(data.quat$dim2.wild[data.quat$biogroup=="pannonian"],
                                   data.quat$dim2.wild[data.quat$culture!="pannonian"])


#testing for difference
#test for normality
shapiro.test(data.quat$dim1.crops[data.quat$sample.cat=="1"])
shapiro.test(data.quat$dim2.crops[data.quat$sample.cat=="1"])
shapiro.test(data.quat$dim1.wild[data.quat$sample.cat=="1"])
shapiro.test(data.quat$dim2.wild[data.quat$sample.cat=="1"])
shapiro.test(data.quat$dim1.crops[data.quat$sample.cat=="2"])
shapiro.test(data.quat$dim2.crops[data.quat$sample.cat=="2"])
shapiro.test(data.quat$dim1.wild[data.quat$sample.cat=="2"])
shapiro.test(data.quat$dim2.wild[data.quat$sample.cat=="2"])
shapiro.test(data.quat$dim1.crops[data.quat$sample.cat=="3"])
shapiro.test(data.quat$dim2.crops[data.quat$sample.cat=="3"])
shapiro.test(data.quat$dim1.wild[data.quat$sample.cat=="3"])
shapiro.test(data.quat$dim2.wild[data.quat$sample.cat=="3"])
shapiro.test(data.quat$dim1.crops[data.quat$sample.cat=="4"])
shapiro.test(data.quat$dim2.crops[data.quat$sample.cat=="4"])
shapiro.test(data.quat$dim1.wild[data.quat$sample.cat=="4"])
shapiro.test(data.quat$dim2.wild[data.quat$sample.cat=="4"])
shapiro.test(data.quat$dim1.crops[data.quat$sample.cat=="5"])
shapiro.test(data.quat$dim2.crops[data.quat$sample.cat=="5"])
shapiro.test(data.quat$dim1.wild[data.quat$sample.cat=="5"])
shapiro.test(data.quat$dim2.wild[data.quat$sample.cat=="5"])

#testing for difference
ca.dim1.crops.1.2<-wilcox.test(data.quat$dim1.crops[data.quat$sample.cat=="1"],
                               data.quat$dim1.crops[data.quat$sample.cat=="2"])
ca.dim1.crops.1.3<-wilcox.test(data.quat$dim1.crops[data.quat$sample.cat=="1"],
                               data.quat$dim1.crops[data.quat$sample.cat=="3"])
ca.dim1.crops.1.4<-wilcox.test(data.quat$dim1.crops[data.quat$sample.cat=="1"],
                               data.quat$dim1.crops[data.quat$sample.cat=="4"])
ca.dim1.crops.1.5<-wilcox.test(data.quat$dim1.crops[data.quat$sample.cat=="1"],
                               data.quat$dim1.crops[data.quat$sample.cat=="5"])
ca.dim1.crops.2.3<-t.test(data.quat$dim1.crops[data.quat$sample.cat=="2"],
                          data.quat$dim1.crops[data.quat$sample.cat=="3"])
ca.dim1.crops.2.4<-t.test(data.quat$dim1.crops[data.quat$sample.cat=="2"],
                          data.quat$dim1.crops[data.quat$sample.cat=="4"])
ca.dim1.crops.2.5<-wilcox.test(data.quat$dim1.crops[data.quat$sample.cat=="2"],
                               data.quat$dim1.crops[data.quat$sample.cat=="5"])
ca.dim1.crops.3.4<-t.test(data.quat$dim1.crops[data.quat$sample.cat=="3"],
                          data.quat$dim1.crops[data.quat$sample.cat=="4"])
ca.dim1.crops.3.5<-wilcox.test(data.quat$dim1.crops[data.quat$sample.cat=="3"],
                               data.quat$dim1.crops[data.quat$sample.cat=="5"])
ca.dim1.crops.4.5<-wilcox.test(data.quat$dim1.crops[data.quat$sample.cat=="4"],
                               data.quat$dim1.crops[data.quat$sample.cat=="5"])
ca.dim2.crops.1.2<-wilcox.test(data.quat$dim2.crops[data.quat$sample.cat=="1"],
                               data.quat$dim2.crops[data.quat$sample.cat=="2"])
ca.dim2.crops.1.3<-wilcox.test(data.quat$dim2.crops[data.quat$sample.cat=="1"],
                               data.quat$dim2.crops[data.quat$sample.cat=="3"])
ca.dim2.crops.1.4<-wilcox.test(data.quat$dim2.crops[data.quat$sample.cat=="1"],
                               data.quat$dim2.crops[data.quat$sample.cat=="4"])
ca.dim2.crops.1.5<-wilcox.test(data.quat$dim2.crops[data.quat$sample.cat=="1"],
                               data.quat$dim2.crops[data.quat$sample.cat=="5"])
ca.dim2.crops.2.3<-t.test(data.quat$dim2.crops[data.quat$sample.cat=="2"],
                          data.quat$dim2.crops[data.quat$sample.cat=="3"])
ca.dim2.crops.2.4<-t.test(data.quat$dim2.crops[data.quat$sample.cat=="2"],
                          data.quat$dim2.crops[data.quat$sample.cat=="4"])
ca.dim2.crops.2.5<-t.test(data.quat$dim2.crops[data.quat$sample.cat=="2"],
                          data.quat$dim2.crops[data.quat$sample.cat=="5"])
ca.dim2.crops.3.4<-t.test(data.quat$dim2.crops[data.quat$sample.cat=="3"],
                          data.quat$dim2.crops[data.quat$sample.cat=="4"])
ca.dim2.crops.3.5<-t.test(data.quat$dim2.crops[data.quat$sample.cat=="3"],
                          data.quat$dim2.crops[data.quat$sample.cat=="5"])
ca.dim2.crops.4.5<-t.test(data.quat$dim2.crops[data.quat$sample.cat=="4"],
                          data.quat$dim2.crops[data.quat$sample.cat=="5"])
ca.dim1.wild.1.2<-wilcox.test(data.quat$dim1.wild[data.quat$sample.cat=="1"],
                              data.quat$dim1.wild[data.quat$sample.cat=="2"])
ca.dim1.wild.1.3<-wilcox.test(data.quat$dim1.wild[data.quat$sample.cat=="1"],
                              data.quat$dim1.wild[data.quat$sample.cat=="3"])
ca.dim1.wild.1.4<-wilcox.test(data.quat$dim1.wild[data.quat$sample.cat=="1"],
                              data.quat$dim1.wild[data.quat$sample.cat=="4"])
ca.dim1.wild.1.5<-wilcox.test(data.quat$dim1.wild[data.quat$sample.cat=="1"],
                              data.quat$dim1.wild[data.quat$sample.cat=="5"])
ca.dim1.wild.2.3<-t.test(data.quat$dim1.wild[data.quat$sample.cat=="2"],
                         data.quat$dim1.wild[data.quat$sample.cat=="3"])
ca.dim1.wild.2.4<-t.test(data.quat$dim1.wild[data.quat$sample.cat=="2"],
                         data.quat$dim1.wild[data.quat$sample.cat=="4"])
ca.dim1.wild.2.5<-wilcox.test(data.quat$dim1.wild[data.quat$sample.cat=="2"],
                              data.quat$dim1.wild[data.quat$sample.cat=="5"])
ca.dim1.wild.3.4<-t.test(data.quat$dim1.wild[data.quat$sample.cat=="3"],
                         data.quat$dim1.wild[data.quat$sample.cat=="4"])
ca.dim1.wild.3.5<-wilcox.test(data.quat$dim1.wild[data.quat$sample.cat=="3"],
                              data.quat$dim1.wild[data.quat$sample.cat=="5"])
ca.dim1.wild.4.5<-wilcox.test(data.quat$dim1.wild[data.quat$sample.cat=="4"],
                              data.quat$dim1.wild[data.quat$sample.cat=="5"])
ca.dim2.wild.1.2<-t.test(data.quat$dim2.wild[data.quat$sample.cat=="1"],
                         data.quat$dim2.wild[data.quat$sample.cat=="2"])
ca.dim2.wild.1.3<-t.test(data.quat$dim2.wild[data.quat$sample.cat=="1"],
                         data.quat$dim2.wild[data.quat$sample.cat=="3"])
ca.dim2.wild.1.4<-t.test(data.quat$dim2.wild[data.quat$sample.cat=="1"],
                         data.quat$dim2.wild[data.quat$sample.cat=="4"])
ca.dim2.wild.1.5<-wilcox.test(data.quat$dim2.wild[data.quat$sample.cat=="1"],
                              data.quat$dim2.wild[data.quat$sample.cat=="5"])
ca.dim2.wild.2.3<-t.test(data.quat$dim2.wild[data.quat$sample.cat=="2"],
                         data.quat$dim2.wild[data.quat$sample.cat=="3"])
ca.dim2.wild.2.4<-t.test(data.quat$dim2.wild[data.quat$sample.cat=="2"],
                         data.quat$dim2.wild[data.quat$sample.cat=="4"])
ca.dim2.wild.2.5<-wilcox.test(data.quat$dim2.wild[data.quat$sample.cat=="2"],
                              data.quat$dim2.wild[data.quat$sample.cat=="5"])
ca.dim2.wild.3.4<-t.test(data.quat$dim2.wild[data.quat$sample.cat=="3"],
                         data.quat$dim2.wild[data.quat$sample.cat=="4"])
ca.dim2.wild.3.5<-wilcox.test(data.quat$dim2.wild[data.quat$sample.cat=="3"],
                              data.quat$dim2.wild[data.quat$sample.cat=="5"])
ca.dim2.wild.4.5<-wilcox.test(data.quat$dim2.wild[data.quat$sample.cat=="4"],
                              data.quat$dim2.wild[data.quat$sample.cat=="5"])













