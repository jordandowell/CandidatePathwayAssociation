install.packages("pastecs")
install.packages("ggExtra")
library(ggExtra)
library(pastecs)
library(ggplot2)

#read in the trait data that was used for GWAS
alltraitdata<-read.csv(paste0(getwd(),"/data/",trait_filename))                       
#set colnames -1 as a list 
Trait_names<-colnames(alltraitdata[,-1])
Trait_names
#get all files that have log.txt from assoc tables

GEMMA_Logs<-list.files(path = paste0(getwd(),"/Tables/Assoc_files"), pattern = "log.txt")

#set a mastersheet

i<-1
#get the name of the log files
triallog<-paste0(getwd(),"/Tables/Assoc_files/",GEMMA_Logs[1])

#import log data except for model betas
logdata<-read.table(triallog,skip = 16, nrows=12, comment.char = "", sep = "#")
Mastershet<- as.data.frame(c(gsub("_.*","",GEMMA_Logs[1]),as.numeric(gsub(".*= ","",logdata[,3]))))


for (i in 1:length(GEMMA_Logs)) {
#get the name of the log files
triallog<-paste0(getwd(),"/Tables/Assoc_files/",GEMMA_Logs[i])

#import log data except for model betas
logdata<-read.table(triallog,skip = 16, nrows=12, comment.char = "", sep = "#")
Mastershet[,i+1]<- c(gsub("_.*","",GEMMA_Logs[i]),as.numeric(gsub(".*= ","",logdata[,3])))
#Mastershet<-rbind(Mastersheet,logentry)

}
View(Mastershet)
Mastershet<-t(Mastershet)[-1,]
rownames(Mastershet)<- Mastershet[,1]
colnames(Mastershet)<-c("Trait","GenotypeObs","PhenotypeObs","Covariates","Phenotypes","GenotypeSNPs","PhenotypeSNPs","REMLENULL","MLENULL","PVENull","SEpve","Vkin","Vpop")

Mastershet(Mastershet)

#get descriptive statistices for All traits
ugly<-t(stat.desc(alltraitdata)[,-1])

All_Descriptive_Statistics<- merge(Mastershet, ugly, by=0, all=TRUE)

All_Descriptive_Statistics[,3:ncol(All_Descriptive_Statistics)]<-sapply(All_Descriptive_Statistics[,3:ncol(All_Descriptive_Statistics)],as.character)
All_Descriptive_Statistics[,3:ncol(All_Descriptive_Statistics)]<-sapply(All_Descriptive_Statistics[,3:ncol(All_Descriptive_Statistics)],as.numeric)

#add in h2kin & h2pop

All_Descriptive_Statistics$h2kin<- as.numeric(All_Descriptive_Statistics$Vkin)/(All_Descriptive_Statistics$var)


All_Descriptive_Statistics$h2pop<- as.numeric(All_Descriptive_Statistics$Vpop)/(All_Descriptive_Statistics$var)



write.csv(All_Descriptive_Statistics, paste0(getwd(),"/Tables/Descriptive_Statistics",trait_filename))

#read in trait metadata

TRAITMETA<-read.csv(paste0(getwd(),"/data/AllTraitsforGWAS_2020_METADATA.csv"))

View(TRAITMETA)
#create subsetdata

data<-data.frame(var1=All_Descriptive_Statistics$h2pop,var2=All_Descriptive_Statistics$h2kin,META1=TRAITMETA$META1,META2=TRAITMETA$META2, META3=TRAITMETA$META3)
#create labels for graphs
toplabel<- expression(italic(h)^2~""[POP])
bottomlabel<- expression(italic(h)^2~""[SNP])
 


# classic plot :
p <- ggplot(data, aes(x=var1, y=var2, color=META3)) +
  geom_point() +
  xlab(toplabel)+
  ylab(bottomlabel)+
   #expand_limits(x=c(0,1.1), y=c(0, 1.1))+
   coord_fixed(ratio = 1, clip = "off") +
  theme_minimal()+
theme(legend.position="right")

p
# with marginal histogram
p1 <- ggMarginal(p+theme_classic()+ theme(legend.position = "left"), type="histogram")

p1
p
cor.test(x=data$var1,y=data$var2)
  