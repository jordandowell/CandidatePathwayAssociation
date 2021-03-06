#### dendrogram of traits in GWAS
library(corrr)
library(Hmisc)
library(ggdendro)
# library(stringi)

#### read in preferences
prefs<-read.table("Scripts/### Preferences ###",header=F,sep="=",skip=1)
  SNPset<-as.character(prefs[2,2])
  pheno.name<-as.character(prefs[1,2])
  multcomp<-as.numeric(as.character(prefs[3,2]))


pheno.data<-read.csv(paste("data/",pheno.name,sep=""))


##### data in the colocate plot
traits<-unique(plot.data$trait)
attempt <- try(expr = 
      {
        traits.env<-paste(traits,"_",envs[q],sep="")
        correlate.data<-pheno.data[,match(traits.env,names(pheno.data))]
      },silent=TRUE)
if(inherits(x = attempt,what = "try-error"))
{
  traits.env<-paste(traits,"_",envs[q],sep="")
  correlate.data<-pheno.data[,grep(traits.env,names(pheno.data))]
}


names(correlate.data)<-traits




### make environment panel correlations
Env.corr<-rcorr(as.matrix(correlate.data),type="pearson")

#### clustering dendrogram to use with collocalization plot
if (length(traits)>=2) {
Env.dist<-as.dist(1-Env.corr$r)    
  Env.clust <- hclust(Env.dist, method = "complete", members=NULL)
    Env.dendro<-ggdendrogram(Env.clust,rotate=T,labels=F)+scale_y_reverse()+theme(axis.text.y=element_blank())+scale_x_continuous(limits=c(0.5,length(traits)+0.5),expand=c(0,0))
      Env.label.order<-Env.clust$labels[Env.clust$order]
}
if (length(traits)==1) {
  Env.dendro <- ggplot() + theme_void()
  Env.label.order <- traits
}
      
  