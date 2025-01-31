#### Code to Recreate Figure S5 #### 

# Load Packages

library(reshape2)
library(ggplot2)
library(corto)
library(ggsignif)

# Define Functions

generate_figureS5<-function(mouse_counts,mouse_metadata,human_DEGs,plot_name){
  mouse_t<-data.frame(t(mouse_counts))
  mouse_t$mouseEns<-row.names(mouse_t)
  mouse_t<-merge(mouse_t,human_DEGs[,c("mouseEns","hgnc_symbol")],by="mouseEns")
  row.names(mouse_t)<-mouse_t$hgnc_symbol
  mouse_t$hgnc_symbol<-NULL
  mouse_t$mouseEns<-NULL
  
  translation_consensus_geneset<-list(c("RPL10","RPL10","RPL12","RPL13A","RPL18","RPL18A","RPL19","RPL27","RPL28","RPL3","RPL35","RPL8","RPLP0","RPLP2","RPS10","RPS15","RPS16","RPS18","RPS19","RPS2","RPS21","RPS23","RPS26","RPS28","RPS29","RPS3","RPS5","RPS8","RPSA"))
  nesmat<-ssgsea(mouse_t,translation_consensus_geneset)
  results_ssgsea_mouse<-data.frame(t(nesmat))
  row.names(results_ssgsea_mouse)<-colnames(mouse_t)
  
  results_ssgsea_mouse<-transform(merge(results_ssgsea_mouse,mouse_metadata,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  colnames(results_ssgsea_mouse)[1]<-"ssGSEA"
  temp<-results_ssgsea_mouse[results_ssgsea_mouse$dose != "Uninfected",c("ssGSEA","mtb_strain","mouse_strain")]
  temp$id<- paste(temp$mouse_strain, temp$mtb_strain)
  temp$mtb_strain<-NULL
  temp$mouse_strain <- NULL
  temp<-melt(temp)
  
  desired_order <- c("C57Bl6 H37RV","C3HeBFeJ H37RV","C57Bl6 HN878","C3HeBFeJ HN878")
  temp$id<-factor(temp$id,levels=desired_order)
  setEPS()
  postscript(paste(plot_name))
  FigureS5<- ggplot(temp, aes(x=id, y=value)) + geom_violin(trim=FALSE)+geom_signif(comparisons = list(c("C3HeBFeJ HN878", "C57Bl6 HN878"),c("C3HeBFeJ H37RV","C57Bl6 H37RV")), map_signif_level=TRUE)+theme_classic()+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
  print(FigureS5)
  dev.off()
}

### Define Directory Path & Files to Upload

directory <- "/Users/kpullen/TransCompR/code/forGithub"
files <- list(
  mouse_metadata = "mouse_metadata.csv",
  mouse_counts = "mouse_normalized_orthologs.csv",
  human_metadata = "human_metadata.csv",
  human_counts = "human_normalized_orthologs.csv",
  human_DEGs = "ortholog_dictionary.csv",
  mouse_full = "mouse_normalized.csv",
  human_full = "human_normalized.csv"
)

# Load Data

mouse_metadata <- read.csv(file.path(directory, files$mouse_metadata), header = TRUE)
row.names(mouse_metadata) <- mouse_metadata$Run
mouse_counts <- read.csv(file.path(directory, files$mouse_counts), header = TRUE, row.names = 1)
human_DEGs <- read.csv(file.path(directory, files$human_DEGs), header = TRUE)

### Generate Figure S5

generate_figureS5(mouse_counts,mouse_metadata,human_DEGs,"~/Figure_S5.eps")






