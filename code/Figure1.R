#### Code to Recreate Figure 1B #### 

# Note: the results will vary slightly depending on the version of biomaRt package and the ensembl database use

### Load Packages   

library(biomaRt)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(TBSignatureProfiler)
library(reshape2)

### Define Functions

get_gene_info <- function(mart, dataset, attributes, filters, values) {
  ensembl <- useMart("ensembl", dataset = dataset)
  res <- getBM(attributes = attributes, filters = filters, values = values, mart = ensembl)
  return(res)
}
generate_figure1B<-function(mouse_counts_infected,mouse_metadata_infected,human_counts_full,human_metadata,plot_name){
  
  # First we must convert the gene names in both datasets to gene symbol format 
  human_ensembl <- colnames(human_counts_full)
  mouse_ensembl <- colnames(mouse_counts_full)
  
  resH <- get_gene_info("ensembl", "hsapiens_gene_ensembl", c("ensembl_gene_id", "hgnc_symbol"), "ensembl_gene_id", human_ensembl) 
  colnames(resH)[1]<-"humEns"
  resM <- get_gene_info("ensembl", "mmusculus_gene_ensembl", c("mgi_symbol", "ensembl_gene_id"), "ensembl_gene_id", mouse_ensembl)
  
  human_counts_transpose<-data.frame(t(human_counts_full))
  human_counts_transpose$humEns<-row.names(human_counts_transpose)
  human_counts_transpose<-merge(human_counts_transpose,resH,by="humEns")
  human_counts_transpose <- human_counts_transpose[!duplicated(human_counts_transpose$hgnc_symbol), ]
  row.names(human_counts_transpose)<-human_counts_transpose$hgnc_symbol
  human_counts_transpose[,c("hgnc_symbol","humEns")]<-NULL
  
  mouse_subset_t<-data.frame(t(mouse_counts_infected))
  mouse_subset_t$mouseEns<-row.names(mouse_subset_t)
  colnames(resM)[2]<-"mouseEns"
  mouse_subset_t<-merge( mouse_subset_t,resM,by="mouseEns")
  mouse_subset_t <- mouse_subset_t[!duplicated(mouse_subset_t$mgi_symbol), ]
  mouse_subset_t$mgi_symbol<-toupper(mouse_subset_t$mgi_symbol)
  
  row.names(mouse_subset_t)<-mouse_subset_t$mgi_symbol
  mouse_subset_t$mouseEns<- NULL
  mouse_subset_t$mgi_symbol<- NULL
  
  row_names <- row.names(human_counts_transpose)
  human_counts_transpose <- as.data.frame(lapply(human_counts_transpose, function(x) as.numeric(as.character(x))))
  row.names(human_counts_transpose)<- row_names
  
  row_names <- row.names(mouse_subset_t)
  mouse_subset_t<- as.data.frame(lapply(mouse_subset_t, function(x) as.numeric(as.character(x))))
  row.names(mouse_subset_t)<- row_names
  
  #run single sample GSEA using the published TB gene sets curated by TBSignatureProfiler
  TBSignature_human<- runTBsigProfiler(input = human_counts_transpose,signatures = TBsignatures,algorithm = "ssGSEA",parallel.sz = 1,update_genes = FALSE)
  TBSignature_mouse<- runTBsigProfiler(input = mouse_subset_t,signatures = TBsignatures,algorithm = "ssGSEA",parallel.sz = 1,update_genes = FALSE)
  
  # remove gene sets with less than 10 genes
  too_small<-c("Blankley_5","Chen_HIV_4","Francisco_OD_2","Gjoen_7","Gliddon_OD_3","Gliddon_OD_4","Gliddon_2_OD_4","Gong_OD_4","Hoang_OD_3","Jacobsen_3","Jenum_8","Kaul_3","
LauxdaCosta_OD_3","Lee_4","Maertzdorf_4","PennNich_RISK_6","Rajan_HIV_5","Roe_3","Roe_OD_4","Sloot_HIV_2","Suliman_4","Suliman_RISK_2","Suliman_RISK_4","Sweeney_OD_3","Thompson_RES_5","Thompson_9","Zhao_NANO_6","Zimmer_RES_3","Chen_5","Gliddon_HIV_3","Kulkarni_HIV_2","LauxdaCosta_OD_3","Natarajan_7","Chen_5")
  
  TBSignature_human<-data.frame(t(TBSignature_human))
  TBSignature_human$group<-human_metadata$group
  TBSignature_human<-TBSignature_human[,!(colnames(TBSignature_human) %in% too_small)]
  
  #Aggregate GSEA scores for each phenotype group (ie- ATB vs LTBI) by taking the median score
  TBSignature_human_aggregate<-aggregate(TBSignature_human[,1:(dim(TBSignature_human)[2]-1)],by=list(TBSignature_human$group),FUN=median)
  TBSignature_human_melt<-melt(TBSignature_human_aggregate)
  TBSignature_human_melt$Group.2<-"human"
  
  TBSignature_mouse<-data.frame(t(TBSignature_mouse))
  TBSignature_mouse$group<-paste(mouse_metadata_infected$mouse_strain,mouse_metadata_infected$mtb_strain,sep="_")
  TBSignature_mouse<-TBSignature_mouse[,!(colnames(TBSignature_mouse) %in% too_small)]
  
  #Aggregate GSEA scores for each phenotype group (ie- by mouse and mtb strain) by taking the median score
  TBSignature_mouse_aggregate<-aggregate(TBSignature_mouse[,1:(dim(TBSignature_mouse)[2]-1)],by=list(TBSignature_mouse$group),FUN=median)
  TBSignature_mouse_melt<-melt(TBSignature_mouse_aggregate)
  TBSignature_mouse_melt$Group.2[TBSignature_mouse_melt$Group.1 == "C3HeBFeJ_HN878"]<-"HN878"
  TBSignature_mouse_melt$Group.2[TBSignature_mouse_melt$Group.1 == "C57Bl6_HN878"]<-"HN878"
  TBSignature_mouse_melt$Group.2[TBSignature_mouse_melt$Group.1 == "C3HeBFeJ_H37RV"]<-"H37RV"
  TBSignature_mouse_melt$Group.2[TBSignature_mouse_melt$Group.1 == "C57Bl6_H37RV"]<-"H37RV"
  
  TBSignature_both<-rbind(TBSignature_human_melt,TBSignature_mouse_melt) 
  
  # For plot formatting purposes only, we label the C3HeBFeJ mice as "TB" and C57Bl6 as "LTBI"
  TBSignature_both$Group.1[TBSignature_both$Group.1 == "C3HeBFeJ_H37RV"]<-"TB"
  TBSignature_both$Group.1[TBSignature_both$Group.1 == "C57Bl6_H37RV"]<-"LTBI"
  TBSignature_both$Group.1[TBSignature_both$Group.1 == "C3HeBFeJ_HN878"]<-"TB"
  TBSignature_both$Group.1[TBSignature_both$Group.1 == "C57Bl6_HN878"]<-"LTBI"

  setEPS()
  postscript(paste(plot_name))
  p<-ggpubr::ggpaired(TBSignature_both, x = "Group.1", y = "value",color = "Group.1",id = "variable", palette = "jco",line.color = "gray", line.size = 0.4,facet.by = "Group.2") 
  p<-p + ggpubr::stat_compare_means(label = "p.format", paired = TRUE,method="wilcox.test")
  print(p)
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

# Load data

mouse_counts_full <- read.csv(file.path(directory, files$mouse_full), header = TRUE, row.names = 1)
human_counts_full <- read.csv(file.path(directory, files$human_full), header = TRUE, row.names = 1)
human_metadata <- read.csv(file.path(directory, files$human_metadata), header = TRUE, row.names = 1)
mouse_metadata <- read.csv(file.path(directory, files$mouse_metadata), header = TRUE, row.names = 1)
mouse_metadata_infected<-mouse_metadata[mouse_metadata$dose != "Uninfected",]
mouse_counts_infected<-mouse_counts_full[row.names(mouse_counts_full) %in% mouse_metadata_infected$Run,]

# Generate Figure 1B

generate_figure1B(mouse_counts_infected,mouse_metadata_infected,human_counts_full,human_metadata,"~/Figure1B.eps")


