#### Code to Recreate Figure S2 #### 

### Load Packages 

library(reshape2)
library(pals)
library(ggExtra)
library(ggplot2)

### Define Functions

run_PCA<-function(dataset){
  pca <- prcomp(dataset, center = FALSE, scale. = FALSE)
  eigen <- pca$sdev^2
  scree<- eigen/sum(eigen)
  PCs <- data.frame(pca$x)
  weights <- data.frame(pca$rotation)
  return(list(eigen = eigen,scree = scree,PCs = PCs,weights = weights))
}
generate_figureS2A<-function(mouse_subset,plot_name){

  results<-run_PCA(mouse_subset)
  list2env(results, envir = .GlobalEnv)
  compMus<-transform(merge(PCs,mouse_metadata[, c("mtb_strain", "dose", "mouse_strain", "day")],by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  
  # we manually invert a set of the PC scores because later in the transCompR pipeline we rotate all of the mPCs such that human TB is aligned in the positive direction.
  # this step does not change any of the results of the model but rather hopefully just makes interpretation of the direction of the loadings more straightforward
  compMus[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]<-compMus[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]*-1 #we invert because later in the transCompR pipeline we rotate all of the mPCs such that human TB is aligned in the positive direction.
  
  setEPS()
  postscript(paste(plot_name))
  FigureS2A<-ggplot(compMus, aes(x = PC1, y = PC2, col = as.factor(mtb_strain))) + 
    geom_point(size = 4) +
    xlab(paste('PC1 (', signif(scree[1]*100,4), '%)')) + 
    ylab(paste('PC2 (', signif(scree[2]*100,4), '%)'))+stat_ellipse()+theme_classic()
  FigureS2A_full<-ggMarginal(FigureS2A,type="boxplot",groupColour = TRUE, groupFill = TRUE,alpha=1)
  print(FigureS2A_full)
  dev.off()
  
  PC1_stat <- pairwise.wilcox.test(x = compMus$PC1, g= compMus$mtb_strain,p.adjust.method = "BH")
  print(PC1_stat)
  
  PC2_stat <- pairwise.wilcox.test(x = compMus$PC2, g= compMus$mtb_strain, p.adjust.method = "BH")
  print(PC2_stat)
}
generate_figureS2B<-function(mouse_subset,plot_name){
  
  results<-run_PCA(mouse_subset)
  list2env(results, envir = .GlobalEnv)
  compMus<-transform(merge(PCs,mouse_metadata[, c("mtb_strain", "dose", "mouse_strain", "day")],by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  
  # we manually invert a set of the PC scores because later in the transCompR pipeline we rotate all of the mPCs such that human TB is aligned in the positive direction.
  # this step does not change any of the results of the model but rather hopefully just makes interpretation of the direction of the loadings more straightforward
  compMus[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]<-compMus[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]*-1 #we invert because later in the transCompR pipeline we rotate all of the mPCs such that human TB is aligned in the positive direction.
  compMus$day[compMus$dose == "Uninfected"]<-NA
  
  setEPS()
  postscript(paste(plot_name))
  FigureS2B<-ggplot(compMus, aes(x = PC1, y = PC2, fill = day)) +
    geom_point(data = subset(compMus, !is.na(day)), 
               shape = 21, size = 4, color = "black", stroke = 0.5) +
    geom_point(data = subset(compMus, is.na(day)), 
               fill = "white", shape = 21, size = 4, color = "black", stroke = 0.5) +
    xlab(paste('PC1 (', signif(scree[1] * 100, 4), '%)')) +
    ylab(paste('PC2 (', signif(scree[2] * 100, 4), '%)')) +
    scale_fill_gradient(low = "#FCE0D9", high = "#581059") + 
    theme_bw()
  print(FigureS2B)
  dev.off()
}
generate_figureS2C<-function(mouse_subset,plot_name){
  
  results<-run_PCA(mouse_subset)
  list2env(results, envir = .GlobalEnv)
  compMus<-transform(merge(PCs,mouse_metadata[, c("mtb_strain", "dose", "mouse_strain", "day")],by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  
  # we manually invert a set of the PC scores because later in the transCompR pipeline we rotate all of the mPCs such that human TB is aligned in the positive direction.
  # this step does not change any of the results of the model but rather hopefully just makes interpretation of the direction of the loadings more straightforward
  compMus[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]<-compMus[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]*-1 #we invert because later in the transCompR pipeline we rotate all of the mPCs such that human TB is aligned in the positive direction.
  
  setEPS()
  postscript(paste(plot_name))
  FigureS2C<-ggplot(compMus, aes(x = PC1, y = PC2, col = dose)) + 
    geom_point() +
    xlab(paste('PC1 (', signif(scree[1]*100,4), '%)')) + 
    ylab(paste('PC2 (', signif(scree[2]*100,4), '%)'))+stat_ellipse()+theme_classic()
  FigureS2C_full<-ggMarginal(FigureS2C,type="boxplot",groupColour = TRUE, groupFill = TRUE,alpha=1)
  print(FigureS2C_full)
  dev.off()
  
  PC1_stat <- pairwise.wilcox.test(x = compMus$PC1, g= compMus$dose, p.adjust.method = "BH")
  print(PC1_stat)
  
  PC2_stat <- pairwise.wilcox.test(x = compMus$PC2, g= compMus$dose, p.adjust.method = "BH")
  print(PC2_stat)
  
}

### Define Directory Path & Files

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

### Load Data

mouse_metadata <- read.csv(file.path(directory, files$mouse_metadata), header = TRUE)
row.names(mouse_metadata) <- mouse_metadata$Run
mouse_counts <- read.csv(file.path(directory, files$mouse_counts), header = TRUE, row.names = 1)

### Generate Figure S2A

generate_figureS2A(mouse_counts,"~/Figure_S2A.eps")

### Generate Figure S2B

generate_figureS2B(mouse_counts,"~/Figure_S2B.eps")

### Generate Figure S2C

generate_figureS2C(mouse_counts,"~/Figure_S2C.eps")



