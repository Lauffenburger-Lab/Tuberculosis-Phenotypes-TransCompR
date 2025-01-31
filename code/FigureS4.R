#### Code to Recreate Figure S4 #### 

#Load Packages

library(effsize)
library(caret)
library(ggplot2)
library(ggExtra)

### Define Functions

run_PCA<-function(dataset){
  pca <- prcomp(dataset, center = FALSE, scale. = FALSE)
  eigen <- pca$sdev^2
  scree<- eigen/sum(eigen)
  PCs <- data.frame(pca$x)
  weights <- data.frame(pca$rotation)
  return(list(eigen = eigen,scree = scree,PCs = PCs,weights = weights))
}
generate_figureS4A<-function(human_subset, plot_name){
  results<-run_PCA(human_subset)
  list2env(results, envir = .GlobalEnv)
  compHuman<-transform(merge(PCs,human_metadata[, c("cohort","group")],by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

  # we manually invert a set of the PC scores such that human TB is aligned in the positive direction.
  # this step does not change any of the results of the model but rather hopefully just makes interpretation of the direction of the loadings more straightforward
  compHuman[,c(1:8,11,15:17,19:22)]<-compHuman[,c(1:8,11,15:17,19:22)]*-1
  
  setEPS()
  postscript(paste(plot_name))
  Figure_S4A<-ggplot(compHuman, aes(x = PC1, y = PC2, col = group)) + 
    geom_point() +
    xlab(paste('PC1 (', signif(scree[1]*100,4), '%)')) + 
    ylab(paste('PC2 (', signif(scree[2]*100,4), '%)'))+stat_ellipse()+theme_classic()
  Figure_S4A<-ggMarginal(Figure_S4A,type="boxplot",groupColour = TRUE, groupFill = TRUE,alpha=1)
  print(Figure_S4A)
  dev.off()
  
  PC1_stat <- wilcox.test(PC1 ~ group, data = compHuman)
  print(PC1_stat)
  
  PC2_stat <- wilcox.test(PC2 ~ group, data = compHuman)
  print(PC2_stat)
  
  cumulative_variance<-cumsum(scree)*100
  keep_PCs_variance<-1
  while (cumulative_variance[keep_PCs_variance]<75){
    keep_PCs_variance<-keep_PCs_variance+1
  }
  return(compHuman[,c(1:keep_PCs_variance,dim(compHuman)[2])])
  
}
generate_figureS4B<-function(human_PCs, plot_name){
  
  hum_scores_aTB<-human_PCs[human_PCs$group=="TB",]
  hum_scores_LTBI<-human_PCs[human_PCs$group=="LTBI",]
  cohens_d<- data.frame(matrix(ncol = (dim(human_PCs)[2])-1, nrow = 3)) 
  for (i in 1:(dim(human_PCs)[2]-1)){
    effect_size<-cohen.d(hum_scores_aTB[,i], hum_scores_LTBI[,i])
    cohens_d[1,i]<-effect_size[["estimate"]]
    cohens_d[2,i]<-effect_size[["conf.int"]][["lower"]]
    cohens_d[3,i]<-effect_size[["conf.int"]][["upper"]]
  }
  cohens_d<-t(cohens_d)
  colnames(cohens_d)<-c("coeff","lower","upper")
  cohens_d<-data.frame(cohens_d)
  cohens_d["PCs"]<-c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23")
  
  setEPS()
  postscript(paste(plot_name))
  Figure_S4B<-plotrix::plotCI(x = 1:(dim(human_PCs)[2]-1),               # plotrix plot with confidence intervals
                              y = cohens_d$coeff,
                              li = cohens_d$lower,
                              ui = cohens_d$upper,
                              xlab = "Human Principal Components (hPCs)",
                              ylab = "Effect Size (Cohen's d)")
  print(Figure_S4B)
  abline(h = 1, col = "orange", lty = 2, lwd = 2) 
  abline(h = 0, col = "grey", lty = 2, lwd = 2) 
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

### Load Data

human_metadata <- read.csv(file.path(directory, files$human_metadata), header = TRUE)
row.names(human_metadata) <- human_metadata$Run
human_counts <- read.csv(file.path(directory, files$human_counts), header = TRUE, row.names = 1)
human_DEGs <- read.csv(file.path(directory, files$human_DEGs), header = TRUE)

### Generate Figure S4A

human_PCs<-generate_figureS4A(human_counts,"~/Figure_S4A.eps")

### Generate Figure S4B

generate_figureS4B(human_PCs,"~/Figure_S4B.eps")
