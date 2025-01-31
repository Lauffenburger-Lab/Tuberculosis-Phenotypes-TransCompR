#### Code to Recreate Figure 3 #### 

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
generate_figure3A<-function(mouse_subset, mouse_metadata, plot_name){
  
  # Run PCA
  results<-run_PCA(mouse_subset)
  list2env(results, envir = .GlobalEnv)
  compMus<-transform(merge(PCs,mouse_metadata[, c("mtb_strain", "dose", "mouse_strain", "day")],by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  
  # we manually invert a set of the PC scores because later in the transCompR pipeline we rotate all of the mPCs such that human TB is aligned in the positive direction.
  # this step does not change any of the results of the model but rather hopefully just makes interpretation of the direction of the loadings more straightforward
  compMus[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]<-compMus[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]*-1
  
  # Visualize results
  
  setEPS()
  postscript(paste(plot_name))
  Figure3A<-ggplot(compMus, aes(x = PC1, y = PC2, col = mouse_strain)) + 
    geom_point() +
    xlab(paste('PC1 (', signif(scree[1]*100,4), '%)')) + 
    ylab(paste('PC2 (', signif(scree[2]*100,4), '%)'))+stat_ellipse()+theme_classic()
  Figure3A_full<-ggMarginal(Figure3A,type="boxplot",groupColour = TRUE, groupFill = TRUE,alpha=1)
  print(Figure3A_full)
  dev.off()
  
  PC1_stat <- wilcox.test(PC1 ~ mouse_strain, data = compMus)
  print(PC1_stat)
  
  PC2_stat <- wilcox.test(PC2 ~ mouse_strain, data = compMus)
  print(PC2_stat)
  
}
generate_figure3B<-function(mouse_subset, mouse_metadata, plot_name){
  
  results<-run_PCA(mouse_subset)
  list2env(results, envir = .GlobalEnv)
  compMus<-transform(merge(PCs,mouse_metadata[, c("mtb_strain", "dose", "mouse_strain", "day")],by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  
  # we manually invert a set of the PC scores because later in the transCompR pipeline we rotate all of the mPCs such that human TB is aligned in the positive direction.
  # this step does not change any of the results of the model but rather hopefully just makes interpretation of the direction of the loadings more straightforward
  compMus[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]<-compMus[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]*-1  
  
  compMus$mouse_strain <- as.numeric(compMus$mouse_strain == "C3HeBFeJ")
  compMus$mtb_strain <- as.numeric(compMus$mtb_strain == "HN878")
  compMus$days <- mouse_metadata$day
  compMus$bacterial_burden <- as.numeric(factor(compMus$dose, levels = c("Uninfected", "Low", "High")))
  
  # for the regressions on number of days infected and mtb strain, the scores for uninfected mice are removed 
  compMus_no_control <- compMus[compMus$dose != "Uninfected", ]
  
  # logistic regression to identify if distribution of scores on any given PC can be explained by mouse strain or mtb strain
 
  z_values <- data.frame(matrix(ncol = 20, nrow = 2))
  p_values <- data.frame(matrix(ncol = 20, nrow = 2))
  
  perform_logistic_regression <- function(i) {
    logistic_regression_mouse_strain <- glm(factor(mouse_strain) ~ compMus[, i], family = binomial(link = 'logit'), data = compMus)
    logistic_regression_mtb_strain <- glm(factor(mtb_strain) ~ compMus_no_control[, i], family = binomial(link = 'logit'), data = compMus_no_control)
    
    z_values[2, i] <<- summary(logistic_regression_mouse_strain)$coefficients[2, "z value"]
    z_values[1, i] <<- summary(logistic_regression_mtb_strain)$coefficients[2, "z value"]
    
    p_values[2, i] <<- summary(logistic_regression_mouse_strain)$coefficients[2, 4]
    p_values[1, i] <<- summary(logistic_regression_mtb_strain)$coefficients[2, 4]
  }
  
  lapply(1:20, perform_logistic_regression)
  
  row.names(z_values) <- c("mtb_strain", "mouse_strain")
  row.names(p_values) <- c("mtb_strain", "mouse_strain")
  
  # linear regression to identify if distribution of scores on any given PC can be explained by number of days infected or the bacterial "dose" each animal was exposed to
  
  t_values <- data.frame(matrix(ncol = 20, nrow = 2))
  p_linear_values <- data.frame(matrix(ncol = 20, nrow = 2))
  
  perform_linear_regression <- function(i) {
    linear_regression_bacterial_burden <- lm(as.numeric(bacterial_burden) ~ compMus[, i], data = compMus)
    linear_regression_days <- lm(as.numeric(days) ~ compMus_no_control[, i], data = compMus_no_control)
    
    t_values[1, i] <<- summary(linear_regression_days)$coefficients[2, "t value"]
    t_values[2, i] <<- summary(linear_regression_bacterial_burden)$coefficients[2, "t value"]
    
    p_linear_values[1, i] <<- summary(linear_regression_days)$coefficients[2, 4]
    p_linear_values[2, i] <<- summary(linear_regression_bacterial_burden)$coefficients[2, 4]
  }
  
  lapply(1:20, perform_linear_regression)
  
  row.names(t_values) <- c("days", "bacterial_burden")
  row.names(p_linear_values) <- c("days", "bacterial_burden")
  
  # Combine and adjust p-values for multiple hypothesis testing
  p_all <- rbind(p_linear_values, p_values)
  values_all <- rbind(t_values, z_values)
  p_all$id <- row.names(p_all)
  p_values_BH <- melt(p_all)
  p_values_BH$value <- p.adjust(p_values_BH$value, method = "BH")
  
  values_melted <- values_all
  values_melted$id <- row.names(values_melted)
  values_melted <- melt(values_melted)
  values_melted$pval <- p_values_BH$value
  
  values_melted$stars <- cut(values_melted$pval, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", ".", ""))
  
  # Visualize results
  
  setEPS()
  postscript(paste(plot_name))
  Figure_3B<-ggplot(aes(x = variable, y = as.factor(id), fill = value), data = values_melted) +
    geom_tile(stat = "identity") +
    scale_fill_gradientn(colors = as.vector(ocean.curl(100)), limits = c(-9, 9)) +
    geom_text(aes(label = stars), color = "black", size = 5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(Figure_3B)
  dev.off()
}
generate_figure3C<-function(mouse_subset, mouse_metadata, plot_name){
  
  # Run PCA
  results<-run_PCA(mouse_subset)
  list2env(results, envir = .GlobalEnv)
  cumulative_variance <- cumsum(scree)
  
  # Visualize results
  
  setEPS()
  postscript(paste(plot_name))
  plot(cumulative_variance*100,type='b',col = '#10676E',pch = 21,bg = '#10676E',xlab='mouse principal components (mPCs)',ylab='Variance Explained (%)')
  abline(h = 75, col = "black",lty=2)
  dev.off()
  
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

### Generate Figure 3A
generate_figure3A(mouse_counts, mouse_metadata, "~/Figure_3A.eps")

### Generate Figure 3B
generate_figure3B(mouse_counts, mouse_metadata, "~/Figure_3B.eps")

### Generate Figure 3C
generate_figure3C(mouse_counts, mouse_metadata, "~/Figure_3C.eps")







