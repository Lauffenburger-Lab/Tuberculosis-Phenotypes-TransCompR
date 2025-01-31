#### Code to Recreate Figure S3 #### 

#Load Packages

library(caret)
library(ggplot2)
library(ggExtra)
library(dplyr)
library(reshape2)
library(pals)

### Define Functions

run_PCA<-function(dataset){
  pca <- prcomp(dataset, center = FALSE, scale. = FALSE)
  eigen <- pca$sdev^2
  scree<- eigen/sum(eigen)
  PCs <- data.frame(pca$x)
  weights <- data.frame(pca$rotation)
  return(list(eigen = eigen,scree = scree,PCs = PCs,weights = weights))
}
generate_figureS3A<-function(mouse_subset,human_projected,plot_name){
  
  results<-run_PCA(mouse_subset)
  list2env(results, envir = .GlobalEnv)
  
  #calculate variance of the projected human data
  
  human_variance<-apply(human_projected[,1:20],2,var)
  human_total_variance<-sum(apply(human_projected[,1:88],2,var))
  human_variance_fraction<-human_variance/human_total_variance
  
  #format both species variance data
  
  variance_both_species<-data.frame(cbind(scree[1:20],human_variance_fraction))*100
  colnames(variance_both_species)<-c("mouse","human")
  variance_both_species["id"]<-1:20
  variance_both_species$mouse<-cumsum(variance_both_species$mouse)
  variance_both_species$human<-cumsum(variance_both_species$human)
  zero_row<- data.frame(mouse = 0, human = 0, id = 0)
  variance_both_species<-rbind(zero_row,variance_both_species)
  
  setEPS()
  postscript(paste(plot_name))
  Figure_S3A<-ggplot() +
    geom_line(data = variance_both_species, aes(x = as.numeric(id), y = mouse), color = "blue", size = 1, linetype = "solid") +
    geom_line(data = variance_both_species, aes(x = as.numeric(id), y = human), color = "red", size = 1, linetype = "solid") +
    theme_classic() + geom_hline(yintercept = 75, color = "grey", linetype = "dashed") +
    geom_vline(xintercept = 7, color = "grey", linetype = "dashed") +
    labs(x = "Mouse Principal Components (mPCs)", y = "Cumulative Variance (%)")
  print(Figure_S3A)
  dev.off()
  
}
generate_figureS3B<-function(human_projected,human_metadata,PC_num,plot_name){
  
  #logistic regressions - disease state, biological sex, cohort, BCG Vaccination
  
  human_projected_w_metadata<-transform(merge(human_projected[,1:PC_num],human_metadata[, c("group", "cohort","sex","vaccination","age")],by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  human_projected_w_sex<-human_projected_w_metadata[!(human_projected_w_metadata$sex == ""),]
  human_projected_w_age<-human_projected_w_metadata[!(is.na(human_projected_w_metadata$age)),]
  human_projected_w_vaccination<-human_projected_w_metadata[!(human_projected_w_metadata$vaccination ==""),]
  human_projected_w_vaccination<-human_projected_w_vaccination[!(human_projected_w_vaccination$vaccination =="not_known"),]
  
  human_projected_w_metadata$cohort <- as.numeric(human_projected_w_metadata$cohort == "London")
  human_projected_w_metadata$group <- as.numeric(human_projected_w_metadata$group == "TB")
  human_projected_w_sex$sex<-as.numeric(human_projected_w_sex$sex == "M")
  human_projected_w_vaccination$vaccination <-as.numeric(human_projected_w_vaccination$vaccination == "yes")
  
  z_values <- data.frame(matrix(ncol = PC_num, nrow = 4))
  p_values <- data.frame(matrix(ncol = PC_num, nrow = 4))
  
  perform_logistic_regression <- function(i) {
    logistic_regression_disease <- glm(factor(group) ~ human_projected_w_metadata[, i], family = binomial(link = 'logit'), data = human_projected_w_metadata)
    logistic_regression_cohort <- glm(factor(cohort) ~ human_projected_w_metadata[, i], family = binomial(link = 'logit'), data = human_projected_w_metadata)
    logistic_regression_sex <- glm(factor(sex) ~ human_projected_w_sex[, i], family = binomial(link = 'logit'), data = human_projected_w_sex)
    logistic_regression_vaccination <- glm(factor(vaccination) ~ human_projected_w_vaccination[, i], family = binomial(link = 'logit'), data = human_projected_w_vaccination)
    
    z_values[1, i] <<- summary(logistic_regression_disease)$coefficients[2, "z value"]
    z_values[2, i] <<- summary(logistic_regression_sex)$coefficients[2, "z value"]
    z_values[3, i] <<- summary(logistic_regression_cohort)$coefficients[2, "z value"]
    z_values[4, i] <<- summary(logistic_regression_vaccination)$coefficients[2, "z value"]
    
    p_values[1, i] <<- summary(logistic_regression_disease)$coefficients[2, 4]
    p_values[2, i] <<- summary(logistic_regression_sex)$coefficients[2, 4]
    p_values[3, i] <<- summary(logistic_regression_cohort)$coefficients[2, 4]
    p_values[4, i] <<- summary(logistic_regression_vaccination)$coefficients[2, 4]
  }
  
  lapply(1:PC_num, perform_logistic_regression)
  row.names(z_values) <- c("Disease State", "Biological Sex","Cohort","BCG Vaccination")
  row.names(p_values) <- c("Disease State", "Biological Sex","Cohort","BCG Vaccination")
  
  # linear regression on age
  
  t_values <- data.frame(matrix(ncol = PC_num, nrow = 1))
  p_linear_values <- data.frame(matrix(ncol = PC_num, nrow = 1))
  
  perform_linear_regression <- function(i) {
    linear_regression_age <- lm(as.numeric(age) ~ human_projected_w_age[, i], data = human_projected_w_age)
    t_values[1, i] <<- summary(linear_regression_age)$coefficients[2, "t value"]
    p_linear_values[1, i] <<- summary(linear_regression_age)$coefficients[2, 4]
  }
  
  lapply(1:PC_num, perform_linear_regression)
  row.names(t_values) <- c("age")
  row.names(p_linear_values) <- c("age")
  
  # Combine and adjust p-values for multiple hypothesis correction
  
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
  
  setEPS()
  postscript(paste(plot_name))
  Figure_S3B<-ggplot(aes(x = variable, y = as.factor(id), fill = value), data = values_melted) +
    geom_tile(stat = "identity") +
    scale_fill_gradientn(colors = as.vector(ocean.curl(100)), limits = c(-5, 5)) +
    geom_text(aes(label = stars), color = "black", size = 5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(Figure_S3B)
  dev.off()
  
}
generate_figure4A<-function(mouse_subset,human_subset,plot_name){
  
  results<-run_PCA(mouse_subset)
  list2env(results, envir = .GlobalEnv)
  
  # Arrange and subset genes
  mouse_gene_order <- row.names(weights)
  human_gene_order <- human_DEGs %>%
    arrange(sapply(mouseEns, function(y) which(y == mouse_gene_order)))
  human_subset <- human_subset[, human_gene_order$humEns]
  
  human_projected <- data.frame(as.matrix(human_subset)%*%as.matrix(weights)) # (h human samples)x(g homologs), (g homologs)x(88 mice) -> hx88
  human_projected <-transform(merge(human_projected,human_metadata[, c("cohort", "group")],by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  
  # we manually invert a set of the PC scores such that human TB is aligned in the positive direction.
  # this step does not change any of the results of the model but rather hopefully just makes interpretation of the direction of the loadings more straightforward
  human_projected[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]<-human_projected[,c("PC2","PC5","PC6","PC7","PC11","PC13", "PC14","PC16","PC20")]*-1  
  
  setEPS()
  postscript(paste(plot_name))
  Figure4A <- ggplot(human_projected, aes(x = PC1, y = PC2, col = group)) + geom_point() + stat_ellipse()+theme_classic()
  Figure4A_full<-ggMarginal(Figure4A,type="boxplot",groupColour = TRUE, groupFill = TRUE,alpha=1)
  print(Figure4A_full)
  dev.off()
  
  PC1_stat <- wilcox.test(PC1 ~ group, data = human_projected)
  print(PC1_stat)
  
  PC2_stat <- wilcox.test(PC2 ~ group, data = human_projected)
  print(PC2_stat)
  
  return(human_projected)
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

mouse_metadata <- read.csv(file.path(directory, files$mouse_metadata), header = TRUE)
row.names(mouse_metadata) <- mouse_metadata$Run
mouse_counts <- read.csv(file.path(directory, files$mouse_counts), header = TRUE, row.names = 1)
human_metadata <- read.csv(file.path(directory, files$human_metadata), header = TRUE)
row.names(human_metadata) <- human_metadata$Run
human_counts <- read.csv(file.path(directory, files$human_counts), header = TRUE, row.names = 1)
human_DEGs <- read.csv(file.path(directory, files$human_DEGs), header = TRUE)

### Generate Figure S3A 

# you need to built the TransComp-R model first using generate_figure4a

human_projected<-generate_figure4A(mouse_counts,human_counts,"~/Figure_4A.eps")
generate_figureS3A(mouse_counts,human_projected,"~/Figure_S3A.eps")

### Generate Figure S3B

generate_figureS3B(human_projected,human_metadata,7,"~/Figure_S3B.eps")

