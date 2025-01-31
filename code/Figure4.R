#### Code to Recreate Figure 4 #### 

#Load Packages

library(effsize)
library(caret)
library(ggplot2)
library(ggExtra)
library(ggsignif)
library(pROC)
library(dplyr)
library(reshape2)


### Define Functions

run_PCA<-function(dataset){
  pca <- prcomp(dataset, center = FALSE, scale. = FALSE)
  eigen <- pca$sdev^2
  scree<- eigen/sum(eigen)
  PCs <- data.frame(pca$x)
  weights <- data.frame(pca$rotation)
  return(list(eigen = eigen,scree = scree,PCs = PCs,weights = weights))
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
generate_figure4B<-function(mouse_subset,human_projected,plot_name){
  
  results<-run_PCA(mouse_subset)
  list2env(results, envir = .GlobalEnv)
  
  #calculate variance of the projected human data
  
  human_variance<-apply(human_projected[,1:20],2,var)
  human_total_variance<-sum(apply(human_projected[,1:88],2,var))
  human_variance_fraction<-human_variance/human_total_variance
  
  #format both species variance data
  
  variance_both_species<-data.frame(cbind(scree[1:20],human_variance_fraction))*100
  colnames(variance_both_species)<-c("mouse","human")
  variance_both_species["id"]<-row.names(variance_both_species)
  variance_both_species_melted<-melt(variance_both_species)
  variance_both_species_melted$id <- factor(variance_both_species_melted$id,levels = unique(variance_both_species_melted$id))
  
  
  setEPS()
  postscript(paste(plot_name))
  Figure_4B<-ggplot(variance_both_species_melted, aes(fill=variable, y=value, x=id)) +
    geom_bar(position="dodge", stat="identity")+xlab("Mouse Principal Components (mPCs)")+ylab("Variance Explained (%)")+theme_classic()
  print(Figure_4B)
  dev.off()
  
  return(variance_both_species)
}
generate_figure4C<-function(human_projected,variance_both_species,plot_name){
  
  variance_cumulative<-variance_both_species
  variance_cumulative$mouse <-cumsum(variance_cumulative$mouse)
  variance_cumulative$human <-cumsum(variance_cumulative$human)
  variance_cumulative["id"]<-row.names(variance_cumulative)
  variance_cumulative_melted<-melt(variance_cumulative)
  variance_cumulative_melted$id <- factor(variance_cumulative_melted$id,levels = unique(variance_cumulative_melted$id))
  
  # Filter down to mPCs that describe 75% of projected human variance
  keep_PCs_variance<-1
  while (variance_cumulative$human[keep_PCs_variance]<75){
    keep_PCs_variance<-keep_PCs_variance+1
  }
  
  human_projected_subset<-human_projected[,c(1:keep_PCs_variance,dim(human_projected)[2])]
  human_projected_aTB<-human_projected_subset[human_projected_subset$group=="TB",]
  human_projected_LTBI<-human_projected_subset[human_projected_subset$group=="LTBI",]
  
  #Calculate effect size between ATB and LTBI projected scores using cohen's D 
  
  cohens_d<- data.frame(matrix(ncol = 7, nrow = 3)) 
  for (i in 1:(dim(cohens_d)[2])){
    effect_size<-cohen.d(human_projected_aTB[,i], human_projected_LTBI[,i])
    cohens_d[1,i]<-effect_size[["estimate"]]
    cohens_d[2,i]<-effect_size[["conf.int"]][["lower"]]
    cohens_d[3,i]<-effect_size[["conf.int"]][["upper"]]
  }
  
  cohens_d<-t(cohens_d)
  colnames(cohens_d)<-c("coeff","lower","upper")
  cohens_d<-data.frame(cohens_d)
  
  setEPS()
  postscript(paste(plot_name))
  Figure_4C<-plotrix::plotCI(x = 1:keep_PCs_variance,               # plotrix plot with confidence intervals
                             y = cohens_d$coeff,
                             li = cohens_d$lower,
                             ui = cohens_d$upper,
                             xlab = "Mouse Principal Components (mPCs)",
                             ylab = "Effect Size (Cohen's d)")
  print(Figure_4C)
  abline(h = 1, col = "orange", lty = 2, lwd = 2) 
  abline(h = 0, col = "grey", lty = 2, lwd = 2) 
  dev.off()
  
  # Keep mPCs with an effect size greater than 1 
  
  keep_PCs<-list()
  for (i in 1:dim(cohens_d)[1]){
    if (cohens_d[i,1]>0){
      if (cohens_d[i,2]>1){
        keep_PCs<-append(keep_PCs,i)
      }
    }else{
      if (cohens_d[i,3]<(-1)){
        keep_PCs<-append(keep_PCs,i)
      }
    }
  }
  keep_PCs<-unlist(keep_PCs)
  return(list(keep_PCs = keep_PCs,keep_PCs_variance = keep_PCs_variance))
}
generate_figure4D<-function(human_projected,keep_PCs_variance,plot_name1,plot_name2){
  
  human_projected_subset<-human_projected[,c(1:keep_PCs_variance,dim(human_projected)[2])]
  coefs_cvs <- data.frame(matrix(nrow = keep_PCs_variance, ncol = 5))
  pval <- data.frame(matrix(nrow = keep_PCs_variance, ncol = 5))
  set.seed(123)
  folds <- caret::createFolds(human_projected_subset$group, k = 5)
  for (i  in 1:5){
    for (j  in 1:keep_PCs_variance){
      training <- human_projected_subset[-folds[[i]],c(j,8)] # 67 training samples
      testing  <- human_projected_subset[folds[[i]],c(j,8)]
      y<-as.factor((as.numeric(factor(training$group))-1))
      
      glmFitCVHum <- train(
        x=data.frame(training[,1]), # 96 humans x 5 hPCs matrix of predictor variables
        y=y, # 2-level (LTBI/TB) factor for TB status for 96 humans
        method = "glm",
        family=binomial('logit'), 
        maxit=100
      )
      valuesHum <- summary(glmFitCVHum)$coef
      pval[j,i] <- summary(glmFitCVHum)$coef[2,4]
      coefs_cvs[j,i] <- summary(glmFitCVHum)$coef[-1,1]
    }}
  
  coefs_cvs["id"]<-c("mPC1","mPC2","mPC3","mPC4","mPC5","mPC6","mPC7")
  coefs_cvs_melted<-melt(coefs_cvs)
  
  pval["id"]<-c("mPC1","mPC2","mPC3","mPC4","mPC5","mPC6","mPC7")
  pval_melted<-melt(pval)
  pval_melted$p_adj <- p.adjust(pval_melted$value, method = "BH")
  
  setEPS()
  postscript(paste(plot_name1))
  Figure_4D_1 <- ggplot(coefs_cvs_melted, aes(x=id, y=value)) + geom_boxplot()+theme_classic()+xlab("Mouse Principal Components")+ylab("coefficient")
  print(Figure_4D_1)
  dev.off()
  
  setEPS()
  postscript(paste(plot_name2))
  Figure_4D_2 <- ggplot(pval_melted,aes(x=id, y=-1*log10(p_adj))) + geom_boxplot()+theme_classic()+geom_hline(yintercept =-1*log10(0.05))+xlab("Mouse Principal Components")+ylab("-log(p-value)")
  print(Figure_4D_2)
  dev.off()
}
generate_figure4E <- function(human_projected, PC_indices, plot_name) {
  hum_scores_regression <- human_projected[, c(PC_indices, dim(human_projected)[2])]
  
  cross_validate <- function(hum_scores_regression, feature_index,full) {
    human_preds <- data.frame(matrix(ncol = 100, nrow = 89))
    
    for (k in 1:100) {
      set.seed(k)
      folds <- createFolds(hum_scores_regression$group, k = 5)
      for (i in 1:5) {
        training <- hum_scores_regression[-folds[[i]], ]
        testing <- hum_scores_regression[folds[[i]], ]
        X_sel <- training[, -c(6)]
        y <- as.numeric(factor(training$group, levels = c("LTBI", "TB"))) - 1
        if (full == 0){
          temp_x <- data.frame(feat = as.numeric(X_sel[, feature_index]))
          temp_test <- data.frame(feat = as.numeric(testing[, feature_index]))
        } 
        if (full == 1){
          temp_x <- X_sel
          temp_test <- testing
        } 
        
        glmFitCVHum <- train(
          x = temp_x,
          y = y,
          method = "glm",
          family = binomial('logit'),
          maxit = 100)
        
        classes <- predict(glmFitCVHum, newdata = temp_test)
        human_preds[folds[[i]], k] <- classes
      }
    }
    return(rowMeans(human_preds))
  }
  
  # Function to calculate ROC
  plot_roc <- function(predictions, group) {
    roc(group ~ predictions, plot = TRUE, print.auc = TRUE, legacy.axes = TRUE, col = "black")
  }
  
  roc_results <- list()
  for (index in PC_indices) {
    PC_preds <- cross_validate(hum_scores_regression, paste0("PC", index),0)
    roc_results[[paste0("PC", index)]] <- plot_roc(PC_preds, hum_scores_regression$group)
  }
  
  # Full model analysis
  full_model_preds <- cross_validate(hum_scores_regression, keep_PCs,1)
  roc_results[["full_model"]] <- plot_roc(full_model_preds, hum_scores_regression$group)
  
  setEPS()
  postscript(paste(plot_name))
  Figure3E <- ggroc(roc_results, legacy.axes = TRUE) +
    scale_color_manual(values = c("deepskyblue", "brown1", "chartreuse", "darkorchid", "darkorange", "black")) +
    xlab("FPR") + ylab("TPR") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "darkgrey", linetype = "dashed") +
    theme_classic() +
    theme(legend.position = "none")
  print(Figure3E)
  dev.off()
}
generate_figure4F<-function(human_projected, PC_indices, plot_name){
  
  # Note: results will vary slightly each time due to setting a random seed
  
  ## Full Model Validation
  # This loop with run 100 rounds of 5-fold CV to predict disease phenotype from the model including the 5 mPCs we identified as translatable components (TCs)
  
  hum_scores_regression <- human_projected[, c(PC_indices, dim(human_projected)[2])]
  auc_mouse<-data.frame(matrix(ncol=5,nrow=100))
  
  for (j in 1:100){
    set.seed(sample(1:10000, 1))
    folds <- createFolds(hum_scores_regression$group, k = 5)
    for (i  in 1:5){
      training <- hum_scores_regression[-folds[[i]],] # 67 training samples
      testing  <- hum_scores_regression[folds[[i]],] #21
      X_sel<-training[,-c(6)]
      y<-(as.numeric(factor(training$group, levels = c("LTBI","TB")))-1)
      logistic_reg <- train(
        x=X_sel,
        y=y,
        method = "glm",
        family=binomial('logit'), 
        maxit=100,
        trControl = trainControl(method = "cv", number = 5))
      classes <- predict(logistic_reg, newdata = testing[,-c(6)])
      testing$group<- (as.numeric(factor(testing$group, levels = c("LTBI","TB")))-1)
      auc_score <-  pROC::roc(group ~ as.numeric(classes), data = testing)
      auc_mouse[j,i]<-auc_score[["auc"]]
    }
  }
  
  results_fullmodel_auc<-data.frame(variable="full model",value = rowMeans(auc_mouse))
 
  ## Comparison of full model to a model with five random mPCs
  # This loop with run 100 rounds of 5-fold CV to predict disease phenotype from the model including the 5 randomly selected mPCs
  
  auc_rand_PCs<-data.frame(matrix(ncol=5,nrow=100))
  
  for (j in 1:100){
    set.seed(sample(1:10000, 1))
    folds <- createFolds(hum_scores_regression$group, k = 5)
    rand_PCs<-sample(c(1:88), 5, replace = FALSE)
    rand_PCs[6]<-dim(human_projected)[2]
    for (i  in 1:5){
      training <- human_projected[-folds[[i]],rand_PCs] # 67 training samples
      testing  <- human_projected[folds[[i]],rand_PCs]
      X_sel<-training[,-c(6)]
      y<-(as.numeric(factor(training$group, levels = c("LTBI","TB")))-1)
      logistic_rand<- train(
        x=X_sel,
        y=y, 
        method = "glm",
        family=binomial('logit'), 
        maxit=100,
        trControl = trainControl(method = "cv", number = 5))
      classes <- predict(logistic_rand, newdata = testing[,-c(6)])
      testing$group<- (as.numeric(factor(testing$group, levels = c("LTBI","TB")))-1)
      auc_rand <-  pROC::roc(group ~ as.numeric(classes), data = testing)
      auc_rand_PCs[j,i]<-auc_rand[["auc"]]
    }
  }
  results_random_PCs_auc<-data.frame(variable="random PCs",value = rowMeans(auc_rand_PCs))
  
  ## Comparison of full model to a model with shuffled human phenotype labels
  # This loop with run 100 rounds of 5-fold CV to predict disease phenotype from the model including the 5 mPCs we identified as translatable components (TCs) but with randomly shuffled human phenotype labels
  
  auc_shuffled_labels<-data.frame(matrix(ncol=5,nrow=100))
  hum_scores_shuffled <- human_projected[, c(PC_indices, dim(human_projected)[2])]
  
  for (j in 1:100){
    set.seed(sample(1:10000, 1))
    hum_scores_shuffled$group<-sample(hum_scores_shuffled$group)
    folds <- createFolds(hum_scores_shuffled$group, k = 5)
    for (i  in 1:5){
      training <- hum_scores_shuffled[-folds[[i]],] # 67 training samples
      testing  <- hum_scores_shuffled[folds[[i]],]
      X_sel<-training[,-c(6)]
      y<-(as.numeric(factor(training$group, levels = c("LTBI","TB")))-1)
      logistic_shuffled <- train(
        x=X_sel, 
        y=y, 
        method = "glm",
        family=binomial('logit'), 
        maxit=100,
        trControl = trainControl(method = "cv", number = 5))
      classes <- predict(logistic_shuffled, newdata = testing[,-c(6)])
      testing$group<- (as.numeric(factor(testing$group, levels = c("LTBI","TB")))-1)
      auc_shuffled <-  pROC::roc(group ~ as.numeric(classes), data = testing)
      auc_shuffled_labels[j,i]<-auc_shuffled[["auc"]]
    }
  }
  results_shuffled_labels_auc<-data.frame(variable="shuffled labels",value = rowMeans(auc_shuffled_labels))
  
  full_results_auc<-rbind(results_fullmodel_auc,results_random_PCs_auc)
  full_results_auc<-rbind(full_results_auc,results_shuffled_labels_auc)
  
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }
  
  setEPS()
  postscript(paste(plot_name))
  Figure_4F<-ggplot(full_results_auc, aes(x=variable, y=value)) + geom_violin(trim = FALSE)+stat_summary(fun.data=data_summary)+theme_classic()+geom_signif(comparisons = list(c("full model", "random PCs"),c("full model", "shuffled labels")))
  print(Figure_4F)
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

mouse_metadata <- read.csv(file.path(directory, files$mouse_metadata), header = TRUE)
row.names(mouse_metadata) <- mouse_metadata$Run
mouse_counts <- read.csv(file.path(directory, files$mouse_counts), header = TRUE, row.names = 1)
human_metadata <- read.csv(file.path(directory, files$human_metadata), header = TRUE)
row.names(human_metadata) <- human_metadata$Run
human_counts <- read.csv(file.path(directory, files$human_counts), header = TRUE, row.names = 1)
human_DEGs <- read.csv(file.path(directory, files$human_DEGs), header = TRUE)
  
### Generate Figure 4A
  
human_projected<-generate_figure4A(mouse_counts,human_counts,"~/Figure_4A.eps")

### Generate Figure 4B

variance_both_species<-generate_figure4B(mouse_counts,human_projected,"~/Figure_4B.eps")

### Generate Figure 4C

results<-generate_figure4C(human_projected,variance_both_species,"~/Figure_4C.eps")
list2env(results, envir = .GlobalEnv)

### Generate Figure 4D

generate_figure4D(human_projected,keep_PCs_variance,"~/Figure_4D_Top.eps","~/Figure_4D_Bottom.eps")

### Generate Figure 4E

generate_figure4E(human_projected,keep_PCs,"~/Figure_4E.eps")

### Generate Figure 4F

generate_figure4F(human_projected,keep_PCs,"~/Figure_4F.eps")








