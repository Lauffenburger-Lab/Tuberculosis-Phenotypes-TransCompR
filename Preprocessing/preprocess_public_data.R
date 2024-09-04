
wdir <- # set your working directory

### Load Packages   
library(tximport)
library(DESeq2)
library(GEOquery)
library(biomaRt)
library(homologene)

### Define Functions
getSampData <- function(){
  gse_mouse <- getGEO('GSE137092',GSEMatrix=TRUE) # this is a subseries
  mouse_metadata <- pData(phenoData(gse_mouse[[1]]))[,c(21,22)]
  colnames(mouse_metadata)[2] <- c("GRA")
  strSamp <- strsplit(mouse_metadata$description,"_")
  mouse_metadata$mtb_strain <- unlist(lapply(strSamp, `[[`, 1))
  mouse_metadata$dose <- unlist(lapply(strSamp,"[[",3))
  mouse_metadata$mouse_strain <- unlist(lapply(strSamp,"[[",4))
  
  gse_London <- getGEO('GSE107991',GSEMatrix=TRUE)
  london_metadata <- pData(phenoData(gse_London[[1]]))[, c(1,8,10,11)] # title, full, tissue, group,
  colnames(london_metadata) <- c("title","description","tissue","group")
  
  gse_CapeTown <- getGEO('GSE107992',GSEMatrix=TRUE)
  capetown_metadata <- pData(phenoData(gse_CapeTown[[1]]))[, c(1,8,10,11)] # title, full, tissue, group
  colnames(capetown_metadata) <- c("title","description","tissue","group")
  
  sampList <- list("mouse_metadata"=mouse_metadata, "london_metadata"=london_metadata, "capetown_metadata"=capetown_metadata)
  return(sampList)
}
processData <- function(dir,sampleData,SRRGSM){
  
  # Function that processes raw fastq files and outputs mean-centered, z-scored data
  # Input:
  #   - dir: path for overall directory for dataset (includes fastq and RSEM folders)
  #   - sampleData: dataframe of sample information (e.g. TB status) for specified dataset
  #   - SRRGSM: dataframe of mappings between SRR and GSM ids
  
  fqFiles <- list.files(file.path(dir,"fastq/"))
  rsemFiles <- list.files(file.path(dir, "RSEM/"), pattern="\\.genes.results")
  SRRs <- substr(rsemFiles,1,11)
  names(rsemFiles) <- SRRs
  sampleData <- merge(sampleData, SRRGSM[,c("Sample.Name","Run")], by.x=0, by.y="Sample.Name")
  rownames(sampleData) <- sampleData$Run
  
  # tximport to get counts data from raw fastq
  txi.rsem <- tximport(file.path(file.path(dir,"RSEM/"),rsemFiles), type="rsem", txIn=FALSE, txOut=FALSE)
  txi.rsem$length[txi.rsem$length == 0] <- 1
  colnames(txi.rsem$counts) <- names(rsemFiles)
  
  # DESeq2 to variance-stabilize transform data
  dds <- DESeqDataSetFromTximport(txi.rsem, sampleData, ~1)
  dds <- estimateSizeFactors(dds)
  dds <- dds[rowSums(counts(dds))>=10,]
  
  vsd <- assay(vst(dds))
  vsd_norm <- data.frame(t(apply(vsd, 1, scale))) # center and z-score over all samples (ie. across row)
  colnames(vsd_norm) <- rownames(sampleData)
  data <- list(vsd_norm,sampleData)
  return(data)
  
}
homologFilter <- function(humEnsembl, musEnsembl,host_name) {
  ##### Takes human and mouse ENSEMBL symbols
  ##### Filters out genes that don't have a 1:1 homolog between the species
  ##### Returns filtered human and mouse gene lists with only 1:1 homologs
  
  ### convert ENSEMBL -> HGNC 
  ensemblH <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host=host_name)
  resH <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                filters="ensembl_gene_id",
                values=humEnsembl,
                mart=ensemblH)
  
  ### generate HGNC:MGI homologs, remove duplicates
  homologs <- human2mouse(unique(resH$hgnc_symbol))
  dupsH <- unique(homologs$humanGene[duplicated(homologs$humanGene)])
  dupsM <- unique(homologs$mouseGene[duplicated(homologs$mouseGene)])
  keep <- homologs[which(!((homologs$humanGene %in% dupsH)|(homologs$mouseGene %in% dupsM))),]
  
  ### generate MGI:ENSEMBL dataframe for mouse genes, remove duplicates
  ensemblM <- useMart("ensembl", dataset="mmusculus_gene_ensembl",host=host_name)
  resM <- getBM(attributes=c("mgi_symbol","ensembl_gene_id"), 
                filters="ensembl_gene_id", 
                values=musEnsembl, 
                mart=ensemblM) # generate  MGI:ENSEMBL df using homologs list of mouse genes
  dupsMGI <- unique(resM$mgi_symbol[duplicated(resM$mgi_symbol)])
  keep2 <- keep[!(keep$mouseGene %in% dupsMGI),] # don't include mouse genes with duplicate MGI:ENSEMBL encoding in homologs list
  
  finalList <- merge(keep2, resM, by.x="mouseGene", by.y="mgi_symbol")
  finalList <- merge(finalList, resH, by.x="humanGene", by.y="hgnc_symbol")
  colnames(finalList)[c(5,6)] <- c("mouseEns", "humEns")
  
  return(finalList)
  
}

## Load and format all datasets 

sampInfo <- getSampData()
list2env(sampInfo, envir = .GlobalEnv)

### Mouse data formatting and normalization
SRRGSM_mouse <- read.csv(file.path(wdir,"MT_SraRunTable_bloodSeq.txt"))
fqMus <- list.files(file.path(wdir,"MT_data_bloodSeq/fastq/"))
rsemMus <- list.files(file.path(wdir, "MT_data_bloodSeq/RSEM"), pattern="\\.genes.results")
SRRs <- substr(rsemMus,1,11)
names(rsemMus) <- SRRs
mouse_metadata <- merge(mouse_metadata, SRRGSM_mouse[,c("Sample.Name","Run")], by.x=0, by.y="Sample.Name")
rownames(mouse_metadata) <- mouse_metadata$Run

txi.rsem_mouse<- tximport(file.path(file.path(wdir,"MT_data_bloodSeq/RSEM/"),rsemMus), type="rsem", txIn=FALSE, txOut=FALSE)
txi.rsem_mouse$length[txi.rsem_mouse$length == 0] <- 1
colnames(txi.rsem_mouse$counts) <- names(rsemMus)

raw_mouse_counts<-txi.rsem_mouse[["counts"]]
write.csv(raw_mouse_counts,"~/raw_mouse_counts.csv")

dds_mouse <- DESeqDataSetFromTximport(txi.rsem_mouse, mouse_metadata, ~1)
dds_mouse <- estimateSizeFactors(dds_mouse)
dds_mouse <- dds_mouse[rowSums(counts(dds_mouse))>=10,]
vsd_mouse <- assay(vst(dds_mouse))
vsd_mouse_normalized <- data.frame(apply(vsd_mouse, 1, scale))
rownames(vsd_mouse_normalized) <- rownames(mouse_metadata)
write.csv(vsd_mouse_normalized,"~/normalized_mouse_counts.csv")

### Human data formatting and normalization
SRRGSM_London <- read.csv(file.path(wdir,"SinghaniaLon_SraRunTable_bloodSeq.txt"))
SRRGSM_CapeTown <- read.csv(file.path(wdir,"SinghaniaSA_SraRunTable_bloodSeq.txt"))
SRRGSM_human <- rbind(SRRGSM_London[,c("Run","Sample.Name")], SRRGSM_CapeTown[,c("Run","Sample.Name")])

human_metadata <- rbind(london_metadata,capetown_metadata)

fqHum <- list.files(file.path(wdir,"Singhania_data/fastq/"))
rsemHum <- list.files(file.path(wdir, "Singhania_data/RSEM"), pattern="\\.genes.results")
SRRs <- substr(rsemHum,1,10) 
names(rsemHum) <- SRRs
keep <- SRRGSM_human[SRRGSM_human$Run %in% SRRs,] # can only use SRR/GSMs that I have RSEM files for
human_metadata <- merge(human_metadata, keep, by.x=0, by.y="Sample.Name")
rownames(human_metadata) <- human_metadata$Run
human_metadata$cohort <- sapply(strsplit(human_metadata$description, "_"), function(x) x[2]) 

london_SRRs <- human_metadata$Run[human_metadata$Row.names %in% rownames(london_metadata)]
txi.rsem_london <- tximport(file.path(paste0(wdir,"/Singhania_data/RSEM"),rsemHum[london_SRRs]), type="rsem", txIn=FALSE, txOut=FALSE)
txi.rsem_london$length[txi.rsem_london$length == 0] <- 1
colnames(txi.rsem_london$counts) <- london_SRRs

capetown_SRRs <- human_metadata$Run[human_metadata$Row.names %in% rownames(capetown_metadata)]
txi.rsem_capetown <- tximport(file.path(paste0(wdir,"/Singhania_data/RSEM"),rsemHum[capetown_SRRs]), type="rsem", txIn=FALSE, txOut=FALSE)
txi.rsem_capetown$length[txi.rsem_capetown$length == 0] <- 1
colnames(txi.rsem_capetown$counts) <- capetown_SRRs 

dds_london <- DESeqDataSetFromTximport(txi.rsem_london, human_metadata[london_SRRs,], ~1)
dds_london <- estimateSizeFactors(dds_london)
dds_london <- dds_london[rowSums(counts(dds_london))>=10,]

dds_capetown <- DESeqDataSetFromTximport(txi.rsem_capetown, human_metadata[capetown_SRRs,], ~1)
dds_capetown <- estimateSizeFactors(dds_capetown)
dds_capetown <- dds_capetown[rowSums(counts(dds_capetown))>=10,]

# filter to only include genes present in both cohorts
keepHumGenes <- intersect(rownames(dds_london),rownames(dds_capetown))
dds_london <- dds_london[keepHumGenes,]
dds_capetown <- dds_capetown[keepHumGenes,]

# vst transform, mean center, and scale
vsd_london <- assay(vst(dds_london)) # vst transform
vsd_london_normalized <- data.frame(apply(vsd_london, 1, scale)) # mean-centering and z-scoring across samples (ie over row)
rownames(vsd_london_normalized) <- rownames(human_metadata[london_SRRs,])
vsd_london_normalized$cohort <- 'London'
vsd_london_normalized$group <- human_metadata[london_SRRs,'group']

vsd_capetown <- assay(vst(dds_capetown)) # vst transform
vsd_capetown_normalized <- data.frame(apply(vsd_capetown, 1, scale)) # mean-centering and z-scoring across samples (ie over row)
rownames(vsd_capetown_normalized) <- rownames(human_metadata[capetown_SRRs,])
vsd_capetown_normalized$cohort <- 'CapeTown'
vsd_capetown_normalized$group <- human_metadata[capetown_SRRs,'group']

human_normalized <- rbind(vsd_london_normalized, vsd_capetown_normalized)

### Human DEG Analysis 

dds_londonDE <- DESeqDataSetFromTximport(txi.rsem_london, human_metadata[london_SRRs,], ~group)
dds_londonDE <- estimateSizeFactors(dds_londonDE)
dds_londonDE <- dds_londonDE[rowSums(counts(dds_londonDE))>=10,]
dds_londonDE <- DESeq(dds_londonDE)
londonDE <- results(dds_londonDE)
londonDE <- londonDE[complete.cases(londonDE),]
londonDE_filt <- data.frame(londonDE[abs(londonDE$log2FoldChange)>0.5 & londonDE$padj<0.1,])
DEGs_london <- rownames(londonDE_filt)

dds_capetownDE <- DESeqDataSetFromTximport(txi.rsem_capetown, human_metadata[capetown_SRRs,], ~group)
dds_capetownDE <- estimateSizeFactors(dds_capetownDE)
dds_capetownDE <- dds_capetownDE[rowSums(counts(dds_capetownDE))>=10,]
dds_capetownDE <- DESeq(dds_capetownDE)
capetownDE <- results(dds_capetownDE)
capetownDE <- capetownDE[complete.cases(capetownDE),]
capetownDE_filt <- data.frame(capetownDE[abs(capetownDE$log2FoldChange)>0.5 & capetownDE$padj<0.1,])
DEGs_capetown<- rownames(capetownDE_filt)

DEGs <- union(DEGs_london,DEGs_capetown)


setEPS()
postscript("Figure_S2A.eps")
londonDE$pch <- 19
londonDE$pch[londonDE$pvalue ==0] <- 6
cols <- densCols(londonDE$log2FoldChange, londonDE$padj)
plot(londonDE$log2FoldChange,-log10(londonDE$padj),pch=londonDE$pch, col=cols,cex=0.5,xlim=c(-6, 6))
abline(v=0)
abline(v=c(-0.5,0.5))
abline(h=-log10(0.1))
dev.off()

setEPS()
postscript("Figure_S2B.eps")
capetownDE$pch <- 19
capetownDE$pch[capetownDE$pvalue ==0] <- 6
cols <- densCols(capetownDE$log2FoldChange, capetownDE$padj)
plot(capetownDE$log2FoldChange,-log10(capetownDE$padj),pch=capetownDE$pch, col=cols,cex=0.5,xlim=c(-6, 6))
abline(v=0)
abline(v=c(-0.5,0.5))
abline(h=-log10(0.1))
dev.off()

### Match human DEGs with mouse homologs

homologs <- homologFilter(colnames(human_normalized), colnames(vsd_mouse_normalized),"https://oct2022.archive.ensembl.org")
keepHomologs <- homologs[which(homologs$humEns %in% DEGs),c('mouseEns','humEns')]

human_normalized_homologues<-human_normalized[,colnames(human_normalized) %in% keepHomologs$humEns]
mouse_normalized_homologues<-vsd_mouse_normalized[,colnames(vsd_mouse_normalized) %in% keepHomologs$mouseEns]
write.csv(mouse_normalized_homologues,"mouse_norm.csv")
write.csv(human_normalized_homologues,"human_norm.csv")



