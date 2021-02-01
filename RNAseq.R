# load packages
library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(gridExtra)
library(DESeq2)
library(piano)
library(pheatmap)
library(genefilter)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(cowplot)
library(tidyverse)
library(dplyr)
library(hexbin)
library(reshape2)

# tmsg <- function(text = NULL) { return(paste0(text, " [", Sys.info()[['nodename']], ":", Sys.getpid(), "]")) }
tmsg <- function(text = NULL) { message(paste0(" [", Sys.info()[['nodename']], ":", Sys.getpid(), "] ", text)) }

Sys.setenv(LANGUAGE = "en") # Display error message in English
options(stringsAsFactors = FALSE) # Forbid chr convert to factor

# create a folder
if ( ! file.exists('~/Desktop/Renbo'))
  dir.create('~/Desktop/Renbo')

# set work direction
setwd("~/Desktop/Renbo")
load("rnaseq_rawdata.Rdata")
load("rnaseq_normdata.RData")
ls()

colnames(mat_cpm.norm)
colnames(mat_tpm.norm)
groups

################################################
###### convert from count to rlog matrix #######
######             by DESeq2             #######
################################################
# set put the group info
counts<- express[,groups$sampleID]
colData <- data.frame(row.names=colnames(counts), cell=colnames(counts),group=groups$Group)
# colData <- data.frame(row.names=colnames(counts), group=groups$Group)
dds <- DESeqDataSetFromMatrix(countData = round(counts,0),
                              colData = colData,
                              design = ~ group)

# change group levels
colData(dds)$group <- factor(colData(dds)$group,levels=c("EC","SC"))
# only keep genes that have counts higher than 0 in any sample
keep <- apply(counts(dds), 1, function(x) any(x >= 0)) 
dds <- dds[keep,]

# Remove transcripts do not show variance across samples
dds <- estimateSizeFactors(dds)
# sds <- rowSds(counts(dds, normalized = TRUE))
# sh <- shorth(sds)
# dds <- dds[sds >= sh,]


# PCA
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="group")
vstdata=assay(vsdata)

#variance stabilization
dds.norm <- varianceStabilizingTransformation(dds, blind=TRUE)

par(cex = 0.7)
n.sample=ncol(dds.norm)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(log2(counts(dds)+1), col = cols,main="expression value",las=2)
boxplot(vstdata, col = cols,main="expression value",las=2)

#how many genes left
dim(dds)

DEres <- list()
design(dds) <- ~ group
rnaRaw <- DESeq(dds, betaPrior = FALSE)

# fitting model
plotDispEsts(rnaRaw, main="Dispersion plot")

## group "EC","SC"
DEres[["EvsS"]] <- results(rnaRaw, contrast = c("group","EC","SC"))
DEres[["SvsE"]] <- results(rnaRaw, contrast = c("group","SC","EC"))

DEres[["EvsS.Shr"]] <- lfcShrink(rnaRaw, contrast = c("group","EC","SC"), res = DEres[["EvsS"]])
DEres[["SvsE.Shr"]] <- lfcShrink(rnaRaw, contrast = c("group","SC","EC"), res = DEres[["SvsE"]])

dat<- as.data.frame(DEres[["EvsS.Shr"]])
dat2<- as.data.frame(DEres[["SvsE.Shr"]])

dat<- dat[!is.na(dat$log2FoldChange),]
dat2<- dat2[!is.na(dat2$log2FoldChange),]
ids<- intersect(rownames(dat),rownames(dat2))
dat<- dat[ids,]
dat2<- dat2[ids,]
dim(dat)
dim(dat2)

# correlation
plot(dat$log2FoldChange,dat2$log2FoldChange)
cor.test(dat$log2FoldChange,dat2$log2FoldChange)


DEres_dt <- tibble(ID = rownames(dat), 
                   ECvsSC = 0,
                   SCvsEC = 0)


colnames(DEres_dt)[1:3]

DEres[["ECvsSC"]]<- dat
DEres[["SCvsEC"]]<- dat2

for (i in colnames(DEres_dt)[2:3]) {
  dat<- DEres[[i]]
  DEres_dt[(!is.na(dat$padj) & !is.na(dat$log2FoldChange)) & (dat[,"log2FoldChange"] > 1 & dat[,"padj"] < 0.05),i] = 1
  DEres_dt[(!is.na(dat$padj) & !is.na(dat$log2FoldChange)) & (dat[,"log2FoldChange"] < -1 & dat[,"padj"] < 0.05),i] = -1
  cat(i,"\n")
}

head(DEres_dt)

table(DEres_dt$ECvsSC)
table(DEres_dt$SCvsEC)

# Using MSigDB gene set collections
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")



names(gmts)<- c("hmark2gene","onco2gene","tf2gene","kegg2gene","mir2gene","pathway2gene","go2gene","allpathway2gene","chrom2gene","immune2gene")

library(clusterProfiler)
library(enrichplot)       
library(org.Hs.eg.db)
library(enrichplot)

geneset<- lapply(1:length(names(gmts)),function(x){
  x<- read.gmt(gmts[[names(gmts)[x]]])
  x <- x %>% tidyr::separate(ont, c("term"), "%")
  x <- x %>% dplyr::select(term, gene) #TERM2GENE
})

names(geneset)<- c("hmark2gene","onco2gene","tf2gene","kegg2gene","mir2gene","pathway2gene","go2gene","allpathway2gene","chrom2gene","immune2gene")

genelists<- list()
genelists[["EC_up"]]<- data.frame(DEres_dt[DEres_dt$ECvsSC == 1,"ID"])[,1]
genelists[["EC_down"]]<- data.frame(DEres_dt[DEres_dt$ECvsSC == -1,"ID"])[,1]
genelists[["SC_up"]]<- data.frame(DEres_dt[DEres_dt$SCvsEC == 1,"ID"])[,1]
genelists[["SC_down"]]<- data.frame(DEres_dt[DEres_dt$SCvsEC == -1,"ID"])[,1]

###################################################################################
## ++++++++++ Enrichment function for running all interested genesets ++++++++++ ##
###################################################################################
# Enrichment function
rEnrich <- function(geneList, term2gene, backgroundList){
  # enrichment function
  enrich<- enricher(gene= geneList, 
                    TERM2GENE = term2gene,
                    pvalueCutoff = 1,
                    universe = backgroundList,
                    pAdjustMethod = "BH",
                    minGSSize = 20, 
                    maxGSSize = 2000, 
                    qvalueCutoff = 1)
  enrich@result
}

# run Enrichment function
rEn = function(input,gmt,backgroundList) {
  Res <- list()
  for (i in names(input)) {
    geneList=input[[i]]
    resTab <- rEnrich(geneList=geneList, term2gene=gmt, backgroundList=backgroundList)
    Res[[i]] <- resTab
  }
  Res
}
####################################################################################

## run all geneset
# run Enrichment
# input is the fold change cutoff genelist
# backgroundList<- rownames(DEGres[[1]])
backgroundList<- DEres_dt$ID
test<- lapply(1:length(names(geneset)), function(x){
  dat<- rEn(genelists,geneset[[x]],backgroundList)
  return(dat)
})

names(test)<- names(geneset)

####################################################################################
tmp<- lapply(1:length(names(test)),function(x){
  tab<- test[[x]]
  re = do.call("rbind",tab)
  re$sample=rep(names(tab),sapply(tab, nrow))
  return(re)
})

names(tmp)<- names(test)

####################################################################################
nes_mat<- lapply(1:length(names(tmp)),function(x){
  tab<- tmp[[x]]
  es_mat <- tidyr::spread(tab[,c(2,9,10)],'sample',"Count")
  es_mat[is.na(es_mat)] <- 0
  rownames(es_mat)<- es_mat$Description
  es_mat$Description<- NULL
  
  tab$GeneRatio<- tab$GeneRatio <- sapply(tab$GeneRatio,function(x){as.numeric(as.character(strsplit(x,"/")[[1]][2]))})
  es_mat2 <- tidyr::spread(tab[,c(2,3,10)],'sample',"GeneRatio")
  for(i in 1:ncol(es_mat2)){
    es_mat2[,i][is.na(es_mat2[,i])] <- round(mean(es_mat2[,i], na.rm = TRUE))
  }
  rownames(es_mat2)<- es_mat2$Description
  es_mat2$Description<- NULL
  es_mat<- es_mat/es_mat2
  return(es_mat)
})

names(nes_mat)<- names(tmp)


nes_mat_padj<- lapply(1:length(names(tmp)),function(x){
  tab<- tmp[[x]]
  es_mat <- tidyr::spread(tab[,c(2,6,10)],'sample',"p.adjust")
  es_mat[is.na(es_mat)] <- 1
  rownames(es_mat)<- es_mat$Description
  es_mat$Description<- NULL
  return(es_mat)
})

names(nes_mat_padj)<- names(tmp)


