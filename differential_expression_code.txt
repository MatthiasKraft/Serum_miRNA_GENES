### Code used for differential expression analysis ####

library("edgeR")
library(DESeq2)
library(pheatmap)
library("RColorBrewer")
library(ggplot2)

##—————# Count table cleaning & preparation—————
filename_phenotypes<-“path_to_file”

clean_mirna_data <-function(filepath_counts){
 #this function reads in a file (mirdeep2 quantifier output) and removes duplicate miRNAs (keeping the highest value), and performs stringent pre-filtering
library(data.table)
counts <- read.table(filepath,sep="\t",header=TRUE,row.names=NULL)
no_dups<-counts[which(counts$read_count == ave(counts$read_count,counts$miRNA,FUN=function(x) max(x))),]
no_dups<-no_dups[!duplicated(no_dups$miRNA),]
rownames(no_dups)<-no_dups$miRNA
t1<-no_dups[,!(colnames(no_dups) %like% 'norm')]
t1<-t1[,-(1:4)]
stringent_filtered<-t1[(rowMeans(t1>0)>=0.9),]
#filters out all miRNAs with 0 counts in more than 10% of the samples
return(stringent_filtered)}

#Start Analysis

countData<-clean_mirna_data(“countFile.txt”)
###match columns to phenodata
countData<-countData[,rownames(filenmame_phenotypes)]
colData3<-read.table(path,sep=" ",header=TRUE,row.names=1)
colData$condition<-factor(colData$condition,levels=c("Control","RAO”))
dds<-DESeqDataSetFromMatrix(countData = countData,colData = filename_phenotypes, design=~ABS.414nm+fam1+fam2+condition)
dds<-estimateSizeFactors(dds)
idx<-rowMeans(counts(dds8,normalized=TRUE))>13
dds<-DESeq(dds[idx,])
res<-results(dds,contrast = c("condition","RAO","Control”))
sig<-subset(res,res$padj<0.05)


# plot the count boxplots for the significant genes:
name="Boxplots_counts_"
ending=".pdf"
for( mirna in rownames(sig)){ #for all sign. DE. microRNAs do
  filename<-paste(name,mirna,ending,sep="")
  pdf(filename)
  mirna_counts<-plotCounts(dds, gene=mirna, intgroup="condition",returnData = TRUE)
  bp<-ggplot(mirna_counts, aes(x=condition, y=count)) + geom_boxplot(lwd=0.5,outlier.shape=NA) +labs(title=mirna)+ylab("Normalized read counts")+xlab("Condition")+theme(plot.title = element_text(hjust = 0.5,size=12,face="bold"))+theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold"))
  bp2<-bp+geom_point(position= position_jitter(width=0.2),size=1)
  ggsave(filename, plot=bp2, height=4, width=4.2, units="in", dpi=150)}




##########
# EdgeR

dge<-DGEList(counts=countData,group=colData3$condition)
des<-model.matrix(~condition +ABS.414nm, data = colData)
colnames(des)
dge<-calcNormFactors(dge)

#Estimate dispersions
dge<- estimateDisp(dge, des)

DESeqMeans<-as.data.frame(res$baseMean)
rownames(DESeqMeans)<-rownames(res)
colnames(DESeqMeans)<-"DESeq_baseMean"
#fit the model
fit <- glmFit(dge, des)
for (i in 2:3){
  LRT<-glmLRT(fit, coef = i)
  FDR<-p.adjust(LRT$table$PValue, method = "fdr")
  x<-cbind(LRT$table, FDR,DESeq_baseMean=DESeqMeans[rownames(LRT),])
  #Select significant DEGs
  if(i>2){
    x<-x[FDR<fdr_t,]
  }
  LRT<-x[order(x$FDR),]
  #Save the results
  write.table(LRT, file = paste(colnames(des)[i],"txt",sep="."),sep = "\t", quote = FALSE)
  
}

—#### sessionINfo()##########
#R version 3.3.1 (2016-06-21)
#Platform: x86_64-apple-darwin13.4.0 (64-bit)
#Running under: OS X 10.11.6 (El Capitan)

#locale:
#[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
#[8] methods   base     

#other attached packages:
 #[1] edgeR_3.14.0               limma_3.28.21             
 #[3] DESeq2_1.12.4              SummarizedExperiment_1.2.3
 #[5] Biobase_2.32.0             GenomicRanges_1.24.3      
 #[7] GenomeInfoDb_1.8.7         IRanges_2.6.1             
 #[9] S4Vectors_0.10.3           BiocGenerics_0.18.0       
#
#loaded via a namespace (and not attached):
 #[1] genefilter_1.54.2    locfit_1.5-9.1       splines_3.3.1       
 #[4] lattice_0.20-34      colorspace_1.3-2     htmltools_0.3.5     
 #[7] base64enc_0.1-3      survival_2.40-1      XML_3.98-1.5        
#[10] foreign_0.8-67       DBI_0.6              BiocParallel_1.6.6  
#[13] RColorBrewer_1.1-2   plyr_1.8.4           stringr_1.2.0       
#[16] zlibbioc_1.18.0      munsell_0.4.3        gtable_0.2.0        
#[19] htmlwidgets_0.8      memoise_1.0.0        latticeExtra_0.6-28 
#[22] knitr_1.15.1         geneplotter_1.50.0   AnnotationDbi_1.34.4
#[25] htmlTable_1.9        Rcpp_0.12.9          acepack_1.4.1       
#[28] xtable_1.8-2         backports_1.0.5      scales_0.4.1        
#[31] checkmate_1.8.2      Hmisc_4.0-2          annotate_1.50.1     
#[34] XVector_0.12.1       gridExtra_2.2.1      ggplot2_2.2.1       
#[37] digest_0.6.12        stringi_1.1.2        grid_3.3.1          
#[40] tools_3.3.1          bitops_1.0-6         magrittr_1.5        
#[43] lazyeval_0.2.0       RCurl_1.95-4.8       tibble_1.2          
#[46] RSQLite_1.1-2        Formula_1.2-1        cluster_2.0.5       
#[49] Matrix_1.2-8         data.table_1.10.4    assertthat_0.1      
#[52] rpart_4.1-10         nnet_7.3-12 

