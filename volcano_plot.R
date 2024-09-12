#GUKORA VOLCANO PLOT
#PART II: DUKORESHEJE DATA ZAWE
#volcano plot of differentially abundant genes
library(DESeq2)
library(pasilla)
library(EnhancedVolcano)
#Kubwira system aho data ziherereye
pasCts2 <- system.file("extdata2",
                       "alain_genes.tsv",
                       package="pasilla", mustWork=TRUE)
pasAnno2 <- system.file("extdata2",
                        "sample_par.csv",
                        package="pasilla", mustWork=TRUE)
################################
#Kuloading data muri R
cts2 <- as.matrix(read.csv(pasCts2,sep="\t",row.names="gene_name"))

#Kuloading data descriptions (metadata) muri R
coldata2 <- read.csv(pasAnno2, row.names=1)

#Guhitamo factors ziza kugenderwaho 
coldata2 <- coldata2[,c("Location","Substrate")]
rownames(coldata2)<-rownames(coldata2) 

#Guhuza amazina ya rows muri metadata na data
cts2 <- cts2[, rownames(coldata2)]
#Gukora a DESeq object
dds2 <- DESeqDataSetFromMatrix(countData = cts2,
                               colData = coldata2,
                               design = ~ Substrate)
#################################################################
#Gukuramo genes zifite counts nke cyane
keep2 <- rowSums(counts(dds2)) >= 2
dds2 <- dds2[keep2,]
##################################
#Gukora dataset ya fold changes na p-values
dds2 <- DESeq(dds2)
res2 <- results(dds2)
###############################
#Kureba uko bimeze na MA plot
summary(res2)
plotMA(res2, ylim=c(-2,2))

#Gushyira ku murongo results, no kuzisavinga 
resOrdered2 <- res2[order(res2$pvalue),]
write.csv(as.data.frame(resOrdered2),
          file="alain_results_genes.csv")
######################################################
#Gushushanya volcano plot
head(res2)
volcano_plot<-EnhancedVolcano(res2,
                              lab = rownames(res2),
                              x = 'log2FoldChange',
                              y = 'pvalue')

volcano_plot
