---
title: "DESeq2"
author: "Matthew Oldach"
date: "September 21, 2015"
output: html_document
---

Basically, I'm just following along with the DESeq2 Manual (http://www.bioconductor.org/packages/2.13/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf), the packages Vignette, and from a 2014 MEGA workshop (https://github.com/z0on, https://github.com/rachelwright8/Ahya-White-Syndromes).

## Guide to using the DESeq2 analyis package

### Read in counts and simplify sample names
#### Each gene is designated as an isogroup we only see isogroup0, isogroup1, isogroup10, isogroup100, etc. in the head(). Counts were normalized to median (or mean) counts not control genes due to legit up-regulation of abundant genes

```r
countData <- read.table("allcounts.txt")
names(countData)=sub(".sam.counts","",names(countData))
names(countData)
```

```
##  [1] "F1D3Q" "F1DFM" "F1DNM" "F1N1Q" "F1N3Q" "F1NFM" "F1NNM" "F2D3Q"
##  [9] "F2DFM" "F2DNM" "F2N1Q" "F2N3Q" "F2NFM" "F2NNM" "F3DFM" "F3DNM"
## [17] "F3N1Q" "F3N3Q" "F3NFM" "F3NNM"
```

```r
head(countData)
```

```
##               F1D3Q F1DFM F1DNM F1N1Q F1N3Q F1NFM F1NNM F2D3Q F2DFM F2DNM
## isogroup1         6     9     3     2     0     3     7     7     5     4
## isogroup10        0     0     1     0     0     0     0     1     1     0
## isogroup100       7    18     1     3     9    11    10     8     7     7
## isogroup1000      2    27     6    21     5    32    14     1    17    23
## isogroup10000     3     0    46     0     0     0     0     0     1     2
## isogroup10001    53   164    19    61   170    87    61    27    34    44
##               F2N1Q F2N3Q F2NFM F2NNM F3DFM F3DNM F3N1Q F3N3Q F3NFM F3NNM
## isogroup1         5     3     4     2     2     0     4     3     1     0
## isogroup10        2     0     0     0     0     0     0     0     0     0
## isogroup100       4     8     5     4     6     4     6     1     3     1
## isogroup1000     26     3     5    26     4    14     7     3    13    10
## isogroup10000    60     0     1     0     0     0     1     3     4   143
## isogroup10001    28   116    64    37    17    51     6    21    13    10
```

```r
write.csv(countData, "countData.csv", quote=F)
totalCounts=colSums(countData)
```
#### Total counts per sample 

```r
totalCounts
```

```
##   F1D3Q   F1DFM   F1DNM   F1N1Q   F1N3Q   F1NFM   F1NNM   F2D3Q   F2DFM 
## 1453783 2536648 1155904 1516181 2422117 2164573 1965200  838351 1384846 
##   F2DNM   F2N1Q   F2N3Q   F2NFM   F2NNM   F3DFM   F3DNM   F3N1Q   F3N3Q 
## 1623607 1534203 1765733 1143139 1436923  871014 1586063  765926  995658 
##   F3NFM   F3NNM 
##  678658 1151182
```

```r
min(totalCounts) 
```

```
## [1] 678658
```

```r
max(totalCounts) 
```

```
## [1] 2536648
```

```r
mean(totalCounts) 
```

```
## [1] 1449485
```

### Construct conditions

```r
indiv= c("F1", "F1", "F1", "F1", "F1", "F1", "F1", "F2", "F2", "F2", "F2", "F2", "F2", "F2", "F3", "F3", "F3","F3", "F3")  
lunar=c("3Q", "FM", "NM", "1Q","3Q", "FM", "NM", "3Q", "FM", "NM", "1Q", "3Q", "FM", "NM", "FM", "NM", "1Q", "3Q", "FM")
dayTime=c("day","day","day", "night", "night", "night", "night", "day","day","day", "night", "night", "night", "night", "day","day", "night", "night", "night")
g=data.frame(indiv, lunar, dayTime)
g
```

```
##    indiv lunar dayTime
## 1     F1    3Q     day
## 2     F1    FM     day
## 3     F1    NM     day
## 4     F1    1Q   night
## 5     F1    3Q   night
## 6     F1    FM   night
## 7     F1    NM   night
## 8     F2    3Q     day
## 9     F2    FM     day
## 10    F2    NM     day
## 11    F2    1Q   night
## 12    F2    3Q   night
## 13    F2    FM   night
## 14    F2    NM   night
## 15    F3    FM     day
## 16    F3    NM     day
## 17    F3    1Q   night
## 18    F3    3Q   night
## 19    F3    FM   night
```

```r
colData<- g
```
### Make big dataframe

```r
dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~lunar+dayTime) 
```

```
## Error in eval(expr, envir, enclos): could not find function "DESeqDataSetFromMatrix"
```

## Sample Outlier Detection With arrayQualityMetrics (MEGA2014)

```r
vsd=varianceStabilizingTransformation(dds, blind=TRUE)
```

```
## Error in eval(expr, envir, enclos): could not find function "varianceStabilizingTransformation"
```

```r
head(assay(vsd),10)
```

```
## Error in head(assay(vsd), 10): could not find function "assay"
```

```r
e=ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
```

```
## Error in eval(expr, envir, enclos): could not find function "ExpressionSet"
```

```r
arrayQualityMetrics(e, intgroup=c("indiv", "lunar", "dayTime"), force=T)
```

```
## Error in eval(expr, envir, enclos): could not find function "arrayQualityMetrics"
```
No outliers detected, moving on

# Differential Expression Analysis (DESeq2) 

```r
dds <- DESeq(dds)
```

```
## Error in eval(expr, envir, enclos): could not find function "DESeq"
```

```r
res<-results(dds)
```

```
## Error in eval(expr, envir, enclos): could not find function "results"
```

```r
summary(res)
```

```
## Error in summary(res): object 'res' not found
```

```r
head(res)
```

```
## Error in head(res): object 'res' not found
```

```r
mcols(res,use.names=TRUE)
```

```
## Error in eval(expr, envir, enclos): could not find function "mcols"
```

```r
write.csv(as.data.frame(res),file="results_DESeq2.csv")
```

```
## Error in as.data.frame(res): object 'res' not found
```

# Visulizations 

#### Create a MA plot showing DEGs (MEGA2014)
Lets plot the log2 fold changes against the mean normalised counts, colouring in red those genes that are significatn at 10% FDR

```r
plotMA(res, alpha=0.01)
```

```
## Error in eval(expr, envir, enclos): could not find function "plotMA"
```

#### plotDispEsts(dds) (MEGA2014)

```r
plotDispEsts(dds, main="Dispersion plot", ylim = c(1e-6, 1e1) )
```

```
## Error in eval(expr, envir, enclos): could not find function "plotDispEsts"
```

The level of dispersion is related to the biological variation seen in each treatment. As read count increases dispersion decreases, as we would expect. 

However, a lot of those black dot values below the red fitted line are probably underestimates of dispersion based on the small samples sizes used (only 3 replicates values measured in this example). So, to be conservative DESeq moves all these values below the red line up to the fitted value, BUT keeps all the empirical values above the red line, even if some of these are over estimates of dispersion.

#### Plot a histogram of unadjusted p-values after filtering to look at their distribution (MEGA2014)

```r
hist(res$pvalue, breaks=100, col="skyblue",border="slateblue",main="unadjustedP-value Histogram")
```

```
## Error in hist(res$pvalue, breaks = 100, col = "skyblue", border = "slateblue", : object 'res' not found
```
The pileups in the right hand side of distribution are from very low count genes that are given discrete values.


#### Regularized log transformation for clustering/heatmaps, etc

```r
rld <- rlogTransformation(dds, blind=TRUE)
```

```
## Error in eval(expr, envir, enclos): could not find function "rlogTransformation"
```

```r
# we choose blind so that the initial conditions setting does not influence the outcome, ie we want to see if the conditions cluster based purely on the individual 
head(assay(rld))
```

```
## Error in head(assay(rld)): could not find function "assay"
```

```r
hist(assay(rld))
```

```
## Error in hist(assay(rld)): could not find function "assay"
```
#### 

```r
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
```

```
## Error in is.data.frame(x): could not find function "counts"
```

```r
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1), ylim = c(0,2.5))
```

```
## Error in eval(expr, envir, enclos): could not find function "meanSdPlot"
```

```r
meanSdPlot(assay(rld[notAllZero,]), ylim = c(0,2.5))
```

```
## Error in eval(expr, envir, enclos): could not find function "meanSdPlot"
```

```r
meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,2.5))
```

```
## Error in eval(expr, envir, enclos): could not find function "meanSdPlot"
```

```r
qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
```

```
## Error in quantile(res$baseMean[res$baseMean > 0], 0:7/7): object 'res' not found
```

```r
# "cut" the genes into the bins
bins <- cut( res$baseMean, qs )
```

```
## Error in cut(res$baseMean, qs): object 'res' not found
```

```r
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
```

```
## Error in paste0("~", round(0.5 * qs[-1] + 0.5 * qs[-length(qs)])): object 'qs' not found
```

```r
# calculate the ratio of ?p? values less than .01 for each bin
ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
```

```
## Error in tapply(res$pvalue, bins, function(p) mean(p < 0.01, na.rm = TRUE)): object 'bins' not found
```

```r
# plot these ratios
barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")
```

```
## Error in barplot(ratios, xlab = "mean normalized count", ylab = "ratio of small $p$ values"): object 'ratios' not found
```
## Use contrasts to compare levels of lunar period and the dayTime

```r
#------full moon to new moon
resAH=results(dds, contrast=c("lunar", "FM", "NM"))
```

```
## Error in eval(expr, envir, enclos): could not find function "results"
```

```r
table(resAH$padj<0.1)
```

```
## Error in table(resAH$padj < 0.1): object 'resAH' not found
```

```r
table(resAH$pvalue<0.05)
```

```
## Error in table(resAH$pvalue < 0.05): object 'resAH' not found
```

```r
#------full moon to 1Q
resDH=results(dds, contrast=c("lunar","FM","1Q"))
```

```
## Error in eval(expr, envir, enclos): could not find function "results"
```

```r
table(resDH$padj<0.1)
```

```
## Error in table(resDH$padj < 0.1): object 'resDH' not found
```

```r
table(resDH$pvalue<0.05)
```

```
## Error in table(resDH$pvalue < 0.05): object 'resDH' not found
```

```r
#------full moon to 3Q
resDA=results(dds, contrast=c("lunar","FM","3Q"))
```

```
## Error in eval(expr, envir, enclos): could not find function "results"
```

```r
table(resDA$padj<0.1)
```

```
## Error in table(resDA$padj < 0.1): object 'resDA' not found
```

```r
table(resDA$pvalue<0.05)
```

```
## Error in table(resDA$pvalue < 0.05): object 'resDA' not found
```

```r
#------new moon to 3Q
resNM=results(dds, contrast=c("lunar","NM","3Q"))
```

```
## Error in eval(expr, envir, enclos): could not find function "results"
```

```r
table(resNM$padj<0.1)
```

```
## Error in table(resNM$padj < 0.1): object 'resNM' not found
```

```r
table(resNM$pvalue<0.05)
```

```
## Error in table(resNM$pvalue < 0.05): object 'resNM' not found
```

```r
#------new moon to 1Q
resDD=results(dds, contrast=c("lunar","NM","1Q"))
```

```
## Error in eval(expr, envir, enclos): could not find function "results"
```

```r
table(resDD$padj<0.1)
```

```
## Error in table(resDD$padj < 0.1): object 'resDD' not found
```

```r
table(resDD$pvalue<0.05)
```

```
## Error in table(resDD$pvalue < 0.05): object 'resDD' not found
```

```r
#-----1Q to 3Q
resDV=results(dds, contrast=c("lunar","1Q","3Q"))
```

```
## Error in eval(expr, envir, enclos): could not find function "results"
```

```r
table(resDV$padj<0.1)
```

```
## Error in table(resDV$padj < 0.1): object 'resDV' not found
```

```r
table(resDV$pvalue<0.05)
```

```
## Error in table(resDV$pvalue < 0.05): object 'resDV' not found
```

```r
#------day to night
resDN=results(dds,contrast=c("dayTime","day","night"))
```

```
## Error in eval(expr, envir, enclos): could not find function "results"
```

```r
table(resDN$padj<0.1)
```

```
## Error in table(resDN$padj < 0.1): object 'resDN' not found
```

```r
table(resDN$pvalue<0.05)
```

```
## Error in table(resDN$pvalue < 0.05): object 'resDN' not found
```

```r
save(dds, res, resAH, resDH, resDA, resNM, resDD, resDV, resDN, file="DESeq2Results.Rdata")
```

```
## Error in save(dds, res, resAH, resDH, resDA, resNM, resDD, resDV, resDN, : objects 'dds', 'res', 'resAH', 'resDH', 'resDA', 'resNM', 'resDD', 'resDV', 'resDN' not found
```
## Exporting and log GO

```r
vsd=getVarianceStabilizedData(dds) 
```

```
## Error in eval(expr, envir, enclos): could not find function "getVarianceStabilizedData"
```

```r
head(vsd)
```

```
## Error in head(vsd): object 'vsd' not found
```

```r
colnames(vsd)=paste(g$indiv, g$lunar, g$dayTime, sep="")
```

```
## Error in colnames(vsd) = paste(g$indiv, g$lunar, g$dayTime, sep = ""): object 'vsd' not found
```

```r
###--------------get pvals from contrasts
head(resDH)
```

```
## Error in head(resDH): object 'resDH' not found
```

```r
valsDH=cbind(resDH$pvalue, resDH$padj)
```

```
## Error in cbind(resDH$pvalue, resDH$padj): object 'resDH' not found
```

```r
head(valsDH)
```

```
## Error in head(valsDH): object 'valsDH' not found
```

```r
colnames(valsDH)=c("pvalDH", "padjDH")
```

```
## Error in colnames(valsDH) = c("pvalDH", "padjDH"): object 'valsDH' not found
```

```r
head(resAH)
```

```
## Error in head(resAH): object 'resAH' not found
```

```r
valsAH=cbind(resAH$pvalue, resAH$padj)
```

```
## Error in cbind(resAH$pvalue, resAH$padj): object 'resAH' not found
```

```r
head(valsAH)
```

```
## Error in head(valsAH): object 'valsAH' not found
```

```r
colnames(valsAH)=c("pvalAH", "padjAH")
```

```
## Error in colnames(valsAH) = c("pvalAH", "padjAH"): object 'valsAH' not found
```

```r
head(resDA)
```

```
## Error in head(resDA): object 'resDA' not found
```

```r
valsDA=cbind(resDA$pvalue, resDA$padj)
```

```
## Error in cbind(resDA$pvalue, resDA$padj): object 'resDA' not found
```

```r
head(valsDA)
```

```
## Error in head(valsDA): object 'valsDA' not found
```

```r
colnames(valsDA)=c("pvalDA", "padjDA")
```

```
## Error in colnames(valsDA) = c("pvalDA", "padjDA"): object 'valsDA' not found
```

```r
head(resNM)
```

```
## Error in head(resNM): object 'resNM' not found
```

```r
valsNM=cbind(resNM$pvalue, resNM$padj)
```

```
## Error in cbind(resNM$pvalue, resNM$padj): object 'resNM' not found
```

```r
head(valsNM)
```

```
## Error in head(valsNM): object 'valsNM' not found
```

```r
colnames(valsNM)=c("pvalNM", "padjNM")
```

```
## Error in colnames(valsNM) = c("pvalNM", "padjNM"): object 'valsNM' not found
```

```r
head(resDD)
```

```
## Error in head(resDD): object 'resDD' not found
```

```r
valsDD=cbind(resDD$pvalue, resDD$padj)
```

```
## Error in cbind(resDD$pvalue, resDD$padj): object 'resDD' not found
```

```r
head(valsDD)
```

```
## Error in head(valsDD): object 'valsDD' not found
```

```r
colnames(valsDD)=c("pvalDD", "padjDD")
```

```
## Error in colnames(valsDD) = c("pvalDD", "padjDD"): object 'valsDD' not found
```

```r
head(resDV)
```

```
## Error in head(resDV): object 'resDV' not found
```

```r
valsDV=cbind(resDV$pvalue, resDV$padj)
```

```
## Error in cbind(resDV$pvalue, resDV$padj): object 'resDV' not found
```

```r
head(valsDV)
```

```
## Error in head(valsDV): object 'valsDV' not found
```

```r
colnames(valsDV)=c("pvalDV", "padjDV")
```

```
## Error in colnames(valsDV) = c("pvalDV", "padjDV"): object 'valsDV' not found
```

```r
head(resDN)
```

```
## Error in head(resDN): object 'resDN' not found
```

```r
valsDN=cbind(resDN$pvalue, resDN$padj)
```

```
## Error in cbind(resDN$pvalue, resDN$padj): object 'resDN' not found
```

```r
head(valsDN)
```

```
## Error in head(valsDN): object 'valsDN' not found
```

```r
colnames(valsDN)=c("pvalDN","padjDN")
```

```
## Error in colnames(valsDN) = c("pvalDN", "padjDN"): object 'valsDN' not found
```

```r
results=cbind(valsDH,valsAH,valsDA,valsNM,valsDD,valsDV,valsDN)
```

```
## Error in cbind(valsDH, valsAH, valsDA, valsNM, valsDD, valsDV, valsDN): object 'valsDH' not found
```

```r
head(results)
```

```
## Error in head(results): object 'results' not found
```
#### Explore Contrasts

```r
######-------------make vsd and pvals table
vsdpvals=cbind(vsd,results)
```

```
## Error in cbind(vsd, results): object 'vsd' not found
```

```r
head(vsdpvals)
```

```
## Error in head(vsdpvals): object 'vsdpvals' not found
```

```r
write.csv(vsdpvals, "VSDandPVALS_lunar.csv", quote=F)
```

```
## Error in is.data.frame(x): object 'vsdpvals' not found
```

```r
VSD_PVals=read.csv("VSDandPVALS_lunar.csv")
```

```
## Warning in file(file, "rt"): cannot open file 'VSDandPVALS_lunar.csv': No
## such file or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

####---------------get data for WGCNA

```r
head(vsdpvals)
```

```
## Error in head(vsdpvals): object 'vsdpvals' not found
```

```r
vsdp=as.data.frame(vsdpvals)
```

```
## Error in as.data.frame(vsdpvals): object 'vsdpvals' not found
```

```r
AH01=vsdp[vsdp$pvalAH<0.1 & !is.na(vsdp$pvalAH),]
```

```
## Error in eval(expr, envir, enclos): object 'vsdp' not found
```

```r
length(AH01[,1]) 
```

```
## Error in eval(expr, envir, enclos): object 'AH01' not found
```

```r
DA01=vsdp[vsdp$pvalDA<0.1 & !is.na(vsdp$pvalDA),]
```

```
## Error in eval(expr, envir, enclos): object 'vsdp' not found
```

```r
length(DA01[,1]) 
```

```
## Error in eval(expr, envir, enclos): object 'DA01' not found
```

```r
DH01=vsdp[vsdp$pvalDH<0.1 & !is.na(vsdp$pvalDH),]
```

```
## Error in eval(expr, envir, enclos): object 'vsdp' not found
```

```r
length(DH01[,1]) 
```

```
## Error in eval(expr, envir, enclos): object 'DH01' not found
```

```r
NM01=vsdp[vsdp$pvalNM<0.1 & !is.na(vsdp$pvalNM),]
```

```
## Error in eval(expr, envir, enclos): object 'vsdp' not found
```

```r
length(NM01[,1])
```

```
## Error in eval(expr, envir, enclos): object 'NM01' not found
```

```r
DD01=vsdp[vsdp$pvalDD<0.1 & !is.na(vsdp$pvalDD),]
```

```
## Error in eval(expr, envir, enclos): object 'vsdp' not found
```

```r
length(DD01[,1])  
```

```
## Error in eval(expr, envir, enclos): object 'DD01' not found
```

```r
DV01=vsdp[vsdp$pvalDV<0.1 & !is.na(vsdp$pvalDV),]
```

```
## Error in eval(expr, envir, enclos): object 'vsdp' not found
```

```r
length(DV01[,1]) 
```

```
## Error in eval(expr, envir, enclos): object 'DV01' not found
```

```r
DN01=vsdp[vsdp$pvalDN<0.1 & !is.na(vsdp$pvalDN),]
```

```
## Error in eval(expr, envir, enclos): object 'vsdp' not found
```

```r
length(DN01[,1])
```

```
## Error in eval(expr, envir, enclos): object 'DN01' not found
```

```r
degs01=union(row.names(AH01),row.names(DA01))
```

```
## Error in row.names(AH01): object 'AH01' not found
```

```r
degs01=union(degs01,row.names(DH01))
```

```
## Error in as.vector(x): object 'degs01' not found
```

```r
degs01=union(degs01,row.names(NM01))
```

```
## Error in as.vector(x): object 'degs01' not found
```

```r
degs01=union(degs01,row.names(DD01))
```

```
## Error in as.vector(x): object 'degs01' not found
```

```r
degs01=union(degs01,row.names(DV01))
```

```
## Error in as.vector(x): object 'degs01' not found
```

```r
degs01=union(degs01,row.names(DN01))
```

```
## Error in as.vector(x): object 'degs01' not found
```

```r
length(degs01)
```

```
## Error in eval(expr, envir, enclos): object 'degs01' not found
```

```r
wdegs01=vsd[(row.names(vsd) %in% degs01),]
```

```
## Error in eval(expr, envir, enclos): object 'vsd' not found
```

```r
head(wdegs01)
```

```
## Error in head(wdegs01): object 'wdegs01' not found
```

```r
write.csv(wdegs01, "genes4WGCNA.csv", quote=F)
```

```
## Error in is.data.frame(x): object 'wdegs01' not found
```

```r
######----------log GO
head(resDA)
```

```
## Error in head(resDA): object 'resDA' not found
```

```r
logs=data.frame(cbind("gene"=row.names(resDA),"logP"=round(-log(resDA$pvalue+1e-10,10),1)))
```

```
## Error in row.names(resDA): object 'resDA' not found
```

```r
logs$logP=as.numeric(as.character(logs$logP))
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
sign=rep(1,nrow(logs))
```

```
## Error in nrow(logs): object 'logs' not found
```

```r
sign[resDA$log2FoldChange<0]=-1  ##change to correct model
```

```
## Error in sign[resDA$log2FoldChange < 0] = -1: object 'resDA' not found
```

```r
table(sign)
```

```
## Error in unique.default(x, nmax = nmax): unique() applies only to vectors
```

```r
logs$logP=logs$logP*sign
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
write.table(logs,quote=F,row.names=F,file="GO_DA_logP.csv",sep=",")
```

```
## Error in is.data.frame(x): object 'logs' not found
```

```r
head(resAH)
```

```
## Error in head(resAH): object 'resAH' not found
```

```r
logs=data.frame(cbind("gene"=row.names(resAH),"logP"=round(-log(resAH$pvalue+1e-10,10),1)))
```

```
## Error in row.names(resAH): object 'resAH' not found
```

```r
logs$logP=as.numeric(as.character(logs$logP))
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
sign=rep(1,nrow(logs))
```

```
## Error in nrow(logs): object 'logs' not found
```

```r
sign[resAH$log2FoldChange<0]=-1  ##change to correct model
```

```
## Error in sign[resAH$log2FoldChange < 0] = -1: object 'resAH' not found
```

```r
table(sign)
```

```
## Error in unique.default(x, nmax = nmax): unique() applies only to vectors
```

```r
logs$logP=logs$logP*sign
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
write.table(logs,quote=F,row.names=F,file="GO_AH_logP.csv",sep=",")
```

```
## Error in is.data.frame(x): object 'logs' not found
```

```r
head(resDH)
```

```
## Error in head(resDH): object 'resDH' not found
```

```r
logs=data.frame(cbind("gene"=row.names(resDH),"logP"=round(-log(resDH$pvalue+1e-10,10),1)))
```

```
## Error in row.names(resDH): object 'resDH' not found
```

```r
logs$logP=as.numeric(as.character(logs$logP))
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
sign=rep(1,nrow(logs))
```

```
## Error in nrow(logs): object 'logs' not found
```

```r
sign[resDH$log2FoldChange<0]=-1  ##change to correct model
```

```
## Error in sign[resDH$log2FoldChange < 0] = -1: object 'resDH' not found
```

```r
table(sign)
```

```
## Error in unique.default(x, nmax = nmax): unique() applies only to vectors
```

```r
logs$logP=logs$logP*sign
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
write.table(logs,quote=F,row.names=F,file="GO_DH_logP.csv",sep=",")
```

```
## Error in is.data.frame(x): object 'logs' not found
```

```r
head(resNM)
```

```
## Error in head(resNM): object 'resNM' not found
```

```r
logs=data.frame(cbind("gene"=row.names(resNM),"logP"=round(-log(resNM$pvalue+1e-10,10),1)))
```

```
## Error in row.names(resNM): object 'resNM' not found
```

```r
logs$logP=as.numeric(as.character(logs$logP))
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
sign=rep(1,nrow(logs))
```

```
## Error in nrow(logs): object 'logs' not found
```

```r
sign[resNM$log2FoldChange<0]=-1  ##change to correct model
```

```
## Error in sign[resNM$log2FoldChange < 0] = -1: object 'resNM' not found
```

```r
table(sign)
```

```
## Error in unique.default(x, nmax = nmax): unique() applies only to vectors
```

```r
logs$logP=logs$logP*sign
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
write.table(logs,quote=F,row.names=F,file="GO_NM_logP.csv",sep=",")
```

```
## Error in is.data.frame(x): object 'logs' not found
```

```r
head(resDD)
```

```
## Error in head(resDD): object 'resDD' not found
```

```r
logs=data.frame(cbind("gene"=row.names(resDD),"logP"=round(-log(resDD$pvalue+1e-10,10),1)))
```

```
## Error in row.names(resDD): object 'resDD' not found
```

```r
logs$logP=as.numeric(as.character(logs$logP))
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
sign=rep(1,nrow(logs))
```

```
## Error in nrow(logs): object 'logs' not found
```

```r
sign[resDD$log2FoldChange<0]=-1  ##change to correct model
```

```
## Error in sign[resDD$log2FoldChange < 0] = -1: object 'resDD' not found
```

```r
table(sign)
```

```
## Error in unique.default(x, nmax = nmax): unique() applies only to vectors
```

```r
logs$logP=logs$logP*sign
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
write.table(logs,quote=F,row.names=F,file="GO_DD_logP.csv",sep=",")
```

```
## Error in is.data.frame(x): object 'logs' not found
```

```r
head(resDV)
```

```
## Error in head(resDV): object 'resDV' not found
```

```r
logs=data.frame(cbind("gene"=row.names(resDV),"logP"=round(-log(resDV$pvalue+1e-10,10),1)))
```

```
## Error in row.names(resDV): object 'resDV' not found
```

```r
logs$logP=as.numeric(as.character(logs$logP))
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
sign=rep(1,nrow(logs))
```

```
## Error in nrow(logs): object 'logs' not found
```

```r
sign[resDV$log2FoldChange<0]=-1  ##change to correct model
```

```
## Error in sign[resDV$log2FoldChange < 0] = -1: object 'resDV' not found
```

```r
table(sign)
```

```
## Error in unique.default(x, nmax = nmax): unique() applies only to vectors
```

```r
logs$logP=logs$logP*sign
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
write.table(logs,quote=F,row.names=F,file="GO_DV_logP.csv",sep=",")
```

```
## Error in is.data.frame(x): object 'logs' not found
```

```r
head(resDN)
```

```
## Error in head(resDN): object 'resDN' not found
```

```r
logs=data.frame(cbind("gene"=row.names(resDN),"logP"=round(-log(resDN$pvalue+1e-10,10),1)))
```

```
## Error in row.names(resDN): object 'resDN' not found
```

```r
logs$logP=as.numeric(as.character(logs$logP))
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
sign=rep(1,nrow(logs))
```

```
## Error in nrow(logs): object 'logs' not found
```

```r
sign[resDN$log2FoldChange<0]=-1  ##change to correct model
```

```
## Error in sign[resDN$log2FoldChange < 0] = -1: object 'resDN' not found
```

```r
table(sign)
```

```
## Error in unique.default(x, nmax = nmax): unique() applies only to vectors
```

```r
logs$logP=logs$logP*sign
```

```
## Error in eval(expr, envir, enclos): object 'logs' not found
```

```r
write.table(logs,quote=F,row.names=F,file="GO_DN_logP.csv",sep=",")
```

```
## Error in is.data.frame(x): object 'logs' not found
```

```r
#-------------write results tables; includes log2fc and pvals
write.table(results(dds), file="DESeq.results.txt", quote=FALSE, sep="\t");  
```

```
## Error in is.data.frame(x): could not find function "results"
```

```r
write.table(resDA, file="DESeq.results.DA.txt", quote=F, sep="\t")
```

```
## Error in is.data.frame(x): object 'resDA' not found
```

```r
write.table(resAH, file="DESeq.results.AH.txt", quote=F, sep="\t")
```

```
## Error in is.data.frame(x): object 'resAH' not found
```

```r
write.table(resDH, file="DESeq.results.DH.txt", quote=F, sep="\t")
```

```
## Error in is.data.frame(x): object 'resDH' not found
```

```r
write.table(resNM, file="DESeq.results.NM.txt", quote=F, sep="\t")
```

```
## Error in is.data.frame(x): object 'resNM' not found
```

```r
write.table(resDD, file="DESeq.results.DD.txt", quote=F, sep="\t")
```

```
## Error in is.data.frame(x): object 'resDD' not found
```

```r
write.table(resDV, file="DESeq.results.DV.txt", quote=F, sep="\t")
```

```
## Error in is.data.frame(x): object 'resDV' not found
```

```r
write.table(resDN, file="DESeq.results.DN.txt", quote=F, sep="\t")
```

```
## Error in is.data.frame(x): object 'resDN' not found
```

```r
DESeq_DA=read.table("DESeq.results.DA.txt")
```

```
## Warning in file(file, "rt"): cannot open file 'DESeq.results.DA.txt': No
## such file or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
mean(DESeq_DA[,5])
```

```
## Error in mean(DESeq_DA[, 5]): object 'DESeq_DA' not found
```
you may want omit the NA's from this I'm not sure if I should be getting these

#### Write annotated results tables

```r
gg=read.delim("transcriptome_iso2gene.tab", sep="\t")
head(gg)
```

```
##   isogroup22917
## 1 isogroup49376
## 2 isogroup38329
## 3 isogroup40883
## 4 isogroup23336
## 5    isogroup96
## 6 isogroup58103
##    PREDICTED..F.box.LRR.repeat.protein.13.like.E.blastx..2e.11
## 1                             Protein CBG08605 E(blastx)=2e-13
## 2                                 GULO protein E(blastx)=9e-31
## 3                PB1 domain containing protein E(blastx)=4e-06
## 4 PREDICTED: leucine--tRNA ligase, cytoplasmic E(blastx)=9e-16
## 5   Cleavage stimulation factor 77 kDa subunit E(blastx)=3e-33
## 6    Appr-1-p processing enzyme family protein E(blastx)=2e-08
```

```r
resDA=as.data.frame(resDA)
```

```
## Error in as.data.frame(resDA): object 'resDA' not found
```

```r
resDA$X=row.names(resDA)
```

```
## Error in row.names(resDA): object 'resDA' not found
```

```r
names(resDA)
```

```
## Error in eval(expr, envir, enclos): object 'resDA' not found
```

```r
resDA=resDA[c(7,1:6)]
```

```
## Error in eval(expr, envir, enclos): object 'resDA' not found
```

```r
head(resDA)
```

```
## Error in head(resDA): object 'resDA' not found
```

```r
resDAannot=merge(resDA,gg,by=1)
```

```
## Error in merge(resDA, gg, by = 1): object 'resDA' not found
```

```r
head(resDAannot)
```

```
## Error in head(resDAannot): object 'resDAannot' not found
```

```r
resDAannot <- resDAannot[order(resDAannot$padj),]
```

```
## Error in eval(expr, envir, enclos): object 'resDAannot' not found
```

```r
write.table(resDAannot, "annotated_results_DA.txt", sep="\t", quote=F, row.names=F)
```

```
## Error in is.data.frame(x): object 'resDAannot' not found
```

```r
resDH=as.data.frame(resDH)
```

```
## Error in as.data.frame(resDH): object 'resDH' not found
```

```r
resDH$X=row.names(resDH)
```

```
## Error in row.names(resDH): object 'resDH' not found
```

```r
names(resDH)
```

```
## Error in eval(expr, envir, enclos): object 'resDH' not found
```

```r
resDH=resDH[c(7,1:6)]
```

```
## Error in eval(expr, envir, enclos): object 'resDH' not found
```

```r
resDHannot=merge(resDH,gg,by=1)
```

```
## Error in merge(resDH, gg, by = 1): object 'resDH' not found
```

```r
head(resDHannot)
```

```
## Error in head(resDHannot): object 'resDHannot' not found
```

```r
resDHannot <- resDHannot[order(resDHannot$padj),]
```

```
## Error in eval(expr, envir, enclos): object 'resDHannot' not found
```

```r
write.table(resDHannot, "annotated_results_DH.txt", sep="\t", quote=F, row.names=F)
```

```
## Error in is.data.frame(x): object 'resDHannot' not found
```

```r
resAH=as.data.frame(resAH)
```

```
## Error in as.data.frame(resAH): object 'resAH' not found
```

```r
resAH$X=row.names(resAH)
```

```
## Error in row.names(resAH): object 'resAH' not found
```

```r
names(resAH)
```

```
## Error in eval(expr, envir, enclos): object 'resAH' not found
```

```r
resAH=resAH[c(7,1:6)]
```

```
## Error in eval(expr, envir, enclos): object 'resAH' not found
```

```r
resAHannot=merge(resAH,gg,by=1)
```

```
## Error in merge(resAH, gg, by = 1): object 'resAH' not found
```

```r
head(resAHannot)
```

```
## Error in head(resAHannot): object 'resAHannot' not found
```

```r
resAHannot <- resAHannot[order(resAHannot$padj),]
```

```
## Error in eval(expr, envir, enclos): object 'resAHannot' not found
```

```r
write.table(resAHannot, "annotated_results_AH.txt", sep="\t", quote=F, row.names=F)
```

```
## Error in is.data.frame(x): object 'resAHannot' not found
```

```r
resNM=as.data.frame(resNM)
```

```
## Error in as.data.frame(resNM): object 'resNM' not found
```

```r
resNM$X=row.names(resNM)
```

```
## Error in row.names(resNM): object 'resNM' not found
```

```r
names(resNM)
```

```
## Error in eval(expr, envir, enclos): object 'resNM' not found
```

```r
resNM=resNM[c(7,1:6)]
```

```
## Error in eval(expr, envir, enclos): object 'resNM' not found
```

```r
resNMannot=merge(resNM,gg,by=1)
```

```
## Error in merge(resNM, gg, by = 1): object 'resNM' not found
```

```r
head(resNMannot)
```

```
## Error in head(resNMannot): object 'resNMannot' not found
```

```r
resNMannot <- resNMannot[order(resNMannot$padj),]
```

```
## Error in eval(expr, envir, enclos): object 'resNMannot' not found
```

```r
write.table(resNMannot, "annotated_results_NM.txt", sep="\t", quote=F, row.names=F)
```

```
## Error in is.data.frame(x): object 'resNMannot' not found
```

```r
resDD=as.data.frame(resDD)
```

```
## Error in as.data.frame(resDD): object 'resDD' not found
```

```r
resDD$X=row.names(resDD)
```

```
## Error in row.names(resDD): object 'resDD' not found
```

```r
names(resDD)
```

```
## Error in eval(expr, envir, enclos): object 'resDD' not found
```

```r
resDD=resDD[c(7,1:6)]
```

```
## Error in eval(expr, envir, enclos): object 'resDD' not found
```

```r
resDDannot=merge(resDD,gg,by=1)
```

```
## Error in merge(resDD, gg, by = 1): object 'resDD' not found
```

```r
head(resDDannot)
```

```
## Error in head(resDDannot): object 'resDDannot' not found
```

```r
resDDannot <- resDDannot[order(resDDannot$padj),]
```

```
## Error in eval(expr, envir, enclos): object 'resDDannot' not found
```

```r
write.table(resDDannot, "annotated_results_DD.txt", sep="\t", quote=F, row.names=F)
```

```
## Error in is.data.frame(x): object 'resDDannot' not found
```

```r
resDV=as.data.frame(resDV)
```

```
## Error in as.data.frame(resDV): object 'resDV' not found
```

```r
resDV$X=row.names(resDV)
```

```
## Error in row.names(resDV): object 'resDV' not found
```

```r
names(resDV)
```

```
## Error in eval(expr, envir, enclos): object 'resDV' not found
```

```r
resDV=resDV[c(7,1:6)]
```

```
## Error in eval(expr, envir, enclos): object 'resDV' not found
```

```r
resDVannot=merge(resDV,gg,by=1)
```

```
## Error in merge(resDV, gg, by = 1): object 'resDV' not found
```

```r
head(resDVannot)
```

```
## Error in head(resDVannot): object 'resDVannot' not found
```

```r
resDVannot <- resDVannot[order(resDVannot$padj),]
```

```
## Error in eval(expr, envir, enclos): object 'resDVannot' not found
```

```r
write.table(resDVannot, "annotated_results_DV.txt", sep="\t", quote=F, row.names=F)
```

```
## Error in is.data.frame(x): object 'resDVannot' not found
```

```r
resDN=as.data.frame(resDN)
```

```
## Error in as.data.frame(resDN): object 'resDN' not found
```

```r
resDN$X=row.names(resDN)
```

```
## Error in row.names(resDN): object 'resDN' not found
```

```r
names(resDN)
```

```
## Error in eval(expr, envir, enclos): object 'resDN' not found
```

```r
resDN=resDN[c(7,1:6)]
```

```
## Error in eval(expr, envir, enclos): object 'resDN' not found
```

```r
resDNannot=merge(resDN,gg,by=1)
```

```
## Error in merge(resDN, gg, by = 1): object 'resDN' not found
```

```r
head(resDNannot)
```

```
## Error in head(resDNannot): object 'resDNannot' not found
```

```r
resDNannot <- resDNannot[order(resDNannot$padj),]
```

```
## Error in eval(expr, envir, enclos): object 'resDNannot' not found
```

```r
write.table(resDNannot, "annotated_results_DN.txt", sep="\t", quote=F, row.names=F)
```

```
## Error in is.data.frame(x): object 'resDNannot' not found
```
