# Script to perform analysis on case data only
# 
# Author: ccole
###############################################################################

ver = '1.7'

## function to draw boxplot of gene expression by genotype
boxplotGene = function (gene, counts) {
   list = list()
   list[[1]] = t(counts[gene,colnames(counts.dat) %in% wt])
   list[[2]] = t(counts[gene,colnames(counts.dat) %in% het])
   list[[3]] = t(counts[gene,colnames(counts.dat) %in% chet])
   
   boxplot(list,main=gene,ylab="Gene Expression (Normalised read counts)",xlab="Sample Genotype",names=c("Eczema\nWT","Eczema\nHet","Eczema\nCmpdHet"),col=rainbow(3))
}

## extract commandline args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
   if (args[1] == '--version') {
      cat(sprintf("edgeR.R version\t%s\n", packageVersion('edgeR')))
      cat(R.version.string,"\n")
      cat(sprintf("Script version\t%s\n",ver))
      quit('no')
   }
   stop('ERROR - Not enough arguments')
}
library(edgeR)

outPrefix = 'EdgeR_analysis'
countsFile = args[1]
genotypesFile = args[2]
if (length(args) == 3) {
   outPrefix = args[3]
}

## read counts data. Assumes file to look like this:
## 
##  GeneID  GeneName        eczema105       eczema106 ...
##  ENSG00000000003 TSPAN6  130     214 ...
##  ENSG00000000005 TNMD    10      28 ...
counts.dat = read.delim(countsFile,head=T,row.names=1)

## read genotypes/gender data. Assumes file looks like this:
##  
##  Sample  Genotype        Gender
##  eczema105       het     m
##  eczema106       wt      m
genotypes.dat = read.delim(genotypesFile,head=T)

if (nrow(genotypes.dat) != length(counts.dat)-1) {
   stop(sprintf("ERROR - no. columns in %s doesn't agree with no. rows in %s", countsFile, genotypesFile))
}

## remove gene names column
gene.names = data.frame(counts.dat[,1])
rownames(gene.names) <- rownames(counts.dat)
counts.dat = counts.dat[,-1]

wt = genotypes.dat[genotypes.dat$Genotype == 'wt',1]
het = genotypes.dat[genotypes.dat$Genotype == 'het',1]
chet = genotypes.dat[genotypes.dat$Genotype == 'cmpdhet',1]

genotypes = genotypes.dat[genotypes.dat$Sample == colnames(counts.dat),2]
genders = genotypes.dat[genotypes.dat$Sample == colnames(counts.dat),3]
design = model.matrix(~genders+as.factor(genotypes))

d = DGEList(counts.dat,group = genotypes)
d = calcNormFactors(d)
d = estimateGLMCommonDisp(d)
paste(sprintf("Square root of common dispersion: %.2f", sqrt(d$common.dispersion)))
d = estimateGLMTrendedDisp(d,design)
d = estimateGLMTagwiseDisp(d,design)
fit <- glmFit(d, design)

# this is wt vs cmpdhet
if (packageVersion("edgeR")$major > 2) {
   paste("WARNING - you're running a newer version of edgeR than required. Results may differ from those published.")
   lrt.wtvscmpdhet <- glmLRT(fit,contrast=c(0,0,0,-1))
} else {
   lrt.wtvscmpdhet <- glmLRT(d, fit,contrast=c(0,0,0,-1))
}
#topTags(lrt.wtvscmpdhet)
all.cmpdhet=topTags(lrt.wtvscmpdhet,n=2423432)$table
all.cmpdhet = data.frame(GeneID=rownames(all.cmpdhet),GeneName=gene.names[rownames(all.cmpdhet),],logFC=round(all.cmpdhet$logFC,digits=3),FC=round(2^all.cmpdhet$logFC,digits=3),logCPM=round(all.cmpdhet$logCPM,digits=3),CPM=round(2^all.cmpdhet$logCPM,digits=2),LR=round(all.cmpdhet$LR,digits=2),PValue=signif(all.cmpdhet$PValue,digits=3),FDR=signif(all.cmpdhet$FDR,digits=3))
write.table(all.cmpdhet[all.cmpdhet$PValue < 0.01,],file=sprintf("%s_FLG_wt_cases_vs_FLG_cmpdhet_cases.csv",outPrefix),sep=',',row.names=FALSE)

# this is wt vs het
if (packageVersion("edgeR")$major > 2) {
   lrt.wtvshet = glmLRT(fit,contrast=c(0,0,1,-1))
} else {
   lrt.wtvshet = glmLRT(d,fit,contrast=c(0,0,1,-1))
}
all.het=topTags(lrt.wtvshet,n=2423432)$table
all.het = data.frame(GeneID=rownames(all.het),GeneName=gene.names[rownames(all.het),],logFC=round(all.het$logFC,digits=3),FC=round(2^all.het$logFC,digits=3),logCPM=round(all.het$logCPM,digits=3),CPM=round(2^all.het$logCPM,digits=2),LR=round(all.het$LR,digits=2),PValue=signif(all.het$PValue,digits=3),FDR=signif(all.het$FDR,digits=3))
write.table(all.het[all.het$PValue < 0.01,],file=sprintf("%s_FLG_wt_cases_vs_FLG_het_cases.csv",outPrefix),sep=',',row.names=FALSE)

