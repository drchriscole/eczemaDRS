# Script to perform EdgeR analysis on case-control comparison.
#
# Using GLMs to model all the genotypes and to control for gender biases.
# 
# Author: ccole
###############################################################################

ver = '1.5'

## function to draw boxplot of gene expression by genotype
boxplotGene = function (gene, counts) {
   list = list()
   list[[1]] = t(counts[gene,colnames(counts) %in% ctrl.wt])
   list[[2]] = t(counts[gene,colnames(counts) %in% case.wt])
   list[[3]] = t(counts[gene,colnames(counts) %in% case.het])
   list[[4]] = t(counts[gene,colnames(counts) %in% case.chet])

   par(mgp=c(3,1.5,0))
   boxplot(list,main=gene,ylab="Gene Expression (Normalised read counts)",xlab="Sample Genotype",names=c("Ctrl\nWT","Eczema\nWT","Eczema\nHet","Eczema\nCmpdHet"),col=rainbow(length(list)), varwidth=T)
}

## extract commandline args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
   if (args[1] == '--version') {
      cat(sprintf("edgeR.R version\t%s\n", packageVersion('edgeR')))
      cat(R.version.string,"\n")
      cat(sprintf("Script version\t%s\n",ver))
      quit('no')
   }
   stop('ERROR - Not enough arguments')
}
library(edgeR)

paste(sprintf("Script version: %s",ver))
paste(sprintf("R version: %s.%s",R.Version()$major,R.Version()$minor))
paste(sprintf("EdgeR version: %s", packageVersion('edgeR')))

outPrefix = 'edgeR'
countsFile = args[1]
ctrlGenotypesFile = args[2]
caseGenotypesFile = args[3]
if (length(args) == 4) {
   outPrefix = args[4]
}

paste("Reading in data...")
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
genotypes.ctrl = read.delim(ctrlGenotypesFile,head=T)
genotypes.case = read.delim(caseGenotypesFile,head=T)
genotypes.ctrl[,2] <- paste(genotypes.ctrl[,2],"ctrl",sep='_')
genotypes.dat = rbind(genotypes.ctrl,genotypes.case)

if ((nrow(genotypes.case) + nrow(genotypes.ctrl)) != length(counts.dat)-1) {
   stop(sprintf("ERROR - no. columns in %s doesn't agree with no. rows in %s", countsFile, genotypesFile))
}

## remove gene names column
gene.names = data.frame(counts.dat[,1])
rownames(gene.names) <- rownames(counts.dat)
counts.dat = counts.dat[,-1]

## separate controls from patients

ctrl.wt = genotypes.ctrl[genotypes.ctrl$Genotype == 'wt_ctrl',1]
ctrl.het = genotypes.ctrl[genotypes.ctrl$Genotype == 'het_ctrl',1]
case.wt = genotypes.case[genotypes.case$Genotype == 'wt',1]
case.het = genotypes.case[genotypes.case$Genotype == 'het',1]
case.chet = genotypes.case[genotypes.case$Genotype == 'cmpdhet',1]

## drop ctrl het samples from count data as too few to be meaningful
counts.dat = subset(counts.dat, select = -c(ctrl.het))

## ...and from genotype data
genotypes.dat = genotypes.dat[genotypes.dat$Genotype != 'het_ctrl',]

paste("Removing low count genes...")
# remove low counts. Require at least 2 reads per sample (on average) for each gene
counts.dat = counts.dat[rowSums(counts.dat)>= (length(counts.dat)*2),]

paste("Preparing data for edgeR...")
genotypes = genotypes.dat[genotypes.dat$Sample == colnames(counts.dat),2]
genders = genotypes.dat[genotypes.dat$Sample == colnames(counts.dat),3]
design = model.matrix(~genders+as.factor(genotypes))

paste("Performing edgeR GLM fit...")
d = DGEList(counts.dat,group = genotypes)
d = calcNormFactors(d)
counts.norm = data.frame(cpm(d, normalized.lib.sizes=T))
d = estimateGLMCommonDisp(d)
paste(sprintf("Square root of common dispersion: %.2f", sqrt(d$common.dispersion)))
d = estimateGLMTrendedDisp(d,design)
d = estimateGLMTagwiseDisp(d,design)
#pdf(file="BCV_plot.pdf",width=6,height=6)
#plotBCV(d)
#dev.off()
fit <- glmFit(d, design)

paste("Outputting significant (FDR < 0.05) gene lists....")
# this is ctrl wt vs cmpdhet
if (packageVersion("edgeR")$major > 2) {
   paste("WARNING - you're running a newer version of edgeR than required. Results may differ from those published.")
   lrt.ctrlvscmpdhet <- glmLRT( fit,contrast=c(0,0,0,0,-1))
} else {
   lrt.ctrlvscmpdhet <- glmLRT(d, fit,contrast=c(0,0,0,0,-1))
}
#topTags(lrt.wtvscmpdhet)
all.cmpdhet=topTags(lrt.ctrlvscmpdhet,n=2423432)$table
all.cmpdhet = data.frame(GeneID=rownames(all.cmpdhet),GeneName=gene.names[rownames(all.cmpdhet),],logFC=round(all.cmpdhet$logFC,digits=3),FC=round(2^all.cmpdhet$logFC,digits=3),logCPM=round(all.cmpdhet$logCPM,digits=3),CPM=round(2^all.cmpdhet$logCPM,digits=2),LR=round(all.cmpdhet$LR,digits=2),PValue=signif(all.cmpdhet$PValue,digits=2),FDR=signif(all.cmpdhet$FDR,digits=2))
write.table(all.cmpdhet[all.cmpdhet$FDR < 0.05,],file=sprintf("%s_ctrl_wt_vs_cmpdhet.csv",outPrefix),sep=',',row.names=FALSE)

# this is ctrl wt vs het
if (packageVersion("edgeR")$major > 2) {
   lrt.ctrlwtvshet = glmLRT(fit,contrast=c(0,0,1,0,-1))
} else {
   lrt.ctrlwtvshet = glmLRT(d,fit,contrast=c(0,0,1,0,-1))
}
all.het=topTags(lrt.ctrlwtvshet,n=2423432)$table
all.het = data.frame(GeneID=rownames(all.het),GeneName=gene.names[rownames(all.het),],logFC=round(all.het$logFC,digits=3),FC=round(2^all.het$logFC,digits=3),logCPM=round(all.het$logCPM,digits=3),CPM=round(2^all.het$logCPM,digits=2),LR=round(all.het$LR,digits=2),PValue=signif(all.het$PValue,digits=2),FDR=signif(all.het$FDR,digits=2))
write.table(all.het[all.het$FDR < 0.05,],file=sprintf("%s_ctrl_wt_vs_het.csv",outPrefix),sep=',',row.names=FALSE)

# ctrl wt vs eczema wt
if (packageVersion("edgeR")$major > 2) {
   ltr.ctrlwtvseczemawt = glmLRT(fit,contrast=c(0,0,0,1,-1))
} else {
   ltr.ctrlwtvseczemawt = glmLRT(d,fit,contrast=c(0,0,0,1,-1))
}
all.wt=topTags(ltr.ctrlwtvseczemawt,n=2423432)$table
all.wt = data.frame(GeneID=rownames(all.wt),GeneName=gene.names[rownames(all.wt),],logFC=round(all.wt$logFC,digits=3),FC=round(2^all.wt$logFC,digits=3),logCPM=round(all.wt$logCPM,digits=3),CPM=round(2^all.wt$logCPM,digits=2),LR=round(all.wt$LR,digits=2),PValue=signif(all.wt$PValue,digits=2),FDR=signif(all.wt$FDR,digits=2))
write.table(all.wt[all.wt$FDR < 0.05,],file=sprintf("%s_ctrl_wt_vs_wt.csv",outPrefix),sep=',',row.names=FALSE)
