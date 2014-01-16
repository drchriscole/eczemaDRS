# Script to perform EdgeR analysis on case-control comparison.
#
# This uses the exact test only.
# 
# Author: ccole
###############################################################################

ver = '0.5'

## function to draw boxplot of gene expression by genotype
boxplotGene = function (gene, counts) {
   list = list()
   list[[1]] = t(counts[gene,grep('cntrl',names(counts))])
   list[[2]] = t(counts[gene,grep('eczema',names(counts))])
   
   boxplot(list,main=gene,ylab="Gene Expression (Normalised read counts)",xlab="Sample",names=c("Control","Eczema"),col=rainbow(length(list)))
}

## extract commandline args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
   if (args[1] == '--version') {
      cat(sprintf("edgeR version\t%s\n", packageVersion('edgeR')))
      cat(R.version.string,"\n")
      cat(sprintf("Script version\t%s\n",ver))
      quit('no')
   }
   stop('ERROR - Not enough arguments')
}

library(edgeR)

outPrefix = 'EdgeR_analysis_all_cases_vs_all_controls'
countsFile = args[1]
if (length(args) == 2) {
   outPrefix = args[2]
}

paste(sprintf("Script version: %s",ver))
paste(sprintf("R version: %s.%s",R.Version()$major,R.Version()$minor))
paste(sprintf("EdgeR version: %s", packageVersion('edgeR')))

paste("Reading in data...")
## read counts data. Assumes file to look like this:
## 
##  GeneID  GeneName        eczema105       eczema106 ...
##  ENSG00000000003 TSPAN6  130     214 ...
##  ENSG00000000005 TNMD    10      28 ...
counts.dat = read.delim(countsFile,head=T,row.names=1)

## remove gene names column
gene.names = data.frame(counts.dat[,1])
rownames(gene.names) <- rownames(counts.dat)
counts.dat = counts.dat[,-1]

paste("Removing low count genes...")
# remove low counts. Require at least 2 reads per sample (on average) for each gene
counts.dat = counts.dat[rowSums(counts.dat)>= (length(counts.dat)*2),]

paste("Preparing data for edgeR...")
groups = c(rep('cntrl',10),rep('eczema',26))

paste("Performing edgeR exact test...")
d = DGEList(counts.dat,group = groups)
d = calcNormFactors(d)
counts.norm = data.frame(cpm(d, normalized.lib.sizes=T))
d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)
et <- exactTest(d)

# output data
res.all = topTags(et,n=2394742329)$table
res.all = data.frame(GeneID=rownames(res.all),GeneName=gene.names[rownames(res.all),],logFC=round(res.all$logFC,digits=3),FC=round(2^res.all$logFC,digits=3),logCPM=round(res.all$logCPM,digits=3),CPM=round(2^res.all$logCPM,digits=2),PValue=signif(res.all$PValue,digits=3),FDR=signif(res.all$FDR,digits=3))
write.table(res.all[res.all$FDR < 0.05,],file=sprintf("%s.csv",outPrefix),sep=',',row.names=FALSE)
