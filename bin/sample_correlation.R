# Script to generate correlation heatmap and clustering dendrogram from 
# sample gene expression 
# 
# Author: Chris Cole
###############################################################################

ver='0.5'

## extract commandline args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
   if (args[1] == '--version') {
      cat(sprintf("gplots version\t%s\n", packageVersion('gplots')))
      cat(R.version.string,"\n")
      cat(sprintf("Script version\t%s\n",ver))
      quit('no')
   }
   stop('ERROR - Not enough arguments')
}

library(gplots)

paste(sprintf("Script version: %s",ver))
paste(sprintf("R version: %s.%s",R.Version()$major,R.Version()$minor))
paste(sprintf("gplots version: %s", packageVersion('gplots')))

outPrefix = 'sample_correlation'
countsFile = args[1]
ctrlGenotypesFile = args[2]
caseGenotypesFile = args[3]
if (length(args) == 4) {
   outPrefix = args[4]
}

## read counts data. Assumes file to look like this:
## 
##  GeneID  GeneName        eczema105       eczema106 ...
##  ENSG00000000003 TSPAN6  130     214 ...
##  ENSG00000000005 TNMD    10      28 ...
counts.dat = read.delim(countsFile,head=T,row.names=1)

## read genotypes/gender data. Assumes file looks like this:
##  
##  Sample  Genotype        Gender  Age
##  eczema105       het     m  16
##  eczema106       wt      m  15
genotypes.ctrl = read.delim(ctrlGenotypesFile,head=T)
genotypes.case = read.delim(caseGenotypesFile,head=T)
genotypes.ctrl[,2] <- paste(genotypes.ctrl[,2],"ctrl",sep='_')
genotypes.dat = rbind(genotypes.ctrl,genotypes.case)

if (nrow(genotypes.dat) != length(counts.dat)-1) {
   stop(sprintf("ERROR - no. columns in %s doesn't agree with no. rows in %s", countsFile, genotypesFile))
}

## remove gene names column and store ID vs name lookup
gene.names = data.frame(id=rownames(counts.dat), name=counts.dat[,1], stringsAsFactors=FALSE)
#rownames(gene.names) <- counts.dat[,1]
counts.dat = counts.dat[,-1]

## define genotypes
ctrl.wt = genotypes.ctrl[genotypes.ctrl$Genotype == 'wt_ctrl',1]
ctrl.het = genotypes.ctrl[genotypes.ctrl$Genotype == 'het_ctrl',1]
case.wt = genotypes.case[genotypes.case$Genotype == 'wt',1]
case.het = genotypes.case[genotypes.case$Genotype == 'het',1]
case.chet = genotypes.case[genotypes.case$Genotype == 'cmpdhet',1]

## replace sample names with phenotype/genotype labels
names = vector()
names[which(colnames(counts.dat) %in% ctrl.wt)] <- 'Control (Wild-type)'
names[which(colnames(counts.dat) %in% ctrl.het)] <- 'Control (heterozygous)'
names[which(colnames(counts.dat) %in% case.wt)] <- 'Case (Wild-type)'
names[which(colnames(counts.dat) %in% case.het)] <- 'Case (heterozygous)'
names[which(colnames(counts.dat) %in% case.chet)] <- 'Case (compound het)'

## plot correlation matrix of all-against-all
pdf(file=sprintf("%s_heatmap.pdf", outPrefix),height=7,width=8)
heatmap.2(cor(counts.dat),Colv=F,Rowv=F,trace='none',dendrogram="column",margins=c(7.5,7.5),labRow=names, labCol=names,keysize=1,lhei = c(2, 8),density.info='none',lwid=c(1,5))
dev.off()

# append gender and age to sample names for dendrogram plot
genders = genotypes.dat[genotypes.dat$Sample == colnames(counts.dat),3]
ages = genotypes.dat[genotypes.dat$Sample == colnames(counts.dat),4]
names=paste(names,genders,sep=" ")
names=paste(names,ages,sep=" ")

## plot dendrogram
colnames(counts.dat) <- names
d = dist(t(cor(counts.dat)))
clust  = hclust(d)
pdf(file=sprintf("%s_dendrogram.pdf", outPrefix),height=7,width=8)
par(mai=c(0.6,0.6,0.4,2.1))
plot(as.dendrogram(clust),axes=FALSE,main="",xlab="",ylab="",sub="",horiz=TRUE)
dev.off()
