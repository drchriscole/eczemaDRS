########################################
# Plot gene expression across genotypes
# 
# Author: ccole
########################################

ver = '0.3'

## function to draw boxplot of gene expression by genotype
boxplotGene = function (gene, counts) {
   list = list()
   list[[1]] = t(counts[gene,colnames(counts) %in% ctrl.wt])
   list[[2]] = t(counts[gene,colnames(counts) %in% case.wt])
   list[[3]] = t(counts[gene,colnames(counts) %in% case.het])
   list[[4]] = t(counts[gene,colnames(counts) %in% case.chet])
   
   
   # constuct plot bit by bit to ensure axis labels are ok
   boxplot(list,cex.axis=0.8,varwidth=TRUE, col=rainbow(4),axes=FALSE)
   box()
   axis(2)
   mtext("Gene Expression (Normalised read counts)",side=2,line=3,cex=1.2)
   par(mgp=c(3.5,1.5,0))
   axis(1,at=seq(1:4),cex.axis=0.9,labels=c("Control\nWild-type","Eczema case\nwild-type","Eczema case\nheterozygote","Eczema case\ncompound het"))
   mtext("Phenotype and FLG genotype",side=1,line=3,cex=1.2)
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

gene = 'FLG'
countsFile = args[1]
ctrlGenotypesFile = args[2]
caseGenotypesFile = args[3]
if (length(args) > 3) {
   gene = args[4]
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

paste("Removing low count genes...")
# remove low counts. Require at least 2 reads per sample (on average) for each gene
counts.dat = counts.dat[rowSums(counts.dat)>= (length(counts.dat)*2),]

ctrl.wt = genotypes.ctrl[genotypes.ctrl$Genotype == 'wt_ctrl',1]
ctrl.het = genotypes.ctrl[genotypes.ctrl$Genotype == 'het_ctrl',1]
case.wt = genotypes.case[genotypes.case$Genotype == 'wt',1]
case.het = genotypes.case[genotypes.case$Genotype == 'het',1]
case.chet = genotypes.case[genotypes.case$Genotype == 'cmpdhet',1]

## drop ctrl het samples from count data as too few to be meaningful
counts.dat = subset(counts.dat, select = -c(ctrl.het))

## ...and from genotype data
genotypes.dat = genotypes.dat[genotypes.dat$Genotype != 'het_ctrl',]


genotypes = genotypes.dat[genotypes.dat$Sample == colnames(counts.dat),2]
genders = genotypes.dat[genotypes.dat$Sample == colnames(counts.dat),3]
design = model.matrix(~genders+genotypes)

d = DGEList(counts.dat,group = genotypes)
d = calcNormFactors(d)
counts.norm = data.frame(cpm(d, normalized.lib.sizes=T))

## plot data
pdf(file=sprintf("%s_boxplot.pdf",gene))
boxplotGene(gene.names[gene.names$name==gene,1],counts.norm)
dev.off()

