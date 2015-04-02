#
#
## Script to look at gene expression correlation with FLG for eczema patient
## samples only. As there are not many significantly diffex, genes thresholding with 
## logFC and read counts is used instead.
#
#

ver = '1.8'

# calc std error function
stderr <- function(x) {
   s=sd(x)
   se=s/sqrt(length(x))
   return(se)
}

# function to calculate correlation with FLG expression
flgCor = function (x, geneID) {
   val = cor(x,mean.all[geneID,],method='pearson')
   return(val)
}

# function to plot errorbars
error.bar <- function(x, y, upper, lower=upper, length=0.05,...){
   if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
      stop("vectors must be same length")
   arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# given a correlation coefficient, r, and the number of datapoints, n,
# calculate the t statistic
findt = function (r,n) {
   df = n-2
   den = 1 - r^2
   return(r * sqrt(df/den))
}

# given a correlation coefficient, r, assume 3 datapoints and calculate the p-value for the correlation
calcPval = function (r) {
   n = 3
   t = findt(r,n)
   if (t < 0) {
      return(pt(t,n-2))
   } else {
      return(1-pt(t,n-2))
   }
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
} else if (length(args) > 4) {
  stop('ERROR - too many arguments')
}

library(edgeR)
library(sqldf)

outPrefix = 'FLG_correlation_eczema'
geneID = 'ENSG00000143631'
countsFile = args[1]
genotypesFile = args[2]
if (length(args) >= 3) {
   outPrefix = args[3]
} else if (length(args) == 4) {
  geneID = args[4]
} 

paste(sprintf("Performing correlations against geneID %s",geneID))

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

## use edgeR to normalise counts
d = DGEList(counts=counts.dat)
d = calcNormFactors(d)
counts.norm = data.frame(cpm(d, normalized.lib.sizes=T))

## save sample names for each genotype
wt = genotypes.dat[genotypes.dat$Genotype == 'wt',1]
het = genotypes.dat[genotypes.dat$Genotype == 'het',1]
chet = genotypes.dat[genotypes.dat$Genotype == 'cmpdhet',1]

# classify samples by genotype
counts.cmpdhet = counts.norm[,colnames(counts.norm) %in% chet]
counts.het = counts.norm[,colnames(counts.norm) %in% het]
counts.wt = counts.norm[,colnames(counts.norm) %in% wt]

# calc means
mean.het = apply(counts.het,1,mean)
mean.cmpdhet = apply(counts.cmpdhet,1,mean)
mean.wt = apply(counts.wt,1,mean)
mean.all = cbind("WT"=mean.wt,"Het"=mean.het,"CmpdHet"=mean.cmpdhet)

# calc std deviation
sd.het = apply(counts.het,1,sd)
sd.cmpdhet = apply(counts.cmpdhet,1,sd)
sd.wt = apply(counts.wt,1,sd)
sd.all = cbind("WT-sd"=sd.wt,"Het-sd"=sd.het,"CmpdHet-sd"=sd.cmpdhet)

# calc std errors
se.het = apply(counts.het,1,stderr)
se.cmpdhet = apply(counts.cmpdhet,1,stderr)
se.wt = apply(counts.wt,1,stderr)
se.all = cbind("WT"=se.wt,"Het"=se.het,"CmpdHet"=se.cmpdhet)

# calc correlation of expression with FLG 
if (geneID %in% rownames(mean.all)) {
  flg.cor = apply(mean.all[rowSums(mean.all) > 100,],1,flgCor,geneID)
} else {
  stop(sprintf("GeneID '%s' not found in data. Either it is not a valid ID or gene is not expressed.", geneID))
}


# calc significant of correlation
flg.cor.pval = sapply(flg.cor,calcPval)

# create data frame with gene expression (Eczema WT, Het & Cmpdhet), FLG correlation and Fold-change
ratio = mean.het/mean.wt
tmp = cbind(mean.all[rowSums(mean.all) > 100,],sd.all[rowSums(mean.all) > 100,],"cor"=flg.cor,"FC"=ratio[rowSums(mean.all) > 100],"logFC"=log2(ratio[rowSums(mean.all) > 100]),pval=flg.cor.pval)
dat.het = data.frame(tmp)
ratio = mean.cmpdhet/mean.wt
tmp = cbind(mean.all[rowSums(mean.all) > 100,],sd.all[rowSums(mean.all) > 100,],"cor"=flg.cor,"FC"=ratio[rowSums(mean.all) > 100],"logFC"=log2(ratio[rowSums(mean.all) > 100]),pval=flg.cor.pval)
dat.cmpdhet = data.frame(tmp)

# plot scatter of logFC vs correlation
pdf(file=sprintf("%s.pdf",outPrefix),width=9,height=6)
par(mfrow=c(1,2),mai=c(0.8,0.9,0.8,0.2),mgp=c(2.1,0.8,0))
plot(dat.het$cor,dat.het$FC,pch=19,cex=0.3,main="Heterozygous vs Wild-type",xlab="",ylab="",ylim=c(0.1,6),log='y')
points(dat.het["ENSG00000143631",7],dat.het["ENSG00000143631",8],col="orange",pch=19,cex=0.8)
abline(h=0.5,lty=2,col="grey")
abline(h=2,lty=2,col="grey")
abline(v=0.8,lty=2,col="grey")
abline(v=-0.8,lty=2,col="grey")
legend("topright",legend=c("FLG"),pch=19,cex=0.8,col=c("orange"),bg="white")
plot(dat.cmpdhet$cor,dat.cmpdhet$FC,pch=19,cex=0.3,main="Compound Het. vs Wild-type",xlab="",ylab="",ylim=c(0.1,6),log='y')
points(dat.cmpdhet["ENSG00000143631",7],dat.cmpdhet["ENSG00000143631",8],col="orange",pch=19,cex=0.8)
abline(h=0.5,lty=2,col="grey")
abline(h=2,lty=2,col="grey")
abline(v=0.8,lty=2,col="grey")
abline(v=-0.8,lty=2,col="grey")
legend("topright",legend=c("FLG"),pch=19,cex=0.8,col=c("orange"),bg="white")
mtext("Correlation with FLG Expression",side=1,outer=T,line=-1)
mtext(expression("Fold-change "(log[2]~scale)),side=2,outer=T,line=-2,las=3)
dev.off()

# select genes that full within certain log2 and r criteria
#dat.het.select = sqldf("SELECT * from 'dat.het' where pval < 0.05 and abs(logFC) > 0.5",row.names=T)
dat.cmpdhet.select = sqldf("SELECT * from 'dat.cmpdhet' where pval < 0.05 and abs(logFC) > 0.5",row.names=T)

output = cbind('GeneID'=rownames(dat.cmpdhet.select[order(dat.cmpdhet.select$pval),]),GeneName=gene.names[rownames(dat.cmpdhet.select[order(dat.cmpdhet.select$pval),]),],signif(dat.cmpdhet.select[order(dat.cmpdhet.select$pval),],digit=2))
write.table(output,file=sprintf("%s.csv",outPrefix),sep=",",row.names=FALSE)
