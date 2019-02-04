########################################
# Plot gene expression across genotypes
# 
# Author: ccole
########################################

ver = '0.1'
library(edgeR)

paste(sprintf("Script version: %s",ver))
paste(sprintf("R version: %s.%s",R.Version()$major,R.Version()$minor))
paste(sprintf("EdgeR version: %s", packageVersion('edgeR')))

genes = c('C11orf30', 'LRRC32')
countsFile = 'data/all_gene_expression.tsv'
ctrlGenotypesFile = 'data/control_FLG_genotypes.dat'
caseGenotypesFile = 'data/eczema_FLG_genotypes.dat'


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
if (!is.element(genes,gene.names$name)) {
  stop(sprintf("Gene '%s' not found in data. Either it is not a valid gene or it is not expressed.", gene))
}
#rownames(gene.names) <- counts.dat[,1]
counts.dat = counts.dat[,-1]

paste("Removing low count genes...")
# remove low counts. Require at least 2 reads per sample (on average) for each gene
counts.dat = counts.dat[rowSums(counts.dat)>= (length(counts.dat)*2),]

# ignoring genotypes as only interested in cases vs ctrls
# and we know first 10 are the controls
genotypes = c(rep('ctrl', 10), rep('case',26))
genders = genotypes.dat[genotypes.dat$Sample == colnames(counts.dat),3]
design = model.matrix(~genders+genotypes)

d = DGEList(counts.dat,group = genotypes)
d = calcNormFactors(d)
counts.norm = data.frame(cpm(d, normalized.lib.sizes=T))

## plot data

list = list()
list[[1]] = t(counts.norm[gene.names[gene.names$name==genes[1],1],1:10])
list[[2]] = t(counts.norm[gene.names[gene.names$name==genes[1],1],11:36])
list[[3]] = t(counts.norm[gene.names[gene.names$name==genes[2],1],1:10])
list[[4]] = t(counts.norm[gene.names[gene.names$name==genes[2],1],11:36])


# constuct plot bit by bit to ensure axis labels are ok
pdf("EMSY_LRRC32_boxplot.pdf")
boxplot(list,cex.axis=0.8,varwidth=FALSE, col=c('#fff8c7','#d6ca79', '#c9dfff', '#7ea0d4'),axes=FALSE)
box()
axis(2)
mtext("Gene Expression (Normalised read counts)",side=2,line=3,cex=1.2)
par(mgp=c(3.5,2.5,0))
axis(1,at=seq(1:4),cex.axis=0.9,
     labels=c(sprintf("Control\nSkin\nn = 10"),
              sprintf("Atopic\nSkin\nn = 26"),
              sprintf("Control\nSkin\nn = 10"),
              sprintf("Atopic\nSkin\nn = 26")
     ))
mtext("EMSY", side=3, line=1, cex=1, adj=0)
mtext("LRRC32", side=3, line=1, cex=1, adj=1)
dev.off()

