#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

version = "0.2"
countsFile = '../../data/all_gene_expression.tsv'
ctrlGenotypesFile = '../../data/control_FLG_genotypes.dat'
caseGenotypesFile = '../../data/eczema_FLG_genotypes.dat'

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

# ignoring genotypes as only interested in cases vs ctrls
# and we know first 10 are the controls
genotypes = c(rep('ctrl', 10), rep('case',26))
genders = genotypes.dat[genotypes.dat$Sample == colnames(counts.dat),3]
design = model.matrix(~genders+genotypes)

d = DGEList(counts.dat,group = genotypes)
d = calcNormFactors(d)
counts.norm = data.frame(cpm(d, normalized.lib.sizes=T))


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Skin DRS Data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         textInput("gene", "Gene Symbol", "FLG")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("gene_boxplot"),
         downloadLink("downloadPlot", "Download Plot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   

  output$gene_boxplot <- renderPlot({
      if (!is.element(input$gene,gene.names$name)) {
        stop(sprintf("Gene '%s' not found in data. Either it is not a valid gene or it is not expressed.", input$gene))
      }
    
    list = list()
    list[[1]] = t(counts.norm[gene.names[gene.names$name==input$gene,1],1:10])
    list[[2]] = t(counts.norm[gene.names[gene.names$name==input$gene,1],11:36])
    boxplot(list,cex.axis=0.8,varwidth=FALSE, col=c('#fff8c7','#d6ca79'),axes=FALSE)
    box()
    axis(2)
    mtext("Gene Expression (Normalised read counts)",side=2,line=3,cex=1.2)
    par(mgp=c(3.5,2.5,0))
    axis(1,at=seq(1:2),cex.axis=0.9,
         labels=c(sprintf("Control\nSkin\nn = 10"),
                  sprintf("Atopic\nSkin\nn = 26")
         ))
    title(input$gene)
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function(){paste(input$gene, '.pdf', sep = '')},
    
    content = function(file){
      cairo_pdf(filename = file,
                width = 18, height = 10, pointsize = 12, family = "sans", bg = "transparent",
                antialias = "subpixel",fallback_resolution = 300)
      
      list = list()
      list[[1]] = t(counts.norm[gene.names[gene.names$name==input$gene,1],1:10])
      list[[2]] = t(counts.norm[gene.names[gene.names$name==input$gene,1],11:36])
      boxplot(list,cex.axis=0.8,varwidth=FALSE, col=c('#fff8c7','#d6ca79'),axes=FALSE)
      box()
      axis(2)
      mtext("Gene Expression (Normalised read counts)",side=2,line=3,cex=1.2)
      par(mgp=c(3.5,2.5,0))
      axis(1,at=seq(1:2),cex.axis=0.9,
           labels=c(sprintf("Control\nSkin\nn = 10"),
                    sprintf("Atopic\nSkin\nn = 26")
           ))
      title(input$gene)
      dev.off()
    },
    
    contentType = "application/pdf"
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)
