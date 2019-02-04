#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Skin DRS Data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         selectInput("gene",
                     label = "Select gene symbol:",
                     choices = c("FLG",
                                 "C11orf30",
                                 "LRRC32"),
                     selected = "FLG")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         textOutput("selected_txt")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
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

  output$selected_txt <- renderText({
    paste("You selected", input$gene)
  })  
  
    output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

