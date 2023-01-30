library(shiny)
library(citeseqApp)
data(sce)
data(all.sce)
data(clusters.adt)
inlist = all.sce
adtcls = clusters.adt

 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText("Explore CITE-seq subclusters"),
    selectInput("clpick", "ADT cluster", choices=names(inlist), selected="3"),
    selectInput("baseADT", "Base for smooths", choices=rowData(altExp(sce))$ID, selected="CD127",
       multiple=FALSE),
    uiOutput("feats"), width=2
    ),
   mainPanel(
    tabsetPanel(
     tabPanel("tsne", helpText("Guide to ADT-based clusters"), plotOutput("tsne")),
     tabPanel("heatmap", helpText("Guide to protein abundance profiles"), plotOutput("heatmap")),
     tabPanel("boxplots", plotOutput("boxplots")),
     tabPanel("smooths", plotOutput("smooths")),
     tabPanel("stats", dataTableOutput("stats")),
     tabPanel("about", helpText("This app helps to explore RNA-based subclusters of ADT-based clusters
 formed according to ch 12.6.1 of the OSCA book.  Inputs are a basic SingleCellExperiment with
logcounts for RNA and ADT features, a list of SCE formed using scran::quickSubCluster, and
the vector of assignments from cells to ADT subclusters."),
   helpText(" "),
   helpText("The TSNE map of ADT clusters is
given, along with a heatmap for ADT abundances, as guides."),  
   helpText(" "),
   helpText("Choose an ADT-based cluster using the labeling on the
TSNE map and F tests will be performed (using limma) to identify genes whose mean abundances vary strongly
across RNA-based subclusters.  Boxplot tab presents marginal distributions of expression
of selected genes in RNA-based subclusters.  Smooths tab depicts association between RNA abundance and protein
abundance for selected genes and a given protein in the ADT panel."))
     )
    )
   )
  )
