library(shiny)
library(citeseqApp)
data(sce)
data(all.sce)
data(clusters.adt)
inlist = all.sce
adtcls = clusters.adt

 server = function(input, output, session) {
  output$tsne = renderPlot({
   plotTSNE(altExp(sce), colour_by="label", text_by = "label", text_color="red")
   })
  output$heatmap = renderPlot({
   se.avg = sumCountsAcrossCells(altExp(sce), adtcls, exprs_values = "logcounts", average=TRUE)
   avg = assay(se.avg)
   pheatmap::pheatmap(avg - rowMeans(avg), breaks=seq(-3, 3, length.out=101))
   })
  featdata = reactive({
     get_subclustering_features(inlist, input$clpick, n=10) 
     })
  output$feats = renderUI({
    scl = featdata()$feat
    checkboxGroupInput("genes", "genes for boxplots and smooths", choices=scl$Symbol, selected=scl$Symbol[1:3])
    })
  output$stats = renderDataTable({
    tab = featdata()$stats
    cl = which(sapply(tab, is.numeric))
    for (j in cl) tab[[j]] = round(tab[[j]], 4)
    tab
    })
  output$boxplots = renderPlot({
    featdata()$feat
    plotExpression(inlist[[input$clpick]], x="subcluster", features=input$genes,
         swap_rownames = "Symbol", ncol=length(input$genes))
    })
  output$smooths = renderPlot({
    featdata()$feat
    plotExpression(inlist[[input$clpick]], x=input$baseADT, features=input$genes,
         show_smooth=TRUE, show_se=FALSE,
         swap_rownames = "Symbol", ncol=length(input$genes))
    })
 }
