#' get lmFit for heterogeneity across subclusters
#' @param inlist list of SingleCellExperiments (SCEs) formed by scran::quickSubCluster
#' @param clname character(1) name of cluster SCE to assess
#' @note It is assumed that 'logcounts' is an assay element,
#' and that 'subcluster' is a colData element of each SCE in inlist
#' @examples
#' data(all.sce)
#' lm3 = get_subcl_LM(all.sce, "3")
#' names(lm3)
#' @export
get_subcl_LM = function(inlist, clname) {
  se = inlist[[clname]]
  x = assay(se, "logcounts")
  mm = model.matrix(~subcluster, data=as.data.frame(colData(se)))
  lmFit(x, mm)
}

#' get lmFit F-stat based collection of n genes most varying in mean across subclusters
#' @param inlist list of SingleCellExperiments (SCEs) formed by scran::quickSubCluster
#' @param clname character(1) name of cluster SCE to assess
#' @param n numeric(1) number to preserve
#' @examples
#' data(all.sce)
#' scl = get_subclustering_features(all.sce, "3", 10)
#' names(scl)
#' @export
get_subclustering_features = function(inlist, clname, n=20) {
  lm1 = get_subcl_LM( inlist, clname )
  p = seq(2, ncol(lm1$coef)) # to get F stats
  suppressWarnings({
    lm1 = eBayes(lm1) # lots of zeroes
  })
  tt = topTable(lm1, p, n=n)
  list(feat=rowData(inlist[[clname]][rownames(tt),]), stats=tt)
}
#zz = lmFit(assay(all.sce[["3"]], "logcounts"), model.matrix(~subcluster, data=as.data.frame(colData(all.sce[["3"]]))#, data=all.sce[[3]])
#)
#library(citeseqApp)
#data(all.sce)
#plotExpression(all.sce[["3"]], x="subcluster", features=c("GZMH", "IL7R", "KLRB1"), swap_rownames="Symbol", ncol=3)
#options(bitmapType="cairo")
#plotExpression(all.sce[["3"]], x="subcluster", features=c("GZMH", "IL7R", "KLRB1"), swap_rownames="Symbol", ncol=3)
#dev.off()
#plotExpression(all.sce[["3"]], x="subcluster", features=c("GZMH", "IL7R", "KLRB1"), swap_rownames="Symbol", ncol=3)
#all.sce[[3]]
#zz = lmFit(logcounts~subcluster, data=all.sce[[3]])
#?lmFit
#zz = lmFit(assay(all.sce[["3"]], "logcounts"), model.matrix(~subcluster, data=as.data.frame(colData(all.sce[["3"]]))#, data=all.sce[[3]])
#)
#)
#ezz = eBayes(zz)
#ezz[[1]][1,]
#topTable(ezz, 2:3)
#rowData(all.sce[[3]])[rownames(.Last.value),]
#plotExpression(all.sce[["3"]], x="subcluster", features=c("GZMH", "IL7R", "KLRB1", "LTB", "NKG7"), swap_rownames="Symbol", ncol=5)
