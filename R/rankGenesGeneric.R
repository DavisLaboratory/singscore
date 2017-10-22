#' @include simpleScoring.R CoreFuns.R
NULL

#'@title  Rank genes by the gene expression intensities
#'@description The \code{rankExpr} function is a generic function that takes
#'mutilple combinations of inputs.If input is S4 object of \code{DGEList, or
#'ExpressionSet}, the generic function will extract the gene epxression matrix
#'genes on the rows and samples on columns from the object and rank the genes.
#'Given a matrix of gene expression that has samples in columns, genes in rows,
#'and values being gene expression intensity. \code{rankExpr} ranks gene
#'expression in each sample.
#'
#'@param expreMatrix A gene expression matrix or S4 object
#'  (ExpressionSet,DGEList,data.frame)
#'@param tiesMethod A character indicating what method to use when dealing with
#'  ties
#'@name rankExpr
#'
#'@return The ranked gene expression matrix that has samples in columns and
#'  genes in rows
#' @examples
#' rankExpr(toy_expr) # toy_expr is a gene expression matrix, tiesMethod = 'min'
#' e <- Biobase::ExpressionSet(assayData = as.matrix(toy_expr))
#' rankExpr(e) # or it can be a ExpressionSet object
#'@seealso \code{\link{rank}} \linkS4class{GeneSet}
#'
#'@export
setGeneric("rankExpr",
           function(expreMatrix, tiesMethod) standardGeneric("rankExpr"))

#' @rdname rankExpr
setMethod("rankExpr", signature('matrix'), function(expreMatrix){
  rankGenes(expreMatrix)
})
#' @rdname rankExpr
setMethod("rankExpr", signature('matrix', 'character'),
          function(expreMatrix, tiesMethod){
            rankGenes(expreMatrix,tiesMethod)
          })
#' @rdname rankExpr
setMethod("rankExpr", signature('data.frame'), function(expreMatrix){
  rankGenes(as.matrix(expreMatrix))
})
#' @rdname rankExpr
setMethod("rankExpr", signature('data.frame', 'character'),
          function(expreMatrix, tiesMethod){
            rankGenes(expreMatrix,tiesMethod)
          })

#' @rdname rankExpr
setMethod("rankExpr", signature('DGEList', 'character'),
          function(expreMatrix, tiesMethod){
            rankGenes(expreMatrix$counts,tiesMethod)
          })
#' @rdname rankExpr
setMethod("rankExpr", signature('DGEList'), function(expreMatrix){
  rankGenes(expreMatrix$counts)
})
#' @rdname rankExpr
setMethod("rankExpr", signature('ExpressionSet'), function(expreMatrix){
  rankGenes(exprs(expreMatrix))
})
#' @rdname rankExpr
setMethod("rankExpr", signature('ExpressionSet', 'character'),
          function(expreMatrix,tiesMethod){
            rankGenes(exprs(expreMatrix),tiesMethod)
          })
