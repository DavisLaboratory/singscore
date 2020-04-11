#' @include singscore.R rankAndScoring.R
NULL

#'@title  Rank genes by the gene expression intensities
#'@description The \code{rankGenes} function is a generic function that can deal
#'  with mutilple types of inputs. Given a matrix of gene expression that has
#'  samples in columns, genes in rows, and values being gene expression
#'  intensity,\code{rankGenes} ranks gene expression intensities in each sample.
#'
#'  It can also work with S4 objects that have gene expression matrix as a
#'  component (i.e ExpressionSet, DGEList,SummarizedExperiment). It calls the
#'  \code{rank} function in the base package which ranks the gene expression
#'  matrix by its absolute expression level. If the input is S4 object of
#'  \code{DGEList, ExpressionSet, or SummarizedExperiment}, it will extract the
#'  gene expression matrix from the object and rank the genes. The default
#'  'tiesMethod' is set to 'min'.
#'
#'@param expreMatrix matrix, data.frame, ExpressionSet, DGEList or
#'  SummarizedExperiment storing gene expression measurements
#'@param tiesMethod character, indicating what method to use when dealing with
#'  ties
#'@param stableGenes character, containing a list of stable genes to be used to
#'  rank genes using expression of stable genes. This is required when using the
#'  stable genes dependent version of singscore
#'
#'@name rankGenes
#'
#'@return The ranked gene expression matrix that has samples in columns and
#'  genes in rows
#'@examples
#' rankGenes(toy_expr_se) # toy_expr_se is a gene expression dataset
#' 
#' # ExpressionSet object
#' emat <- SummarizedExperiment::assay(toy_expr_se)
#' e <- Biobase::ExpressionSet(assayData = as.matrix(emat))
#' rankGenes(e)
#' 
#' #scoring using the stable version of singscore
#' rankGenes(e, stableGenes = c('2', '20', '25'))
#' 
#' \dontrun{
#' #for real cancer or blood datasets, use getStableGenes()
#' rankGenes(cancer_expr, stableGenes = getStableGenes(5))
#' rankGenes(blood_expr, stableGenes = getStableGenes(5, type = 'blood'))
#' }
#'@seealso \code{\link{rank}} [ExpressionSet][ExpressionSet-class]
#'  [SummarizedExperiment][SummarizedExperiment-class] [DGEList][DGEList-class]
#'
#'
#'@export
setGeneric("rankGenes",
           function(expreMatrix, 
                    tiesMethod = 'min',
                    stableGenes = NULL) standardGeneric("rankGenes"))

#' @rdname rankGenes
setMethod("rankGenes",
          signature('matrix','ANY', 'ANY'), 
          function(expreMatrix, tiesMethod = 'min', stableGenes = NULL){
            stopifnot(tiesMethod %in% c("max", "average","min"))
            if (is.null(stableGenes)) {
              rankMat = rankExpr(expreMatrix, tiesMethod)
            } else {
              stopifnot(is.character(stableGenes))
              stopifnot(length(stableGenes) > 0)
              rankMat = rankExprStable(expreMatrix, tiesMethod, stableGenes)
            }
            return(rankMat)
})

#' @rdname rankGenes
setMethod("rankGenes", 
          signature('data.frame', 'ANY', 'ANY'), 
          function(expreMatrix, tiesMethod = 'min', stableGenes = NULL){
            rankGenes(as.matrix(expreMatrix), tiesMethod, stableGenes)
})

#' @rdname rankGenes
setMethod("rankGenes", 
          signature('DGEList','ANY', 'ANY'),
          function(expreMatrix, tiesMethod = 'min', stableGenes = NULL){
            rankGenes(expreMatrix$counts, tiesMethod, stableGenes)
})
#' @rdname rankGenes
setMethod("rankGenes", 
          signature('ExpressionSet', 'ANY', 'ANY'),
          function(expreMatrix, tiesMethod = 'min', stableGenes = NULL){
            rankGenes(Biobase::exprs(expreMatrix), tiesMethod, stableGenes)
})

#' @rdname rankGenes
setMethod("rankGenes",
          signature('SummarizedExperiment', 'ANY', 'ANY'),
          function(expreMatrix, tiesMethod = 'min', stableGenes = NULL){
            rankGenes(SummarizedExperiment::assay(expreMatrix), tiesMethod, stableGenes)
          })
