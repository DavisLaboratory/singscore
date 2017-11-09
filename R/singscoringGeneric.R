#' @include singscore.R CoreFuns.R
NULL

#'@title single-sample gene-set scoring method
#'@description This function takes a ranked gene expression matrix obtained from
#'  \code{rankExpr} (or \code{rankGenes}) function and gene sets as input and
#'  calulate the scores for each individual sample against gene set. It returns
#'  a data.frame consists of scores and dispersions for each sample. The gene
#'  sets can be in vector format or GeneSet S4 object. Down set can be null if 
#'  no down-regulated gene sets are available.
#'
#' @param rankData A matrix-like object, ranked gene expression matrix data
#' @param subSamples A vector of sample labels/indices that will be
#'   used to subset the rankData matrix. All samples will be scored by default.
#' @param upSet A GeneSet object or vector of gene ids of up regulated gene set
#' @param downSet A GeneSet object or vector of gene ids of down regulated gene 
#'   set
#' @param centerScore A Boolean, specifying whether scores should be centred, 
#'   default as TRUE
#' @param dispersionFun A character, dispersion function with default as 'mad'
#' @return A data.frame consists of scores and dispersions for all samples
#'
#' @examples
#' ranked <- rankExpr(toy_expr)
#' scoredf <- singscoring(ranked, upSet = toy_up, downSet = toy_dn)
#'
#'@seealso \code{\link{rank}} \linkS4class{GeneSet}
#'
#'@export
setGeneric("singscoring",
           function(rankData,
                    upSet,
                    downSet = NULL,
                    subSamples = NULL,
                    centerScore = TRUE,
                    dispersionFun = 'mad')
             standardGeneric("singscoring"))

#' @rdname singscoring
setMethod("singscoring", signature(
  rankData = 'ANY',
  upSet = 'vector',
  downSet = 'missing'
),
function(rankData,
         upSet,
         downSet = NULL,
         subSamples = NULL,
         centerScore = TRUE,
         dispersionFun = 'mad') {
  upSet <- GSEABase::GeneSet(as.character(upSet))
  df <- simpleScore( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = TRUE,
                     dispersionFun = 'mad')
  return(df)
})

#' @rdname singscoring
setMethod("singscoring", signature(
  rankData = 'ANY',
  upSet = 'GeneSet',
  downSet = 'missing'
),
function(rankData,
         upSet,
         downSet = NULL,
         subSamples = NULL,
         centerScore = TRUE,
         dispersionFun = 'mad') {
  df <- simpleScore( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = TRUE,
                     dispersionFun = 'mad')
  return(df)
})
#' @rdname singscoring
setMethod("singscoring", signature(
  rankData = 'ANY',
  upSet = 'vector',
  downSet = 'vector'
),
function(rankData,
         upSet,
         downSet = NULL,
         subSamples = NULL,
         centerScore = TRUE,
         dispersionFun = 'mad') {
  upSet <- GSEABase::GeneSet(as.character(upSet))
  downSet <- GSEABase::GeneSet(as.character(downSet))
  if( ! is.null(downSet)){
  }
  df <- simpleScore( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = TRUE,
                     dispersionFun = 'mad')
  return(df)
})

#' @rdname singscoring
setMethod("singscoring", signature(
  rankData = 'ANY',
  upSet = 'GeneSet',
  downSet = 'GeneSet'
),
function(rankData,
         upSet,
         downSet = NULL,
         subSamples = NULL,
         centerScore = TRUE,
         dispersionFun = 'mad') {
 
  df <- simpleScore( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = TRUE,
                     dispersionFun = 'mad')
  return(df)
})

#' @rdname singscoring
setMethod("singscoring", signature(
  rankData = 'ANY',
  upSet = 'GeneSet',
  downSet = 'vector'
),
function(rankData,
         upSet,
         downSet = NULL,
         subSamples = NULL,
         centerScore = TRUE,
         dispersionFun = 'mad') {
  downSet <- GSEABase::GeneSet(as.character(downSet))
  
  df <- simpleScore( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = TRUE,
                     dispersionFun = 'mad')
  return(df)
})

#' @rdname singscoring
setMethod("singscoring", signature(
  rankData = 'ANY',
  upSet = 'vector',
  downSet = 'GeneSet'
),
function(rankData,
         upSet,
         downSet,
         subSamples = NULL,
         centerScore = TRUE,
         dispersionFun = 'mad') {
  upSet <- GSEABase::GeneSet(as.character(upSet))
  
  df <- simpleScore( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = TRUE,
                     dispersionFun = 'mad')
  return(df)
})
