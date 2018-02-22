#' @include singscore.R CoreFuns.R
NULL

#'@title single-sample gene-set scoring method
#'
#'@description This function takes a ranked gene expression matrix obtained from
#'  \code{rankGenes()} function and gene sets as inputs and it calculates the 
#'  scores for each individual sample against gene set. It returns a data.frame 
#'  consists of scores and dispersions for each sample. The gene sets can be in 
#'  vector format or GeneSet S4 object (GSEABase package).
#'  Down set can be null if down-regulated gene sets are unavailable.
#'
#' @param rankData A matrix-like object, ranked gene expression matrix data
#' @param subSamples A vector of sample labels/indices that will be
#'   used to subset the rankData matrix. All samples will be scored by default.
#' @param upSet A GeneSet object or vector of gene Ids of up-regulated gene set
#' @param downSet A GeneSet object or vector of gene Ids of down-regulated gene 
#'   set
#' @param centerScore A Boolean, specifying whether scores should be centered, 
#'   default as TRUE
#' @param dispersionFun A character, dispersion function with default as 'mad'
#' @return A data.frame consists of scores and dispersions for all samples
#'
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#' # toy_gs_up is a GeneSet object, alternatively a vector of gene ids may also 
#' # be supplied.
#'@seealso 
#'\code{\link{rank}} 
#'\code{"\linkS4class{GeneSet}"}
#'
#'@export
setGeneric("simpleScore",
            function(rankData,
                    upSet,
                    downSet = NULL,
                    subSamples = NULL,
                    centerScore = TRUE,
                    dispersionFun = 'mad')
             standardGeneric("simpleScore"))

#' @rdname simpleScore
setMethod("simpleScore", signature(
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
  df <- singscoring( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = centerScore,
                     dispersionFun = dispersionFun)
  return(df)
})

#' @rdname simpleScore
setMethod("simpleScore", signature(
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
  df <- singscoring( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = centerScore,
                     dispersionFun = dispersionFun)
  return(df)
})
#' @rdname simpleScore
setMethod("simpleScore", signature(
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
  df <- singscoring( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = centerScore,
                     dispersionFun = dispersionFun)
  return(df)
})

#' @rdname simpleScore
setMethod("simpleScore", signature(
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
 
  df <- singscoring( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = centerScore,
                     dispersionFun = dispersionFun)
  return(df)
})

#' @rdname simpleScore
setMethod("simpleScore", signature(
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
  
  df <- singscoring( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = centerScore,
                     dispersionFun = dispersionFun)
  return(df)
})

#' @rdname simpleScore
setMethod("simpleScore", signature(
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
  
  df <- singscoring( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = centerScore,
                     dispersionFun = dispersionFun)
  return(df)
})
