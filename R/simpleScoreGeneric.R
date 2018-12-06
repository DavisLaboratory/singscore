#' @include singscore.R rankAndScoring.R
NULL

#'@title single-sample gene-set scoring method
#'
#'@description This function computes 'singscores' using a ranked gene
#'  expression matrix obtained from the [rankGenes()] function and a gene set or
#'  a pair of up-regulated and down-regulated gene sets. It returns a data.frame
#'  of scores and dispersions for each sample. The gene sets can be in vector
#'  format or as GeneSet objects (from GSEABase packages). If samples need to be
#'  scored against a single gene set, the \code{upSet} argument should be used
#'  to pass the gene set while the \code{downSet} argument is set to
#'  \code{NULL}. This setting is ideal for gene sets representing gene
#'  ontologies where the nature of the genes is unknown (up- or down-regulated).
#'  
#'@param rankData A matrix object, ranked gene expression matrix data generated
#'  using the [rankGenes()] function
#'@param subSamples A vector of sample labels/indices that will be used to
#'  subset the rankData matrix. All samples will be scored if not provided
#'@param upSet A GeneSet object or character vector of gene IDs of up-regulated
#'  gene set or a gene set where the nature of genes is not known
#'@param downSet A GeneSet object or character vector of gene IDs of
#'  down-regulated gene set or NULL where only a single gene set is provided
#'@param centerScore A Boolean, specifying whether scores should be centered
#'  around 0, default as TRUE
#'@param dispersionFun A function, dispersion function with default being `mad`
#'@param knownDirection A boolean flag, it deterimines whether the scoring
#'  method should derive the scores in a directional mannar when the gene
#'  signature only contains one set of gene set (passing the gene set via
#'  upSet). It is default as TRUE but one can set the argument to be FALSE to
#'  derive the score for a single gene set in a undirectional way. This
#'  parameter becomes irrelevant when both upSet and downSet are provided.
#'@return A data.frame consists of singscores and dispersions for all samples
#'
#' @examples
#' ranked <- rankGenes(toy_expr_se)
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
                    dispersionFun = mad,
                    knownDirection = TRUE)
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
         dispersionFun = mad,
         knownDirection = TRUE) {
  stopifnot(is.logical(centerScore), is.logical(knownDirection))
  upSet <- GSEABase::GeneSet(as.character(upSet))
  if(knownDirection){
    df <- singscoring( rankData,
                       upSet = upSet,
                       downSet = downSet,
                       subSamples = subSamples,
                       centerScore = centerScore,
                       dispersionFun = dispersionFun)
  } else {
    df <- singscoringOneGS(rankData,
                           upSet = upSet,
                           subSamples = subSamples,
                           centerScore = centerScore,
                           dispersionFun = dispersionFun)
  }
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
         dispersionFun = mad,
         knownDirection = TRUE) {
  stopifnot(is.logical(centerScore), is.logical(knownDirection))
  if(knownDirection){
    df <- singscoring(
      rankData,
      upSet = upSet,
      subSamples = subSamples,
      centerScore = centerScore,
      dispersionFun = dispersionFun
    )
  } else {
    df <- singscoringOneGS(
      rankData,
      upSet = upSet,
      subSamples = subSamples,
      dispersionFun = dispersionFun,
      centerScore = centerScore
    )
  }
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
         dispersionFun = mad,
         knownDirection = TRUE) {
  stopifnot(is.logical(centerScore), is.logical(knownDirection))
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
         dispersionFun = mad,
         knownDirection = TRUE) {
  stopifnot(is.logical(centerScore), is.logical(knownDirection))
  df <- singscoring( rankData,
                     upSet = upSet,
                     downSet = downSet,
                     subSamples = subSamples,
                     centerScore = centerScore)
  return(df)
})


