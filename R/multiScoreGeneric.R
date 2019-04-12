#' @include singscore.R rankAndScoring.R
NULL

#'@title single-sample gene-set scoring method for multiple signatures
#'
#'@description This function computes 'singscores' using a ranked gene
#'  expression matrix obtained from the [rankGenes()] function and a
#'  GeneSetCollection object or a list of GeneSet objects. It returns a list of
#'  two matrices containing the scores and dispersions. This function should be
#'  used when scoring needs to be performed for multiple signatures. It is
#'  faster than applying [simpleScore()] across the different signatures
#'  independently.
#'
#'@inheritParams simpleScore
#'@param upSetColc A GeneSetCollection object or a list of GeneSet objects of
#'  up-regulated (or mixed, see \code{\link{simpleScore}}) gene sets.
#'@param downSetColc A GeneSetCollection object or a list of GeneSet objects of
#'  down-regulated gene sets. NULL otherwise.
#'@return A list of two matrices containing the scores and dispersions
#'
#' @examples
#' ranked <- rankGenes(toy_expr_se)
#' setName(toy_gs_up)  = "toy_gs_up"
#' setName(toy_gs_dn)  = "toy_gs_dn"
#' gslist <- list(toy_gs_up, toy_gs_dn)
#' 
#' gscolc <- GSEABase::GeneSetCollection(gslist)
#' scoredf <- multiScore(ranked, upSetColc = gscolc)
#'@seealso \code{\link{rank}} \code{"\linkS4class{GeneSet}"}
#'
#'@export
setGeneric("multiScore",
           function(rankData,
                    upSetColc,
                    downSetColc,
                    subSamples = NULL,
                    centerScore = TRUE,
                    dispersionFun = mad,
                    knownDirection = TRUE)
             standardGeneric("multiScore"))

#' @rdname multiScore
setMethod("multiScore", signature(
  rankData = 'matrix',
  upSet = 'GeneSetCollection',
  downSet = 'missing'
),
function(rankData,
         upSetColc,
         downSetColc,
         subSamples = NULL,
         centerScore = TRUE,
         dispersionFun = mad,
         knownDirection = TRUE) {
  
  stopifnot(is.logical(centerScore), is.logical(knownDirection))
  matlist = multiSingscore(
    rankData,
    upSetColc,
    downSetColc = NULL,
    subSamples = subSamples,
    centerScore = centerScore,
    dispersionFun = mad,
    knownDirection = knownDirection
  )
  return(matlist)
})

#' @rdname multiScore
setMethod("multiScore", signature(
  rankData = 'matrix',
  upSet = 'GeneSetCollection',
  downSet = 'GeneSetCollection'
),
function(rankData,
         upSetColc,
         downSetColc,
         subSamples = NULL,
         centerScore = TRUE,
         dispersionFun = mad,
         knownDirection = TRUE) {
  
  stopifnot(is.logical(centerScore), is.logical(knownDirection))
  matlist = multiSingscore(
    rankData,
    upSetColc,
    downSetColc,
    subSamples = subSamples,
    centerScore = centerScore,
    dispersionFun = mad,
    knownDirection = knownDirection
  )
  return(matlist)
})

#' @rdname multiScore
setMethod("multiScore", signature(
  rankData = 'matrix',
  upSet = 'list',
  downSet = 'missing'
),
function(rankData,
         upSetColc,
         downSetColc,
         subSamples = NULL,
         centerScore = TRUE,
         dispersionFun = mad,
         knownDirection = TRUE) {
  
  stopifnot(is.logical(centerScore), is.logical(knownDirection))
  stopifnot(all(lapply(upSetColc, class) %in% 'GeneSet'))
  
  upSetColc = GSEABase::GeneSetCollection(upSetColc)
  
  matlist = multiSingscore(
    rankData,
    upSetColc,
    downSetColc = NULL,
    subSamples = subSamples,
    centerScore = centerScore,
    dispersionFun = mad,
    knownDirection = knownDirection
  )
  return(matlist)
})

#' @rdname multiScore
setMethod("multiScore", signature(
  rankData = 'matrix',
  upSet = 'list',
  downSet = 'list'
),
function(rankData,
         upSetColc,
         downSetColc,
         subSamples = NULL,
         centerScore = TRUE,
         dispersionFun = mad,
         knownDirection = TRUE) {
  
  stopifnot(is.logical(centerScore), is.logical(knownDirection))
  stopifnot(all(lapply(upSetColc, class) %in% 'GeneSet'))
  stopifnot(all(lapply(downSetColc, class) %in% 'GeneSet'))
  
  upSetColc = GSEABase::GeneSetCollection(upSetColc)
  downSetColc = GSEABase::GeneSetCollection(downSetColc)
  
  matlist = multiSingscore(
    rankData,
    upSetColc,
    downSetColc,
    subSamples = subSamples,
    centerScore = centerScore,
    dispersionFun = mad,
    knownDirection = knownDirection
  )
  return(matlist)
})
