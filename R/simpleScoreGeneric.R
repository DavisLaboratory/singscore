#' @include singscore.R rankAndScoring.R
NULL

#'@title single-sample gene-set scoring method
#'
#'@description This function computes 'singscores' using an **unmodified**ranked
#'  gene expression matrix obtained from the [rankGenes()] function and a gene
#'  set or a pair of up-regulated and down-regulated gene sets. It returns a
#'  data.frame of scores and dispersions for each sample. The gene sets can be
#'  in vector format or as GeneSet objects (from GSEABase packages). If samples
#'  need to be scored against a single gene set, the \code{upSet} argument
#'  should be used to pass the gene set while the \code{downSet} argument is set
#'  to \code{NULL}. This setting is ideal for gene sets representing gene
#'  ontologies where the nature of the genes is unknown (up- or down-regulated).
#'
#'@param rankData A matrix object, ranked gene expression matrix data generated
#'  using the [rankGenes()] function (make sure this matrix is not modified)
#'@param subSamples A vector of sample labels/indices that will be used to
#'  subset the rankData matrix. All samples will be scored if not provided
#'@param upSet A GeneSet object or character vector of gene IDs of up-regulated
#'  gene set or a gene set where the nature of genes is not known
#'@param downSet A GeneSet object or character vector of gene IDs of
#'  down-regulated gene set or NULL where only a single gene set is provided
#'@param centerScore A Boolean, specifying whether scores should be centered
#'  around 0, default as TRUE. Note: scores never centered if `knownDirection =
#'  FALSE`
#'@param dispersionFun A function, dispersion function with default being `mad`
#'@param knownDirection A boolean, determining whether the gene set should be
#'  considered to be directional or not. A gene set is directional if the type
#'  of genes in it are known i.e. up- or down-regulated. This should be set to
#'  TRUE if the gene set is composed of both up- AND down-regulated genes.
#'  Defaults to TRUE. This parameter becomes irrelevant when both upSet(Colc)
#'  and downSet(Colc) are provided.
#'@return A data.frame consists of singscores and dispersions for all samples
#'
#' @examples
#' ranked <- rankGenes(toy_expr_se)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#' # toy_gs_up is a GeneSet object, alternatively a vector of gene ids may also
#' # be supplied.
#'@seealso \code{\link{rank}} \code{"\linkS4class{GeneSet}"}
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
  rankData = 'matrix',
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
  upSet = GSEABase::GeneSet(as.character(upSet))
  df = singleSingscore(
    rankData,
    upSet,
    downSet = downSet,
    subSamples = subSamples,
    centerScore = centerScore,
    dispersionFun = mad,
    knownDirection = knownDirection
  )
  return(df)
})

#' @rdname simpleScore
setMethod("simpleScore", signature(
  rankData = 'matrix',
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
  df = singleSingscore(
    rankData,
    upSet,
    downSet = downSet,
    subSamples = subSamples,
    centerScore = centerScore,
    dispersionFun = mad,
    knownDirection = knownDirection
  )
  return(df)
})

#' @rdname simpleScore
setMethod("simpleScore", signature(
  rankData = 'matrix',
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
  upSet = GSEABase::GeneSet(as.character(upSet))
  downSet = GSEABase::GeneSet(as.character(downSet))
  df = singleSingscore(
    rankData,
    upSet,
    downSet = downSet,
    subSamples = subSamples,
    centerScore = centerScore,
    dispersionFun = mad,
    knownDirection = knownDirection
  )
  return(df)
})

#' @rdname simpleScore
setMethod("simpleScore", signature(
  rankData = 'matrix',
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
  df = singleSingscore(
    rankData,
    upSet,
    downSet = downSet,
    subSamples = subSamples,
    centerScore = centerScore,
    dispersionFun = mad,
    knownDirection = knownDirection
  )

  return(df)
})
