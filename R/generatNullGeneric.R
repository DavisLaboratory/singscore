# generateNull generic function
#' @include singscore.R permuTest.R
NULL
#' @title Permutation test for the derived scores of each sample
#'
#' @description This function generates a number of random gene sets that
#'   have the same number of genes as the scored gene set. It scores each random
#'   gene set and returns a matrix of scores for all samples.
#'   The empirical scores are used to calculate the empirical p-values and plot
#'   the null distribution. The implementation uses [BiocParallel::bplapply()]
#'   for easy access to parallel backends. Note that one should pass the same
#'   values to the `upSet`, `downSet`, `centerScore` and `bidirectional`
#'   arguments as what they provide for the `simpleScore()` function to generate
#'   a proper null distribution.
#' 
#' @inheritParams simpleScore
#' @param B integer, the number of permutation repeats or the number of random
#' gene sets to be generated, default as 1000
#' @param ncores, integer, the number of CPU cores the function can use
#' @param seed integer, set the seed for randomisation
#' @param useBPPARAM, the backend the function uses, if NULL is provided, the
#' function uses the default parallel backend which is the first on the list
#' returned by \code{BiocParallel::registered()} i.e
#' \code{BiocParallel::registered()[[1]]} for your machine. It can be changed
#' explicitly by passing a BPPARAM
#'
#' @return A matrix of empirical scores for all samples
#' @seealso
#' [Post about BiocParallel](http://lcolladotor.github.io/2016/03/07/BiocParallel/#.WgXMF61L28U)
#' `browseVignettes("BiocParallel")`
#' @author Ruqian Lyu
#' @export
#' @examples
#' ranked <- rankGenes(toy_expr_se)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#'
#' # find out what backends can be registered on your machine
#' BiocParallel::registered()
#' # the first one is the default backend
#' # ncores = ncores <- parallel::detectCores() - 2
#' permuteResult = generateNull(upSet = toy_gs_up, downSet = toy_gs_dn, ranked,
#' centerScore = TRUE, B =10, seed = 1, ncores = 1 )
setGeneric("generateNull",
           function(upSet,
                    downSet = NULL,
                    rankData,
                    subSamples = NULL,
                    centerScore = TRUE,
                    knownDirection = TRUE,
                    B = 1000,
                    ncores = 1,
                    seed = sample.int(1E6, 1),
                    useBPPARAM = NULL)
             standardGeneric("generateNull"))

#' @rdname generateNull
setMethod("generateNull", signature(
  upSet = 'vector',
  downSet = 'missing'
),
function(upSet,
         downSet = NULL,
         rankData,
         subSamples = NULL,
         centerScore = TRUE,
         knownDirection = TRUE,
         B = 1000,
         ncores = 1,
         seed = sample.int(1E6, 1),
         useBPPARAM = NULL) {

  upSet <- GSEABase::GeneSet(as.character(upSet))
  plt <- generateNull(
    upSet = upSet,
    downSet = downSet,
    rankData = rankData,
    subSamples = subSamples,
    centerScore = centerScore,
    knownDirection = knownDirection,
    B = B,
    ncores = ncores,
    seed = seed,
    useBPPARAM = useBPPARAM
  )
  return(plt)
})

#' @rdname generateNull
setMethod("generateNull", signature(
  upSet = 'GeneSet',
  downSet = 'missing'
),
function(upSet,
         downSet = NULL,
         rankData,
         subSamples = NULL,
         centerScore = TRUE,
         knownDirection = TRUE,
         B = 1000,
         ncores = 1,
         seed = sample.int(1E6, 1),
         useBPPARAM = NULL) {

  stopifnot(is.logical(centerScore), is.logical(knownDirection), B%%1==0,
            ncores%%1==0)
  plt <- generateNull_intl(
    upSet = upSet,
    downSet = downSet,
    rankData = rankData,
    subSamples = subSamples,
    centerScore = centerScore,
    knownDirection = knownDirection,
    B = B,
    ncores = ncores,
    seed = seed,
    useBPPARAM = useBPPARAM
  )
  return(plt)
})

#' @rdname generateNull
setMethod("generateNull", signature(
  upSet = 'vector',
  downSet = 'vector'
),
function(upSet,
         downSet = NULL,
         rankData,
         subSamples = NULL,
         centerScore = TRUE,
         knownDirection = TRUE,
         B = 1000,
         ncores = 1,
         seed = sample.int(1E6, 1),
         useBPPARAM = NULL) {

  upSet <- GSEABase::GeneSet(as.character(upSet))
  downSet <- GSEABase::GeneSet(as.character(downSet))
  plt <- generateNull(
    upSet = upSet,
    downSet = downSet,
    rankData = rankData,
    subSamples = subSamples,
    centerScore = centerScore,
    knownDirection = knownDirection,
    B = B,
    ncores = ncores,
    seed = seed,
    useBPPARAM = useBPPARAM
  )
  return(plt)
})

#' @rdname generateNull
setMethod("generateNull", signature(
  upSet = 'GeneSet',
  downSet = 'GeneSet'
),
function(upSet,
         downSet = NULL,
         rankData,
         subSamples = NULL,
         centerScore = TRUE,
         knownDirection = TRUE,
         B = 1000,
         ncores = 1,
         seed = sample.int(1E6, 1),
         useBPPARAM = NULL) {

  stopifnot(is.logical(centerScore), is.logical(knownDirection), B%%1==0,
            ncores%%1==0)
  plt <- generateNull_intl(
    upSet = upSet,
    downSet = downSet,
    rankData = rankData,
    subSamples = subSamples,
    centerScore = centerScore,
    knownDirection = knownDirection,
    B = B,
    ncores = ncores,
    seed = seed,
    useBPPARAM = useBPPARAM
  )
  return(plt)
})
