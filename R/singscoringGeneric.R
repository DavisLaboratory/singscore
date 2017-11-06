#' @include singscore.R CoreFuns.R
NULL

#'@title Calculate scores for the ranked gene expression matrix against gene set
#'@description This function takes a ranked gene expression matrix obtained
#'   from \code{rankExpr} (or \code{rankGenes}) function and gene sets as input.
#'   It returns a data.frame consists of scores and dispersions for each sample.
#'   The gene sets can be in vector format or GeneSet S4 object.
#'
#' @param rankData A matrix-like object, ranked gene expression matrix data
#' @param subSamples A character or vector of sample labels/indices that will be
#'   used to subset the rankData matrix
#' @param upSet A GeneSet object or vector of gene ids, up regulated gene set
#' @param downSet A GeneSet object or vector of gene ids, down regulated gene set
#' @param centerScore A Boolean, specifying whether scores should be centred
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
           function(rankData,subSamples,upSet,downSet,centerScore,dispersionFun) standardGeneric("singscoring"))

#' @rdname singscoring
setMethod("singscoring", signature(rankData = 'matrix',subSamples='ANY',upSet ='vector',downSet = 'vector', centerScore='ANY',dispersionFun='ANY'),
          function(rankData,subSamples,upSet,downSet,centerScore,dispersionFun){
            upSet <- GSEABase::GeneSet(as.character(upSet))
            downSet <- GSEABase::GeneSet(as.character(downSet))

            df <- simpleScore(rankData,upSet = upSet, downSet = downSet)
            return(df)
            })
