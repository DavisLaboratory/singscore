#' @include singscore.R CoreFuns.R
NULL

#' Plot the densities of ranks for one sample 
#' @description This function plots the density and the rugs of the ranks for 
#' genes in the gene set in a sample. It takes a single column data frame, 
#' which is a one column subset of the ranked data obtained from [rankGenes()]
#' function, as well as the gene sets as inputs.
#'
#' @param rankData one column of the ranked gene expression matrix obtained from
#'   the [rankGenes()] function, use `drop = FALSE` when subsetting the ranked 
#'   gene expression matrix, see examples.
#' @param isInteractive Boolean, determine whether the returned plot is
#'   interactive
#' @param textSize numberic, set the size of text on the plot
#' @param upSet GeneSet object or a vector of gene Ids, up-regulated gene set
#' @param downSet GeneSet object or a vector of gene Ids, down-regulated gene 
#' set
#'
#' @return A ggplot object (optionally interactive) demonstrating the rank
#'   density along with rug plot
#' @examples
#' ranked <- rankGenes(toy_expr)
#' plotRankDensity(ranked[,2,drop = FALSE], upSet = toy_gs_up)
#' 
#'@export
setGeneric("plotRankDensity",
           function (rankData,
                        upSet,
                        downSet = NULL,
                        isInteractive = FALSE,
                        textSize = 1.5)
             standardGeneric("plotRankDensity"))

#' @rdname plotRankDensity
setMethod("plotRankDensity", signature(
  rankData = 'ANY',
  upSet = 'vector',
  downSet = 'missing'
),
function(rankData,
         upSet,
         downSet = NULL,
         isInteractive = FALSE,
         textSize = 1.5) {
  upSet <- GSEABase::GeneSet(as.character(upSet))
  plt <- plotRankDensity_intl(rankData,
                              upSet = upSet,
                              downSet = downSet,
                              isInteractive = isInteractive,
                              textSize = textSize)
  return(plt)
})

#' @rdname plotRankDensity
setMethod("plotRankDensity", signature(
  rankData = 'ANY',
  upSet = 'GeneSet',
  downSet = 'missing'
),
function(rankData,
         upSet,
         downSet = NULL,
         isInteractive = FALSE,
         textSize = 1.5) {
  plt <- plotRankDensity_intl(rankData,
                              upSet = upSet,
                              downSet = downSet,
                              isInteractive = isInteractive,
                              textSize = textSize)
  return(plt)
})
#' @rdname plotRankDensity
setMethod("plotRankDensity", signature(
  rankData = 'ANY',
  upSet = 'vector',
  downSet = 'vector'
),
function(rankData,
         upSet,
         downSet = NULL,
         isInteractive = FALSE,
         textSize = 1.5) {
  upSet <- GSEABase::GeneSet(as.character(upSet))
  downSet <- GSEABase::GeneSet(as.character(downSet))
  plt <- plotRankDensity_intl(rankData,
                              upSet = upSet,
                              downSet = downSet,
                              isInteractive = isInteractive,
                              textSize = textSize)
  return(plt)
})

#' @rdname plotRankDensity
setMethod("plotRankDensity", signature(
  rankData = 'ANY',
  upSet = 'GeneSet',
  downSet = 'GeneSet'
),
function(rankData,
         upSet,
         downSet = NULL,
         isInteractive = FALSE,
         textSize = 1.5) {
  
  plt <- plotRankDensity_intl(rankData,
                              upSet = upSet,
                              downSet = downSet,
                              isInteractive = isInteractive,
                              textSize = textSize)
  return(plt)
})

#' @rdname plotRankDensity
setMethod("plotRankDensity", signature(
  rankData = 'ANY',
  upSet = 'GeneSet',
  downSet = 'vector'
),
function(rankData,
         upSet,
         downSet = NULL,
         isInteractive = FALSE,
         textSize = 1.5) {
  downSet <- GSEABase::GeneSet(as.character(downSet))
  
  plt <- plotRankDensity_intl(rankData,
                              upSet = upSet,
                              downSet = downSet,
                              isInteractive = isInteractive,
                              textSize = textSize)
  return(plt)
})

#' @rdname plotRankDensity
setMethod("plotRankDensity", signature(
  rankData = 'ANY',
  upSet = 'vector',
  downSet = 'GeneSet'
),
function(rankData,
         upSet,
         downSet = NULL,
         isInteractive = FALSE,
         textSize = 1.5) {
  upSet <- GSEABase::GeneSet(as.character(upSet))
  
  plt <- plotRankDensity_intl(rankData,
                              upSet = upSet,
                              downSet = downSet,
                              isInteractive = isInteractive,
                              textSize = textSize)
  return(plt)
})
