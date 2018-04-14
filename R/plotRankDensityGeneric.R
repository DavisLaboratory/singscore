#' @include singscore.R plot.R
NULL

#' Plot the densities of ranks for one sample
#' @description This function takes a single-column data frame, which is a
#'   single-column subset of the ranked matrix data generated using
#'   [rankGenes()] function, and the gene sets of interest as inputs. It plots
#'   the density of ranks for genes in the gene set and overlays a barcode plot
#'   of these ranks. Ranks are normalized by dividing them by the maximum rank.
#'   Densities are estimated using KDE.
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
#' @return  A ggplot object (or a plotly object) with a rank density plot
#'  overlayed with a barcode plot
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
