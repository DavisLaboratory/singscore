#' @include singscore.R
NULL

#' Rank gene expression matrix
#' @description Given a gene expression matrix and a tiesMethod (character),
#'   this fucntion calls the 'rank' function from 'base' package which ranks the
#'   gene expression matrix by gene's expression level. The default tiesMethod
#'   is 'min'. There is a generic version of this function, details can be found
#'   via the link in the see also section down at the bottom. It is suggested to
#'   use the generic function 'rankExpr' which can accept multiple data formats
#'   as input.
#' @param exprsM A matrix, gene expression matrix
#' @param tiesMethod A character, default as 'min'
#' @return A matrix that has samples in colunm and genes in rows. Values are the
#'   ranks of each gene in each sample.
#' @keywords internal
#' @seealso
#' [rankGenes()]
#' @examples
#' \dontrun{ranked <- rankExpr(toy_expr_se)}
rankExpr <- function(exprsM, tiesMethod = "min") {
  rname<- rownames(exprsM)
  cname <- colnames(exprsM)
  rankedData <- matrixStats::colRanks(as.matrix(exprsM),
                                      ties.method = tiesMethod,
                                      preserveShape = TRUE)
  rownames(rankedData) <- rname
  colnames(rankedData) <- cname
  return (rankedData)
}

#define the main helper functions
#helper: check if all genes in the gene set are in the data
checkGenes <- function (geneset, background) {
  #check if there are some missing genes in the geneset
  missingGenes = setdiff(geneset, background)
  if (length(missingGenes) > 0) {
    warningMsg = paste(length(missingGenes), "genes missing:", sep = ' ')
    warningMsg = paste(warningMsg, paste(missingGenes, collapse = ', '), sep = ' ')
    warning(warningMsg)
  }

  return(setdiff(geneset, missingGenes))
}

#helper: compute the theoretical boundaries of scores when direction is KNOWN.
# Should be used after filtering out missing genes
calcBounds <- function(gsSize, bgSize) {
  lowBound = (gsSize + 1) / 2
  upBound = (2 * bgSize - gsSize + 1) / 2

  return(list('lowBound' = lowBound, 'upBound' = upBound))
}

#helper: compute the theoretical boundaries of scores when direction is UNKNOWN.
# Should be used after filtering out missing genes
calcBoundsUnknownDir <- function(gsSize, bgSize) {
  gsSize = floor(gsSize / 2) #number of unique ranks in the geneset
  bgSize = ceiling(bgSize / 2) #number of unique ranks in the matrix

  return(calcBounds(gsSize, bgSize))
}

#helper: empty data frame of results
emptyScoreDf <- function(onegs = TRUE) {
  if (onegs) {
    df = data.frame("TotalScore" = numeric(),
                    "TotalDispersion" = numeric())
  } else{
    df = data.frame(
      "TotalScore" = numeric(),
      "TotalDispersion" = numeric(),
      "UpScore" = numeric(),
      "UpDispersion" = numeric(),
      "DownScore" = numeric(),
      "DownDispersion" = numeric()
    )
  }

  return(df)
}

#helper: compute scores for one gene set given a rank matrix and bounds
# Assume genes are filtered
calcScores <- function(ranks, geneset, dispersionFun, bounds, scoffset = 0) {
  #compute raw score
  gsranks = ranks[geneIds(geneset), , drop = FALSE]
  gsscore = colMeans(gsranks)
  gsscore = (gsscore - bounds$lowBound) / (bounds$upBound - bounds$lowBound)

  #calculate dispersion
  dispersion = apply(gsranks, 2, dispersionFun)

  return(list('TotalScore' = gsscore + scoffset, 'TotalDispersion' = dispersion))
}

#helper: compute scores for two gene sets given a rank matrix and bounds
# Assume genes are filtered
calcScoresUpDn <- function(ranks, upset, downset, dispersionFun, upbounds, downbounds, scoffset = 0) {
  #compute independent scores
  upscores = calcScores(ranks, upset, dispersionFun, upbounds, scoffset)
  dnscores = calcScores(ranks, downset, dispersionFun, downbounds, 0)
  dnscores$TotalScore = 1 - dnscores$TotalScore + scoffset

  #compute total scores
  totalscores = list(
    'TotalScore' = upscores$TotalScore + dnscores$TotalScore,
    'TotalDispersion' = upscores$TotalDispersion + dnscores$TotalDispersion,
    'UpScore' = upscores$TotalScore,
    'UpDispersion' = upscores$TotalDispersion,
    'DownScore' = dnscores$TotalScore,
    'DownDispersion' = dnscores$TotalDispersion
  )

  return(totalscores)
}

singleSingscore <-
  function (rankData,
            upSet,
            downSet = NULL,
            subSamples = NULL,
            centerScore = TRUE,
            dispersionFun = mad,
            knownDirection = TRUE) {

    #1. subset the data for samples whose calculation is to be performed
    if (!is.null(subSamples)) {
      rankData <- rankData[, subSamples, drop = FALSE]
    }

    #2. center scores?
    center_offset = ifelse(centerScore, -0.5, 0)

    #3.1 filter genes - return empty data.frame if no scores
    geneIds(upSet) = checkGenes(geneIds(upSet), rownames(rankData))
    if (length(geneIds(upSet)) == 0)
      return(emptyScoreDf(onegs = is.null(downSet)))

    #known direction
    if (!knownDirection) {
      #2.3 do not center scores
      warning('\'centerScore\' is disabled for this setting')
      center_offset = 0

      #modify rank matrix
      rankData = apply(rankData, 2, function(x) {
        x = abs(x - ceiling(median(x)))
        return(x)
      })

      #4.3 calculate bounds
      upset_bounds = calcBoundsUnknownDir(length(geneIds(upSet)), nrow(rankData))
    } else {
      #4.1 calculate bounds
      upset_bounds = calcBounds(length(geneIds(upSet)), nrow(rankData))
    }

    #for two gene sets
    if (!is.null(downSet)) {
      #3.2 filter genes - return empty data.frame if no scores
      geneIds(downSet) = checkGenes(geneIds(downSet), rownames(rankData))
      if (length(geneIds(downSet)) == 0)
        return(emptyScoreDf(onegs = FALSE))

      #4.2 calculate bounds
      downset_bounds = calcBounds(length(geneIds(downSet)), nrow(rankData))

      #5.2 compute scores
      scores = calcScoresUpDn(
        rankData,
        upSet,
        downSet,
        dispersionFun,
        upset_bounds,
        downset_bounds,
        center_offset
      )
    } else{
      #5.1 & 5.3 compute scores
      scores = calcScores(rankData,
                          upSet,
                          dispersionFun,
                          upset_bounds,
                          center_offset)
    }

    #6 process the results
    scores = data.frame(scores)
    rownames(scores) = colnames(rankData)

    return(scores)
  }
