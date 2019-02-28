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

################################################################################
####========================= singscoring() function ===========================
################################################################################

#' Single-sample Gene-set scoring method
#'
#' @description This function takes a ranked gene expression matrix obtained
#'   from \code{rankGenes()}  function and a GeneSet object as input parameters.
#'   It returns a data.frame consists of calculated scores and dispersions for
#'   each sample. It is suggested to use the generic version of of this method
#'   which can work with gene set stored in both vector or GeneSet.
#'
#' @param rankData A matrix, ranked gene expression matrix data
#' @param upSet A GeneSet object, up regulated gene set
#' @param downSet A GeneSet object, down regulated gene set
#' @param subSamples A character or vector of sample labels/indices that will be
#'   used to subset the rankData matrix.All samples will be scored by default.
#' @param centerScore A Boolean, specifying whether scores should be centred
#' @param dispersionFun A function, dispersion function with default as 'mad'
#' @return A data.frame consists of scores and dispersions for all samples
#' @keywords internal
#' @seealso
#' [simpleScore()]
#' \code{"\linkS4class{GeneSet}"}

singscoring <- function (rankData, upSet, downSet = NULL, subSamples = NULL,
                         centerScore = TRUE, dispersionFun = mad) {

  #subset the data for samples whose calculation is to be performed
  if (!is.null(subSamples)) {
    rankData <- rankData[, subSamples, drop = FALSE]
  }

  # #check genes in the gene set - return empty results if 0 genes are measured
  # upSet = checkGenes(geneIds(upSet), rownames(rankData))
  # if (length(upSet) == 0)
  #   return(emptyScoreDf(onegs = is.null(downSet)))

  #values needed for calculating the boundaries
  upSigSize <- length(geneIds(upSet))
  nTotalGenes <- nrow(rankData)

  #check if there are some missing genes in the geneset
  missingGenes <- setdiff(geneIds(upSet), rownames(rankData))
  if (length(missingGenes) > 0) {
    warningMsg <-
      paste(length(missingGenes), "genes missing:", sep = ' ')
    warningMsg <- paste(warningMsg,
                        paste(missingGenes, collapse = ', '), sep = ' ')
    warning(warningMsg)
  }

  #remove missing genes from signature for further analysis
  geneIds(upSet) <- setdiff(geneIds(upSet), missingGenes)
  upRanks <- rankData[geneIds(upSet), , drop = FALSE]
  upScore <- colMeans(upRanks)
  lowBound <- (upSigSize + 1) / 2
  upBound  <- (2 * nTotalGenes - upSigSize + 1) / 2
  normUpScore <- (upScore - lowBound) / (upBound - lowBound)

  #calculate dispersion
  upDispersion <- apply(upRanks, 2, dispersionFun)

  #if just the up score is provided
  if (is.null(downSet)) {
    #if centering
    if (centerScore) {
      normUpScore <- normUpScore - 0.5
    }

    scoredf <-
      data.frame("TotalScore" = normUpScore, "TotalDispersion" = upDispersion)
    rownames(scoredf) <- colnames(rankData)
    return(scoredf)
  } else {
    scoreDown <-
      singscoring(rankData, downSet, NULL, subSamples, FALSE, dispersionFun)
    normDownScore <- 1 - scoreDown$TotalScore
    downDispersion <- scoreDown$TotalDispersion

    #if centering
    if (centerScore) {
      normUpScore <- normUpScore - 0.5
      normDownScore <- normDownScore - 0.5
    }

    scoredf <-
      data.frame(
        "TotalScore" = normUpScore + normDownScore,
        "TotalDispersion" = upDispersion + downDispersion,
        "UpScore" = normUpScore,
        "UpDispersion" = upDispersion,
        "DownScore" = normDownScore,
        "DownDispersion" = downDispersion
      )
    rownames(scoredf) <- colnames(rankData)
    return(scoredf)
  }
}

################################################################################
#### =============================== singscoringOneGS()=========================
################################################################################
#' scoring single gene set signature when direction of gene set is unknown
#'
#' @param rankData A matrix, ranked gene expression matrix data
#' @param upSet A GeneSet object, up regulated gene set
#' @param downSet A GeneSet object, down regulated gene set
#' @param subSamples A character or vector of sample labels/indices that will be
#'   used to subset the rankData matrix.All samples will be scored by default.
#' @param centerScore A Boolean, specifying whether scores should be centred,
#'   default as FALSE
#' @param dispersionFun A character, dispersion function with default as 'mad'
#' @return A data.frame consists of scores and dispersions for all samples
#' @keywords internal
#'
singscoringOneGS <- function (rankData, upSet, subSamples = NULL,
                              centerScore = FALSE, dispersionFun = mad) {

  #subset the data for samples whose calculation is to be performed
  if (!is.null(subSamples)) {
    rankData <- rankData[, subSamples, drop = FALSE]
  }

  #median center ranks
  rankData = apply(rankData, 2, function(x) {
    x = abs(x - ceiling(median(x)))
    return(x)
  })

  #values needed for calculating the boundaries
  upSigSize <- floor(length(geneIds(upSet))/2)
  #number of unique ranks in the geneset
  nTotalGenes <- ceiling(nrow(rankData)/2)
  #number of unique ranks in the matrix

  #check if there are some missing genes in the geneset
  missingGenes <- setdiff(geneIds(upSet), rownames(rankData))
  if (length(missingGenes) > 0) {
    warningMsg <-
      paste(length(missingGenes), "genes missing:", sep = ' ')
    warningMsg <- paste(warningMsg,
                        paste(missingGenes, collapse = ', '), sep = ' ')
    warning(warningMsg)
  }

  #remove missing genes from signature for further analysis
  geneIds(upSet) <- setdiff(geneIds(upSet), missingGenes)
  upRanks <- rankData[geneIds(upSet), , drop = FALSE]
  upScore <- colMeans(upRanks)
  lowBound <- (upSigSize + 1) / 2
  upBound <- (2 * nTotalGenes - upSigSize + 1) / 2
  normUpScore <- (upScore - lowBound) / (upBound - lowBound)

  #calculate dispersion
  upDispersion <- apply(upRanks, 2, dispersionFun)

  #if centering
  if (centerScore) {
    normUpScore <- normUpScore - 0.5
  }

  scoredf <-
    data.frame("TotalScore" = normUpScore, "TotalDispersion" = upDispersion)
  rownames(scoredf) <- colnames(rankData)
  return(scoredf)
}
