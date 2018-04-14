#' @include singscore.R
NULL

################################################################################
####=========================== rankExpr() function ===========================
################################################################################

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
#' \dontrun{ranked <- rankExpr(toy_expr)}
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
#' @param centerScore A Boolean, specifying whether scores should be centred
#' @param dispersionFun A character, dispersion function with default as 'mad'
#' @return A data.frame consists of scores and dispersions for all samples
#' @keywords internal
#' 
singscoringOneGS <- function (rankData, upSet, subSamples = NULL,
                              centerScore = TRUE, dispersionFun = mad) {
  
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
