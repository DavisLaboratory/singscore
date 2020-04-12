#' @include singscore.R
NULL

rankExpr <- function(exprsM, tiesMethod = "min") {
  rname= rownames(exprsM)
  cname = colnames(exprsM)
  rankedData = matrixStats::colRanks(as.matrix(exprsM),
                                      ties.method = tiesMethod,
                                      preserveShape = TRUE)
  rownames(rankedData) = rname
  colnames(rankedData) = cname
  
  #indicator of the type of ranks
  attr(rankedData, 'stable') = FALSE
  return (rankedData)
}

rankExprStable <- function(exprsM, tiesMethod = "min", stgenes) {
  stgenes = intersect(stgenes, rownames(exprsM))
  stopifnot(length(stgenes) > 0)
  
  rname = rownames(exprsM)
  cname = colnames(exprsM)
  
  #compute ranks
  rankedData = apply(exprsM, 2, function(x) {
    rowSums(outer(x, x[stgenes], '>')) + 1
  })
  
  #normlise ranks
  rankedData = rankedData / (length(stgenes) + 1)
  
  rownames(rankedData) = rname
  colnames(rankedData) = cname
  
  #indicator of the type of ranks
  attr(rankedData, 'stable') = TRUE
  return(rankedData)
}

#define the main helper functions
#helper: check if all genes in the gene set are in the data
checkGenes <- function (geneset, background) {
  #check if there are some missing genes in the geneset
  missingGenes = setdiff(geneIds(geneset), background)
  if (length(missingGenes) > 0) {
    warningMsg = paste(length(missingGenes), "genes missing:", sep = ' ')
    warningMsg = paste(warningMsg, paste(missingGenes, collapse = ', '), sep = ' ')
    warning(warningMsg)
  }
  geneIds(geneset) = setdiff(geneIds(geneset), missingGenes)

  return(geneset)
}

checkGenesMulti <- function(geneset_colc, background) {
  geneset_colc = endoapply(geneset_colc, function(geneset) {
    geneIds(geneset) = intersect(geneIds(geneset), background)
    return(geneset)
  })

  #remove genesets with no genes
  geneset_colc = geneset_colc[sapply(geneset_colc, function(x) length(geneIds(x))) > 0]

  return(geneset_colc)
}

#helper: compute the theoretical boundaries of scores when direction is KNOWN.
# Should be used after filtering out missing genes
calcBounds <- function(gsSize, bgSize, stableSc = FALSE) {
  if (stableSc) {
    lowBound = 0
    upBound = 1
  } else {
    lowBound = (gsSize + 1) / 2
    upBound = (2 * bgSize - gsSize + 1) / 2
  }

  return(list('lowBound' = lowBound, 'upBound' = upBound))
}

calcBoundsMulti <- function(gsSizes, bgSize, stableSc = FALSE) {
  bounds = lapply(gsSizes, calcBounds, bgSize, stableSc)

  return(bounds)
}

#helper: compute the theoretical boundaries of scores when direction is UNKNOWN.
# Should be used after filtering out missing genes
calcBoundsUnknownDir <- function(gsSize, bgSize, stableSc = FALSE) {
  gsSize = floor(gsSize / 2) #number of unique ranks in the geneset
  bgSize = ceiling(bgSize / 2) #number of unique ranks in the matrix

  return(calcBounds(gsSize, bgSize, stableSc))
}

calcBoundsUnknownDirMulti <- function(gsSizes, bgSize, stableSc = FALSE) {
  gsSizes = floor(gsSizes / 2) #number of unique ranks in the geneset
  bgSize = ceiling(bgSize / 2) #number of unique ranks in the matrix

  return(calcBoundsMulti(gsSizes, bgSize, stableSc))
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
    'TotalDispersion' = (upscores$TotalDispersion + dnscores$TotalDispersion) / 2,
    'UpScore' = upscores$TotalScore,
    'UpDispersion' = upscores$TotalDispersion,
    'DownScore' = dnscores$TotalScore,
    'DownDispersion' = dnscores$TotalDispersion
  )

  return(totalscores)
}

#singscore function for a single geneset/signature
singleSingscore <-
  function (rankData,
            upSet,
            downSet = NULL,
            subSamples = NULL,
            centerScore = TRUE,
            dispersionFun = mad,
            knownDirection = TRUE) {
    #0. determine type of singscore being used
    stableSc = attr(rankData, 'stable')

    #1. subset the data for samples whose calculation is to be performed
    if (!is.null(subSamples)) {
      rankData <- rankData[, subSamples, drop = FALSE]
    }

    #2. center scores?
    center_offset = ifelse(centerScore, -0.5, 0)

    #3.1 filter genes - return empty data.frame if no scores
    upSet = checkGenes(upSet, rownames(rankData))
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
      upset_bounds = calcBoundsUnknownDir(length(geneIds(upSet)), nrow(rankData), stableSc)
    } else {
      #4.1 calculate bounds
      upset_bounds = calcBounds(length(geneIds(upSet)), nrow(rankData), stableSc)
    }

    #for two gene sets
    if (!is.null(downSet)) {
      #3.2 filter genes - return empty data.frame if no scores
      downSet = checkGenes(downSet, rownames(rankData))
      if (length(geneIds(downSet)) == 0)
        return(emptyScoreDf(onegs = FALSE))

      #4.2 calculate bounds
      downset_bounds = calcBounds(length(geneIds(downSet)), nrow(rankData), stableSc)

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

#singscore function for a multiple geneset/signature
multiSingscore <-
  function (rankData,
            upSetColc,
            downSetColc = NULL,
            subSamples = NULL,
            centerScore = TRUE,
            dispersionFun = mad,
            knownDirection = TRUE) {
    #0. determine type of singscore being used
    stableSc = attr(rankData, 'stable')
    
    #1. subset the data for samples whose calculation is to be performed
    if (!is.null(subSamples)) {
      rankData <- rankData[, subSamples, drop = FALSE]
    }

    #2. center scores?
    center_offset = ifelse(centerScore, -0.5, 0)

    #3.1 filter genes - empty genesets will be removed
    upNames = names(upSetColc)
    upSetColc = checkGenesMulti(upSetColc, rownames(rankData))
    if (length(upSetColc) == 0)
      return(NULL)

    #known direction
    upSizes = sapply(upSetColc, function(x) length(geneIds(x)))
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
      upset_bounds = calcBoundsUnknownDirMulti(upSizes, nrow(rankData), stableSc)
    } else {
      #4.1 calculate bounds
      upset_bounds = calcBoundsMulti(upSizes, nrow(rankData), stableSc)
    }

    #for two gene sets
    if (!is.null(downSetColc)) {
      #3.2 filter genes - empty genesets will be removed
      downNames = names(downSetColc)
      downSetColc = checkGenesMulti(downSetColc, rownames(rankData))
      if (length(downSetColc) == 0)
        return(NULL)

      #drop those not matching with upSets
      stopifnot(all(upNames == downNames))
      commonNames = intersect(upNames, downNames)
      downSetColc = downSetColc[names(downSetColc) %in% commonNames]

      #4.2 calculate bounds
      downSizes = sapply(downSetColc, function(x) length(geneIds(x)))
      downset_bounds = calcBoundsMulti(downSizes, nrow(rankData), stableSc)

      #check that names match
      upselect = names(upSetColc) %in% commonNames
      upSetColc = upSetColc[upselect]
      upset_bounds = lapply(upset_bounds, function(x) x[upselect])

      #5.2 compute scores
      scores = mapply(
        calcScoresUpDn,
        list(rankData),
        upSetColc,
        downSetColc,
        list(dispersionFun),
        upset_bounds,
        downset_bounds,
        center_offset
      )
    } else{
      #5.1 & 5.3 compute scores
      scores = mapply(
        calcScores,
        list(rankData),
        upSetColc,
        list(dispersionFun),
        upset_bounds,
        center_offset
      )
    }

    #6 process the results
    dispersions = scores['TotalDispersion', ]
    dispersions = t(mapply(c, dispersions))
    scores = scores['TotalScore', ]
    scores = t(mapply(c, scores))
    colnames(scores) = colnames(dispersions) = colnames(rankData)
    rownames(scores) = rownames(dispersions) = upNames

    return(list('Scores' = scores, 'Dispersions' = dispersions))
  }
