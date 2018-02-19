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
#' [rankExpr()]
#' @examples
#' \dontrun{ranked <- rankExpr(toy_expr)}
rankExpr <- function(exprsM, tiesMethod = "min") {
    rankedData <- apply(exprsM, 2, rank, ties.method = tiesMethod)
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
#' @param dispersionFun A character, dispersion function with default as 'mad'
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
  upSigSize <- length(GSEABase::geneIds(upSet))
  nTotalGenes <- nrow(rankData)
  
  #check if there are some missing genes in the geneset
  missingGenes <- setdiff(GSEABase::geneIds(upSet), rownames(rankData))
  if (length(missingGenes) > 0) {
    warningMsg <-
      paste(length(missingGenes), "genes missing:", sep = ' ')
    warningMsg <- paste(warningMsg,
                        paste(missingGenes, collapse = ', '), sep = ' ')
    warning(warningMsg)
  }
  
  #remove missing genes from signature for further analysis
  GSEABase::geneIds(upSet) <- setdiff(GSEABase::geneIds(upSet), missingGenes)
  upRanks <- rankData[GSEABase::geneIds(upSet), , drop = FALSE]
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
#### =============================== plotDispersion() ==========================
################################################################################

#' Plot the score v.s. despersion for all samples
#' @description This function takes the output from the simpleScore() function
#'   and plots scatter plots for the score vs the dispersion for the total
#'   score, the up score and the down score of samples. The plots 
#' @param scoredf data.frame, results of the simpleScore() function
#' @param annot any annotation provided by the user that needs to be plot
#' annot must be ordered in the same was as the scores
#' @param alpha numeric, set the transparency of points
#' @param size numeric, Set the size of each point
#' @param textSize numeric, relative text sizes for title, labels and axis
#' values
#' @param isInteractive Boolean, Determine whether the plot is interactive
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#' plotDispersion(scoredf)
#' plotDispersion(scoredf, isInteractive = TRUE)
#' @return A ggplot object
#' @export
plotDispersion <- function(scoredf, annot = NULL, alpha = 1, size = 1,
                           textSize = 1.5, isInteractive=FALSE){
  if (is.null(annot)) {
    annot = rep('', nrow(scoredf))
  }
  #name annots
  annot = as.factor(annot)
  names(annot) = rownames(scoredf)

  #transform data for plot
  plotdf = scoredf
  plotdf['SampleID'] = rownames(plotdf)
  plotdf['Class'] = annot

  if (ncol(scoredf) > 2) {
    total = cbind(plotdf[, c(1:2, 7:8)], 'Total Score')
    up = cbind(plotdf[, c(3:4, 7:8)], 'Up Score')
    down = cbind(plotdf[, c(5:6, 7:8)], 'Down Score')
    colnames(total) = colnames(up) = colnames(down) = 1:ncol(total)
    plotdf = rbind(total, up, down)
  }
  colnames(plotdf)[1:4] = c('Score', 'Dispersion', 'SampleID', 'Annotation')
  Annotation <- NULL
  Score <- NULL
  Dispersion <- NULL
  SampleID <- NULL
  
  #Scatter plot
  p = with(plotdf,{ ggplot(plotdf, aes(Score, Dispersion, text = SampleID))})
  #colour by classification
  if (is.null(annot)) {
    p = p + geom_point(alpha = alpha, size = size)
  } else{
    p = p + geom_point(aes(colour = Annotation), alpha = alpha, size = size)
  }
  #up/down?
  if (ncol(scoredf) > 2) {
    colnames(plotdf)[5] = 'Type'
    p = p + facet_wrap( ~ plotdf$Type, scales = 'free')
  }
  n_color = length(unique(plotdf$Annotation))
  #plot properties
  if(n_color == 1) {
    p = p + 
      scale_color_manual(values = RColorBrewer::brewer.pal(8,'Set1')[4])
  } else if(n_color > 9){
    p = p + 
      scale_color_manual(values = grDevices::terrain.colors(n_color))
  }else {
    p = p + 
      scale_color_brewer(palette = 'RdYlBu',direction = 1)
  }
  
    p = p +
    ggtitle('Score vs Dispersion') +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = rel(textSize)),
      axis.text.x = element_text(angle = 0, size = rel(textSize)),
      axis.text.y = element_text(angle = 0, size = rel(textSize)),
      strip.background = element_rect(colour = "#f0f0f0",
                                      fill = "#f0f0f0"),
      strip.text = element_text(size = rel(textSize)),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.margin = margin(unit(0, "cm")),
      legend.title = element_text(face = "italic"),
      plot.title = element_text(
        face = "bold",
        size = rel(textSize),
        hjust = 0.5
      )
    )

  #if no annotation provided
  if(all(annot %in% '')){
    p = p + theme(legend.position="none")
  }
  if (isInteractive) {

    # replace params as ggplot objects are mutable
    oldparams = p$layers[[1]]$aes_params
    p$layers[[1]]$aes_params = NULL
    ply = plotly::ggplotly(p)
    p$layers[[1]]$aes_params = oldparams

    return(ply)
  }
  else{
    return(p)
  }
}

################################################################################
#### =============================== plotScoreLandscape() ======================
################################################################################

#' Plot landscape of two gene signatures scores
#' @description This function takes two data frames which are outputs from the
#'   simpleScore() function and plots the relationship between the two gene set
#'   scores for samples in the gene expression matrix.
#'
#' @param scoredf1 data.frame, result of the simpleScore() function which scores
#'   the gene expression matrix against a gene set of interest
#' @param scoredf2 data.frame, result of the simpleScore() function which scores
#'   the gene expression matrix against another gene set of interest
#' @param scorenames character vector of length 2, names for the two scored gene
#'   set/signatures stored in scoredf1 and scoredf2
#' @param isInteractive boolean, whether the plot is interactive default as
#'   FALSE
#' @param textSize numeric, set the text size for the plot, default as 1.5
#' @param hexMin integer, the threshold which decides whether hex bin plot or
#'   scatter plot is displayed, default as 100
#' @return A ggplot object, a scatter plot, demostrating the relationship
#'   between scores from two signatures on the same set of samples.
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#' scoredf2 <- simpleScore(ranked, upSet = toy_gs_up)
#' plotScoreLandscape(scoredf, scoredf2)
#' @export
plotScoreLandscape <- function(scoredf1, scoredf2, scorenames = c(),
                               textSize = 1.5, isInteractive = FALSE, 
                               hexMin = 100){
  if (length(scorenames) == 0){
    scorenames = c('Signature 1', 'Signature 2')
  }
  plotdf = data.frame(scoredf1$TotalScore, scoredf2$TotalScore)
  colnames(plotdf) = scorenames

  # generate labels
  pxlab = paste0('`', scorenames[1], '`')
  pylab = paste0('`', scorenames[2], '`')
  if(nrow(scoredf1) < hexMin){
    p = ggplot(plotdf, aes_string(pxlab, pylab)) +
      geom_point(colour = 'blue') +
      scale_fill_distiller(palette = 'RdPu', direction = 1)
    p = p +
      ggtitle('Signature landscape')
  }else{
    p = ggplot(plotdf, aes_string(pxlab, pylab)) +
      geom_hex(colour = 'white') +
      scale_fill_distiller(palette = 'RdPu', direction = 1)
    p = p +
      ggtitle('Signature landscape')
  }
  
    p = p+ 
      theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = rel(textSize)),
      axis.text.x = element_text(angle = 0, size = rel(textSize)),
      axis.text.y = element_text(angle = 0, size = rel(textSize)),
      strip.background = element_rect(colour = "#f0f0f0",
                                      fill = "#f0f0f0"),
      strip.text = element_text(size = rel(textSize)),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      legend.margin = margin(unit(0, "cm")),
      legend.title = element_text(face = "italic",
                                  size = rel(max(1, textSize * 0.55))),
      legend.text = element_text(size = rel(max(1, textSize * 0.5))),
      plot.title = element_text(
        face = "bold",
        size = rel(textSize),
        hjust = 0.5
      )
    )

  if (isInteractive) {
    #replace params as ggplot objects are mutable
    oldparams = p$layers[[1]]$aes_params
    p$layers[[1]]$aes_params = NULL
    ply = plotly::ggplotly(p)
    p$layers[[1]]$aes_params = oldparams

    return(ply)
  } else{
    return(p)
  }
}

################################################################################
####============================ projectScoreLandscape() =======================
################################################################################

#' Project data on the landscape plot obtained from \code{plotScoreLandscape()}
#' 
#'@description This function takes the output (ggplot object) of the
#'  \code{plotScoreLandscape()}  and a dataset. It projects the data onto the 
#'  ggplot object and returns a ggplot object with projected data points.
#' @param plotObj a dataframe, resulted from [plotScoreLandscape()]
#' @param subSamples vector of character or indices for subsetting the scoredfs,
#'  default as NULL and all samples in scoredfs will be plotted
#' @param sampleLabels vector of character, sample names to display, ordered in
#'   the same way as samples are ordered in the 'scoredf1' data matrix, default 
#'   as NULL which means points are labelled by sample names.
#' @param annot vector of characters, annotations to colour the data and should 
#'  have the sample of length with number of samples in scoredfs
#' @inheritParams plotScoreLandscape
#' 
#' @return New data points on the already plotted ggplot object from
#'   plotScoreLanscape()
#' @seealso 
#' [plotScoreLandscape()]
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf1 <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#' scoredf2 <- simpleScore(ranked, upSet = toy_gs_up)
#' psl <- plotScoreLandscape(scoredf1, scoredf2)
#' projectScoreLandscape(psl,scoredf1, scoredf2)
#' @export
projectScoreLandscape <- function(plotObj = NULL,
                                  scoredf1,
                                  scoredf2,
                                  subSamples = NULL,
                                  sampleLabels = NULL,
                                  annot = NULL,
                                  isInteractive = FALSE){
  #create data frame with the new data
  #subsetting the two data frames, scoredfs
  if(! is.null(subSamples)){
    scoredf1 <- scoredf1[subSamples,]
    scoredf2 <- scoredf2[subSamples,]
    if(anyNA(scoredf1)){
      message('some selected samples not exist in provided scoredf1')
      scoredf1 <- na.omit(scoredf1)
    }
    if(anyNA(scoredf2)){
      message('some selected samples not exist in provided scoredf2')
      scoredf2 <- na.omit(scoredf2)
      
    }
  }
  #if no sample labels are provided, the plot labels the points using sample 
  #names
  if (is.null(sampleLabels)) {
    sampleLabels <- rownames(scoredf1)
  }else{
    if(length(sampleLabels) != nrow(scoredf1))
    stop("sampleLabels must have same number of labels with samples")
  }
  if (is.null(annot)) {
    annot = ''
  }
  newdata = data.frame(scoredf1$TotalScore, scoredf2$TotalScore, sampleLabels)
  
  if (! is.ggplot(plotObj)) {
    stop('Please provide a ggplot object (',
         class(plotObj)[1], ' object given)')
  }

  plabs = c(plotObj$labels$x, plotObj$labels$y)
  Annotation <- NULL
  SampleLabel <- NULL
  colnames(newdata) = c(plabs, 'SampleLabel')
  newdata[, 'Annotation'] = as.factor(annot) #need to make it work for factor
  
  #need to deal with legends in both interactive and non-interactive
  if (!isInteractive) {
    #add layer with new data
    pproj = plotObj + geom_point(
      data = newdata,
      aes(text = SampleLabel, colour = Annotation),
      shape = 21,
      fill = 'white',
      size = 2,
      stroke = 2
    ) + ggsci::scale_color_npg()
    
    
    #label samples
    pproj = pproj +
      ggrepel::geom_label_repel(
        data = newdata,
        aes(label = SampleLabel, colour = Annotation),
        show.legend = FALSE
      ) 
  } else if(isInteractive) {
    #replace params as ggplot objects are mutable
    oldparams = plotObj$layers[[1]]$aes_params
    plotObj$layers[[1]]$aes_params = NULL
    ply = plotly::ggplotly(plotObj)
    plotObj$layers[[1]]$aes_params = oldparams

    #add layer with new data
    npgpal =ggsci::pal_npg('nrc')(length(levels(newdata$Annotation)))
    ply = ply %>%
      plotly::add_trace(data = newdata,
                color = ~Annotation,
                colors = npgpal,
                type = 'scatter',
                mode = 'markers',
                marker = list(
                  size = 10,
                  line = list(color = 'white', width = 2)
                ),
                text = paste('Cell line:', newdata$SampleLabel)) %>%
     plotly::layout(showlegend = TRUE,
             legend = list(
               orientation = 'h',
               xanchor = 'center',
               x = 0.5,
               yanchor = 'top',
               y = -0.2
             ))

    return(ply)
  }

  return(pproj)
}

################################################################################
#### =============================== plotRankDensity_intl() ====================
################################################################################

#' Plot the densities of ranks for one sample
#' 
#' @description This function takes a single column data frame, which is a
#' subset of the ranked data obtained from [rankExpr()]function, and gene sets,
#' and it returns plots visualising the density and the rugs of the ran ks.
#'
#' @param rankData one column of the ranked gene expression matrix obtained from
#' the [rankGenes()] function, use drop = FALSE when subsetting the ranked gene 
#' expression matrix, see examples.
#' @param isInteractive Boolean, determin whether the returned plot is
#'   interactive
#' @param textSize numberic, set the size of text on the plot
#' @param upSet GeneSet object, up regulated gene set
#' @param downSet GeneSet object, down regulated gene set
#'
#' @return A ggplot object (optionally interactive) demonstrating the rank
#'   density along with rug plot

#' @seealso 
#' \code{"\linkS4class{GeneSet}"}
plotRankDensity_intl <- function (rankData,
                             upSet,
                             downSet = NULL,
                             isInteractive = FALSE,
                             textSize = 1.5) {
  #values needed for calculating the boundaries
  upSigSize = length(GSEABase::geneIds(upSet))
  nTotalGenes = nrow(rankData)
  #browser()
  #check if there are some missing genes in the geneset
  missingGenes = setdiff(GSEABase::geneIds(upSet), rownames(rankData))
  if (length(missingGenes) > 0) {
    warningMsg = paste(length(missingGenes), 'genes missing:', sep = ' ')
    warningMsg = paste(warningMsg, paste(missingGenes, collapse = ', '),
                       sep = ' ')
    warning(warningMsg)
  }

  #remove missing genes from signature for further analysis
  GSEABase::geneIds(upSet) = setdiff(GSEABase::geneIds(upSet), missingGenes)
  upRanks = rankData[GSEABase::geneIds(upSet), , drop = FALSE] / nrow(rankData)
  upRank = data.frame(upRanks, type = "Up Gene-set")
  allRanks = upRank

  if (!is.null(downSet)) {
    #check if there are some missing genes in the geneset
    missingGenes = setdiff(GSEABase::geneIds(downSet), rownames(rankData))
    if (length(missingGenes) > 0) {
      warningMsg = paste(length(missingGenes), 'genes missing:',
                         sep = ' ')
      warningMsg = paste(warningMsg,
                         paste(missingGenes, collapse = ', '), sep = ' ')
      warning(warningMsg)
    }

    #remove missing genes from signature for further analysis
    GSEABase::geneIds(downSet) = setdiff(GSEABase::geneIds(downSet), missingGenes)
    downRanks = rankData[GSEABase::geneIds(downSet), , drop = FALSE] / nrow(rankData)
    downRank = data.frame(downRanks, type =  "Down Gene-set")
    allRanks = rbind(upRank, downRank)
  }
  Ranks <- NULL
  upDown <- NULL
  EntrezID <- NULL
  colnames(allRanks) <- c("Ranks", "upDown")
  allRanks$EntrezID <- row.names(allRanks)


  #bar plot preparations
  ymap = c(0, 0)
  yendmap = ymap + 0.3
  colmap = c(RColorBrewer::brewer.pal(8, "Set1")[c(1, 2)])
  typemap = c('Up-regulated gene', 'Down-regulated gene')
  names(colmap) = names(ymap)  = c('Up Gene-set', 'Down Gene-set')
  names(yendmap) = names(typemap) = c('Up Gene-set', 'Down Gene-set')
  ..density.. <- NULL
  
  #plot density and calculate max density and barcode line heights and
  #positions
  p =with(allRanks,{
    ggplot(allRanks, aes(x = Ranks, col = upDown)) +
      stat_density(aes(y = ..density..), geom = 'line', position = 'identity')
  })

  dens = ggplot_build(p)$data[[1]]$density
  ymap[1] = round(max(dens), digits = 1) + 0.1
  ymap[2] = round(min(dens), digits = 1) - 0.1
  bcheight = (max(dens) - min(dens))
  bcheight = bcheight/ifelse(is.null(downSet), 4, 3)
  yendmap = ymap + c(1, -1) * bcheight
  
  #plot barcode plot
  #text aes useful for the plotly plot, so supress the warnings
  #
  p = suppressWarnings( p + geom_segment(aes(
    y = ymap[upDown],
    xend = Ranks,
    yend = yendmap[upDown],
    text = paste0(typemap[upDown], '\nGene symbol: ', EntrezID)
  ),alpha = 0.8) +
    scale_colour_manual(values = colmap,
                        guide = guide_legend(title = "Type")))

  #publication quality plot
  p = p + ggtitle('Rank density') +
    xlab('Normalised Ranks') +
    ylab('Density') +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = rel(textSize)),
      axis.text.x = element_text(angle = 0, size = rel(textSize)),
      axis.text.y = element_text(angle = 0, size = rel(textSize)),
      strip.background = element_rect(colour = "#f0f0f0",
                                      fill = "#f0f0f0"),
      strip.text = element_text(size = rel(textSize)),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.margin = margin(unit(0, "cm")),
      legend.title = element_text(size = rel(textSize * 0.8),
                                  face="italic"),
      legend.text = element_text(size = rel(textSize * 0.8)),
      plot.title = element_text(
        face = "bold",
        size = rel(textSize),
        hjust = 0.5
      )
    )

  #if single geneset, remove legend
  if (is.null(downSet)) {
    p = p + theme(legend.position = 'none')
  }

  if (isInteractive) {
    #Horizontal legend not supported by plotly yet so re-orient after
    #creating plotly object
    ply = suppressWarnings(plotly::ggplotly(p, tooltip = c('text', 'x')))
    ply = ply %>% plotly::layout(
      legend = list(
        orientation = 'h',
        xanchor = 'center',
        x = 0.5,
        yanchor = 'top',
        y = -0.25
      ), yaxis =list(
        fixedrange = TRUE
      ))
    return(ply)
  } else{
    return(p)
  }
}

#' @title Permutation test for the derived scores of each sample
#'
#' @description This function randomly generates a number of gene sets which
#'   have the same number of genes as the scored gene set. It scores each random
#'   gene set and return a matrix with scores for each permutation across all
#'   samples. The empirical scores are used to calculate the empirical p value
#'   and plot the null distribution. The implementation uses
#'   [BiocParallel::bplapply()] for easy access to parallel backends.
#' @param n_up integer,  size of up set
#' @param n_down integer, size of down set
#' @param rankData matrix, outcome of function [rankExpr()]
#' @param B integer, the number of permutation repeats default as 1000
#' @param seed integer, set the seed for randomisation
#'
#' @return A matrix of empircal scores for each sample
#' @seealso 
#' [Post about BiocParallel](http://lcolladotor.github.io/2016/03/07/BiocParallel/#.WgXMF61L28U)
#' `browseVignettes("BiocParallel")`
#' @author Ruqian Lyu
#' @export
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#' n_up = length(GSEABase::geneIds(toy_gs_up))
#' n_down = length(GSEABase::geneIds(toy_gs_dn))
#' # find out what backends can be registered on your machine
#' BiocParallel::registered()
#' # the first one is the default backend, and it can be changed explicitly.
#' permuteResult = generateNull(n_up = n_up, n_down = n_down, ranked, B =10,
#' seed = 1) # call the permutation function to generate the empirical scores 
#' #for B times.

generateNull <- function(n_up, n_down, rankData, B = 1000, seed = 1){
  # A copy of singscoring method for usage in the Bpapply function  
  singscoring <- function (rankData, upSet, downSet = NULL, subSamples = NULL,
                             centerScore = TRUE, dispersionFun = mad) {
      #subset the data for samples whose calculation is to be performed
      if (!is.null(subSamples)) {
        rankData <- rankData[, subSamples, drop = FALSE]
      }
      #values needed for calculating the boundaries
      upSigSize <- length(GSEABase::geneIds(upSet))
      nTotalGenes <- nrow(rankData)
      
      #check if there are some missing genes in the geneset
      missingGenes <- setdiff(GSEABase::geneIds(upSet), rownames(rankData))
      if (length(missingGenes) > 0) {
        warningMsg <-
          paste(length(missingGenes), "genes missing:", sep = ' ')
        warningMsg <- paste(warningMsg,
                            paste(missingGenes, collapse = ', '), sep = ' ')
        warning(warningMsg)
      }
      
      #remove missing genes from signature for further analysis
      GSEABase::geneIds(upSet) <- setdiff(GSEABase::geneIds(upSet), missingGenes)
      upRanks <- rankData[GSEABase::geneIds(upSet), , drop = FALSE]
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
    
    set.seed(seed)
    all_genes <- rownames(rankData)
    totalNo <- n_up + n_down
    temSets <- BiocParallel::bplapply(1:B, function(i) {
      tms <-  sample(all_genes, size = totalNo, replace = FALSE)
      tms
    })
    r <- BiocParallel::bplapply(1:B, function(i) {
      tms <- temSets[[i]]
      if (n_down > 0) {
        upSet <- GSEABase::GeneSet(as.character(tms[1:n_up]))
        downSet <-  GSEABase::GeneSet(as.character(tms[-(1:n_up)]))
        ss <-  singscoring(rankData, upSet = upSet, downSet = downSet)
      } else {
        #else all the random generated genes are in upSet
        ss <- singscoring(rankData, upSet = GeneSet(as.character(tms)))
      }
      ss[, 1]
    })
    r <- plyr::ldply(r)
    colnames(r) <- colnames(rankData)
    return(r)
}

#' Calculate the empirical p values
#'
#' @description This function takes the permutation results, which is the
#'   empirical scores, and the calculated sample scores using [simpleScore()] as
#'   input. It calculates the empirical p values of the simple sample scoring
#'   test using formula p = (r+1)/(m+1) where r is the number of empirical
#'   scores that are larger than the obtained score and m is the total number of
#'   permutation run which is the B parameter in [generateNull()]
#'
#' @param permuResult A matrix, result from [generateNull()] function
#' @param scoredf A dataframe, result from [simpleScore()] function
#'
#' @return pvals for each sample, the calculated empirical p values for all
#'   empirical sample scores null distribution
#' @author Ruqian Lyu
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#' n_up = length(GSEABase::geneIds(toy_gs_up))
#' n_down = length(GSEABase::geneIds(toy_gs_dn))
#' # find out what backends can be registered on your machine
#' BiocParallel::registered()
#' # the first one is the default backend, and it can be changed explicitly.
#' permuteResult = generateNull(n_up = n_up, n_down = n_down, ranked, B =10,
#' seed = 1) # call the permutation function to generate the empirical scores 
#' #for B times.
#' pvals <- getPvals(permuteResult,scoredf)
#' @export
getPvals <- function(permuResult,scoredf){
  resultSc <- t(scoredf[, 1, drop = FALSE])
# combine the permutation with the result score for the computation of P values
# p = (r+1)/(m+1)
  empirScore_re <- rbind(permuResult, as.character(resultSc))
  
  # x[length(x)] is the calculated score
  
  pvals <- apply(empirScore_re,2,function(x){
    length(x[x >= x[length(x)]]) / length(x)
  })
  return(pvals)
}

#' Plot the empirical null distribution using the permutation result
#' 
#' @description This function takes the results from function [generateNull()] 
#' and plots the density curves of empirical scores for the provided samples via
#' \code{sampleNames} parameter. It can plot null distribution for a single 
#' sample or multiple samples.
#' 
#' @param permuResult A matrix, outcome from function [generateNull()]
#' @param scoredf A dataframe, outcome from function [simpleScore()]
#' @param pvals A vector, outcome of function [getPvals()]
#' @param sampleNames A vector of character, sample names or multiple sample 
#' labels
#' @param alpha numeric,ggplot theme element
#' @param size numeric,ggplot theme element
#' @param textSize numeric,ggplot theme element
#' @param labelSize numeric, size of label texts.
#' @param cutoff double, the cutoff value for determining significance
#' @return a ggplot object
#' @author Ruqian Lyu
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#' n_up = length(GSEABase::geneIds(toy_gs_up))
#' n_down = length(GSEABase::geneIds(toy_gs_dn))
#' # find out what backends can be registered on your machine
#' BiocParallel::registered()
#' # the first one is the default backend, and it can be changed explicitly.
#' permuteResult = generateNull(n_up = n_up, n_down = n_down, ranked, B =10,
#' seed = 1) # call the permutation function to generate the empirical scores 
#' #for B times.
#' pvals <- getPvals(permuteResult,scoredf)
#' plotNull(permuteResult,scoredf,pvals,sampleNames = names(pvals))
#' plotNull(permuteResult,scoredf,pvals,sampleNames = names(pvals)[1])
#' @export
plotNull <- function(permuResult, scoredf, pvals, sampleNames = NULL,
                      cutoff = 0.01,alpha = 1, size = 1, 
                     textSize = 2,labelSize = 5){
  quantile_title <- as.character((1 - cutoff)*100)
  if(!is.null(sampleNames)){
    pvals <- pvals[sampleNames, drop = FALSE]
  }
  pval_r <- as.character(format(pvals[sampleNames],scientific = TRUE,
                                digits = 3))
  pvalTitle <- paste0(' p-value = ',pval_r)
  names(pvalTitle) <- names(pvals[sampleNames])
  cutoff_score <- c()
  for(i in 1:length(sampleNames)){
    cutoff_score[i] <- quantile(permuResult[,sampleNames[i]],(1-cutoff))
  }
  names(cutoff_score) = sampleNames
  cutoff_annot = data.frame(sampleNames = sampleNames, 
                            cutoff_score = cutoff_score)
  #pDt <-  as.data.frame(pvals)
  if( !is.null(sampleNames) ){
    if(length(sampleNames)>1){
      dt <- as.data.frame(permuResult[,sampleNames])
      longDt <- reshape::melt(dt,variable_name = "sampleNames")
      resultScs <-  scoredf[,1,drop = FALSE]
      resultScs$sampleNames <-  rownames(resultScs)
      #pDt$sampleNames <- names(pvals)
      sampleLSc <-  merge(longDt, resultScs, by.x = "sampleNames", 
                          by.y = "sampleNames")
      #plotDt  <-  merge(sampleLSc,pDt, by.x = 'sampleNames', 
      #by.y = 'sampleNames')
      sampleLSc <- merge(sampleLSc,cutoff_annot)
      sampleLSc <- merge(sampleLSc, pvalTitle)
      value <- NULL
      TotalScore <- NULL
      #browser()
      plotObj <-  ggplot(data = sampleLSc)+
        geom_density(mapping = aes( x = value), size =1)+
        coord_cartesian(xlim = c(min(permuResult[,sampleNames]) - 0.03,
                                 max(permuResult[,sampleNames])+0.05))+
        facet_grid(sampleNames~.)+
        geom_segment(mapping =  aes(x  = cutoff_score, y = 11, 
                                    xend = cutoff_score, yend = 0), 
                     linetype="dashed", colour = 'blue',size = 1)+
        geom_segment(mapping = aes(x  = TotalScore, y = 11, xend = TotalScore, 
                                   yend = 0),colour = 'red',size = 2)+
        geom_text(mapping = aes(x  = TotalScore-0.01, y = 12, 
                                label = pvalTitle[sampleNames]), 
                  colour = 'red',size = labelSize)+
        geom_text(mapping = aes(x  = cutoff_score,y = 12, 
                                label = paste0(quantile_title,
                                               '%-ile threshold')), 
                  colour = 'blue',size = labelSize)+
        xlab("Scores")+
        ggtitle("Null distribution")
    }else{
      plotDt <- data.frame(sampleNames = sampleNames, 
                           value = permuResult[,sampleNames],
                           TotalScore = scoredf[sampleNames,]$TotalScore)
      plotObj <-  ggplot(data = plotDt)+
        geom_density(mapping = aes( x = value),size = 1)+
        coord_cartesian(xlim = c(min(permuResult[,sampleNames]) - 0.03,
                                 max(permuResult[,sampleNames])+0.05))+
        geom_segment(mapping =  aes(x = cutoff_score, y = 11, 
                                    xend = cutoff_score, yend =0), 
                     linetype="dashed", colour = 'blue',size = 1)+
        geom_segment(mapping = aes(x = TotalScore, y = 11, 
                                   xend = TotalScore, yend =0),
                     colour = 'red',size = 2)+
        geom_text(mapping = aes(x = TotalScore, y = 12, 
                                label = pvalTitle[sampleNames]), 
                  colour = 'red',size = labelSize)+
        geom_text(mapping = aes(x = cutoff_score, y = 12, 
                                label = paste0(quantile_title,
                                               '%-ile threshold')), 
                  colour = 'blue',size = labelSize)+
        xlab("Scores")+
        ggtitle( paste0(sampleNames," null distribution"))
    }
  }else{
    warning("Please provide which sample's null distribution to
            plot by specifying the sampleNames")
  }
  plotObj+
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = rel(textSize)),
      axis.text.x = element_text(angle = 0, size = rel(textSize)),
      axis.text.y = element_text(angle = 0, size = rel(textSize)),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(size = rel(textSize)),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.margin = margin(unit(0, "cm")),
      legend.title = element_text(face = "italic"),
      plot.title = element_text(
        face = "bold",
        size = rel(textSize),
        hjust = 0.5))
}
