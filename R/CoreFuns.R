#' @include simpleScoring.R
NULL


################################################################################
####=========================== rankGenes() function ===========================
################################################################################

#' Rank gene expression matrix
#' @description Given a gene expression matrix and a tiesMethod
#'   character, this fucntion calls the 'rank' function in the base package
#'   The default ties.Method is set to 'min'. There is a generic version of
#'   this function, details can be found in the see also section down at the
#'   bottom. It is suggested to use the generic function 'rankExpr'
#' @param exprsM A matrix, gene expression matrix
#' @param tiesMethod A character, default as 'min'
#'
#' @return A matrix that has samples in colunm and genes in rows. Values are the
#'   ranks of each gene in each sample.
#' @seealso
#' \code{ \link{rankExpr}}
#' @examples
#' ranked <- rankGenes(toy_expr)
#' @export
rankGenes<- function(exprsM, tiesMethod = "min"){
  rankedData<- apply(exprsM, 2, rank, ties.method= tiesMethod)
  return (rankedData)
}

################################################################################
####========================= simpleScore() function ===========================
################################################################################

#' Calculate scores for the ranked gene expression matrix against a gene set
#' @description This function takes a ranked gene expression matrix obtained
#'   from \code{rankExpr} (or \code{rankGenes}) function and a GeneSet object as
#'   input parameters. It returns a data.frame consists of scores and
#'   dispersions for each sample.
#'
#' @param rankData A matrix, ranked gene expression matrix data
#' @param subSamples A character or vector of sample labels/indices that will be
#'   used to subset the rankData matrix
#' @param upSet A GeneSet object, up regulated gene set
#' @param downSet A GeneSet object, down regulated gene set
#' @param centerScore A Boolean, specifying whether scores should be centred
#' @param dispersionFun A character, dispersion function with default as 'mad'
#' @return A data.frame consists of scores and dispersions for all samples
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_up, downSet = toy_dn)
#' @export
simpleScore <-
  function (rankData,
            subSamples = NULL,
            upSet,
            downSet = NULL,
            centerScore = TRUE,
            dispersionFun = mad) {
    #subset the data for samples whose calculation is to be performed
    if (!is.null(subSamples)) {
      rankData <- rankData[, subSamples, drop = F]
    }

    #values needed for calculating the boundaries
    upSigSize <- length(geneIds(upSet))
    nTotalGenes <- nrow(rankData)

    #check if there are some missing genes in the geneset
    missingGenes <- setdiff(geneIds(upSet), rownames(rankData))
    if (length(missingGenes) > 0) {
      warningMsg <- paste(length(missingGenes), "genes missing:", sep = ' ')
      warningMsg <-
        paste(warningMsg, paste(missingGenes, collapse = ', '), sep = ' ')
      warning(warningMsg)
    }

    #remove missing genes from signature for further analysis
    geneIds(upSet) <- setdiff(geneIds(upSet), missingGenes)
    upRanks <- rankData[geneIds(upSet), , drop = F]
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
    } else{
      scoreDown <-
        simpleScore(rankData, subSamples, downSet, NULL, F, dispersionFun)
      normDownScore <- 1 - scoreDown$TotalScore
      downDispersion <- scoreDown$TotalDispersion

      #if centering
      if (centerScore) {
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
#'   score, the up score and the down score
#' @param scoredf data.frame, results of the simpleScore() function
#' @param annot any annotation provided by the user that needs to be plot
# 		annot must be ordered in the same was as the scores
#' @param alpha numeric, set the transparency of points
#' @param size numeric, Set the size of each point
#' @param textSize numeric, relative text sizes for title, labels and axis
#'   values
#' @param isInteractive Boolean, determin whether the plot is interactive
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_up, downSet = toy_dn)
#' plotDispersion(scoredf)
#' plotDispersion(scoredf, isInteractive = TRUE)
#' @return A ggplot object
#' @export
plotDispersion <- function(scoredf, annot = NULL, alpha = 1, size = 1,
                           textSize = 1.5, isInteractive=F){
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

  #Scatter plot
  p = ggplot(plotdf, aes(Score, Dispersion, text = SampleID))
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

  #plot properties
  p = p + scale_colour_manual(values= RColorBrewer::brewer.pal(8,"Set1")[4])+
    # scale_color_npg() +
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
#' @description  This function takes two data frames which are the output from
#' the simpleScore() function and plots the relationship between the two scores.
#'
#' @param scoredf1 data.frame, results of the simpleScore() function
#'
#' @param scoredf2 data.frame, results of the simpleScore() function
#' @param scorenames character, names for the two gene signatures scored stored
#' in scoredf1 and scoredf2
#' @param isInteractive boolean, whether the plot is interactive
#' @inheritParams plotDispersion
#' @return A ggplot object, a scatter plot, demostrating the relationship
#' between scores from two signatures on the same set of samples.
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_up, downSet = toy_dn)
#' scoredf2 <- simpleScore(ranked, upSet = toy_up)
#' plotScoreLandscape(scoredf, scoredf2)
#' @export
plotScoreLandscape <- function(scoredf1, scoredf2, scorenames = c(),
                               textSize = 1.5, isInteractive = F){
  if (length(scorenames) == 0){
    scorenames = c('Signature 1', 'Signature 2')
  }
  plotdf = data.frame(scoredf1$TotalScore, scoredf2$TotalScore)
  colnames(plotdf) = scorenames

  # generate labels
  pxlab = paste0('`', scorenames[1], '`')
  pylab = paste0('`', scorenames[2], '`')

  p = ggplot(plotdf, aes_string(pxlab, pylab)) +
    geom_hex(colour = 'white') +
    scale_fill_distiller(palette = 'RdPu', direction = 1)
  p = p +
    ggtitle('Signature landscape') +
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

#' Project the landscape score obtained from plotScoreLanscape
#' @description This function takes the output of the plotScoreLandscape() and
#'   anew data, and plots the new data onto the ggplot object.
#' @param plotObj a dataframe, resulted from plotScoreLanscape()
#'
#' @param sampleIDs character, sample names to display, ordered in the same way
#'   as samples are ordered in the rank matrix
#' @param annot Annotations to colour the data.
#' @param labels integer, number of samples labeled on the plot in the new data
#' @inheritParams plotScoreLandscape
#' @return New data points on the already plotted ggplot object from
#'   plotScoreLanscape()
#'
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_up, downSet = toy_dn)
#' scoredf2 <- simpleScore(ranked, upSet = toy_up)
#' psl <- plotScoreLandscape(scoredf, scoredf2)
#' projectScoreLandscape(psl,scoredf, scoredf2)
#' @export
projectScoreLandscape <- function(plotObj = NULL,
                                  scoredf1,
                                  scoredf2,
                                  sampleIDs = NULL,
                                  annot = NULL,
                                  labels = 20,
                                  isInteractive = F){
  #create data frame with the new data
  if (is.null(sampleIDs)) {
    sampleIDs = 1:nrow(scoredf1)
  }
  if (is.null(annot)) {
    annot = ''
  }

  newdata = data.frame(scoredf1$TotalScore, scoredf2$TotalScore, sampleIDs)

  if (! is.ggplot(plotObj)) {
    stop('Please provide a ggplot object (',
         class(plotObj)[1], ' object given)')
  }

  plabs = c(plotObj$labels$x, plotObj$labels$y)
  colnames(newdata) = c(plabs, 'SampleID')
  newdata[, 'Annotation'] = as.factor(annot) #need to make it work for factor

  #need to deal with legends in both interactive and non-interactive
  if (!isInteractive) {
    #add layer with new data
    pproj = plotObj + geom_point(
      data = newdata,
      aes(text = SampleID, colour = Annotation),
      shape = 21,
      fill = 'white',
      size = 2,
      stroke = 2
    ) + scale_color_npg()

    #label samples
    if (nrow(newdata) <= labels) {
      pproj = pproj +
        geom_label_repel(
          data = newdata,
          aes(label = SampleID,
              colour = Annotation),
          show.legend = F
        )
    }
  } else if(isInteractive) {
    #replace params as ggplot objects are mutable
    oldparams = plotObj$layers[[1]]$aes_params
    plotObj$layers[[1]]$aes_params = NULL
    ply = plotly::ggplotly(plotObj)
    plotObj$layers[[1]]$aes_params = oldparams

    #add layer with new data
    npgpal = pal_npg('nrc')(length(levels(newdata$Annotation)))
    ply = ply %>%
      add_trace(data = newdata,
                color = ~Annotation,
                colors = npgpal,
                type = 'scatter',
                mode = 'markers',
                marker = list(
                  size = 10,
                  line = list(color = 'white', width = 2)
                ),
                text = paste('Cell line:', newdata$SampleID)) %>%
      layout(showlegend = TRUE,
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
#### =============================== plotRankDensity() =========================
################################################################################

#' Plot the density of ranks for one sample
#' @description This function takes a single column data frame, which is a
#' subset of the ranked data obtained from rankGenes() function, and gene sets,
#' and it returns plots visualising the density and the rugs of the ranks.
#'
#' @param rankData one column of the data.frame obtained from the
#'   rankExpr()/rankGenes() function
#' @param isInteractive Boolean, determin whether the returned plot is
#'   interactive
#' @param textSize numberic, set the size of text on the plot
#' @param upSet GeneSet object, up regulated gene set
#' @param downSet GeneSet object, down regulated gene set
#'
#' @return A ggplot object (optionally interactive) demonstrating the rank
#'   density along with rug plot
#' @examples
#' ranked <- rankExpr(toy_expr)
#'
#' @export
plotRankDensity <- function (rankData,
                             upSet,
                             downSet = NULL,
                             isInteractive = F,
                             textSize = 1.5) {
  #values needed for calculating the boundaries
  upSigSize = length(geneIds(upSet))
  nTotalGenes = nrow(rankData)

  #check if there are some missing genes in the geneset
  missingGenes = setdiff(geneIds(upSet), rownames(rankData))
  if (length(missingGenes) > 0) {
    warningMsg = paste(length(missingGenes), 'genes missing:', sep = ' ')
    warningMsg = paste(warningMsg, paste(missingGenes, collapse = ', '),
                       sep = ' ')
    warning(warningMsg)
  }

  #remove missing genes from signature for further analysis
  geneIds(upSet) = setdiff(geneIds(upSet), missingGenes)
  upRanks = rankData[geneIds(upSet), , drop = F] / nrow(rankData)
  upRank = data.frame(upRanks, type = "Up Gene-set")
  allRanks = upRank

  if (!is.null(downSet)) {
    #check if there are some missing genes in the geneset
    missingGenes = setdiff(geneIds(downSet), rownames(rankData))
    if (length(missingGenes) > 0) {
      warningMsg = paste(length(missingGenes), 'genes missing:',
                         sep = ' ')
      warningMsg = paste(warningMsg,
                         paste(missingGenes, collapse = ', '), sep = ' ')
      warning(warningMsg)
    }

    #remove missing genes from signature for further analysis
    geneIds(downSet) = setdiff(geneIds(downSet), missingGenes)
    downRanks = rankData[geneIds(downSet), , drop = F] / nrow(rankData)
    downRank = data.frame(downRanks, type =  "Down Gene-set")
    allRanks = rbind(upRank, downRank)
  }

  colnames(allRanks) <- c("Ranks", "upDown")
  allRanks$EntrezID <- row.names(allRanks)


  #bar plot preparations
  ymap = c(0, 0)
  yendmap = ymap + 0.3
  colmap = c(RColorBrewer::brewer.pal(8, "Set1")[c(1, 2)])
  typemap = c('Up-regulated gene', 'Down-regulated gene')
  names(colmap) = names(ymap)  = c('Up Gene-set', 'Down Gene-set')
  names(yendmap) = names(typemap) = c('Up Gene-set', 'Down Gene-set')

  #plot density and calculate max density and barcode line heights and
  #positions
  p = ggplot(allRanks, aes(x = Ranks, col = upDown)) +
    stat_density(aes(y = ..density..), geom = 'line', position = 'identity')

  dens =  ggplot_build(p)$data[[1]]$density
  ymap[1] = round(max(dens), digits = 1) + 0.1
  ymap[2] = round(min(dens), digits = 1) - 0.1
  bcheight = (max(dens) - min(dens))
  bcheight = bcheight/ifelse(is.null(downSet), 4, 3)
  yendmap = ymap + c(1, -1) * bcheight

  #plot barcode plot
  p = p + geom_segment(aes(
    y = ymap[upDown],
    xend = Ranks,
    yend = yendmap[upDown],
    text = paste0(typemap[upDown], '\nGene symbol: ', EntrezID)
  ), alpha = 0.8) +
    scale_colour_manual(values = colmap,
                        guide = guide_legend(title = "Type"))

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
    ply = ply %>% layout(
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
