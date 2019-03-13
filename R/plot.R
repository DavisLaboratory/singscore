#' @include singscore.R
#' @import ggplot2
NULL

#default theme
getTheme <- function(rl = 1.2) {
  current_theme = ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = element_rect(colour = 'black', fill = NA),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = rel(rl) * 1.1),
      axis.text = element_text(size = rel(rl)),
      plot.title = element_text(size = rel(rl ^ 2)),
      strip.background = element_rect(fill = NA, colour = 'black'),
      strip.text = element_text(size = rel(rl)),
      legend.text = element_text(size = rel(rl)),
      legend.title = element_text(size = rel(rl), face = 'italic'),
      legend.position = 'bottom',
      legend.direction = 'horizontal'
    )

  return(current_theme)
}

isScoreCol <- function(c) {
  cnames = paste0(rep(c('Total', 'Up', 'Down'), each = 2), c('Score', 'Dispersion'))
  return(c %in% cnames)
}

#process annotation field for plots. annot could be a vector of annotations, or
#a column name from the score dataframe
processAnnotation <- function(df, annot) {
  #return vector of empty strings if nothing is specified
  if (is.null(annot))
    annot = rep('', nrow(df))

  #characters will either be column names or character annotations
  if (is.character(annot) && length(annot) == 1 && annot %in% colnames(df)) {
  	#a column name has been specified, extract the annotation
  	annot = df[[annot]]
  }

  if (is.character(annot)) {
    #convert char annot to factor
    annot = as.factor(annot)
  }

  #check length of annotation matches number of observations
  stopifnot(length(annot) == nrow(df))

  #do nothing if numeric
  return(annot)
}

getColorScale <- function(annot) {
  #specify a discrete scale for categoricals
  if (is.factor(annot)) {
    #specify a discrete scale for categoricals
    if (length(levels(annot)) > 8) {
      warning('Too many levels of the annotation, using default ggplot2 colours')
      return(NULL)
    } else{
      return(ggplot2::scale_colour_brewer(palette = 'Dark2'))
    }
  }

  #specify a continous scale for numerics
  if (is.numeric(annot)) {
    return(ggplot2::scale_colour_viridis_c())
  }

  return(NULL)
}

################################################################################
#### =============================== plotDispersion() ==========================
################################################################################

#' Plot the score v.s. despersion for all samples
#' @description This function takes the output from the simpleScore() function
#'   and generates scatter plots of score vs. dispersion for the total
#'   score, the up score and the down score of samples. If you wish to use the
#'   plotting function but with some customized inputs (instead of outputs from
#'   `simpleScore` function), you need to make sure the formats are the same.
#'   To be specific, you need to have columns names "TotalScore"
#'   "TotalDispersion" "UpScore" "UpDispersion" "DownScore" "DownDispersion"
#'   and rows names as samples.
#' @param scoredf data.frame, generated using the [simpleScore()] function
#' @param annot any numeric, character or factor annotation provided by the user that
#' needs to be plot. Alternatively, this can be a character specifying the
#' column of scoredf holding the annotation. Annotations must be ordered in the
#' same way as the scores
#' @param annot_name character, legend title for the annotation
#' @param alpha numeric, set the transparency of points
#' @param size numeric, set the size of each point
#' @param textSize numeric, relative text sizes for title, labels, and axis
#' values
#' @param isInteractive Boolean, determine whether the plot is interactive
#' @examples
#' ranked <- rankGenes(toy_expr_se)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#' plotDispersion(scoredf)
#' plotDispersion(scoredf, isInteractive = TRUE)
#' @return A ggplot object
#' @export
plotDispersion <- function(scoredf, annot = NULL, annot_name = '', alpha = 1,
                           size = 1, textSize = 1.2, isInteractive=FALSE){
  #parameter type checks
  stopifnot(is.numeric(alpha), is.numeric(size), is.numeric(textSize),
            is.logical(isInteractive))

  #process annotations - set annot_name to annot if its a column reference
  if (is.null(annot_name) &&
      is.character(annot) &&
      length(annot) == 1 &&
      annot %in% colnames(scoredf)) {
    annot_name = annot
  }
  annot = processAnnotation(scoredf, annot)

  #transform data for plot
  plotdf = scoredf
  plotdf['SampleID'] = rownames(plotdf)
  plotdf['Class'] = annot
  plotdf['Type'] = 'Total'

  #for up-down signatures, melt the dataframe
  score_cols = isScoreCol(colnames(plotdf))
  if (sum(score_cols) > 2) {
    nsamp = nrow(plotdf)
    idvars = colnames(plotdf)[!score_cols]
    plotdf = reshape2::melt(plotdf, id.vars = idvars)
    plotdf$Type = rep(c('Total', 'Up', 'Down'), each = nsamp * 2)
    plotdf$variable = rep(c('Score', 'Dispersion'), each = nsamp)
    df_form = stats::as.formula(paste0(paste(idvars, collapse = ' + '), ' ~ variable'))
    plotdf = reshape2::dcast(plotdf, df_form, value.var = 'value')
  } else{
    sc_col = isScoreCol(colnames(plotdf))
    colnames(plotdf)[sc_col] = substring(colnames(plotdf)[sc_col], first = 6)
  }

  #convert Type to factors
  plotdf$Type = factor(plotdf$Type, levels = c('Total', 'Up', 'Down'))

  #setup plot
  p1 = ggplot(plotdf, aes(Score, Dispersion, colour = Class, text = SampleID)) +
    ggtitle('Score vs Dispersion') +
    geom_point(size = size, alpha = alpha) + #add scatter layer
    getColorScale(annot) + #add colour layer
    labs(colour = annot_name) + #annotation title
    getTheme(textSize) #specify the theme

  #facet up/down pair
  if (sum(score_cols) > 2) {
    p1 = p1 + facet_wrap( ~ Type, scales = 'free')
  }

  #remove legend if no annotation provided
  if (all(annot %in% '')) {
    p1 = p1 + theme(legend.position = "none")
  }

  if (isInteractive) {
    ply = plotly::ggplotly(p1)
    return(ply)
  } else{
    return(p1)
  }
}

################################################################################
#### =============================== plotScoreLandscape() ======================
################################################################################

#' Plot landscape of two gene signatures scores
#' @description This function takes two data frames which are outputs from the
#'   simpleScore() function and plots the relationship between the two gene set
#'   scores for samples in the gene expression matrix.Scoredf1 and Scoredf2 are
#'   two scoring results of the same set of samples against two different gene
#'   signatures. If you wish to use the plotting function but with some
#'   customized inputs (instead of outputs from the `simpleScore` function), you
#'   need to make sure the formats are the same To be specific, you need to have
#'   column names "TotalScore" "TotalDispersion" "UpScore" "UpDispersion"
#'   "DownScore" "DownDispersion" and rows names as samples.
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
#' @return A ggplot object, a scatter plot, demonstrating the relationship
#'   between scores from two signatures on the same set of samples.
#' @examples
#' ranked <- rankGenes(toy_expr_se)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#' scoredf2 <- simpleScore(ranked, upSet = toy_gs_up)
#' plotScoreLandscape(scoredf, scoredf2)
#' @export
plotScoreLandscape <- function(scoredf1, scoredf2, scorenames = c(),
                               textSize = 1.2, isInteractive = FALSE,
                               hexMin = 100){
	stopifnot(nrow(scoredf1) == nrow(scoredf2),
			  rownames(scoredf1) == rownames(scoredf2),
			  length(scorenames) %in% c(0, 2))

  #default axes labels
  if (length(scorenames) == 0) {
    scorenames = c('Signature 1', 'Signature 2')
  }

  #create data structure for plot
  plotdf = data.frame(scoredf1$TotalScore, scoredf2$TotalScore)
  colnames(plotdf) = c('sc1', 'sc2')

  #setup plot
  p1 = ggplot(plotdf, aes(sc1, sc2)) +
    ggtitle('Signature landscape') +
    xlab(scorenames[1]) + #label the x-axis
    ylab(scorenames[2]) + #label the y-axis
    getTheme(textSize) #specify the theme

  # choose plot according to the number of observations
  if (nrow(scoredf1) < hexMin) {
    #use scatter plot for fewer points
    p1 = p1 + geom_point(colour = 'blue')
  } else{
  	#number of bins to use for the hexbin plot
  	numbins = round(sqrt(nrow(plotdf)))
    p1 = p1 + geom_hex(colour = 'white', bins = numbins) +
      scale_fill_distiller(palette = 'RdPu', direction = 1)
  }

  if (isInteractive) {
    #replace params as ggplot objects are mutable
    #white colour for bins looks bad on interactive and static looks bad without
    oldparams = p1$layers[[1]]$aes_params
    p1$layers[[1]]$aes_params = NULL
    ply = plotly::ggplotly(p1)
    p1$layers[[1]]$aes_params = oldparams
    return(ply)
  } else{
    return(p1)
  }
}

################################################################################
####============================ projectScoreLandscape() =======================
################################################################################

#' Project data on the landscape plot obtained from \code{plotScoreLandscape()}
#'
#' @description This function takes the output (ggplot object) of the function
#'   \code{plotScoreLandscape()} and a new dataset. It projects the new data
#'   points onto the landscape plot and returns a new ggplot object with
#'   projected data points.
#'
#' @param plotObj a ggplot object, resulted from [plotScoreLandscape()]
#' @param scoredf1 data.frame, result of the simpleScore() function which scores
#'   the gene expression matrix against a gene set of interest
#' @param scoredf2 data.frame, result of the simpleScore() function which scores
#'   the gene expression matrix against another gene set of interest. Scores in
#'   scoredf1 and scoredf2 consist of the new data points that will be projected
#'   on the `plotObj` landscape plot.
#' @param subSamples vector of character or indices for subsetting the scoredfs,
#'   default as NULL and all samples in scoredfs will be plotted. The subsetted
#'   samples are projected onto the landscape plot of `plotObj`.
#' @param sampleLabels vector of character, sample names to display, ordered in
#'   the same way as samples are ordered in the 'scoredfs' data matrix, default
#'   as NULL which means the projected points are not labelled.
#' @param annot any numeric, character or factor annotation provided by the user
#'   that needs to be plot. Alternatively, this can be a character specifying
#'   the column of scoredf1 holding the annotation. Annotations must be ordered
#'   in the same way as the scores
#' @param annot_name character, legend title for the annotation
#' @param isInteractive boolean, whether the plot is interactive default as
#'   FALSE
#'
#' @return New data points on the already plotted ggplot object from
#'   plotScoreLanscape()
#' @seealso [plotScoreLandscape()]
#'  @examples
#'  ranked <- rankGenes(toy_expr_se)
#'  scoredf1 <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#'  scoredf2 <- simpleScore(ranked, upSet = toy_gs_up)
#'  psl <- plotScoreLandscape(scoredf1, scoredf2)
#'  projectScoreLandscape(psl,scoredf1, scoredf2)
#' @export
projectScoreLandscape <- function(plotObj = NULL,
                                  scoredf1,
                                  scoredf2,
                                  annot = NULL,
                                  annot_name = NULL,
                                  subSamples = NULL,
                                  sampleLabels = NULL,
                                  isInteractive = FALSE){
  #stop and print message
  if (!is.ggplot(plotObj)) {
    stop(
      'Please provide a ggplot object generated using plotScoreLandscape() (',
      class(plotObj)[1],
      ' object given)'
    )
  }

  #check that the score dataframe have the same number of samples and names are
  #are the same too
  stopifnot(dim(scoredf1)[1] == dim(scoredf2)[1],
            rownames(scoredf1) == rownames(scoredf2))

  #create data frame with the new data
  #subsetting the two data frames, scoredfs then checks the data dimensions
  if (!is.null(subSamples)) {
    scoredf1 = scoredf1[subSamples, ]
    scoredf2 = scoredf2[subSamples, ]
    if (anyNA(scoredf1)) {
      message('some selected samples not exist in provided scoredf1')
      scoredf1 = na.omit(scoredf1)
    }
    if (anyNA(scoredf2)) {
      message('some selected samples not exist in provided scoredf2')
      scoredf2 = na.omit(scoredf2)
    }
  }

  #if no sample labels are provided
  if (is.null(sampleLabels)) {
    sampleLabels = ""
  } else{
    if (length(sampleLabels) != nrow(scoredf1))
      stop(
        "sampleLabels must contain the same number of labels with the number of samples in scoredf"
      )
  }

  #process annotations - set annot_name to annot if its a column reference
  if (is.null(annot_name) &&
      is.character(annot) &&
      length(annot) == 1 &&
      annot %in% colnames(scoredf1)) {
    annot_name = annot
  }
  annot = processAnnotation(scoredf1, annot)
  newdata = data.frame(
    'sc1' = scoredf1$TotalScore,
    'sc2' = scoredf2$TotalScore,
    'SampleLabel' = sampleLabels,
    'Class' = annot
  )

  #setup plot
  p1 = plotObj +
    geom_point(
      aes(text = SampleLabel, colour = Class),
      shape = 21,
      fill = 'white',
      size = 2,
      stroke = 2,
      data = newdata
    ) + getColorScale(annot) +
    labs(colour = annot_name)

  #remove legend if no annotation provided
  if (all(annot %in% '')) {
    p1 = p1 + guides(colour = FALSE)
  }

  #need to deal with legends in both interactive and non-interactive
  if (!isInteractive) {
    #label samples
    p1 = p1 +
      ggrepel::geom_label_repel(
        data = newdata,
        aes(label = SampleLabel, colour = Class),
        show.legend = FALSE
      )
  } else if(isInteractive) {
    #replace params as ggplot objects are mutable
    oldparams1 = p1$layers[[1]]$aes_params
    oldparams2 = p1$layers[[2]]$aes_params
    p1$layers[[1]]$aes_params = NULL
    p1$layers[[2]]$aes_params = NULL
    ply = plotly::ggplotly(p1)
    p1$layers[[1]]$aes_params = oldparams1
    p1$layers[[2]]$aes_params = oldparams2
    return(ply)
  }

  return(p1)
}

################################################################################
#### =============================== plotRankDensity_intl() ====================
################################################################################

#' Plot the densities of ranks for one sample
#'
#' @description This function takes a single column data frame, which is a
#' subset of the ranked data obtained from [rankGenes()]function and gene sets,
#' and it returns plots visualising the density and the rugs of the ran ks.
#'
#' @param rankData one column of the ranked gene expression matrix obtained from
#' the [rankGenes()] function, use drop = FALSE when subsetting the ranked gene
#' expression matrix, see examples.
#' @param isInteractive Boolean, determin whether the returned plot is
#'   interactive
#' @param textSize numeric, set the size of text on the plot
#' @param upSet GeneSet object, up regulated gene set
#' @param downSet GeneSet object, down regulated gene set
#' @keywords internal
#'
#' @return A ggplot object (optionally interactive) demonstrating the rank
#'   density along with rug plot

#' @seealso
#' \code{"\linkS4class{GeneSet}"}
plotRankDensity_intl <- function (rankData,
                                  upSet,
                                  downSet = NULL,
                                  isInteractive = FALSE,
                                  textSize = 1.2) {
  stopifnot(is.logical(isInteractive), is.numeric(textSize))
  #values needed for calculating the boundaries
  upSigSize = length(geneIds(upSet))
  nTotalGenes = nrow(rankData)

  #remove missing genes from signature for further analysis
  missingGenes = setdiff(geneIds(upSet), rownames(rankData))
  geneIds(upSet) = setdiff(geneIds(upSet), missingGenes)
  upRanks = rankData[geneIds(upSet), , drop = FALSE] / nrow(rankData)
  upRank = data.frame(upRanks, type = "Up Gene-set")
  allRanks = upRank

  if (!is.null(downSet)) {
  	#remove missing genes from signature for further analysis
    missingGenes = setdiff(geneIds(downSet), rownames(rankData))
    geneIds(downSet) = setdiff(geneIds(downSet), missingGenes)
    downRanks = rankData[geneIds(downSet), , drop = FALSE] / nrow(rankData)
    downRank = data.frame(downRanks, type =  "Down Gene-set")
    allRanks = rbind(upRank, downRank)
  }
  colnames(allRanks) = c("Ranks", "upDown")
  allRanks$EntrezID = row.names(allRanks)

  #barcode plot preparations
  ymap = c(0, 0)
  yendmap = ymap + 0.3
  colmap = c(RColorBrewer::brewer.pal(8, "Accent")[c(1, 2)])
  typemap = c('Up-regulated gene', 'Down-regulated gene')
  names(colmap) = names(ymap)  = c('Up Gene-set', 'Down Gene-set')
  names(yendmap) = names(typemap) = c('Up Gene-set', 'Down Gene-set')

  #plot density and calculate max density and barcode line heights and
  #positions
  p1 = ggplot(allRanks, aes(x = Ranks, col = upDown)) +
  	stat_density(aes(y = ..density..), geom = 'line', position = 'identity')

  #update y-positions for barcode of up-set
  dens = ggplot_build(p1)$data[[1]]$density
  ymap[1] = round(max(dens), digits = 1) + 0.1
  ymap[2] = round(min(dens), digits = 1) - 0.1
  bcheight = (max(dens) - min(dens))
  bcheight = bcheight / ifelse(is.null(downSet), 4, 3)
  yendmap = ymap + c(1,-1) * bcheight

  #plot barcode plot
  #text aes useful for the plotly plot, so supress the warnings
  p1 = p1 +
  	geom_segment(aes(
  		y = ymap[upDown],
  		xend = Ranks,
  		yend = yendmap[upDown],
  		text = paste0(typemap[upDown], '\nGene symbol: ', EntrezID)
  	), alpha = 0.75) +
  	scale_colour_manual(values = colmap) +
  	labs(colour = 'Type') +
  	ggtitle('Rank density') +
  	xlab('Normalised Ranks') +
  	ylab('Density') +
  	getTheme(textSize)

  #if single geneset, remove legend
  if (is.null(downSet)) {
    p1 = p1 + theme(legend.position = 'none')
  }

  if (isInteractive) {
  	#Horizontal legend not supported by plotly yet so re-orient after
  	#creating plotly object
  	ply = plotly::ggplotly(p1, tooltip = c('text', 'x'))
  	ply = ply %>% plotly::layout(
  		legend = list(
  			orientation = 'h',
  			xanchor = 'center',
  			x = 0.5,
  			yanchor = 'top',
  			y = -0.25
  		),
  		yaxis = list(fixedrange = TRUE)
  	)
  	return(ply)
  } else{
  	return(p1)
  }
}

#' Plot the empirically estimated null distribution and associated p-values
#'
#' @description This function takes the results from function [generateNull()]
#' and plots the density curves of permuted scores for the provided samples via
#' \code{sampleNames} parameter. It can plot null distribution(s) for a single
#' sample or multiple samples.
#'
#' @param permuteResult A matrix, null distributions for each sample generated
#'   using the [generateNull()] function
#' @param scoredf A dataframe, singscores generated using the [simpleScore()]
#'   function
#' @param pvals A vector, estimated p-values using the [getPvals()] function
#' `permuteResult`,`scoredf` and `pvals` are the results for the same samples.
#'
#' @param sampleNames A character vector, sample IDs for which null
#'   distributions will be plotted
#' @param textSize numeric, size of axes labels, axes values and title
#' @param labelSize numeric, size of label texts
#' @param cutoff numeric, the cutoff value for determining significance
#' @return a ggplot object
#' @author Ruqian Lyu
#' @examples
#' ranked <- rankGenes(toy_expr_se)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)

#' # find out what backends can be registered on your machine
#' BiocParallel::registered()
#' # the first one is the default backend, and it can be changed explicitly.
#' permuteResult = generateNull(upSet = toy_gs_up, downSet = toy_gs_dn, ranked,
#' B =10, seed = 1,useBPPARAM = NULL)
#' # call the permutation function to generate the empirical scores
#' #for B times.
#' pvals <- getPvals(permuteResult,scoredf)
#' # plot for all samples
#' plotNull(permuteResult,scoredf,pvals,sampleNames = names(pvals))
#' #plot for the first sample
#' plotNull(permuteResult,scoredf,pvals,sampleNames = names(pvals)[1])
#' @export
plotNull <- function(permuteResult,
                     scoredf,
                     pvals,
                     sampleNames = NULL,
                     cutoff = 0.01,
                     textSize = 2,
                     labelSize = 5) {

  stopifnot(!is.null(sampleNames))

  quantile_title <- as.character((1 - cutoff)*100)
  if(is.null(sampleNames)){
    warning("Please provide which sample's null distribution to
            plot by specifying the sampleNames argument.")
  } else {
    pvals <- pvals[sampleNames, drop = FALSE]
    pval_r <- as.character(format(pvals[sampleNames],scientific = TRUE,
                                  digits = 3))
    pvalTitle <- paste0(' p-value = ',pval_r)
    names(pvalTitle) <- names(pvals[sampleNames])
    cutoff_score <- c()
    for(i in 1:length(sampleNames)){
      cutoff_score[i] <- quantile(permuteResult[,sampleNames[i]],(1-cutoff))
    }
    names(cutoff_score) <-  sampleNames
    cutoff_annot  <-  data.frame(sampleNames = sampleNames,
                                 cutoff_score = cutoff_score)
    #pDt <-  as.data.frame(pvals)
    if(length(sampleNames)>1){
      dt <- as.data.frame(permuteResult[,sampleNames])
      longDt <- reshape::melt(dt,variable_name = "sampleNames")
      resultScs <- scoredf[,1,drop = FALSE]
      resultScs$sampleNames <-  rownames(resultScs)
      #pDt$sampleNames <- names(pvals)
      sampleLSc <-  merge(longDt, resultScs, by.x = "sampleNames",
                          by.y = "sampleNames")
      #plotDt  <-  merge(sampleLSc,pDt, by.x = 'sampleNames',
      #by.y = 'sampleNames')
      sampleLSc <- merge(sampleLSc,cutoff_annot)
      sampleLSc <- merge(sampleLSc, pvalTitle)

      xlimStart <- min(dt, scoredf[,1]) - 0.01
      xlimEnd <- max(dt, scoredf[,1]) + 0.02
      value <- NULL
      TotalScore <- NULL
      #browser()
      plotObj <-  ggplot(data = sampleLSc)+
        geom_density(mapping = aes( x = value), size =1)+
        coord_cartesian(xlim = c(xlimStart,xlimEnd))+
        facet_grid(sampleNames~.)+
        geom_segment(mapping =  aes(x  = cutoff_score, y = 11,
                                    xend = cutoff_score, yend = 0),
                     linetype="dashed", colour = 'blue',size = 1)+
        geom_segment(mapping = aes(x  = TotalScore, y = 11, xend = TotalScore,
                                   yend = 0),colour = 'red',size = 2)+
        geom_text(mapping = aes(x  = TotalScore, y = 12,
                                label = pvalTitle[sampleNames]),
                  colour = 'red',size = labelSize)+
        geom_text(mapping = aes(x  = cutoff_score, y = 12,
                                label = paste0(quantile_title,
                                               '%-ile threshold')),
                  colour = 'blue',size = labelSize)+
        xlab("Scores")+
        ggtitle("Null distribution")
    } else {

      xlimStart <- min(permuteResult[,sampleNames],
                       scoredf[sampleNames,]$TotalScore) - 0.01
      xlimEnd <- max(permuteResult[,sampleNames],
                     scoredf[sampleNames,]$TotalScore) + 0.02
      plotDt <- data.frame(sampleNames = sampleNames,
                           value = permuteResult[,sampleNames],
                           TotalScore = scoredf[sampleNames,]$TotalScore)
      plotObj <-  ggplot(data = plotDt)+
        geom_density(mapping = aes( x = value),size = 1)+
        coord_cartesian(xlim = c(xlimStart,xlimEnd))+
        geom_segment(mapping =  aes(x = cutoff_score, y = 11,
                                    xend = cutoff_score, yend =0),
                     linetype = "dashed", colour = 'blue',size = 1)+
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

}
