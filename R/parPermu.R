#' @include simpleScoring.R CoreFuns.R
NULL
#' @title permutation test for the scored samples with gene set size and
#' @description
#'  This function randomize the genes in the gene set and calculate the sample's score against the
#'  random gene sets which form the null distribution for the simple score test. The resulting
#'  permutation score results are used to calculate the empirical P value.
#' @param n_up numeric  size of up set
#' @param n_down numeric size of down set
#' @param ranked_sample matrix, outcome of function rankGenes()
#' @param B numeric number of permuting times
#' @param seed nunmeric set the seed for randomization
#'
#' @return a list consists of empircal scores, p values and the tested result scores
#' @export
#'
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_up, downSet = toy_dn)
#' n_up = length(GSEABase::geneIds(toy_up))
#' n_down = length(GSEABase::geneIds(toy_dn))
#' #no_cores <- detectCores() - 1
#' #This permutation function can be run using parallel scripts, refer to the
#' #vignette for examples
#' #cl <- makeCluster(no_cores)
#' #registerDoParallel(cl)
#' permuteResult = permute_null(n_up = n_up, n_down = n_down, ranked, B =10,
#' seed = 1)
#' #stopCluster(cl)
#' #registerDoSEQ()
permute_null <- function(n_up, n_down, ranked_sample, B = 100, seed = 1){
  set.seed(seed)
  all_genes  <-  rownames(ranked_sample)
  totalNo   <-    n_up + n_down
  temSets <- sapply(1:B, function(i) {
    sample(all_genes, size = totalNo, replace = FALSE)
  })
  r <- foreach(i = 1:B,
               .combine = rbind,
               .packages = "GSEABase",
               .export = c("GeneSet", "simpleScore", "geneIds")
  ) %dopar% {
    tms <- temSets[, i]

    if (n_down > 0) {
      upSet <-  GeneSet(as.character(tms[1:n_up]))
      downSet <-  GeneSet(as.character(tms[-(1:n_up)]))
      ss = simpleScore(ranked_sample, upSet = upSet, downSet = downSet)
    } else {
      #else all the random generated genes are in upSet
      ss = simpleScore(ranked_sample, upSet = GeneSet(as.character(tms)))
    }
    ss[, 1]
  }
  colnames(r) <- colnames(ranked_sample)
  return(r)
}
#' calculate the p values for tested sample scores
#' @description This function takes the permutation result, which is the
#'   empirical scores, and the test sample scores as input. It calculates the
#'   empirical p values of the simple sample scoring test.
#'
#' @param permuResult A matrix, result from permuteNull() function
#' @param scoredf A dataframe, result from simplescore() function
#'
#' @return pvals, the calculated empirical p values for all empirical null
#'   sample scores distribution
#'
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_up, downSet = toy_dn)
#' n_up = length(GSEABase::geneIds(toy_up))
#' n_down = length(GSEABase::geneIds(toy_dn))
#' #no_cores <- detectCores() - 1
#' #This permutation function can be run using parallel scripts, refer to the
#' #vignette for examples
#' #cl <- makeCluster(no_cores)
#' #registerDoParallel(cl)
#' permuResult = permute_null(n_up = n_up, n_down = n_down, ranked, B = 10,
#'  seed = 1)
#' #stopCluster(cl)
#' #registerDoSEQ()
#' pvals <- get_pval(permuResult,scoredf)
#' @export
get_pval <- function(permuResult,scoredf){
  resultSc <- t(scoredf[, 1, drop = F])
  #combine the permutation with the result score for the computation of P values
  #p = (r+1)/(m+1)
  empirScore_re <- rbind(permuResult, as.character(resultSc))
  pvals <- apply(empirScore_re,2,function(x){
    length(x[x >= x[length(x)]]) / length(x)
  })
  return(pvals)
}

#' Plot the empirical null distribution using the permutation result
#' @description
#' This function takes the result from function permuteNull(), plot the density
#' curve for the given samples and intercept the sample score on the x axis
#' @param permuResult A matrix, outcome from function permuteNull()
#' @param scoredf A dataframe, outcome from function simplescore()
#' @param pvals A vector, outcome from function get_pval()
#' @param sampleLabel character, one sample label or multiple sample labels
#' @param alpha ggplot theme
#' @param size ggplot theme
#' @param textSize ggplot theme
#' @param cutoff double, the cutoff value for significant p values
#'
#' @examples
#' ranked <- rankGenes(toy_expr)
#' scoredf <- simpleScore(ranked, upSet = toy_up, downSet = toy_dn)
#' n_up = length(GSEABase::geneIds(toy_up))
#' n_down = length(GSEABase::geneIds(toy_dn))
#' #This permutation function can be run using parallel scripts, refer to the
#' #vignette for examples
#' #no_cores <- detectCores() - 1
#' #using parallel scripts
#' #cl <- makeCluster(no_cores)
#' #registerDoParallel(cl)
#' permuResult = permute_null(n_up = n_up, n_down = n_down, ranked, B = 10,
#'  seed = 1)
#' #stopCluster(cl)
#' #registerDoSEQ()
#' pvals <- get_pval(permuResult,scoredf)
#' plot_null(permuResult,scoredf,pvals,sampleLabel = names(pvals))
#' plot_null(permuResult,scoredf,pvals,sampleLabel = names(pvals)[1])
#' @export
plot_null <- function(permuResult, scoredf, pvals, sampleLabel = NULL,
                      alpha = 1, size = 1, textSize = 1, cutoff = 0.01){
  if(!is.null(sampleLabel)){
    pvals <- pvals[sampleLabel,drop=F]
  }
  pval_r <- as.character(round(pvals[sampleLabel],3))
  pvalTitle <- paste0(sampleLabel,' pv = ',pval_r)
  names(pvalTitle) <- names(pvals[sampleLabel])
  sig_colr <- c()
  for(i in 1:length(pvals)){
    if(pvals[i] < cutoff){
      sig_colr[i] = 'significant'
    }else
    {
      sig_colr[i] = 'not significant'
    }
  }
  colr_annot = data.frame(sampleLabel = sampleLabel, color = sig_colr)
  #pDt <-  as.data.frame(pvals)
  if( !is.null(sampleLabel) ){
    if(length(sampleLabel)>1){
      dt <- as.data.frame(permuResult[,sampleLabel])
      longDt <- reshape::melt(dt,variable_name = "sampleLabel")
      resultScs <-  scoredf[,1,drop = F]
      resultScs$sampleLabel <-  rownames(resultScs)
      #pDt$sampleLabel <- names(pvals)
      sampleLSc <-  merge(longDt, resultScs, by.x = "sampleLabel", by.y = "sampleLabel")
      #plotDt  <-  merge(sampleLSc,pDt, by.x = 'sampleLabel', by.y = 'sampleLabel')
      sampleLSc <- merge(sampleLSc,colr_annot)
      plotObj <-  ggplot2::ggplot(data = sampleLSc)+
        geom_density(mapping = ggplot2::aes( x = value,fill = color))+
        facet_grid(sampleLabel~., labeller = labeller(sampleLabel= as_labeller(pvalTitle)))+
        geom_vline(mapping =  ggplot2::aes(xintercept  = TotalScore), linetype="dashed", colour = 'red',size = 1)+
        xlab("Null Distribution & Calculated score")+
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
            hjust = 0.5)
        )
    }else{
      plotDt <- data.frame(sampleLabel = sampleLabel, value = permuResult[,sampleLabel])
      plotObj <-  ggplot2::ggplot(data = plotDt)+
        geom_density(mapping = ggplot2::aes(x = value))+
        geom_vline(mapping = ggplot2::aes(xintercept  = scoredf[sampleLabel,1]))+
        xlab(paste0(sampleLabel,"-Null Distribution & Calculated score"))+
        ggtitle(paste0("pval =", round(pvals,3))) +
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
            hjust = 0.5)
        )
    }
  }else{
    warning("Please provide which sample's null distribution to plot")
  }
  plotObj
}

