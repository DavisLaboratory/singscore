#' @include singscore.R rankAndScoring.R
NULL

#' @title Permutation test for the derived scores of each sample
#'
#' @description This function generates a number of random gene sets that
#'   have the same number of genes as the scored gene set. It scores each random
#'   gene set and returns a matrix of scores for all samples. 
#'   The empirical scores are used to calculate the empirical p-values and plot
#'   the null distribution. The implementation uses [BiocParallel::bplapply()] 
#'   for easy access to parallel backends. Note that one should pass the same 
#'   values to the `upSet`, `downSet`, `centerScore` and `knownDirection` 
#'   arguments as what they provide for the `simpleScore()` function to generate
#'   a proper null distribution.
#' @param downSet A GeneSet object, down regulated gene set
#' @param rankData matrix, outcome of function [rankGenes()]
#' @param centerScore A Boolean, specifying whether scores should be centered
#'  around 0, default as TRUE
#' @param knownDirection A boolean flag, it deterimines whether the scoring
#'  method should derive the scores in a directional mannar when the gene
#'  signature only contains one set of gene set (passing the gene set via
#'  upSet). It is default as TRUE but one can set the argument to be FALSE to
#'  derive the score for a single gene set in a undirectional way. This
#'  parameter becomes irrelevant when both upSet and downSet are provided.
#' @param B integer, the number of permutation repeats or the number of random 
#' gene sets to be generated, default as 1,000
#' @param ncores, integer, the number of CPU cores the function can use
#' @param seed integer, set the seed for randomisation
#' @param useBPPARAM, the backend the function uses, if NULL is provided, the
#' function uses the default parallel backend which is the first on the list 
#' returned by \code{BiocParallel::registered()} i.e 
#' \code{BiocParallel::registered()[[1]]} for your machine. It can be changed 
#' explicitly by passing a BPPARAM 
#'
#' @return A matrix of empirical scores for all samples
#' @seealso 
#' [Post about BiocParallel](http://lcolladotor.github.io/2016/03/07/BiocParallel/#.WgXMF61L28U)
#' `browseVignettes("BiocParallel")`
#' @keywords internal
#' @author Ruqian Lyu

generateNull_intl <- function(upSet, downSet = NULL, rankData, 
                              centerScore = TRUE,
                              knownDirection = TRUE, 
                              B = 1000, ncores = 1,
                              seed = sample.int(1E6, 1), useBPPARAM = NULL){
  n_up <- length(GSEABase::geneIds(upSet))
  if(is.null(downSet)){
    n_down <- 0
  } else {
    n_down <- length(GSEABase::geneIds(downSet))
  }
  all_genes <- rownames(rankData)
  totalNo <- n_up + n_down
  
  # If user does not supply a preferred BPPARAM , the function goes with
  # the default one for the system
  
  if(is.null(useBPPARAM)){
    useBPPARAM <- BiocParallel::registered()[[1]]
  }
  
  useBPPARAM$RNGseed <- seed
  useBPPARAM$workers <- ncores
  r <- BiocParallel::bplapply(1:B, function(i) {
    #load libraries on workers
    requireNamespace("GSEABase")
    requireNamespace("singscore")
    
    #execute code
    tms <-  sample(all_genes, size = totalNo, replace = FALSE)
    if (n_down > 0) {
      upSet <- GeneSet(as.character(tms[1:n_up]))
      downSet <-  GeneSet(as.character(tms[-(1:n_up)]))
      ss <- singscoring(rankData, upSet = upSet, downSet = downSet, 
                        centerScore = centerScore)
    } else {
      #else all the random generated genes are in upSet
      if(!knownDirection){
        ss <- singscoringOneGS(rankData, upSet = GeneSet(as.character(tms)), 
                               centerScore = centerScore)
      } else {
        ss <- singscoring(rankData, upSet = GeneSet(as.character(tms)), 
                          centerScore = centerScore)
      }
    }
    ss[, 1]
  }, BPPARAM = useBPPARAM)
  
  r <- plyr::ldply(r)
  colnames(r) <- colnames(rankData)
  return(r)
}
#' Estimate the empirical p-values
#'
#' @description With null distributions estimated using the [generateNull()]
#'   function, p-values are estimated using a one-tailed test. A minimum p-value
#'   of 1/B can be achieved with B permutations.
#'
#' @param permuteResult A matrix, null distributions for each sample generated
#'   using the [generateNull()] function
#' @param scoredf A dataframe, the scored results of samples under test
#'   generated using the [simpleScore()] function
#'
#' @return Estimated p-values for enrichment of the signature in each sample. A
#'   p-value of 1/B indicates that the estimated p-value is less
#'   than or equal to 1/B.

#' @examples
#' ranked <- rankGenes(toy_expr_se)
#' scoredf <- simpleScore(ranked, upSet = toy_gs_up, downSet = toy_gs_dn)
#' # find out what backends can be registered on your machine
#' BiocParallel::registered()
#' # the first one is the default backend, and it can be changed explicitly.
#' # See vignette for more details
#' permuteResult = generateNull(upSet = toy_gs_up, downSet = toy_gs_dn, ranked, 
#' B =10, seed = 1, useBPPARAM = NULL) 
#' 
#' # call the permutation function to generate the empirical scores 
#' # for B times.
#' pvals <- getPvals(permuteResult,scoredf)
#' @export
getPvals <- function(permuteResult,scoredf){
  
  stopifnot(dim(permuteResult)[2] == dim(scoredf)[1])
  
  resultSc <- t(scoredf[, 1, drop = FALSE])
  # combine the permutation with the result score for the computation of P values
  # p = (r+1)/(m+1)
 
  B <- nrow(permuteResult)
  # x[length(x)] is the calculated score
  pvals <- colSums(permuteResult > matrix(1, nrow = nrow(permuteResult), 
                                  ncol = 1) %*% resultSc) / B
  pvals <- sapply(pvals, max, 1/B)
  return(pvals)
}

