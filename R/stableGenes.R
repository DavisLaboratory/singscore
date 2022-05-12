#' @title Get a list of stably expressed genes
#'
#' @description Get a list of genes that are stably expressed in cancer and
#'   normal solid tissue.
#'
#' @param n_stable numeric, number of stable genes to retrieve
#' @param type character, type of stable genes requested, stable genes in
#'   "carcinoma" or stable genes in "blood"
#' @param id character, gene identifier required. This can be either "geneid"
#'   for symbols or "ensembl" ensembl id)
#'
#' @return a character vector with gene IDs sorted by their expected expression
#'   levels in the requested tissue
#' @export
#'
#' @examples
#' getStableGenes(5)
#' getStableGenes(5, id = 'ensembl')
#' getStableGenes(5, type = 'blood')
#' 
getStableGenes <- function(n_stable, type = c('carcinoma', 'blood', 'protein'), id = c('geneid', 'ensembl')) {
  #get params
  type = match.arg(type)
  id = match.arg(id)
  
  #select datasets to use
  dsnames = switch (type,
    'carcinoma' = c('TCGA carcinomas', 'CCLE carcinoma-derived'),
    'blood' = c('TCGA other', 'CCLE other'),
    'protein' = c('CCLE TMT', 'CPTAC Colon LFQ', 'CPTAC Breast Ratios', 'CPTAC3 GBM Ratios')
  )
  
  #compute combined ranks
  rankmat = st_ranks[, c(id, 'type')]
  rankmat$prodrank = rowSums(log(st_ranks[, dsnames]))
  rankmat = reshape2::acast(rankmat,
                            stats::formula(paste(id, 'type', sep = '~')),
                            value.var = 'prodrank',
                            fun.aggregate = mean,
                            na.rm = TRUE)
  rankmat = apply(rankmat, 2, rank, na.last = 'keep')
  rankmat = rankmat[!is.na(rankmat[, 'strank']), , drop = FALSE]
  rankmat = rankmat[rankmat[, 'strank'] <= n_stable, ]
  
  #order genes by expression and return
  genes = rownames(rankmat)[order(rankmat[, 'exprank'])]
  if (length(genes) < n_stable)
    warning('number of stable genes requested exceeds the number analysed')
  
  return(genes)
}
