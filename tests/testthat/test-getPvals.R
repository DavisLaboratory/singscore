context("getPvals")

test_that("getPvals works", {
  ranked <- rankGenes(toy_expr)
  
  n_down = length(GSEABase::geneIds(toy_gs_dn))
  # find out what backends can be registered on your machine
  # the first one is the default backend, and it can be changed explicitly.
  BiocParallel::registered()
  
  
  # call the permutation function to generate the empirical scores for B times.
  permuteResult = generateNull(upSet = toy_gs_up, downSet = toy_gs_dn, ranked, 
                               B = 10, ncores = 1, seed = 1, 
                               useBPPARAM = BiocParallel::registered()[[1]]) 
  scoredf = simpleScore(ranked,upSet = toy_gs_up, downSet = toy_gs_dn)
  pvals = getPvals(permuteResult, scoredf)
  expect_equivalent(length(pvals),2)
})
test_that("getPvals input checking works", {
  ranked <- rankGenes(toy_expr)
  
  n_down = length(GSEABase::geneIds(toy_gs_dn))
  # find out what backends can be registered on your machine
  # the first one is the default backend, and it can be changed explicitly.
  BiocParallel::registered()
  
  
  # call the permutation function to generate the empirical scores for B times.
  permuteResult = generateNull(upSet = toy_gs_up, downSet = toy_gs_dn, ranked, 
                               B = 10, ncores = 1, seed = 1, 
                               useBPPARAM = BiocParallel::registered()[[1]]) 
  scoredf = simpleScore(ranked,upSet = toy_gs_up, downSet = toy_gs_dn)
  
  expect_error(getPvals(permuteResult, scoredf[1,,drop=FALSE]))
})
