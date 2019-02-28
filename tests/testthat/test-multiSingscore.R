library(GSEABase)
context("multiSingscore")

test_that("checkGenesMulti works", {
  #set seed for reproducibility
  set.seed(36)

  #generate geneset collection
  gsl = lapply(1:10, function(x)
    GeneSet(as.character(sample.int(100, 5)), setName = LETTERS[x]))
  gsc = GeneSetCollection(gsl)

  gsc_filt = checkGenesMulti(gsc, as.character(1:20))
  expect_lt(length(gsc_filt), length(gsc))
  expect_equal(all(sapply(gsc_filt, function(x) length(geneIds(x))) < 5), TRUE)
})

test_that("multiSingscore works", {
  #set seed for reproducibility
  set.seed(36)

  #generate geneset collection
  gsl = lapply(1:20, function(x)
    GeneSet(as.character(sample.int(100, 5)), setName = LETTERS[ceiling(x / 2)]))
  gsc_up = GeneSetCollection(gsl[(1:20) %% 2 == 0])
  gsc_dn = GeneSetCollection(gsl[(1:20) %% 2 == 1])

  #generate data
  emat = matrix(rnorm(6 * 80), ncol = 6)
  rownames(emat) = as.character(21:100)
  eranks = rankGenes(emat)

  expect_equal(length(multiSingscore(eranks, gsc_up)), 2)
  expect_equal(ncol(multiSingscore(eranks, gsc_up)$Score), 6)
  expect_equal(ncol(multiSingscore(eranks, gsc_up, gsc_dn)$Score), 6)
  expect_equal(ncol(multiSingscore(eranks, gsc_up, knownDirection = FALSE)$Score), 6)
})
