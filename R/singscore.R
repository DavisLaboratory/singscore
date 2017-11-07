#' Singscore: A package for deriving gene-set score at a single sample level.
#'
#' The packge provides functions of calculating gene-set expression scores at a 
#' single-sample level and it is able to deal with single-sample input. It also 
#' includes efficient visulisation function for exploring the results as well as 
#' permutation test for significance test and plot null distribution.
#' @import methods
#' @import stats
#' @import graphics
#' @import ggplot2
#' @import GSEABase
#' @importClassesFrom edgeR DGEList
#' @importClassesFrom Biobase ExpressionSet
#' @docType package
#' @name singscore
NULL

#' Gene expression dataset of 100 genes in ten samples
#'
#' A gene by sample matrix-like dataset with value representing the gene's
#' expression intensity in a sample
#'
#' @format A data frame with 100 rows and 10 columns:
#' \describe{
#' \item{D_Ctrl_R1}{a control sample label}
#' \item{D_TGFb_R1}{a tgfb sample label}
#'  ...
#'  }
#' @docType data
#' @source \url{http://www.somelinktotgfbdata.info/}
"toy_expr"

#' Gene set of 5 genes
#'
#' A GeneSet object that has ten genes
#'
#'
#' @format A GeneSet obeject with 5 genes
#' @docType data
#' @source \url{http://www.somelinktotgfbdata.info/}
"toy_dn"

#' Gene set of 5 genes
#'
#' A GeneSet object that has ten genes
#'
#'
#' @format A GeneSet obeject with 5 genes
#' @docType data
#' @source \url{http://www.somelinktotgfbdata.info/}
"toy_up"

#' Gene expression data.frame
#'
#' A data.frame object that has 10 samples selected from a TFGb full dataset
#'
#'
#' @format A data.frame obeject
#' @docType data
#' @source \url{http://www.somelinktotgfbdata.info/}
"tgfb_expr_10"

#' Up-set for TGFb gene expression signature
#'
#' A GeneSet object
#'
#'
#' @format A GeneSet obeject
#' @docType data
#' @source \url{http://www.somelinktotgfbdata.info/}
"tgfb_gs_up"

#' Down-set for TGFb gene expression signature
#'
#' A GeneSet object
#'
#'
#' @format A GeneSet obeject
#' @docType data
#' @source \url{http://www.somelinktotgfbdata.info/}
"tgfb_gs_dn"

