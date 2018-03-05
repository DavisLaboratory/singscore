#' singscore: A package for deriving gene-set scores at a single sample level
#'
#' The package provides functions for calculating gene-set enrichment scores at
#' a single-sample level using gene expression data. It includes functions to
#' perform hypothesis testing and provides visualisations to enable diagnosis of
#' scores and gene sets along with visualisations to enable exploration of
#' results.
#' 
#' @import stats
#' @import graphics
#' @import ggplot2
#' @import GSEABase
#' @import tidyr
#' @import methods
#' @import grDevices
#' @importFrom magrittr "%>%"
#' @importClassesFrom edgeR DGEList
#' @importClassesFrom Biobase ExpressionSet
#' @docType package
#' @name singscore
NULL

#' A toy gene expression dataset of two samples
#'
#' A toy dataset consisting of 2 samples with the expression values of 20 genes.
#' The data was created by sampling 2 samples and 20 genes from the dataset by
#' Foroutan et al, 2017.
#'
#' @format A data frame of 2 samples each with 20 genes \describe{
#'   \item{D_Ctrl_R1}{a control sample} \item{D_TGFb_R1}{a TGFb-treated sample}
#'   }
#' @docType data
#' @references Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, Erik W
#'   Thompson, and Melissa J Davis. 2017. “A Transcriptional Program for
#'   Detecting Tgfbeta-Induced Emt in Cancer.” Molecular Cancer Research.
#'   American Association for Cancer Research.
#'   doi:10.1158/1541-7786.MCR-16-0313.
#' @source [Foroutan et
#'   al.,2017](http://mcr.aacrjournals.org/content/early/2017/01/21/1541-7786.MCR-16-0313)
#'
#'
"toy_expr"

#' A gene set object of down-regulated genes for the toy dataset
#'
#' A GeneSet object with 5 genes randomly selected from the toy dataset. These
#' genes are independent of those in [toy_gs_up]
#'
#' @format A GSEABase::GeneSet object with 5 genes
#' @docType data
#' @seealso \code{"\linkS4class{GeneSet}"},[toy_expr],[toy_gs_up]
"toy_gs_dn"

#' A gene set object of up-regulated genes for the toy dataset
#'
#' A GeneSet object with 5 genes randomly selected from the toy dataset. These
#' genes are independent of those in [toy_gs_dn]
#'
#' @format A GeneSet object with 5 genes
#' @docType data
#' @seealso 
#' \code{"\linkS4class{GeneSet}"},[toy_expr],[toy_gs_dn]
"toy_gs_up"

#' A sample gene expression data.frame
#'
#' A microarray gene expression dataset that was originally obtained from the
#' integrated TGFb-EMT data published by (Foroutan et al, 2017). (ComBat
#' corrected values). \code{tgfb_expr_10} is a subset of the integrated TGFb-EMT
#' data consisting of 10 samples (4 TGFb treated and 6 controls) each with
#' expression values for 11900 genes.
#'
#' @references Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, Erik W
#'   Thompson, and Melissa J Davis. 2017. “A Transcriptional Program for
#'   Detecting Tgfbeta-Induced Emt in Cancer.” Molecular Cancer Research.
#'   American Association for Cancer Research.
#'   doi:10.1158/1541-7786.MCR-16-0313.
#' @source [Foroutan et
#'   al,2017](http://mcr.aacrjournals.org/content/early/2017/01/21/1541-7786.MCR-16-0313)
#'
#' @format A data.frame object
#' @docType data
#'   
"tgfb_expr_10"

#' Gene set of up-regulated genes for the TGFb-induced EMT gene signature
#'
#' A GeneSet object that contains the up-regulated genes of the TGFb-induced EMT
#' gene signature that was derived by (Foroutan et al.,2017), using two
#' meta-analysis techniques. The gene signature contains an up-regulated gene
#' set (up-set) and a down-regulated gene set (down-set). Please refer to the
#' vignettes for the steps to acquire the exact data object.
#'
#' @format A GeneSet object
#' @docType data
#' @references Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, Erik W
#' Thompson, and Melissa J Davis. 2017. “A Transcriptional Program for Detecting
#' Tgfbeta-Induced Emt in Cancer.” Molecular Cancer Research. American
#' Association for Cancer Research. doi:10.1158/1541-7786.MCR-16-0313.
#' @source [Foroutan et
#' al,2017](http://mcr.aacrjournals.org/content/early/2017/01/21/1541-7786.MCR-16-0313)
#' @seealso \code{"\linkS4class{GeneSet}"},[tgfb_gs_dn]
"tgfb_gs_up"

#' Gene set of down-regulated genes for the TGFb-induced EMT gene signature
#'
#' A GeneSet object that contains the down-regulated genes of the TGFb-induced
#' EMT gene signature that was derived by (Foroutan et al,2017), using two
#' meta-analysis techniques. The gene signature contains an up-regulated gene
#' set (up-set) and a down-regulated gene set (down-set). Please refer to the
#' vignettes for the steps to acquire the exact data object.
#' @format A GeneSet object
#' @docType data
#' @references Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, Erik W
#'   Thompson, and Melissa J Davis. 2017. “A Transcriptional Program for
#'   Detecting Tgfbeta-Induced Emt in Cancer.” Molecular Cancer Research.
#'   American Association for Cancer Research.
#'   doi:10.1158/1541-7786.MCR-16-0313.
#' @source [Foroutan et
#'   al,2017](http://mcr.aacrjournals.org/content/early/2017/01/21/1541-7786.MCR-16-0313)
#'
#' @seealso \code{"\linkS4class{GeneSet}"},[tgfb_gs_up]
"tgfb_gs_dn"

#' Pre-computed scores of the CCLE dataset against an epithelial gene
#' signature
#'
#' This data.frame stores pre-computed scores of the CCLE dataset [Barretina et
#' al](https://www.nature.com/articles/nature11003) calculated using the
#' [simpleScore()] function against the epithelial gene signature from [Tan,
#' Tuan Zea et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287932/). The
#' data.frame has scores for 55 samples. Please refer to the vignettes for
#' instructions on how to obtain the full datasets.
#'
#' @seealso [scoredf_ccle_mes]
#' @references Barretina, Jordi, Giordano Caponigro, Nicolas Stransky, Kavitha
#'   Venkatesan, Adam A Margolin, Sungjoon Kim, Christopher J Wilson, et al.
#'   2012. “The Cancer Cell Line Encyclopedia Enables Predictive Modelling of
#'   Anticancer Drug Sensitivity.” Nature 483 (7391): 603–7.
#'
#'   Tan, Tuan Zea, Qing Hao Miow, Yoshio Miki, Tetsuo Noda, Seiichi Mori, Ruby
#'   Yun-Ju Huang, and Jean Paul Thiery. 2014–10AD. “Epithelial-Mesenchymal
#'   Transition Spectrum Quantification and Its Efficacy in Deciphering Survival
#'   and Drug Responses of Cancer Patients.” EMBO Molecular Medicine 6 (10).
#'   Oxford, UK: BlackWell Publishing Ltd: 1279–93. doi:10.15252/emmm.201404208.
#'   
"scoredf_ccle_epi"


#' Pre-computed scores of the CCLE dataset against a mesenchymal gene
#' signature
#'
#' This data.frame stores pre-computed scores of the CCLE dataset [Barretina et
#' al](https://www.nature.com/articles/nature11003) calculated using the
#' [simpleScore()] function against the mesenchymal gene signature from [Tan,
#' Tuan Zea et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287932/). The
#' data.frame has scores for 55 samples. Please refer to the vignettes for
#' instructions on how to obtain the full datasets.
#'
#' @seealso [scoredf_ccle_epi]
#' @references Barretina, Jordi, Giordano Caponigro, Nicolas Stransky, Kavitha
#' Venkatesan, Adam A Margolin, Sungjoon Kim, Christopher J Wilson, et al. 2012.
#' “The Cancer Cell Line Encyclopedia Enables Predictive Modelling of Anticancer
#' Drug Sensitivity.” Nature 483 (7391): 603–7.
#'
#' Tan, Tuan Zea, Qing Hao Miow, Yoshio Miki, Tetsuo Noda, Seiichi Mori, Ruby
#' Yun-Ju Huang, and Jean Paul Thiery. 2014–10AD. “Epithelial-Mesenchymal
#' Transition Spectrum Quantification and Its Efficacy in Deciphering Survival
#' and Drug Responses of Cancer Patients.” EMBO Molecular Medicine 6 (10).
#' Oxford, UK: BlackWell Publishing Ltd: 1279–93 doi:10.15252/emmm.201404208.
#' 
"scoredf_ccle_mes"

#' Pre-computed scores of the TCGA breast cancer gene expression matrix
#' against an epithelial signature
#'
#' This data.frame stores pre-computed scores of the
#' [TCGA](https://cancergenome.nih.gov) dataset calculated using the
#' [simpleScore()] function against the epithelial gene signature from [Tan,
#' Tuan Zea et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287932/).
#' Please refer to the vignettes for instructions on how to obtain the full
#' datasets.
#'
#' @seealso [scoredf_tcga_mes]
#' @references Tan, Tuan Zea, Qing Hao Miow, Yoshio Miki, Tetsuo Noda, Seiichi
#'   Mori, Ruby Yun-Ju Huang, and Jean Paul Thiery. 2014–10AD.
#'   “Epithelial-Mesenchymal Transition Spectrum Quantification and Its Efficacy
#'   in Deciphering Survival and Drug Responses of Cancer Patients.” EMBO
#'   Molecular Medicine 6 (10). Oxford, UK: BlackWell Publishing Ltd: 1279–93
#'   doi:10.15252/emmm.201404208.
#'   
"scoredf_tcga_epi"


#' Pre-computed scores of the TCGA breast cancer gene expression matrix
#' against a mesenchymal signature
#'
#' This data.frame stores pre-computed scores of the
#' [TCGA](https://cancergenome.nih.gov) dataset calculated using the
#' [simpleScore()] function against the mesenchymal gene signature from [Tan,
#' Tuan Zea et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287932/).
#' Please refer to the vignettes for instructions on how to obtain the full
#' datasets.
#'
#' @seealso [scoredf_tcga_epi]
#' @references Tan, Tuan Zea, Qing Hao Miow, Yoshio Miki, Tetsuo Noda, Seiichi
#'   Mori, Ruby Yun-Ju Huang, and Jean Paul Thiery. 2014–10AD.
#'   “Epithelial-Mesenchymal Transition Spectrum Quantification and Its Efficacy
#'   in Deciphering Survival and Drug Responses of Cancer Patients.” EMBO
#'   Molecular Medicine 6 (10). Oxford, UK: BlackWell Publishing Ltd: 1279–93
#'   doi:10.15252/emmm.201404208.
#'   
"scoredf_tcga_mes"
