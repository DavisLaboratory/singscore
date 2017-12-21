#' singscore: A package for deriving gene-set score at a single sample level.
#'
#' The packge provides functions of calculating gene-set expression scores at a 
#' single-sample level and it is able to deal with single-sample input. It also 
#' includes efficient visulisation function for exploring the results as well 
#' as permutation test for significance test and plot null distribution.
#' @import methods
#' @import stats
#' @import graphics
#' @import ggplot2
#' @import GSEABase
#' @import tidyr
#' @importFrom magrittr "%>%"
#' @importClassesFrom edgeR DGEList
#' @importClassesFrom Biobase ExpressionSet
#' @docType package
#' @name singscore
NULL

#' Gene expression dataset of two samples
#'
#' A microarray gene expression data.frame dataset that was originially obtained
#' from the integrated TGFb-EMT data published by (Foroutan et al, 2017).
#' (ComBat corrected values).\code{toy_expr} is a subset of the integrated
#' TGFb-EMT data consisting of 2 samples each with 20 genes.
#' 
#' @format A data frame of 2 samples each with 20 genes
#' \describe{
#' \item{D_Ctrl_R1}{a control sample}
#' \item{D_TGFb_R1}{a TGFb-treated sample}
#'  }
#' @docType data
#' @references 
#' Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, 
#' Erik W Thompson, and Melissa J Davis. 2017. 
#' “A Transcriptional Program for Detecting Tgfbeta-Induced Emt in Cancer.” 
#' Molecular Cancer Research. American Association for Cancer Research. 
#' doi:10.1158/1541-7786.MCR-16-0313.
#' @source 
#' [Foroutan et al.,2017](http://mcr.aacrjournals.org/content/early/2017/01/21/1541-7786.MCR-16-0313)
"toy_expr"

#' Gene set object
#'
#' A GeneSet object that represents the down-regulated gene set in a gene
#' signature and it consists of five genes which are randomly selected from the
#' gene list in \code{toy_expr} but has no overlap with \code{toy_gs_up}.
#' 
#' @format A GSEABase::GeneSet obeject with 5 genes
#' @docType data
#' @seealso 
#' \code{"\linkS4class{GeneSet}"},[toy_expr],[toy_gs_up]
"toy_gs_dn"

#' Gene set of 5 genes
#'
#' A GeneSet object that represents the up-regulated gene set in a gene
#' signature and it consists of five genes which are randomly selected from the
#' gene list in \code{toy_expr} but has no overlap with \code{toy_gs_dn}.
#'
#' @format A GeneSet obeject with 5 genes
#' @docType data
#' @seealso 
#' \code{"\linkS4class{GeneSet}"},[toy_expr],[toy_gs_dn]
"toy_gs_up"

#' Gene expression data.frame
#'
#' A microarray gene expression dataset that was originially obtained from the
#' integrated TGFb-EMT data published by (Foroutan et al, 2017). (ComBat
#' corrected values).\code{tgfb_expr_10} is a subset of the integrated TGFb-EMT
#' data consisting of ten samples each with 11900 genes.
#' @references 
#' Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, 
#' Erik W Thompson, and Melissa J Davis. 2017. 
#' “A Transcriptional Program for Detecting Tgfbeta-Induced Emt in Cancer.” 
#' Molecular Cancer Research. American Association for Cancer Research. 
#' doi:10.1158/1541-7786.MCR-16-0313.
#' @source
#' [Foroutan et al,2017](http://mcr.aacrjournals.org/content/early/2017/01/21/1541-7786.MCR-16-0313)
#' @format A data.frame obeject
#' @docType data
"tgfb_expr_10"

#' Gene expression data.frame
#'
#' A microarray gene expression dataset that was originially obtained
#' from the integrated TGFb-EMT data published by (Foroutan et al, 2017). (ComBat
#' corrected values).\code{tgfb_expr_10} is a subset of the integrated TGFb-EMT 
#' data consisting of ten samples each with 11900 genes.
#' @references 
#' Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, 
#' Erik W Thompson, and Melissa J Davis. 2017. 
#' “A Transcriptional Program for Detecting Tgfbeta-Induced Emt in Cancer.” 
#' Molecular Cancer Research. American Association for Cancer Research. 
#' doi:10.1158/1541-7786.MCR-16-0313.
#' @source
#' [Foroutan et al,2017](http://mcr.aacrjournals.org/content/early/2017/01/21/1541-7786.MCR-16-0313)
#' @format A data.frame obeject
#' @docType data
"tgfb_expr_20"

#' Up-set for TGFb-induced EMT gene signature
#'
#' A GeneSet object that contains the up-regulated genes of a TGFb-induced
#' EMT gene signature that was derived by (Foroutan et al.,2017), using two
#' meta-analysis techniques. The gene signature contains a up-regulated gene set
#' (up-set) and a down-regulated gene set (down-set). Please refer to the
#' vignettes for the steps to aquaire the exact data object.
#'
#' @format A GeneSet obeject
#' @docType data
#' @references
#' Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, 
#' Erik W Thompson, and Melissa J Davis. 2017. 
#' “A Transcriptional Program for Detecting Tgfbeta-Induced Emt in Cancer.” 
#' Molecular Cancer Research. American Association for Cancer Research. 
#' doi:10.1158/1541-7786.MCR-16-0313.
#' @source
#' [Foroutan et al,2017](http://mcr.aacrjournals.org/content/early/2017/01/21/1541-7786.MCR-16-0313)
#' @seealso 
#' \code{"\linkS4class{GeneSet}"},[tgfb_gs_dn]
"tgfb_gs_up"

#' Down-set for TGFb gene expression signature
#'
#' A GeneSet object that contains the up regulated gene set of a TGFb-induced
#' EMT gene signature that was derived by (Foroutan et al,2017), using two
#' meta-analysis techniques. The gene signature contains a up-regulated gene set
#' (up-set) and a down-regulated gene set (down-set). Please refer to the
#' vignettes for the steps to aquaire the exact data object.
#' @format A GeneSet obeject
#' @docType data
#' @references 
#' Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, 
#' Erik W Thompson, and Melissa J Davis. 2017. 
#' “A Transcriptional Program for Detecting Tgfbeta-Induced Emt in Cancer.” 
#' Molecular Cancer Research. American Association for Cancer Research. 
#' doi:10.1158/1541-7786.MCR-16-0313.
#' @source 
#' [Foroutan et al,2017](http://mcr.aacrjournals.org/content/early/2017/01/21/1541-7786.MCR-16-0313)
#' @seealso 
#' \code{"\linkS4class{GeneSet}"},[tgfb_gs_up]
"tgfb_gs_dn"

#' Scoring results of a CCLE dataset against an epithelial gene signature
#'
#' This data.frame is the returned results of function [singscoring()] on a
#' CCLE dataset [Barretina et al](https://www.nature.com/articles/nature11003) 
#' consists of 55 samples agaist an epithelial gene signature obtained from 
#' [Tan, Tuan Zea et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287932/) 
#' Please refer to the vignettes for the steps to aquaire the exact data object.
#' @seealso 
#' [scoredf_ccle_mes]
#' @references 
#' Barretina, Jordi, Giordano Caponigro, Nicolas Stransky, Kavitha Venkatesan,
#' Adam A Margolin, Sungjoon Kim, Christopher J Wilson, et al. 2012. 
#' “The Cancer Cell Line Encyclopedia Enables Predictive Modelling of Anticancer 
#' Drug Sensitivity.” Nature 483 (7391): 603–7.
#' 
#' Tan, Tuan Zea, Qing Hao Miow, Yoshio Miki, Tetsuo Noda, Seiichi Mori, 
#' Ruby Yun-Ju Huang, and Jean Paul Thiery. 2014–10AD. 
#' “Epithelial-Mesenchymal Transition Spectrum Quantification and Its Efficacy 
#' in Deciphering Survival and Drug Responses of Cancer Patients.” 
#' EMBO Molecular Medicine 6 (10). Oxford, UK: BlackWell Publishing Ltd: 1279–93.
#' doi:10.15252/emmm.201404208.
#' 
#' Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, 
#' Erik W Thompson, and Melissa J Davis. 2017. 
#' “A Transcriptional Program for Detecting Tgfbeta-Induced Emt in Cancer.” 
#' Molecular Cancer Research. American Association for Cancer Research. 
#' doi:10.1158/1541-7786.MCR-16-0313.
#' 
"scoredf_ccle_epi"


#' Scoring results of a CCLE dataset against an mesenchymal gene signature
#'
#' This data.frame is the returned results of function 
#' [singscoring()] on a CCLE dataset 
#' [Barretina et al](https://www.nature.com/articles/nature11003)
#' consists of 55 samples agaist an mesenchymal gene signature obtained from
#' [Tan, Tuan Zea et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287932/)
#' Please refer to the vignettes for the steps to aquaire the exact data object.
#'
#' @seealso 
#' \code{[scoredf_ccle_epi]}
#' @references 
#' Barretina, Jordi, Giordano Caponigro, Nicolas Stransky, Kavitha Venkatesan,
#' Adam A Margolin, Sungjoon Kim, Christopher J Wilson, et al. 2012. 
#' “The Cancer Cell Line Encyclopedia Enables Predictive Modelling of Anticancer 
#' Drug Sensitivity.” Nature 483 (7391): 603–7.
#' 
#' Tan, Tuan Zea, Qing Hao Miow, Yoshio Miki, Tetsuo Noda, Seiichi Mori, 
#' Ruby Yun-Ju Huang, and Jean Paul Thiery. 2014–10AD. 
#' “Epithelial-Mesenchymal Transition Spectrum Quantification and Its Efficacy 
#' in Deciphering Survival and Drug Responses of Cancer Patients.” 
#' EMBO Molecular Medicine 6 (10). Oxford, UK: BlackWell Publishing Ltd: 1279–93
#' doi:10.15252/emmm.201404208.
#' 
#' Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, 
#' Erik W Thompson, and Melissa J Davis. 2017. 
#' “A Transcriptional Program for Detecting Tgfbeta-Induced Emt in Cancer.” 
#' Molecular Cancer Research. American Association for Cancer Research. 
#' doi:10.1158/1541-7786.MCR-16-0313.
"scoredf_ccle_mes"

#' Scoring results of TCGA tumour gene expression matrix
#'
#' This data.frame is the returned results of function [singscoring()] on
#' tumour samples from [TCGA](https://cancergenome.nih.gov) database against 
#' an mesenchymal gene signature obtained from 
#' [Tan, Tuan Zea et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287932/).
#' Please refer to the vignettes for the steps to aquaire the exact data object.
#'
#' @seealso 
#' [scoredf_tcga_mes]
#' @references 
#' Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, 
#' Erik W Thompson, and Melissa J Davis. 2017. 
#' “A Transcriptional Program for Detecting Tgfbeta-Induced Emt in Cancer.” 
#' Molecular Cancer Research. American Association for Cancer Research. 
#' doi:10.1158/1541-7786.MCR-16-0313.
"scoredf_tcga_epi"

#' Scoring results of TCGA tumour gene expression matrix
#'
#' This data.frame is the returned results of function [singscoring()] on
#' tumour samples from [TCGA](https://cancergenome.nih.gov) database against 
#' an epithelial gene signature obtained from 
#' [Tan, Tuan Zea et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287932/).
#' Please refer to the
#' vignettes for the steps to aquaire the exact data object.
#'
#' @seealso 
#' [scoredf_tcga_epi]
#' @references 
#' Foroutan, Momeneh, Joseph Cursons, Soroor Hediyeh-Zadeh, 
#' Erik W Thompson, and Melissa J Davis. 2017. 
#' “A Transcriptional Program for Detecting Tgfbeta-Induced Emt in Cancer.” 
#' Molecular Cancer Research. American Association for Cancer Research. 
#' doi:10.1158/1541-7786.MCR-16-0313.
"scoredf_tcga_mes"
