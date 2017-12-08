# Overview

'singscore' is a gene-set scoring R package which implements a simple single-sample gene-set (or gene-signature) scoring method. It uses rank-based statistics to analyze the each sample's gene expression profile and scores the expression activities of gene sets at a single-sample level.


## Getting Started

These instructions will get you to install the package up and running on your local machine. Please be aware that the package is still under development and we aim to submit the release version to [Bioconductor](https://www.bioconductor.org). If you experience any issues, please let us know via email 
<ul>
 <li> Ruqian Lyu <lyu.r@wehi.edu.au> </li>
 <li> Momeneh Foroutan <foroutan.m@wehi.edu.au></li> 
 <li> Dharmesh Bhuva <bhuva.d@wehi.edu.au></li>
 <li> Joseph Cursons <cursons.j@wehi.edu.au></li>
 <li> Melissa Davis <davis.m@wehi.edu.au></li>
</ul>

### Prerequisites
To install the package from git hub, install the package 'devtools' first and then use the function `install_github` to install 'singscore' by running the following script


### Installing

```
install.packages('devtools')
devtools::install_github('DavisLaboratory/singscore')
library(singscore)
```


## Running singscoring

```
#load the example data sets
# expression matrix
data('tgfb_expr_10')

# up gene set
data('tgfb_gs_up')
# down gene set
data('tgfb_gs_dn')

# rank the expression matrix first
rankedData <- rankExpr('tgfb_expr_10')

# Call singscoring to score each individual sample
scoredf <- singscoring(rankedData,tgfb_gs_up,tgfb_gs_dn)

scoredf
##             TotalScore TotalDispersion   UpScore UpDispersion    DownScore
## D_Ctrl_R1    0.4119020        5734.697 0.5609641     3119.390 -0.149062139
## D_TGFb_R1    0.7869942        4435.939 0.7493157     2352.886  0.037678558
## D_Ctrl_R2    0.4010359        5722.836 0.5684124     3129.769 -0.167376501
## D_TGFb_R2    0.7707220        4757.663 0.7503566     2470.012  0.020365345
## Hes_Ctrl_R1  0.4979152        5492.292 0.5804649     3134.216 -0.082549688
## Hes_TGFb_R1  0.6761228        5195.030 0.7289403     2416.638 -0.052817510
## Hes_Ctrl_R2  0.5168839        5401.112 0.5881783     3138.664 -0.071294412
## Hes_TGFb_R2  0.6884670        4910.371 0.7389547     2324.717 -0.050487776
## Hil_Ctrl_R1  0.4380088        6078.660 0.5831425     3553.792 -0.145133706
## Hil_Ctrl_R2  0.4350626        5918.539 0.5743386     3396.637 -0.139276000
## Hil_Ctrl_R3  0.4373703        6060.127 0.5823587     3478.180 -0.144988442
## Hil_TGFb_R1  0.7560978        4642.762 0.7339124     2300.995  0.022185474
## Hil_TGFb_R2  0.7315808        4579.751 0.7253359     2342.508  0.006244818
## Hil_TGFb_R3  0.7666968        4613.110 0.7418824     2167.561  0.024814375
##             DownDispersion
## D_Ctrl_R1         2615.306
## D_TGFb_R1         2083.053
## D_Ctrl_R2         2593.067
## D_TGFb_R2         2287.652
## Hes_Ctrl_R1       2358.075
## Hes_TGFb_R1       2778.392
## Hes_Ctrl_R2       2262.448
## Hes_TGFb_R2       2585.654
## Hil_Ctrl_R1       2524.868
## Hil_Ctrl_R2       2521.903
## Hil_Ctrl_R3       2581.948

```



## Produced rank density plot
`plotRankDensity()` takes a single column data frame, which is a subset of the ranked data obtained from `rankExpr()`function, and gene sets, and it returns plots visualising the density and the rugs of the ranks.
```
plotRankDensity(rankData[,1,drop = FALSE], upSet = tgfb_gs_up, 
                downSet = tgfb_gs_dn, isInteractive = FALSE)

```
Which can produce 

![](https://user-images.githubusercontent.com/12887308/33762071-52e5c1fe-dc5f-11e7-86db-b5edebeb5e0b.png)

## Plot Score Landscape

plotScoreLandscape plots the scores of the samples against two different gene signatures in a landscape for exploring the replationships of gene signatures.

```
data("scoredf_tumour_ep")
data("scoredf_tumour_mes")
plotScoreLandscape(scoredf_tumour_ep, scoredf_tumour_mes, 
                   scorenames = c('tcga-EP','tcga-MES'), isInteractive = FALSE)

```
![](https://user-images.githubusercontent.com/12887308/33762072-531c544e-dc5f-11e7-92c8-1791b897f038.png)

**For more examples, please refer to the vigntte included in the package by**
`browseVignettes('singscore')`

### Linked Article
TBA