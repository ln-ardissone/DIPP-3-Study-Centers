# DIPP-3-Study-Centers
16S rRNA Microbiota Analysis of the DIPP cohort

## Introduction

The Diabetes Prediction and Prevention (DIPP) study (http://dipp.utu.fi/index.php?mid=2&language=en)....

## Explanation of Contents

### Inputs

This directory contains intermediate files that are called in ```D3C_Statistical_Analysis_Microbiome.Rmd``` which take a bit of time to generate. It is not necessary to call these files and they can be regenerated from the .Rmd file by switching the corresponding chunk options to ```eval=TRUE```. They are provided to speed up execution of the .Rmd file. A brief description of each file is provided:

- **WUF.dist.all.Rdata**
- **WUF.dist.oulu.Rdata**
- **WUF.dist.tampere.Rdata**
- **WUF.dist.turku.Rdata**
- **dipp_NB.Rdata**
- **dipp_UUF.Rdata**
- **dipp_WUF.Rdata**
- **gg_13_8.97_otus.tree**

### D3C2.0_TechReps.pdf/.Rmd

The Rmarkdown and the resulting .pdf when *knitr* was executed. This document explains the processing and analysis of technical replicates included in the sequencing experimental design. Breifly, we observed that the greatest source of technical variation was across different extractions. Furthermore, we applied a method in which OTUs with a high ratio of technical to biological variation were removed from further analyses. One would not be able to distinguish whether the variation in these OTUs were due to technical processing or true biological differences and were therefore excluding from further consdieration.

### D3C_Statistical_Analysis_Microbiome.pdf/.Rmd

The Rmarkdown and the resulting pdf when *knitr* was executed for the microbiota analyses. Code for analytical results leading to the conclusions presented in the manuscript are provided in these documents as well as supplemental analyses. 

### D3C_function_5.0.r

Contains functions called and used in ```D3C_Statistical_Analysis_Microbiome.Rmd```.

### README.md

This document.

### dipp.Rdata

The input data for ```D3C2.0_TechReps.Rmd```. It contains the complete OTU table which resulted from sequence classifaction, prior to applying any type of OTU filter methods.

### dipp.sam3.genfilt.Rdata

The input data for ```3C_Statistical_Analysis_Microbiome.Rmd```. It contains the OTU table which resulted from ```D3C2.0_TechReps.Rmd``` and has OTUs with high technical variation removed.

## R Session Info

All RMarkdown files were executed from Rstudio with the following session information:

```R version 3.3.2 (2016-10-31)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: macOS Sierra 10.12.5

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gdata_2.18.0    vegan_2.4-2     lattice_0.20-35 permute_0.9-4   geepack_1.2-1   gridExtra_2.2.1 reshape2_1.4.2  plyr_1.8.4     
 [9] lubridate_1.6.0 phyloseq_1.19.1 ggthemr_1.0.2   ggplot2_2.2.1   knitr_1.15.1   

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.10        highr_0.6           XVector_0.14.1      iterators_1.0.8     tools_3.3.2         zlibbioc_1.20.0    
 [7] digest_0.6.12       jsonlite_1.3        evaluate_0.10       tibble_1.3.0        gtable_0.2.0        nlme_3.1-131       
[13] rhdf5_2.18.0        mgcv_1.8-17         Matrix_1.2-8        foreach_1.4.3       igraph_1.0.1        yaml_2.1.14        
[19] parallel_3.3.2      stringr_1.2.0       cluster_2.0.6       gtools_3.5.0        Biostrings_2.42.1   S4Vectors_0.12.2   
[25] IRanges_2.8.2       rprojroot_1.2       stats4_3.3.2        ade4_1.7-6          multtest_2.30.0     grid_3.3.2         
[31] Biobase_2.34.0      data.table_1.10.4   survival_2.41-3     rmarkdown_1.4       magrittr_1.5        backports_1.0.5    
[37] htmltools_0.3.5     MASS_7.3-45         scales_0.4.1        codetools_0.2-15    BiocGenerics_0.20.0 splines_3.3.2      
[43] biomformat_1.2.0    colorspace_1.3-2    ape_4.1             labeling_0.3        stringi_1.1.3       lazyeval_0.2.0     
[49] munsell_0.4.3      ```
