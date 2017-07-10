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
