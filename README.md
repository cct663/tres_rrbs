# Summary

This repository contains all the code needed to reproduce analyses, figures, and the manuscript for 'Corticosterone exposure causes long-term changes in DNA methylation, physiology, and breeding decisions in a wild bird.' Metadata for the paper and data for some analyses are also included in the repository. The sequencing data used for methylation analyses and for genome assembly can be accessed from NCBI and details are included in the manuscript. 

All of the code for analyses is included in the '1_main_rrbs_script.R' file with extensive annotation throughout the file. We provide a brief description here of the key files, but please see the script comments for more details.

The file '2_rrbs_manuscript.rmd' reproduces the full manuscript and supplementary materials along with figures and tables derived from the main code analyses.

A version of this analysis is permanently archived at Zenodo.

[![DOI](https://zenodo.org/badge/315731095.svg)](https://zenodo.org/badge/latestdoi/315731095)


Brief description of folder and file contents:

- 1_Raw_Data folder. Includes data text files.

*cross_year_cort.txt*: Data for the cross year phenotypic effects of corticosterone treatment on year n+1 corticosterone and breeding timing. Each row is one female and each column is one measurement.

*gene.key.txt*: Key of genes in our annotated genome. Used in matching gene locations to names.

*genes_TRES_new.txt*: Key of gene locations in our annotated genome. Used in matching CpGs to genes.

*rrbs_sample_metadata.txt*: Metadata for each RRBS sample. Used for correlation between methylation and corticosterone and for within and between year effects of treatment on corticosterone. Each row is one sample (females may have multiple rows) and each column is one measurement.

- 2_r_scripts. Scripts used in analysis and manuscript production along with some output of figures produced by the script to be embedded in the rendered manuscript.

*1_main_rrbs_script.R*: All analyses in R are included here.

*2_rrbs_manuscript.Rmd*: Produces the manuscript

*references.bib*: List of reference information used in the manuscript.

*figures*: various figure files produced by the analysis script

- 5_temporary files. These can be ignored. Intermediate files produced in parts of the analysis script to make loading back in faster.

- 8_output_for_DAVID. These are output text files formatted to be used with DAVID online software based on the comparisons described in the text.