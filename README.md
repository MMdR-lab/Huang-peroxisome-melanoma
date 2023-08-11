# Mapping melanoma persister cell states in the process of MAPK-targeted therapy with single-cell data
This resource provides the code developed in the study of Huang et al. It reproduces the key results of the study and can be applied to other single-cell cohorts to explore the characteristics of phenotype switching in melanoma.

## Requirements
R (tested in R version 4.0.3 (2020-10-10)).
R libraries: Seurat, dplyr, patchwork, AUCell, GSEABase
## Raw Data
The data is provided in the RAWDATA folder (GSE116237_scRNAseq_expressionMatrix.txt).

In the outputs folder you will also find the processed single-cell gene expression (melanoma.rds).

## Quick start
To reproduce the results reported in Huang et al., download GSE116237_scRNAseq_expressionMatrix.txt from the RAWDATA folder. In R run the Code file single cell v1.R. It will walk you through the different stages of the study.

First, we analyzed the single-cell data to generate various gene signatures that characterize different cell states, including pigmentation, SMC-like, NCSC-like and invasion. We referred to the article Toward Minimal Residual Disease-Directed Therapy in Melanoma for the gene signatures. For more information, see the file signature.csv.

Next, to measure the activity of the gene signatures in each cell, we used the AUCell algorithm (Aibar et al., 2017). The activity of each of the signatures was visualized by i) projecting all 674 cells into a two-dimensional space using t-distributed stochastic neighbor embedding (t-SNE, perplexity = 30, initial_dims = 10) based on the expression of all genes in the signatures and ii) coloring cells according to their binary AUCell score.

Lastly, we identified cell states (pigmentation, SMC-like, NCSC-like and invasion) in melanoma patient-derived xenografts with different therapy stages according to AUCell results and visualized the expression of UGCG, GBA, PEX3 and AGPS with violin plot.

## General notes
The code provided in single cell v1.R reproduces the key results of the study. It also generates the study figures in the images directory. The code refers to the analyses that were performed in the article Toward Minimal Residual Disease-Directed Therapy in Melanoma.

As the results are already provided in the Results directory, it is possible to run only some parts of the code and focus on specific analyses, or apply the approach to other datasets.

## Citation
Huang F et al. Peroxisome disruption alters lipid metabolism and potentiates anti-tumor response with MAPK-targeted therapy in melanoma.
