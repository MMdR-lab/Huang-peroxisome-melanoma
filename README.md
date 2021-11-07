# Mapping phenotype switching in the process of melanoma therapy with single-cell data
This resource provides the code developed in the study of Fan et al. It reproduces the key results of the study and can be applied to other single-cell cohorts to explore the characteristics of phenotype switching in melanoma.

## Requirements
R (tested in R version 4.0.3 (2020-10-10)).
R libraries: Seurat, dplyr, patchwork, AUCell, GSEABase
## Raw Data
The data is provided in the RAW DATA folder (ImmRes_Rfiles.zip).

In the Portal you will also find the processed single-cell gene expression along with interactive views.

Quick start
To reproduce the results reported in Jerby-Arnon et al. download ImmRes_Rfiles.zip from the Single Cell Portal. Unzip the file and move the resulting Data directory to the ImmuneResistance directory.

In R go to the Code directory and run master.code() which is provided in ImmRes_master.R. The master.code() will walk you through the different stages of the study, divided into six main modules:

(1-2) First, analyzing the single-cell data to generate various gene signatures that characterize different cell subtypes and immune resistant cell states. For more information see Mapping immune resistance in melanoma.

(3-5) Next, analyzing independent cohorts obtained from bulk melanoma tumors to explore and test the immune resistance program. For more information see Predicting immunotherapy resistance.

(6) Lastly, performing a pan-cancer analysis to identify drugs that could repress the immune resistance program in cancer cells. For more information see Repressing the immune resistance program.

General notes
The code provided in ImmRes_master.R reproduces the key results of the study. It also generates the study figures and table in the Output directory. The code follows the analyses that were performed in the study in a sequential manner.

As the results are already provided in the Results directory, it is possible to run only some parts of the code and focus on specific analyses, or apply the approach to other datasets.

Citation
Jerby-Arnon L et al. Single-cell RNA-seq of melanoma ecosystems reveals sources of T cell exclusion linked to immunotherapy clinical outcomes.
