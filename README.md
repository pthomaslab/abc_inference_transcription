# Bayesian inference of transcription regulation from time-resolved scRNA-seq data
This repository contains code produced for modelling and data analysis in the manuscript **Global transcription regulation revealed from dynamical correlations in time-resolved single-cell RNA-sequencing**.

You can follow the guidance below to obtain the data and reproduce the results of the manuscript:

## Data
The original scEU-seq data that we analyse were generated and published by the authors of <https://www.science.org/doi/10.1126/science.aax3072>. You can downloaded the following data as CSV files from the GEO repository with accession number [GSE128365](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128365):

|                             Name                               |        Description        |           
| ---------------------------------------------------------------|:--------------------------:
| GSE128365_SupplementaryData_RPE1_unlabeled_unsplicedUMI.csv.gz | unlabelled unspliced UMIs | 
| GSE128365_SupplementaryData_RPE1_unlabeled_splicedUMI.csv.gz   |  unlabelled spliced UMIs  |    
| GSE128365_SupplementaryData_RPE1_labeled_unsplicedUMI.csv.gz   |  labelled unspliced UMIs  |  
| GSE128365_SupplementaryData_RPE1_labeled_splicedUMI.csv.gz     |   labelled spliced UMIs   |  
| GSE128365_SupplementaryData_RPE1_metadata.csv.gz               |   cell-specific metadata  |


### Load required packages and data
After you have obtained the CSV files, you can load all the data for pre-processing in Julia by running the `load_data.jl` script. This script also loads all the required packages for the following downstream tasks.



