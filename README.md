# Bayesian inference of transcription regulation from time-resolved scRNA-seq data
This repository contains code produced for modelling and data analysis in the manuscript **Global transcription regulation revealed from dynamical correlations in time-resolved single-cell RNA-sequencing**.

You can follow the guidance below to obtain the data and reproduce main results of the manuscript by running the script `wrapper.jl`. The subsections of the README file correspond to the different sections of code in the `wrapper.jl` file. It is suggested that the folders `data` and `scripts` and the file `wrapper.jl` of this repository are downloaded and placed in your local directory. Additionally, it is recommended that the data files contained in the Zenodo repository <https://doi.org/10.5281/zenodo.10724157> are downloaded and placed in the same `data` folder.

### 1. Data
The original scEU-seq data that we analyse were generated and published by the authors of <https://www.science.org/doi/10.1126/science.aax3072>. The data can be downloaded as CSV files from the GEO repository with accession number [GSE128365](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128365). It is recommended that the downloaded data are moved to the folder `data/raw_data/`. 

|                             Name                               |        Description        |           
| ---------------------------------------------------------------|---------------------------|
| GSE128365_SupplementaryData_RPE1_unlabeled_unsplicedUMI.csv.gz | unlabelled unspliced UMIs | 
| GSE128365_SupplementaryData_RPE1_unlabeled_splicedUMI.csv.gz   |  unlabelled spliced UMIs  |    
| GSE128365_SupplementaryData_RPE1_labeled_unsplicedUMI.csv.gz   |  labelled unspliced UMIs  |  
| GSE128365_SupplementaryData_RPE1_labeled_splicedUMI.csv.gz     |   labelled spliced UMIs   |  
| GSE128365_SupplementaryData_RPE1_metadata.csv.gz               |   cell-specific metadata  |


#### Load required packages and raw data
The script first loads all the required packages for all downstream tasks. If you are using Julia for the first time, you can add these packages in Julia by typing `] add "package_name"` in the command line. 

Next, we load all the raw data including: the unlabelled/labelled-unspliced/spliced mRNA count matrices `uu_data`,`us_data`,`lu_data`,`ls_data`, the cell-specific coordinates along cell cycle progression `theta`, the single-cell fluorescence levels of the cell cycle reporters Cdt1-RFP and Geminin-GFP `rfp`,`gfp`, as well as the pulse/chase labelling conditions of each cell `experiment`. We also obtain the total mRNA count matrix `total_data` by summing up the 4 mRNA count matrices.

#### Pre-process data
The cells are clustered into 5 cell cycle stages (age groups), encoded by the vector `age`. In addition, pulse-treated and chase-treated cells are separated into `pulse_idx` and `chase_idx`. Next, all useful empirical frequencies/distributions of cells across cell cycle stages and labelling conditions are computed. 

We next proceed with gene filtering and we load the list of selected genes `genes` that we use in all downstream analyses. 

#### Estimate cell-specific capture efficiencies
We next estimate cell-specific capture efficiencies from the data, using each cell's total count and cell cycle stage. 

#### Compute summary statistics
After pre-processing the data, we proceed with computing all summary statistics that we will need for inference. We obtain the matrices `pulse_mean`,`pulse_ff`,`chase_mean`,`chase_ff`,`ratio_data`,`mean_corr_data`,`corr_mean_data`, as well as the estimates for the standard error of these statistics, `pulse_mean_se`,`pulse_ff_se`,`chase_mean_se`,`chase_ff_se`,`ratio_data_se`,`mean_corr_data_se`,`corr_mean_data_se`, computed using bootstrapping.

### 2. Modelling and ABC simulations
We define the model and we run example realisations of the 5 gene models (constant scaling, constant non-scaling, burst frequency, burst size and decay rate). This produces plots of the summary statistics and the characteristic kinetic rate of each simulated gene.

We next set-up the ABC framework. We sample $M = 5 \cdot 10^6$ parameter sets from our prior distribution and numerically simulate our 5 models. We use our estimated capture efficiencies to add technical noise to the output of our models and then we generate summary statistics. We determine a model index `m` $=1, \dots ,5$, corresponding to the $5$ models we want to simulate, the number of simulations (sampled parameter sets) `n_trials` that we want to run, as well as a `submit` index that is used to distinguish multiple parellel simulation runs (e.g. in a high-performance computing cluster). As this is the most computationally intensive step, it can be skipped and our simulation outputs can be obtained in the folder `simulations/` that was downloaded from the Zenodo repository. 

### 3. Parameter inference and model selection

#### Applying the ABC rejection scheme
The complete output of the ABC simulations can be obtained from the Zenodo repository and added in the `data` folder. Next, we proceed with comparing the summary statistics generated by the model with the summary statistics of the data, through the ABC rejection scheme. Computing errors between model and data gives us data files of very large size, so we process them to get a compressed format. Next, for each of the 5 models and each gene, we collect all the particles that achieved an error smaller than our acceptance threshold. The compressed errors and the accepted particles per gene and per model can be found in the `errors/` and `posteriors/` folders that were downloaded from the Zenodo repository. 

#### Computing model probabilities
We proceed with model selection, which involves 2 steps: First, for each gene, we assign a probability to the constant hypothesis and the non-constant hypothesis. Second, we compute the conditional model probabilites over constant models and over non-constant models separately. For all these probabilites, we also compute $95$% confidence intervals by bootstrapping the accepted particles of each gene. Next, we perform model selection and classify genes to the different models. 

### 4. Transcription kinetics
Next, we can inspect the inferred posterior distributions of kinetic rate parameters. We can pick a model index `m` $=1, \dots ,5$ and a gene index `g` and inspect the marginal posterior densities of its parameters or pairs of parameters. We also plot the distributions of point estimates of the kinetic rates across genes and models. Additionally, we compare the inferred gene-specific labelling efficiencies with gene lengths.  

### 5. Statistics recovery and noise decomposition
We can use the inferred parameters and the model to recover statistics of gene expression: we obtain model predictions for the mean and noise levels of gene expression in the absence of technical (downsampling) noise. All the statistics computed in this step can be found in the `recovered_statistics/` folder. 

### 6. Constant genes: Scaling properties of gene expression with cell size
We analyse the 2 classes of constant genes (scaling and non-scaling) to reveal scaling properties of gene expression with cell size. We compare the 2 groups of genes with respect to their mean expression fold change, mean expression rate of change, levels of cell cycle variation and mRNA degradation rates.

### 7. Non-constant genes: Properties of cell cycle-dependent gene expression
We cluster the 3 classes of non-constant genes with respect to their cell cycle-dependent kinetic rates to reveal a general pattern of transcription regulation along the cell cycle. 
