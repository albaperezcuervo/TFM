# A Predictive Framework for Alzheimer´s Disease Status Based on Gut microbiome Data and Machine Learning
These scripts support the analysis of gut microbiome data for the classification of Alzheimer’s disease (AD) status. The workflow begins with data pre-processing to generate .qza artefacts from raw sequencing reads, which are then analysed in R for diversity and taxonomic patterns. Machine learning models (Random Forest, CNN, and MLPNN) are subsequently applied to classify AD status. The code reflects the analytical pipeline used in the study and facilitates reproducibility of the results

### 📑 Table of Contents
- `1) process_data_from_ena` – Scripts for downloading and processing raw microbiome data from ENA.
- `2) R_scripts` – Statistical analyses and diversity metrics calculated using R.
- `3) ANNs_and_RF` – Machine learning implementation in Python, including ANN (CNN, MLPNN) and Random Forest models.
- `data` – Input data used in the pipeline (e.g., abundance matrices, metadata).


Below is a summary of the repository structure and the role of each directory.

### 📁 `data` – Input Resources and Pre-Processed Files
This directory contains essential input files used throughout the pipeline, including scripts for downloading raw sequencing reads and `.qza` outputs derived from the pre-processing workflows.

| File name                        | Description                                                                 |
|----------------------------------|-----------------------------------------------------------------------------|
| `ena-file-PRJNA798058-WGS.sh`   | Shell script to download WGS `.fastq.gz` files from ENA via `wget`.        |
| `ena-file-PRJNA801673-16s.sh`   | Shell script to download 16S `.fastq.gz` files from ENA via `wget`.        |
| `16s.gg2.biom.qza`              | BIOM table from 16S reads, generated using the `pre_process_16s_data` script. |
| `16s.gg2.taxonomy.qza`          | Taxonomic classification of 16S data, also from `pre_process_16s_data`.    |
| `woltka_gg2_WGS.biom.qza`       | Taxonomic and functional table derived from WGS reads via `pre_process_WGS_data`. |
| `woltka_gg2_WGS.taxonomy.qza`   | Taxonomic classification table derived from WGS reads via `pre_process_WGS_data`. |
| `metadata_16s.csv`              | Metadata for 16S samples, including `sampleid`, `diagnosis` (Control or Alzheimer), and `sex` among other clinical and demographic details. |
| `metadata_wgs.csv`              | Metadata for WGS samples, with `sampleid`, `diagnosis` (Control or Alzheimer), `sex`, and additional study-related information. |

### 📁 `1) process_data_from_ena` – Pre-processing Scripts for 16S and WGS Data
This folder contains scripts for processing raw 16S and whole-genome shotgun (WGS) sequencing data. The workflows cover key steps including quality control, taxonomic classification, and output formatting for machine learning applications.

- `16s_manifest.py` – Generates a QIIME2-compatible `manifest.csv` by scanning local 16S `.fastq.gz` files and assigning sample IDs and read directions.
- `bowtie2.pl` – Aligns WGS paired-end reads using Bowtie2 and pipes the output through a formatting command to generate compressed `.sam.gz` files.
- `pre_process_16s_data.sh` – Automates the QIIME2 pipeline for 16S data: importing, denoising with DADA2, taxonomic classification using Greengenes2, and exporting `.qza` artefacts.
- `pre_process_wgs_data.sh` – Executes WGS pre-processing: quality control via MetaWRAP, alignment with Bowtie2, and classification using Woltka and QIIME2. Outputs are saved in `.qza` format.

### 📁 `2) R_scripts` – Statistical Analysis and Data Preparation in R
This directory contains R scripts for the statistical analysis and preparation of microbiome data, prior to and during downstream machine learning. Scripts are grouped according to their purpose in the workflow.
#### 🔹 Pre-analysis Scripts
These scripts handle the integration, rarefaction, and batch correction of phyloseq objects prior to diversity and classification analysis.
- `1. Merge_phyloseq_objects.R` – Combines individual phyloseq objects from 16S and WGS datasets into a unified object. Performs taxonomic cleaning and assigns study identifiers before merging.
- `2. phyloseq_rarefaction_analysis.R` – Performs rarefaction analysis to assess sequencing depth across samples. Outputs rarefied phyloseq objects to normalise read counts.
- `3. Batch_correction.R` – Applies batch effect correction using the MMUPHin package. Produces adjusted phyloseq objects and ordination plots based on UniFrac distances.
#### 🔹 Main Statistical Analysis
These scripts perform core statistical analyses, group comparisons, and data preparation for external tools and machine learning models.
- `alpha_and_beta_diversity.R` – Computes alpha diversity indices (Chao1, Shannon, Simpson) and beta diversity distances (Bray–Curtis, Jaccard, UniFrac), including PERMANOVA tests and ordination visualisations.
- `taxa_abundance_by_diagnosis_and_sex.R` – Assesses differences in taxonomic composition between AD and Control groups, and between sexes. Includes non-parametric testing and visualisation of significant taxa.
- `lefse_and_nn_input.R` – Prepares genus-level abundance tables and metadata for use with LEfSe and as input to neural network models.
- `Lefse_analysis.R` – Processes and visualises the results of LEfSe analysis, highlighting taxa with significant differential abundance via LDA scores.

### 📁 `3) ANNs_and_RF` – Machine Learning Modelling and Evaluation
This folder contains all Python scripts used for machine learning model development, optimisation, and evaluation. The process begins with data formatting and input generation, followed by initial training and tuning of three models (CNN, MLPNN, and Random Forest). 

- `generate_inputs.py` – Prepares input tensors for machine learning. Transforms raw taxonomic abundance and metadata into structured `.h5` and `.npy` formats suitable for CNN, MLP, and RF models. Includes cumulative abundance encoding across taxonomic levels.

- `models_cnn_mlp_rf.py` – Implements CNN, MLPNN, and Random Forest models using raw abundance data without data augmentation. Includes genetic feature selection, hyperparameter tuning via `GridSearchCV`, and evaluation using accuracy, recall, and classification reports.

- `Optimal_number_of_simulated_samples.py` – Systematically evaluates different oversampling sizes using SMOTE to determine the optimal number of synthetic samples per class. Stores performance results for visual inspection and comparison.

- `models_cnn_mlp_rf_smote.py` – Retrains all models using the optimal number of synthetic samples identified above. Applies the same architecture and tuning strategies to assess improvements due to SMOTE-based augmentation.

### 📌 Reproducibility and Requirements
This repository assumes access to QIIME2, Python 3.8+, and R (≥ 4.0) with packages specified in the corresponding script headers. It is recommended to execute the scripts sequentially following the directory order.
