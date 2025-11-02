# AFLP Classification Pipeline

This repository provides a pipeline for extracting **AFLP distributions** from genome datasets and classifying them using **Random Forest**, **MLP**, and **XGBoost** models.

---

## STEP 1: Install Miniconda (or Anaconda)

### Install Miniconda on Linux

```bash
# Download the latest Miniconda installer (for Linux x86_64)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer
bash Miniconda3-latest-Linux-x86_64.sh

# Follow on-screen instructions (press ENTER to accept defaults)
# Then activate conda
source ~/.bashrc

# Check installation
conda --version

## STEP 2: Create the Conda Environment
conda env create -f environment.yml

Then activate it:
conda activate RDD_env

## STEP 3: Run the Example
./example.sh

===============
Script Descriptions

- download_ncbi_dataset.py
(Run python download_ncbi_dataset.py -h for help.)

Downloads the genomes listed in a .tsv file — usually the summary file obtained from the NCBI Datasets
 search page. (e.g. for bacteria --> https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2).

example: python download_ncbi_dataset.py -i ./datasets_lists/dataset_example.tsv

Output:
- Genomes are saved in data/dataset_example/, divided by species.
- The zip subfolder contains the raw downloaded archives (can be deleted after extraction).
Note: NCBI downloads can fail intermittently. The script automatically retries until all genomes are successfully downloaded.

- main.py
(Run python main.py -h for help.)

Performs:
- AFLP distribution feature extraction
- Classification using Random Forest, MLP, and XGBoost

Results are saved in the results/ folder, containing:
- aflp/: AFLP distributions for each k-mer pair
- classifiers/: classification performance reports (Accuracy, F1-score)

Classifier output files:
- classifiers_results.csv — using raw AFLP occurrences
- classifiers_results_binarized.csv — using binarized AFLP presence/absence
- classifiers_results_genus.csv — same as (1) but using genus labels (--genus option)
- classifiers_results_binarized_genus.csv — same as (2) but using genus labels

Folder Structure Overview

methods/
 ├── classifiers/
 │    ├── mlp.py
 │    ├── random_forest.py
 │    ├── xgboost_model.py
 │    └── model_params.py
 └── (other utility scripts)

- New classifiers can be added under methods/classifiers/ as Python classes, then imported in run_classifiers.py.
- New parameter presets can be added as dictionaries in model_params.py.


