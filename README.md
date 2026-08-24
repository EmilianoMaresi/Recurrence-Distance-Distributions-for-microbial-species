# Description

Bioinformatics pipeline that simulates Amplified Fragment Length Polymorphism (AFLP) analysis in silico.

Given a dataset of bacterial genomic FASTAs (e.g., from NCBI), the tool identifies the optimal k-mer pairs and computes the distribution of fragment lengths occurring between them. 
It then evaluates the discriminative power of these fragment distributions through a dual downstream pipeline:

- **Machine Learning Classification:** Benchmarks species separation using 3 classifiers (MLP, Random Forest, and XGBoost) with cross-validation (reporting accuracy and F1-score).

- **Taxonomic Reconstruction:** Performs average-linkage hierarchical clustering to reconstruct bacterial taxonomies, evaluated via ARI, NMI, and silhouette scores.

Wrapped in a custom Docker image for complete reproducibility and ease of deployment.

---

## Installation

### STEP 1: Install Docker

You need **Docker** installed and running on your system.

- **Ubuntu / Debian:**

```bash
sudo apt update
sudo apt install -y docker.io
sudo systemctl start docker
sudo systemctl enable docker
```

Check Docker installation:
```bash
docker --version
docker info
```

### STEP 2: Build the Docker image
Build the Docker image containing all dependencies:
```bash
./install.sh
```
- If the image already exists, the script will ask if you want to rebuild it.
- The Docker image includes all Python packages from environment.yml and makes the pipeline portable and reproducible.

---

## Uninstall
Uninstall the docker image with:
```bash
./uninstall.sh
```
---
## Usage
After installation, all commands are run using `launch_rdd.sh`

### Run Example
```bash
./example.sh
```
---
## Results folder
- Results are saved in `results/`:

```
results/
 ├── aflp/                                  # AFLP distributions for each k-mer pair
 ├── classifiers/                           # Classification performance reports
 ├── clusters/                              # Clustering performance reports
 ├── kmer_pairs.csv                         # optimal kmer-pairs ranked by score
 ├── aflp_validity.csv                      # Valid AFLP distributions
 ├── kmer_k6_multiplicities.csv             # Matrix of kmer occurrences
 ├── kmer_k6_multiplicities_normalized.csv  # Normalized matrix of kmer occurrences
 └── genomes_length.csv                     # Meta information on the fasta genomes ```

### Help

Launcher help:
```bash
./launch_rdd.sh
```

To see the help for the individual tools:
```bash
./launch_rdd.sh download -h
./launch_rdd.sh run -h
```
---
## Scripts Descriptions

`download_ncbi_dataset.py`
Downloads the genomes listed in a .tsv file (usually the NCBI summary file obtained from the NCBI Datasets search page https://www.ncbi.nlm.nih.gov/datasets/genome/. 
- e.g. search for bacterial genomes: https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2

Example:
```bash
./launch_rdd.sh download -i datasets_lists/dataset_example.tsv
```
Output:

- Genomes saved in `data/dataset_example/`, divided by species
- The `zip/` subfolder contains raw downloaded archives (can be deleted after extraction)

`main.py`

Performs:

- AFLP distribution feature extraction
- Classification using Random Forest, MLP, and XGBoost
- Clustering using hierarchical clustering with average-linkage
 
Example:
```bash
./launch_rdd.sh run -i data/dataset_example -s 6 -e 6 -top 2 -cpu 2
```

Results are saved in `results/` folder containing:
- `aflp/` — AFLP distributions
- `classifiers/` — classification reports
- `clustering/` — clustering reports

Classifier output files:
- `classifiers_results.csv` (default) — classifiers performances using AFLP fragment distributions (identified as kmer pairs, e.g. ACTGCC-TTCTCC)
- `classifiers_results_genus.csv` (only with `--genus` option) — classifiers results when using genera as labels
- `classifiers_summary.csv` — Summary of the classifiers results for the three classifiers used to validate

Clustering output files:
- `clustering_results.csv` — clustering performances using AFLP fragment distributions (identified as kmer pairs, e.g. ACTGCC-TTCTCC)
- `clustering_summary.csv` — Summary resuming the performance of the MLP, Random Forest and XGBoost classifiers

---
## Folder Structure Overview
```
methods/
 ├── classifiers/           # Definition of the classifiers and their methods
 │    ├── mlp.py            # MultiLayer Perceptron Classifier (MLP)
 │    ├── random_forest.py  # Random Forest Classifier
 │    ├── xgboost_model.py  # XGBoost Classifier
 │    └── model_params.py   # Model parameters for MLP, Random Forest and XGBoost classifiers
 ├── kmers_distributions.py # computation of optimal kmer-pairs and AFLP fragment distribution extraction
 ├── run_classifiers.py     # classifiers pipeline
 ├── run_clustering.py      # clustering pipeline and methods
 └── pipeline.py            # general pipeline
```

- New classifiers can be added under `methods/classifiers/` as Python classes and imported in `run_classifiers.py`.

- New parameter presets can be added as dictionaries in `model_params.py`.
---

# Install without Docker using conda environments (Ubuntu only)

### STEP 1: Install Miniconda (or Anaconda)

Install Miniconda on Linux

Download the latest Miniconda installer (for Linux x86_64)
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
Run the installer
```
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow on-screen instructions (press ENTER to accept defaults)
Then activate conda
```
source ~/.bashrc
```

Check installation
```
conda --version
```

### STEP 2: Create the Conda Environment
```
conda env create -f environment.yml
```

Then activate it:
```
conda activate RDD_env
```

### STEP 3: Run the Example
```
./example.sh
```

### Scripts Descriptions

- `download_ncbi_dataset.py`
(Run python`download_ncbi_dataset.py -h` for help.)

Downloads the genomes listed in a .tsv file — usually the summary file obtained from the NCBI Datasets
 search page. (e.g. for bacteria --> https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2).

example:
```
python download_ncbi_dataset.py -i ./datasets_lists/dataset_example.tsv
```

Output:
- Genomes are saved in `data/dataset_example/`, divided by species.
- The zip subfolder contains the raw downloaded archives (can be deleted after extraction).
Note: NCBI downloads can fail intermittently. The script automatically retries until all genomes are successfully downloaded.

- main.py
(Run `python main.py -h` for help.)

Performs:
- AFLP distribution feature extraction
- Classification using Random Forest, MLP, and XGBoost

Results are saved in the results/ folder, containing:
- aflp/: AFLP distributions for each k-mer pair
- classifiers/: classification performance reports (Accuracy, F1-score)

Classifier output files:
- `classifiers_results.csv` — using raw AFLP occurrences
- `classifiers_results_binarized.csv` — using binarized AFLP presence/absence
- `classifiers_results_genus.csv` — same as (1) but using genus labels (--genus option)
- `classifiers_results_binarized_genus.csv` — same as (2) but using genus labels

Folder Structure Overview

```
methods/
 ├── classifiers/
 │    ├── mlp.py
 │    ├── random_forest.py
 │    ├── xgboost_model.py
 │    └── model_params.py
 └── (other utility scripts)
```

- New classifiers can be added under methods/classifiers/ as Python classes, then imported in run_classifiers.py.
- New parameter presets can be added as dictionaries in model_params.py.


