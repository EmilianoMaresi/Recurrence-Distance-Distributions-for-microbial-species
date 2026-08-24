# AFLP Classification Pipeline

Bioinformatics pipeline designed to evaluate bacterial species discrimination using an Amplified Fragment Length Polymorphism (AFLP)-inspired approach.
Given a dataset of bacterial genomes (e.g., from NCBI), the tool identifies the optimal k-mer pairs to compute fragment length distributions. It then validates the discriminative power of these distributions through a dual pipeline featuring machine learning classification (MLP, Random Forest, XGBoost) and hierarchical clustering for taxonomy reconstruction.

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

### Download Example Dataset
```bash
./launch_rdd.sh download -i datasets_lists/dataset_example.tsv
```
- Downloads genomes listed in `dataset_example.tsv` into `data/dataset_example/`.
- The script automatically retries failed downloads until all genomes are obtained.

### Run Example
```bash
./example.sh
```
- Extracts AFLP distribution features and performs classification.
- Results are saved in `results/`:

```
results/
 ├── aflp/           # AFLP distributions for each k-mer pair
 └── classifiers/    # Classification performance reports
```

### Help
To see the full list of arguments for each script:
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

Example:
```bash
./launch_rdd.sh run -i data/dataset_example -s 6 -e 6 -top 2 -cpu 2
```

Results:

- Saved in `results/` folder containing:
- `aflp/` — AFLP distributions
- `classifiers/` — classification performance reports

Classifier output files:
- `classifiers_results.csv` — raw AFLP occurrences
- `classifiers_results_binarized.csv` — binarized AFLP presence/absence
- `classifiers_results_genus.csv` — using genus labels (--genus option)
- `classifiers_results_binarized_genus.csv` — binarized with genus labels
---
## Folder Structure Overview
```
methods/
 ├── classifiers/
 │    ├── mlp.py
 │    ├── random_forest.py
 │    ├── xgboost_model.py
 │    └── model_params.py
 └── (other utility scripts)
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


