#!/bin/bash

echo "Example using genome list: dataset_example.tsv"

echo "Download genomes from NCBI"
#python download_ncbi_dataset.py -i ./datasets_lists/dataset_example.tsv

echo "Running main script: Extract AFLP fragment distributions and species classification"
echo "Example with kmer length = 6 for both start and end kmers; select only the top 5 kmer pairs ; 2 cpus (same as default cpus: 2)"
python main.py -i data/dataset_example/ -s 6 -e 6 -top 5 -cpu 10
