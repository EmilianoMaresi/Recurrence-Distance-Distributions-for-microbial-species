#!/bin/bash
set -e

echo "[[[ Example using genome list: dataset_example.tsv ]]]"

./launch_rdd.sh download -i datasets_lists/dataset_example.tsv

./launch_rdd.sh run -i data/dataset_example -s 4 -e 4 -top 2 -cpu 2

