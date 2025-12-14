#!/bin/bash
set -e

echo "[[[ Example using genome list: dataset_example.tsv ]]]"

./launch_rdd.sh download -i datasets_lists/dataset_example.tsv

./launch_rdd.sh run -i data/dataset_example -s 6 -e 6 -top 2 -cpu 2

