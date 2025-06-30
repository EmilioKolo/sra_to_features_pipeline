#!/bin/bash

# Remove Windows carriage returns
sed -i 's/\r//' sra_table_test.txt
# Run code with parallel
cat sra_table_selected.txt | parallel -j2 'docker run -v ../output/:/content/data/output -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) features-pipeline {} /content/data/output'
