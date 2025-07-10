#!/bin/bash

# Build the docker container
docker build -t features-pipeline .
# Remove Windows carriage returns from the sra table
sed -i 's/\r//' sra_table_selected.txt
# Run code with parallel
cat sra_table_selected.txt | parallel 'docker run -v ../output/:/content/data/output -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) features-pipeline {} /content/data/output'
