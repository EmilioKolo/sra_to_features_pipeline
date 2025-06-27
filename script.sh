#!/bin/bash

cat sra_table_selected.txt | parallel -j2 'docker run -v ../output/:/content/data/output -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) features-pipeline {} /content/data/output'
