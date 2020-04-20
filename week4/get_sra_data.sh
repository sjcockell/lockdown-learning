#!/bin/bash

# activate the right conda env
# conda activate msa

# get input data (uniprot ids)
VAR=$(tail -n +2 SraRunTable.txt | cut -d ',' -f 1)

# This is a loop for downloading the data

for i in ${VAR}
    do
        if [ -f ${i}.fastq.gz ]
	    then
		echo "${i} already downloaded"
	else
		echo "(o) Downloading SRA entry: ${i}" 
		# downloading SRA entry
        	fastq-dump --gzip --defline-qual '+' ${i}
		echo "(o) Done downloading ${i}"
	fi
    done


