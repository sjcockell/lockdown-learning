#!/bin/bash

# activate the right conda env
# conda activate msa

# get input data (uniprot ids)
SRA=$(tail -n +2 SraRunTable.txt | cut -d ',' -f 1)

# This is a loop for downloading the data

for i in ${SRA}
    do
	prefetch ${i}
	if [ -f ${i}.fastq.gz ]
            then
                echo "${i} already finished"
        else
                echo "(o) Convert SRA entry: ${i}" 
                # downloading SRA entry
                fastq-dump --gzip --defline-qual '+' ${i}/${i}.sra
                echo "(o) Done convert  ${i}"
        fi
    done


