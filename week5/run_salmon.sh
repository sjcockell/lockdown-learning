#!/bin/bash

# read input
GEO=$(cat geo_accessions.txt)

for i in ${GEO}
    do
	SRR=$(grep ${i} SraRunTable.txt | cut -d ',' -f 1)
	SRR=$(echo ${SRR} | sed 's/ /.fastq.gz /g')
	SRR=${SRR}.fastq.gz
	salmon quant -i gencode_v33_index -l A -o ${i} -r ${SRR}
    done
