#!/bin/bash

# activate the right conda env
conda activate msa

# get input data (uniprot ids)
VAR=$(cat uniprot_ids.txt)

URL="https://www.uniprot.org/uniprot/"

# This is a loop for downloading the data

for i in ${VAR}
    do
        echo "(o) Downloading Uniprot entry: ${i}" 
	# downloading uniprot entry
        wget ${URL}${i}.fasta
	echo "(o) Done downloading ${i}"
	# construct muscle input
	cat ${i}.fasta >> muscle_input.fasta
	rm ${i}.fasta
    done

# run muscle to build MSA
muscle -clw -in muscle_input.fasta -out muscle_output.aln

