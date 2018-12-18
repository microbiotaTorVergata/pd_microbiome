#!/bin/bash
# work dir
PATH_WORK="/../../work_dir"

source activate picrust

mkdir -p Picrust/

normalize_by_copy_number.py \
  -i Picrust/boh.biom \
  -o Picrust/normalized_otus_boh.biom

predict_metagenomes.py \
        -i Picrust/normalized_otus.biom \
        -o Picrust/metagenome_predictions.biom 

predict_metagenomes.py \
        -f \
        -i Picrust/normalized_otus.biom \
        -o Picrust/metagenome_predictions.txt

categorize_by_function.py \
	-i Picrust/metagenome_predictions.biom \
	-c KEGG_Pathways \
	-l 3 \
	-o Picrust/predicted_metagenomes.L3.biom

categorize_by_function.py \
	-f \
	-i Picrust/metagenome_predictions.biom \
	-c KEGG_Pathways \
	-l 3 \
	-o Picrust/predicted_metagenomes.L3.txt


source deactivate
