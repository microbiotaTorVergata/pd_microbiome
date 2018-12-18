#!/bin/bash
# work dir
PATH_WORK="/../../work_dir"
# path to greengenes
PATH_DB_SEQ="/../../gg_13_8_otus/rep_set/97_otus.fasta"
PATH_DB_TAX="/../../gg_13_8_otus/taxonomy/97_otu_taxonomy.fasta"
PATH_DB_ALN="/../../gg_13_8_otus/rep_set_aligned/97_otus.fasta"
# path to bin (usearch61)
PATH_BIN="/../../bin"
export PATH="$PATH_BIN:$PATH"                      		        

# temp dir
mkdir -p tmpdir 

source activate qiime1

export QIIME_HOME=$PATH_WORK                     
export QIIME_CONFIG_FP=$PATH_WORK/qiime_config   

# mapping file
map="mapping.txt"

mkdir -p Phylogeny

# take only otu not discarded by Bokulich filter
python $PATH_WORK/Script_Python/filter_otu_map.py -b Otu_Tables_Closed_NoRep/otu_table_filtered.txt -i Phylogeny/combined_seqs_otus_norep.txt -o Phylogeny/combined_seqs_otus_filtered.txt

# phylogeny
pick_rep_set.py -i Phylogeny/combined_seqs_otus_filtered.txt -f combined_seqs.fna  -o Phylogeny/rep_set.fna -r $PATH_DB_SEQ
align_seqs.py -i Phylogeny/rep_set.fna -m pynast -o Phylogeny -a uclust -t $PATH_DB_ALN
filter_alignment.py -i Phylogeny/rep_set_aligned.fasta -o Phylogeny/
make_phylogeny.py -i Phylogeny/rep_set_aligned_pfiltered.fasta -o Phylogeny/rep_set.tre --root_method midpoint 


source deactivate
