#!/bin/bash
# work dir
PATH_WORK="/../../work_dir"
# path to greengenes
PATH_DB_SEQ="/../../gg_13_8_otus/rep_set/97_otus.fasta"
PATH_DB_TAX="/../../gg_13_8_otus/taxonomy/97_otu_taxonomy.fasta"
# path to bin (usearch61)
PATH_BIN="/../../bin"
export PATH="$PATH_BIN:$PATH" 

# temp dir
mkdir -p tmpdir 

source activate qiime1

export QIIME_HOME=$PATH_WORK                     
export QIIME_CONFIG_FP=$PATH_WORK/qiime_config   


add_qiime_labels.py -m $map -i Chimeras/FilteredData/ -c SampleName -o ./
seqs=combined_seqs.fna

# pick closed parameters file
echo "pick_otus:otu_picking_method usearch61_ref" > usearch_ref_params.txt
echo "pick_otus:similarity 0.97" >> usearch_ref_params.txt
echo "pick_otus:enable_rev_strand_match True" >> usearch_ref_params.txt

# utilizza pick closed
pick_closed_reference_otus.py -i $seqs -o Pick_Closed_Otu -r $PATH_DB_SEQ -t $PATH_DB_TAX -p usearch_ref_params.txt

# lavora sulla OTU TABLE
mkdir -p Otu_Tables_Closed/
cp Pick_Closed_Otu/otu_table.biom Otu_Tables_Closed/otu_table.biom

# bockulich filtering
filter_otus_from_otu_table.py -i Otu_Tables_Closed/otu_table.biom -o Otu_Tables_Closed/otu_table_filtered.biom --min_count_fraction 0.000005

biom convert -i Otu_Tables_Closed/otu_table.biom -o Otu_Tables_Closed/otu_table.txt --to-tsv --header-key taxonomy --output-metadata-id "Consensus Lineage"
biom convert -i Otu_Tables_Closed/otu_table_filtered.biom -o Otu_Tables_Closed/otu_table_filtered.txt --to-tsv --header-key taxonomy --output-metadata-id "Consensus Lineage"

# otu table conversion (phyloseq)
python $WORK/Script_Python/phyloseq_conversion.py -i Otu_Tables_Closed/otu_table.txt -o Otu_Tables_Closed/otu_table_OTU.txt
python $WORK/Script_Python/phyloseq_conversion.py -i Otu_Tables_Closed/otu_table_filtered.txt -o Otu_Tables_Closed/otu_table_filtered_OTU.tx

biom summarize-table -i Otu_Tables_Closed/otu_table_filtered.biom

