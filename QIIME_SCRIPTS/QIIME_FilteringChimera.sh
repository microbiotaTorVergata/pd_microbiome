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

# conversion fastq to fasta and qual
mkdir -p Fasta/ Qual/
for fastq in `ls Trimmed/*.fastq`
do	
	echo $fastq
	name=$(basename $fastq)
	fasta=${name/trimmed.fastq/fasta}
	qual=${name/trimmed.fastq/qual}
	echo $fasta
	echo $qual
	python $PATH_WORK/Script_Python/convert_FASTQ_to_FASTA_QUAL.py -f $fastq -a Fasta/$fasta -q Qual/$qual
done

rm -r Qual  

mkdir -p Chimeras/
mkdir -p Chimeras/IdentifiedChimeras
mkdir -p Chimeras/FilteredData

echo "analisi sequenze chimeriche"
for fileFasta in `ls Fasta/`
do
	ls -lhrt Fasta/$fileFasta
	outdir=${fileFasta/.fasta/}
	identify_chimeric_seqs.py -m usearch61 -i Fasta/$fileFasta -r $PATH_DB_SEQ -o Chimeras/IdentifiedChimeras/$outdir
	outfile=$fileFasta
	filter_fasta.py -f Fasta/$fileFasta -s Chimeras/IdentifiedChimeras/$outdir/chimeras.txt -o Chimeras/FilteredData/$outfile -n
done
echo ""

source deactivate
