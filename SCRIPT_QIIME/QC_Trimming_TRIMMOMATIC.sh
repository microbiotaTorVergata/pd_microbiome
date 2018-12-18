
# path to trimmomatic bin
PATH_TRIMM="/galileo/home/userexternal/dpietruc/bin/Trimmomatic-0.36/"

inp_path="Merged"
out_path="Trimmed"

mkdir -p $out_path

for file in `ls $inp_path/*/*.assembled.fastq*`
	echo $file
	output=$(basename $file)
	output=${output/assembled.fastq/trimmed.fastq}
	output=$out_path/$output
	log=${output/trimmed.fastq/log.txt}
	java -jar $PATH_TRIMM/trimmomatic-0.36.jar SE -phred33 -trimlog $log $file $output LEADING:30 TRAILING:30 AVGQUAL:30 MINLEN:400
	echo ""
done

