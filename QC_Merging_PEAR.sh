
# path to PEAR bin
PATH_BIN="/../../bin/PEAR/bin"
export PATH="$PATH_BIN:$PATH"                      		        

mkdir -p Merged/

for forward in `ls DATA/*_R1_*`
do
	reverse=${forward/R1/R2}
	ls -lhrt $forward
	ls -lhrt $reverse
	name=$(basename $forward)
	name=${name/_L001_R1_001.fastq/}
	outdir=Merged/$name
	mkdir -p $outdir
	pear-0.9.6-bin-64 -f $forward -r $reverse -o $outdir/$name --p-value 0.01 --min-overlap 100 -j 14
done

