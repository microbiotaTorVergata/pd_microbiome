from Bio import SeqIO
import optparse

parser = optparse.OptionParser(description="obtain a FASTA and QUAL file from a FASTA file")
parser.add_option('-f','--fq', dest='fastq', help='fastq')
parser.add_option('-a','--fa', dest='fasta', help='fasta')
parser.add_option('-q','--qual',dest='qual', help='qual')
(options,args)=parser.parse_args()

SeqIO.convert(options.fastq , "fastq" , options.fasta , "fasta")
SeqIO.convert(options.fastq , "fastq" , options.qual , "qual")

