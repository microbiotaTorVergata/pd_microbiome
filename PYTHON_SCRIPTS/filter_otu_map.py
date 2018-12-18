import optparse

parser = optparse.OptionParser(description="Filter a otu mapping file using a biom converted to txt")
parser.add_option('-i', '--inp' , dest='otu_table_inp', help='file input (otu map originale)')
parser.add_option('-b', '--biom', dest='biom', help='file biom da usare per il filtering')
parser.add_option('-o', '--out' , dest='otu_table_out', help='file output (otu map modificata)')
(options,args)=parser.parse_args()

File_bio=open( options.biom , "r" )
ids  = []
for line in File_bio:
	if line[0] != "#":
		ids.append(line.split("\t")[0])
File_bio.close()

print "Filtering otu map"

File_inp=open( options.otu_table_inp , "r" ) 
File_out=open( options.otu_table_out , "w" )      

for line in File_inp:
	elemento_inp = line.split("\t")[0]
	if elemento_inp in ids:
		File_out.write( line )

File_inp.close()
File_out.close()
