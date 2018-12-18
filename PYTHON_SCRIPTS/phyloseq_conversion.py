import optparse


def convertSilvaTaxonomy(stringa_tassonomia):
	if "D_0__" in stringa_tassonomia:
		stringa_tassonomia = stringa_tassonomia.replace("D_0__","k__")
	if "D_1__" in stringa_tassonomia:
		stringa_tassonomia = stringa_tassonomia.replace("D_1__","p__")
	if "D_2__" in stringa_tassonomia:
		stringa_tassonomia = stringa_tassonomia.replace("D_2__","c__")
	if "D_3__" in stringa_tassonomia:
		stringa_tassonomia = stringa_tassonomia.replace("D_3__","o__")
	if "D_4__" in stringa_tassonomia:
		stringa_tassonomia = stringa_tassonomia.replace("D_4__","f__")
	if "D_5__" in stringa_tassonomia:
		stringa_tassonomia = stringa_tassonomia.replace("D_5__","g__")
	if "D_6__" in stringa_tassonomia:
		stringa_tassonomia = stringa_tassonomia.replace("D_6__","s__")
	return(stringa_tassonomia)


parser = optparse.OptionParser(description="")
parser.add_option('-i', '--inp' , dest='otu_table_inp', help='file input (otu table per livello tassonomico generato da summarize taxa)')
parser.add_option('-o', '--out' , dest='otu_table_out', help='file output per phyloseq')
(options,args)=parser.parse_args()

File_inp=open( options.otu_table_inp , "r" )		# open the input file (reading)
File_out=open( options.otu_table_out , "w" )		# open the output file (writing)

File_INP = File_inp.readlines()

# parse the rows in 
for i in range(0,len(File_INP)):
	
	Row = File_INP[i]

	if( i == 0 ):
		File_out.write( Row )
	
	if( i == 1 ):
		File_out.write( Row )
	
	if( i > 1 ):
		Row_splitted = Row.split("\t")		# split the row
		taxonomy = Row_splitted[-1]		# get the taxonomy
		taxonomy = convertSilvaTaxonomy(taxonomy)	# convert the taxonomy
		Row_splitted[-1] = taxonomy			# replace the taxonomy

		#print(Row_splitted[1:len(Row_splitted)-2])
		# convert to integers
		for j in range(1,len(Row_splitted)-1):
			#print Row_splitted[j]
			#print round( float(Row_splitted[j]), 0)
			Row_splitted[j] = str( round( float(Row_splitted[j]), 0) )

		#print(Row_splitted)

		#a = raw_input()
		print(Row_splitted)
		new_Row = "\t".join( Row_splitted )					# create the new row
		print(new_Row)
		File_out.write( new_Row  )							# add the row to the final file


# close the Files
File_inp.close()
File_out.close()
