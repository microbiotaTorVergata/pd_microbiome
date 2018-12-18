library("plyr")
library("phyloseq")
library("ggplot2")
library("vegan")

Save.graph<- function( Graph , OutDir , Name , Name_prefix="" ){
	out_file_name <- paste( Name_prefix,Name,sep="_")
	out_path_name <- paste( OutDir, out_file_name ,  sep="//" )
	print(out_path_name)
	jpeg( paste(out_path_name, ".jpeg") , height=15, width=25, units="cm", res=300)
	print(Graph)
	dev.off()
	return()
}

makeBarplot <- function( dataframe, xvalues, yvalues, xlabel, ylabel, title){
  
  Barplot <- ggplot(dataframe,aes(x = xvalues, y = yvalues, fill = factor(xvalues) )) +
    geom_boxplot(  ) +
    #geom_jitter(aes(x = xvalues, y = yvalues),
    #            position=position_jitter(width=0.1,height=0),
    #            alpha=0.6,
    #            size=3,
    #            show.legend=FALSE) +
    xlab( xlabel ) +
    ylab( ylabel ) +
    ggtitle( title ) #+
    #scale_fill_manual(values = cols) 
  
  return( Barplot )
}


# set work dir
setwd("..")

# qiime files
otu <- "otu_table_filtered_OTU.txt"
map <- "mapping_ad.txt"
tre <- "rep_set.tre"

# create a phyloseq object
qiimedata = import_qiime(otu,map,tre)		      

# rarefaction
set.seed(1); qiimedata = rarefy_even_depth(qiimedata,50000)

map_table <- sample_data(qiimedata)

# tax glom to species level
SPECIES <- tax_glom( qiimedata , "Species", NArm=FALSE)

# diversity indices
DIV <- estimate_richness(SPECIES)


measures_metadata_rare <- merge( map_table , DIV , by="row.names" , all = TRUE )
rownames(measures_metadata_rare) <- measures_metadata_rare$Row.names
measures_metadata_rare$Row.names <- NULL

shannon_PD <- measures_metadata_rare[measures_metadata_rare$Status=="PD",]$Shannon
shannon_CT <- measures_metadata_rare[measures_metadata_rare$Status=="CON",]$Shannon
wilcox.test(shannon_PD,shannon_CT,paired=FALSE)
S <- makeBarplot(dataframe=measures_metadata_rare, xvalues=measures_metadata_rare$Status, yvalues=measures_metadata_rare$Shannon, xlabel="Status", ylabel="Shannon index", title="")
Save.graph(Graph=S , OutDir="..i" , Name="shannon" , Name_prefix="diversity_" )

simpson_PD <- measures_metadata_rare[measures_metadata_rare$Status=="PD",]$Simpson
simpson_CT <- measures_metadata_rare[measures_metadata_rare$Status=="CON",]$Simpson
wilcox.test(simpson_PD,simpson_CT,paired=FALSE)
s<- makeBarplot(dataframe=measures_metadata_rare, xvalues=measures_metadata_rare$Status, yvalues=measures_metadata_rare$Simpson, xlabel="Status", ylabel="Simpson index", title="")
Save.graph(Graph=S , OutDir=".." , Name="shannon" , Name_prefix="diversity_" )


