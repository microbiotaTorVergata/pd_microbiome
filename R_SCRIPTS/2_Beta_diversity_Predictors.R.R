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

# work dir 
setwd("..")

# phyloseq data
otu <- "otu_table_filtered_OTU.txt"
map <- "mapping_ad.txt"
tre <- "rep_set.tre"

# phyloseq object
qiimedata = import_qiime(otu,map,tre)		     

# rarefaction
set.seed(1); qiimedata = rarefy_even_depth(qiimedata,50000)

# metrics
D_bray <- distance(qiimedata,method="bray")
D_Uuni <- distance(qiimedata,method="unifrac")
D_Wuni <- distance(qiimedata,method="wunifrac")

#### PERMANOVA TEST

M<-sample_data(qiimedata)
M$Status <- factor(M$Status)
M$Provincia <- factor(M$Provincia)
M$Gain_5KG <- factor(M$Gain_5KG)
M$Loss_5KG <- factor(M$Loss_5KG)
M$Yogurt <- factor(M$Yogurt)
M$Pasta <- factor(M$Pasta)
M$Pane <- factor(M$Pane)
M$Latticini <- factor(M$Latticini)
M$Frutta_Verdura <- factor(M$Frutta_Verdura)
M$Carne <- factor(M$Carne)
M$Pesce <- factor(M$Pesce)
M$Cereali <- factor(M$Cereali)
M$Legumi <- factor(M$Legumi)
M$Caffe <- factor(M$Caffe)
M$Alcol <- factor(M$Alcol)
M$Pizza <- factor(M$Pizza)
M$Fuma <- factor(M$Fuma)
M$Esercizio_Fisico <- factor(M$Esercizio_Fisico)
M$Parto_Cesareo <- factor(M$Parto_Cesareo)
M$Phenotype <- factor(M$Phenotype)
M$L_Dopa <- factor(M$L_Dopa)
M$DA <- factor(M$DA)
M$iMAo <- factor(M$iMAo)
M$iCOMT <- factor(M$iCOMT)
M$AntiCH <- factor(M$AntiCH)
M$Amantadina <- factor(M$Amantadina)
M$Antiacidi <- factor(M$Antiacidi)

sample_data(qiimedata) <- M

metadata <- as(sample_data(qiimedata), "data.frame")

# only PD-Status
set.seed(1); ad1=adonis2( D_bray ~ Status, data = metadata , permutations = 9999, by = "margin")
print(ad1)
set.seed(1); ad1=adonis2( D_Uuni ~ Status, data = metadata , permutations = 9999, by = "margin")
print(ad1)
set.seed(1); ad1=adonis2( D_Wuni ~ Status, data = metadata , permutations = 9999, by = "margin")
print(ad1)

mergia.dataframe <- function( dataframe1 , dataframe2 ){
	M <- merge( dataframe1 , dataframe2 , by="row.names", all=TRUE)
	rownames(M) <- M$Row.names ; M$Row.names <- NULL
	return(M)
}


# STEPWISE BACKWARD BC
set.seed(1); ad0=adonis2( D_bray ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Pesce + Cereali + Legumi + Caffe + Alcol + Pizza + Fuma + Esercizio_Fisico + Provincia + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad1=adonis2( D_bray ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Pesce + Cereali + Legumi + Caffe + Alcol + Pizza + Fuma + Esercizio_Fisico + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad2=adonis2( D_bray ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Pesce + Cereali + Legumi + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad3=adonis2( D_bray ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Cereali + Legumi + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad4=adonis2( D_bray ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad5=adonis2( D_bray ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad6=adonis2( D_bray ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Fuma + Gain_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad7=adonis2( D_bray ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Fuma + Gain_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad8=adonis2( D_bray ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Frutta_Verdura + Cereali + Caffe + Alcol + Fuma + Gain_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad9=adonis2( D_bray ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Cereali + Caffe + Alcol + Fuma + Gain_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad10=adonis2( D_bray ~ Status + Age + Sex + BMI + Yogurt + Pane + Cereali + Caffe + Alcol + Fuma + Gain_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad11=adonis2( D_bray ~ Status + Age + BMI + Yogurt + Pane + Cereali + Caffe + Alcol + Fuma + Gain_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad12=adonis2( D_bray ~ Status + Age + BMI + Yogurt + Pane + Cereali + Caffe + Fuma + Gain_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad13=adonis2( D_bray ~ Status + Age + BMI + Pane + Cereali + Caffe + Fuma + Gain_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad14=adonis2( D_bray ~ Status + Age + BMI + Pane + Cereali + Caffe + Fuma + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad15=adonis2( D_bray ~ Status + Age + BMI + Cereali + Caffe + Fuma + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad16=adonis2( D_bray ~ Status + Age + BMI + Cereali + Caffe + Fuma, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad17=adonis2( D_bray ~ Status + Age + BMI + Cereali + Fuma, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad18=adonis2( D_bray ~ Status + Age + BMI + Fuma, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad19=adonis2( D_bray ~ Status + Age + BMI, data = metadata , permutations = 9999, by = "margin")

# STEPWISE BACKWARD UU
set.seed(1); ad0=adonis2( D_Uuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Pesce + Cereali + Legumi + Caffe + Alcol + Pizza + Fuma + Esercizio_Fisico + Provincia + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad1=adonis2( D_Uuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Pesce + Cereali + Caffe + Alcol + Pizza + Fuma + Esercizio_Fisico + Provincia + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad2=adonis2( D_Uuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Pesce + Cereali + Caffe + Alcol + Pizza + Fuma + Esercizio_Fisico + Provincia + Gain_5KG + Loss_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad3=adonis2( D_Uuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Pesce + Cereali + Caffe + Alcol + Pizza + Fuma + Esercizio_Fisico + Gain_5KG + Loss_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad4=adonis2( D_Uuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Pizza + Fuma + Esercizio_Fisico + Gain_5KG + Loss_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad5=adonis2( D_Uuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Pizza + Fuma + Esercizio_Fisico + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad6=adonis2( D_Uuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Pizza + Fuma + Esercizio_Fisico + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad7=adonis2( D_Uuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad8=adonis2( D_Uuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Carne + Cereali + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad9=adonis2( D_Uuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Cereali + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad10=adonis2( D_Uuni ~ Status + Age + BMI + Yogurt + Pane + Pasta + Cereali + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad11=adonis2( D_Uuni ~ Status + Age + BMI + Yogurt + Pasta + Cereali + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad12=adonis2( D_Uuni ~ Status + Age + BMI + Yogurt + Cereali + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad13=adonis2( D_Uuni ~ Status + Age + BMI + Yogurt + Cereali + Alcol + Fuma + Esercizio_Fisico + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad14=adonis2( D_Uuni ~ Status + Age + BMI + Yogurt + Cereali + Alcol + Esercizio_Fisico + Gain_5KG, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad15=adonis2( D_Uuni ~ Status + Age + BMI + Yogurt + Cereali + Alcol + Esercizio_Fisico, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad16=adonis2( D_Uuni ~ Status + Age + BMI + Yogurt + Cereali + Esercizio_Fisico, data = metadata , permutations = 9999, by = "margin")

# STEPWISE BACKWARD WU
set.seed(1); ad0=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Pesce + Cereali + Legumi + Caffe + Alcol + Pizza + Fuma + Esercizio_Fisico + Provincia + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad1=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Latticini + Frutta_Verdura + Carne + Pesce + Cereali + Legumi + Caffe + Alcol + Fuma + Esercizio_Fisico + Provincia + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad2=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Frutta_Verdura + Carne + Pesce + Cereali + Legumi + Caffe + Alcol + Fuma + Esercizio_Fisico + Provincia + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad3=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Frutta_Verdura + Carne + Pesce + Cereali + Caffe + Alcol + Fuma + Esercizio_Fisico + Provincia + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad4=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Frutta_Verdura + Carne + Pesce + Cereali + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad5=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Fuma + Esercizio_Fisico + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad6=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Pasta + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Fuma + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad7=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Fuma + Gain_5KG + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad8=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Pane + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Fuma + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad9=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Fuma + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad10=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Frutta_Verdura + Carne + Cereali + Caffe + Alcol + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad11=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Frutta_Verdura + Cereali + Caffe + Alcol + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad12=adonis2( D_Wuni ~ Status + Age + Sex + BMI + Yogurt + Cereali + Caffe + Alcol + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad13=adonis2( D_Wuni ~ Status + Age + BMI + Yogurt + Cereali + Caffe + Alcol + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad14=adonis2( D_Wuni ~ Status + Age + BMI + Yogurt + Cereali + Caffe + Loss_5KG + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad15=adonis2( D_Wuni ~ Status + Age + BMI + Yogurt + Cereali + Caffe + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad16=adonis2( D_Wuni ~ Status + Age + BMI + Yogurt + Caffe + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad17=adonis2( D_Wuni ~ Status + Age + BMI + Caffe + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad18=adonis2( D_Wuni ~ Status + Age + BMI + Parto_Cesareo, data = metadata , permutations = 9999, by = "margin")
set.seed(1); ad19=adonis2( D_Wuni ~ Status + Age + BMI, data = metadata , permutations = 9999, by = "margin")

#### Only patients
otu <- "Otu_Tables_Closed_NoRep/otu_table_filtered_OTU.txt"
map <- "Otu_Tables_Closed_NoRep/new_mapping.txt"
tre <- "Phylogeny/rep_set.tre"


qiimedata = import_qiime(otu,map,tre)		        # crea un oggetto phyloseq

set.seed(1); qiimedata = rarefy_even_depth(qiimedata,50000)

qiimedata = subset_samples(qiimedata, Status == "PD")

D_bray <- distance(qiimedata,method="bray")
D_Uuni <- distance(qiimedata,method="unifrac")
D_Wuni <- distance(qiimedata,method="wunifrac")

ordPCoA_bray <- ordinate(qiimedata, method = "PCoA", distance = D_bray)
ordPCoA_Uuni <- ordinate(qiimedata, method = "PCoA", distance = D_Uuni )
ordPCoA_Wuni <- ordinate(qiimedata, method = "PCoA", distance = D_Wuni )


M<-sample_data(qiimedata)
M$Status <- factor(M$Status)
M$Provincia <- factor(M$Provincia)
M$Gain_5KG <- factor(M$Gain_5KG)
M$Loss_5KG <- factor(M$Loss_5KG)
M$Yogurt <- factor(M$Yogurt)
M$Pasta <- factor(M$Pasta)
M$Pane <- factor(M$Pane)
M$Latticini <- factor(M$Latticini)
M$Frutta_Verdura <- factor(M$Frutta_Verdura)
M$Carne <- factor(M$Carne)
M$Pesce <- factor(M$Pesce)
M$Cereali <- factor(M$Cereali)
M$Legumi <- factor(M$Legumi)
M$Caffe <- factor(M$Caffe)
M$Alcol <- factor(M$Alcol)
M$Pizza <- factor(M$Pizza)
M$Fuma <- factor(M$Fuma)
M$Esercizio_Fisico <- factor(M$Esercizio_Fisico)
M$Phenotype <- factor(M$Phenotype)
M$L_Dopa <- factor(M$L_Dopa)
M$DA <- factor(M$DA)
M$iMAo <- factor(M$iMAo)
M$iCOMT <- factor(M$iCOMT)
M$AntiCH <- factor(M$AntiCH)
M$Amantadina <- factor(M$Amantadina)

sample_data(qiimedata) <- M


metadata <- as(sample_data(qiimedata), "data.frame")


# STEPWISE BACKWARD Bray-Curtis
set.seed(1); ad1 <- adonis2(formula = D_bray ~ Phenotype + Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + DA + iMAo + iCOMT + Amantadina, data = metadata, permutations = 9999) 
set.seed(1); ad2 <- adonis2(formula = D_bray ~ Phenotype + Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + DA + iCOMT + Amantadina, data = metadata, permutations = 9999) 
set.seed(1); ad3 <- adonis2(formula = D_bray ~ Phenotype + Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + iCOMT + Amantadina, data = metadata, permutations = 9999) 
set.seed(1); ad4 <- adonis2(formula = D_bray ~ Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + iCOMT + Amantadina, data = metadata, permutations = 9999) 
set.seed(1); ad5 <- adonis2(formula = D_bray ~ Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + iCOMT, data = metadata, permutations = 9999) 
set.seed(1); ad6 <- adonis2(formula = D_bray ~ Disase_Duration + Staging + Tot_L_Dopa_die + iCOMT, data = metadata, permutations = 9999) 

# STEPWISE BACKWARD UU
set.seed(1); ad1 <- adonis2(formula = D_Uuni ~ Phenotype + Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + DA + iMAo + iCOMT + Amantadina, data = metadata, permutations = 9999) set.seed(1); ad2 <- adonis2(formula = D_Uuni ~ Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + DA + iMAo + iCOMT + Amantadina, data = metadata, permutations = 9999) set.seed(1); ad3 <- adonis2(formula = D_Uuni ~ Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + DA + iCOMT + Amantadina, data = metadata, permutations = 9999)  
set.seed(1); ad5 <- adonis2(formula = D_Uuni ~ Disase_Duration + Staging + Tot_L_Dopa_die + DA + iCOMT + Amantadina, data = metadata, permutations = 9999) 
set.seed(1); ad6 <- adonis2(formula = D_Uuni ~ Disase_Duration + Staging + Tot_L_Dopa_die + iCOMT + Amantadina, data = metadata, permutations = 9999) 
set.seed(1); ad7 <- adonis2(formula = D_Uuni ~ Disase_Duration + Staging + Tot_L_Dopa_die + iCOMT, data = metadata, permutations = 9999) 

# STEPWISE BACKWARD WU
set.seed(1); ad1 <- adonis2(formula = D_Wuni ~ Phenotype + Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + DA + iMAo + iCOMT + Amantadina, data = metadata, permutations = 9999) 
set.seed(1); ad2 <- adonis2(formula = D_Wuni ~ Phenotype + Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + DA + iCOMT + Amantadina, data = metadata, permutations = 9999) 
set.seed(1); ad3 <- adonis2(formula = D_Wuni ~ Phenotype + Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + DA + iCOMT, data = metadata, permutations = 9999) 
set.seed(1); ad4 <- adonis2(formula = D_Wuni ~ Phenotype + Disase_Duration + Staging + Tot_L_Dopa_die + L_Dopa + iCOMT, data = metadata, permutations = 9999) 
set.seed(1); ad5 <- adonis2(formula = D_Wuni ~ Phenotype + Disase_Duration + Staging + Tot_L_Dopa_die + iCOMT, data = metadata, permutations = 9999) 
set.seed(1); ad6 <- adonis2(formula = D_Wuni ~ Disase_Duration + Staging + Tot_L_Dopa_die + iCOMT, data = metadata, permutations = 9999) 