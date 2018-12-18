library("plyr")
library("phyloseq")
library("ggplot2")
library("reshape2")
library("scales")
library("vegan")

# insert work dir
setwd("..")

# path to QIIME data
otu <- "otu_table_filtered_OTU.txt"
map <- "mapping_ad.txt"
tre <- "rep_set.tre"


qiimedata = import_qiime(otu,map)

# rarefaction
set.seed(1)
qiimedata <- rarefy_even_depth(qiimedata,sample.size=50000)

# load sample data
TT<-sample_data(qiimedata)

# logistic regression 

y=ifelse(TT$Status=="PD",1,0)
                                       
summary(glm(y~TT$Sex,family=binomial))$coef[2,4]
summary(glm(y~TT$Age,family=binomial))$coef[2,4]
summary(glm(y~TT$Parto,family=binomial))$coef[2,4]
summary(glm(y~TT$BMI,family=binomial))$coef[2,4]
summary(glm(y~TT$Gain_5KG,family=binomial))$coef[2,4]
summary(glm(y~TT$Loss_5KG,family=binomial))$coef[2,4]
summary(glm(y~TT$Yogurt,family=binomial))$coef[2,4]
summary(glm(y~TT$Pane,family=binomial))$coef[2,4]
summary(glm(y~TT$Pasta,family=binomial))$coef[2,4]
summary(glm(y~TT$Latticini,family=binomial))$coef[2,4]
summary(glm(y~TT$Frutta_Verdura,family=binomial))$coef[2,4]
summary(glm(y~TT$Carne,family=binomial))$coef[2,4]
summary(glm(y~TT$Pesce,family=binomial))$coef[2,4]
summary(glm(y~TT$Cereali,family=binomial))$coef[2,4]
summary(glm(y~TT$Legumi,family=binomial))$coef[2,4]
summary(glm(y~TT$Caffe,family=binomial))$coef[2,4]
summary(glm(y~TT$Alcol,family=binomial))$coef[2,4]
summary(glm(y~TT$Pizza,family=binomial))$coef[2,4]
summary(glm(y~TT$Fuma,family=binomial))$coef[2,4]
summary(glm(y~TT$Esercizio_Fisico,family=binomial))$coef[2,4]