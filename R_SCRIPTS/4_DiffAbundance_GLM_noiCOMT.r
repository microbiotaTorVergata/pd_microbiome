library("phyloseq")
library("MASS")
library("pscl")
library("flexmix")
library("AICcmodavg")
library("lmtest")
library("someMTP")

# mergia due dataframe in base alle rownames
mergia.dataframe <- function( dataframe1 , dataframe2 ){
	M <- merge( dataframe1 , dataframe2 , by="row.names", all=TRUE)
	rownames(M) <- M$Row.names ; M$Row.names <- NULL
	return(M)
}

take_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

# work dir
setwd("C:\\Users\\danie\\Desktop\\Data3")

# inserisce i dati di input per phyloseq
otu <- "otu_table_filtered_OTU.txt"
map <- "mapping_ad.txt"
tre <- "rep_set.tre"

# select taxonomic level ("Family" or "Genus")
taxa="Family"

# out1 = significant taxa ; out2 = fold change ; out3 = frequency
out="AD_Family_sign.txt"
out2="AD_Family_FC.txt"
out3="AD_Family_PC.txt"

# qiimedata
qiimedata = import_qiime(otu,map,tre)

# rarefaction 
set.seed(1)
qiimedata <- rarefy_even_depth(qiimedata,sample.size=50000)		

# only patients which dont take iCOMT
qiimedata <- subset_samples(qiimedata, (Status == "PD" & iCOMT == "0") | (Status=="CON"))       
       

# save original data (used to compute FC and frequency)
qiimedata_original <- qiimedata

# select only PD and HC
pd_qiimedata <- subset_samples(qiimedata , Status == "PD" ) 
ct_qiimedata <- subset_samples(qiimedata , Status == "CON" )

# select only otus in 50% of PD or 50% of HC
pd_qiimedata = filter_taxa(pd_qiimedata, function(x) sum(x > 1) > (0.5*length(x)), TRUE)
ct_qiimedata = filter_taxa(ct_qiimedata, function(x) sum(x > 1) > (0.5*length(x)), TRUE)

# union of otus in 50% of PD or 50% of HC
tax <- union( rownames(tax_table(pd_qiimedata)) , rownames(tax_table(ct_qiimedata)) )

# select otu from qiimedata
qiimedata <- take_taxa(qiimedata,tax)

# tax glom 
qiimedata <- tax_glom( qiimedata , taxa, NArm=FALSE)
qiimedata_original <- tax_glom( qiimedata_original , taxa , NArm=FALSE )


M<-sample_data(qiimedata)
#M$Age <- factor(M$Age)
M$Yogurt <- factor(M$Yogurt)
M$Loss_5KG <- factor(M$Loss_5KG)
M$Esercizio_Fisico <- factor(M$Esercizio_Fisico)
M$Cereali <- factor(M$Cereali)
M$Yogurt <- factor(M$Yogurt)
sample_data(qiimedata) <- M


# extract data
otu_table <- otu_table(qiimedata)
map_table <- sample_data(qiimedata)
tax_table <- tax_table(qiimedata)

# vector of BIC 
names <- c()		# rownames
model_nb <- c()		# negative binomial		
model_ps <- c()		# poisson
model_hd <- c()		# hurdle
model_zi <- c()		# zero inflated nb
model_zi2 <- c()	# zero inflated poi
model_mp <- c()		# misture di poisson
for (r in 1:nrow( otu_table )  ){
	r=r
	name   <- rownames( otu_table)[r]
	names  <- c(names,name)						
	values <- t(otu_table[r,]); colnames(values) <- "taxa_freq"		
	TABLE <- mergia.dataframe( values , map_table )
	
	# negative binomial and poisson
	m2 <- glm( taxa_freq ~ Status + Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE, family=poisson(), control=list(maxit=1e3)); b2 <- BIC(m2)	
	m1 <- try(glm.nb( taxa_freq ~ Status + Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico , data=TABLE, maxit=1e3)); 
		if(inherits(m1,"try-error")){
				b1<-NA;} else {b1<-useBIC(m1)}
	
	# zero inflated model ONLY if there is at least one "0"
	if( min(TABLE$taxa_freq)==0 ){
		#print("zero")
		m3 <- try(hurdle( formula = taxa_freq ~ Status + Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE )); 
			if(inherits(m3,"try-error")){
				b3<-NA;} else {b3<-useBIC(m3)}
		m4 <- try(zeroinfl( formula = taxa_freq ~ Status + Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico , data=TABLE, dist="negbin")); 
			if(inherits(m4,"try-error")){
				b4<-NA;} else {b4<-useBIC(m4)}
		m5 <- try(zeroinfl( formula = taxa_freq ~ Status + Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico , data=TABLE, dist="poisson")); 
			if(inherits(m5,"try-error")){
				b5<-NA;} else {b5<-useBIC(m5)}
	} else { b3 <- NA; b4 <- NA; b5 <- NA;}

	model_nb <- c(model_nb,b1)				
	model_ps <- c(model_ps,b2)		
	model_hd <- c(model_hd,b3)		
	model_zi <- c(model_zi,b4)		
	model_zi2 <- c(model_zi2,b5)

}

matrice_distribuzione <- data.frame(model_nb,model_ps,model_hd,model_zi,model_zi2)
colnames(matrice_distribuzione) <- c("neg_bin","poisson","hurdle","zero_inflated_neg_bin","zero_inflated_poi")
rownames(matrice_distribuzione) <- names


# select the best model for each otu
matrice_distribuzione_final <- c()
for (r in 1:nrow( matrice_distribuzione )  ){
	Row <- matrice_distribuzione[r,]
	min_index <- which(Row == min(Row,na.rm=TRUE))
	print(colnames( matrice_distribuzione )[min_index])
	matrice_distribuzione_final <- c(matrice_distribuzione_final,colnames( matrice_distribuzione )[min_index])
}
matrice_distribuzione_final <- data.frame(matrice_distribuzione_final)
colnames(matrice_distribuzione_final) <- c("dist")
rownames(matrice_distribuzione_final) <- names


names <- c()			# otu names
pvalues <- c()			# pvalues
for (r in 1:nrow( otu_table )  ){
	cat("\n")
	r=r
	name <- rownames( otu_table)[r]					
	taxa <- tax_table[name,]
	print(r)	
	print(name)
	print(taxa)	
	names  <- c(names,name)
	values <- t(otu_table[r,]); colnames(values) <- "taxa_freq"
	TABLE <- mergia.dataframe( values , map_table )
	
	if(matrice_distribuzione_final[c(name),]=="neg_bin"){
		model_nb_1 <- glm.nb( formula = taxa_freq ~ Status + Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE, maxit=1e3)
		model_nb_2 <- try(glm.nb( formula = taxa_freq ~ Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE, maxit=1e3))
		if(inherits(model_nb_2,"try-error")){
				pvalues <- c(pvalues, NA); Residuals <- c(Residuals,NA); } 
			else {
		R1 <- anova( model_nb_1 , model_nb_2 , test="Chisq")
		pvalue_r1  <- R1$`Pr(Chi)`[2] 
		pvalues <- c(pvalues, pvalue_r1)
		Variance_r1 <- sum(resid(model_nb_2)^2)/model_nb_2$df.residual
		Residuals <- c(Residuals,Variance_r1)
		}	
	}
	
	if(matrice_distribuzione_final[c(name),]=="zero_inflated_neg_bin"){
		model_zinb_1 <- zeroinfl( formula = taxa_freq ~ Status + Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE, dist="negbin")
		model_zinb_2 <- try(zeroinfl( formula = taxa_freq ~ Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE, dist="negbin"))
		if(inherits(model_zinb_2,"try-error")){
				pvalues <- c(pvalues, NA); Residuals <- c(Residuals,NA); } 
			else {
		R2 <- lrtest(model_zinb_1,model_zinb_2)
		pvalue_r2  <- R2$`Pr(>Chisq)`[2]
		pvalues <- c(pvalues, pvalue_r2)
		Variance_r2 <- sum(resid(model_zinb_2)^2)/model_zinb_2$df.residual
		Residuals <- c(Residuals,Variance_r2)
		}
	}

	if(matrice_distribuzione_final[c(name),]=="hurdle"){
		model_h_1 <- hurdle( formula = taxa_freq ~ Status + Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE, dist="negbin")
		model_h_2 <- try(hurdle( formula = taxa_freq ~ Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE, dist="negbin"))
		if(inherits(model_h_2,"try-error")){
				pvalues <- c(pvalues, NA); Residuals <- c(Residuals,NA); } 
			else {
			R3 <- lrtest(model_h_1,model_h_2)
			pvalue_r3  <- R3$`Pr(>Chisq)`[2]
			pvalues <- c(pvalues, pvalue_r3)
			Variance_r3 <- sum(resid(model_h_2)^2)/model_h_2$df.residual
			Residuals <- c(Residuals,Variance_r3)
			}
	}
	
	if(matrice_distribuzione_final[c(name),]=="poisson"){
		model_p_1 <- glm( taxa_freq ~ Status + Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE, family=poisson(), control=list(maxit=1e3))
		model_p_2 <- glm( taxa_freq ~ Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE, family=poisson(), control=list(maxit=1e3))
		R4 <- anova( model_p_1 , model_p_2 , test="Chisq")
		pvalue_r4  <- R4$`Pr(>Chi)`[2] 
		pvalues <- c(pvalues, pvalue_r4)
		Variance_r4 <- sum(resid(model_p_2)^2)/model_p_2$df.residual
		Residuals <- c(Residuals,Variance_r4)
	}
	
	if(matrice_distribuzione_final[c(name),]=="zero_inflated_poi"){
		model_zipo_1 <- zeroinfl( formula = taxa_freq ~ Status + Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE, dist="poisson")
		model_zipo_2 <- zeroinfl( formula = taxa_freq ~ Sex + Age + Loss_5KG + BMI + Yogurt + Cereali + Esercizio_Fisico, data=TABLE, dist="poisson")
		R5 <- lrtest( model_zipo_1 , model_zipo_2)
		pvalue_r5  <- R5$`Pr(>Chisq)`[2] 
		pvalues <- c(pvalues, pvalue_r5)
		Variance_r5 <- sum(resid(model_zipo_2)^2)/model_zipo_2$df.residual
		Residuals <- c(Residuals,Variance_r5)
	}

}

# multiple test correction
fdr <- p.adjust(pvalues,n=length(pvalues),method="fdr")
# matrix
matrice <- data.frame(pvalues,fdr); colnames(matrice) <- c("P","Pc"); rownames(matrice) <- names
# taxonomy
sigtab <- mergia.dataframe(matrice, as(tax_table, "matrix"))
# save file
write.table(sigtab,file=out,sep="\t",quote=FALSE)


# FC
otu_table_original <- otu_table(qiimedata_original)
tax_table_original <- tax_table(qiimedata_original)
map_table_original <- sample_data(qiimedata_original)

names <- c()
means_pd <- c()			
means_ct <- c()		
fold_ch <- c()			
for (r in 1:nrow( otu_table_original )  ){

	cat("\n")
	r=r
	name <- rownames( otu_table_original )[r]				
	taxa <- tax_table_original[name,]
	print(r)	
	print(name)
	print(taxa)	
	names  <- c(names,name)
	values <- t(otu_table_original[r,]); colnames(values) <- "taxa_freq"
	TABLE <- mergia.dataframe( values , map_table_original )

	pd_values <- TABLE[TABLE$Status=="PD",]$taxa_freq
	ct_values <- TABLE[TABLE$Status=="CON",]$taxa_freq
	print(mean(pd_values))
	print(mean(ct_values))
	print(median(pd_values))
	print(median(ct_values))
	

	mean_pd <- mean( pd_values )
	means_pd <- c(means_pd,mean_pd)
	mean_ct <- mean( ct_values )
	means_ct <- c(means_ct,mean_ct)
	fc <- log( mean_pd / mean_ct , 2 )
	fold_ch <- c(fold_ch,fc)
}

matrice <- data.frame(fold_ch,means_pd,means_ct); colnames(matrice) <- c("Fc","mean(PD)","mean(CT)"); rownames(matrice) <- names
sigtab <- mergia.dataframe(matrice, as(tax_table_original, "matrix"))
write.table(sigtab,file=out2,sep="\t",quote=FALSE)


qiimedata_original <- transform_sample_counts(qiimedata_original, function(x) {x/sum(x)*100} )

otu_table_original <- otu_table(qiimedata_original)
tax_table_original <- tax_table(qiimedata_original)
map_table_original <- sample_data(qiimedata_original)

# ABUNDANCE
names <- c()
perc_pd <- c()			
perc_ct <- c()			
median_perc_pd <- c()
median_perc_ct <- c()
for (r in 1:nrow( otu_table_original )  ){

	cat("\n")
	r=r
	name <- rownames( otu_table_original )[r]					
	taxa <- tax_table_original[name,]
	print(r)	
	print(name)
	print(taxa)	
	names  <- c(names,name)
	values <- t(otu_table_original[r,]); colnames(values) <- "taxa_freq"		
	TABLE <- mergia.dataframe( values , map_table_original )

	pd_values <- TABLE[TABLE$Status=="PD",]$taxa_freq
	ct_values <- TABLE[TABLE$Status=="CON",]$taxa_freq
	mean_perc_pd <- mean(pd_values) 
	mean_perc_ct <- mean(ct_values)
	perc_pd <- c(perc_pd,mean_perc_pd)
	perc_ct <- c(perc_ct,mean_perc_ct)
	print(mean(pd_values))
	print(mean(ct_values))
	print(median(pd_values))
	print(median(ct_values))
	median_perc_pd <- c(median_perc_pd,median(pd_values))
	median_perc_ct <- c(median_perc_ct,median(ct_values))
}


matrice <- data.frame(perc_pd,perc_ct,median_perc_pd,median_perc_ct)
colnames(matrice) <- c("Mean Perc(PD)","Mean Perc(CT)","Median Perc(PD)","Median Perc(CT)")
rownames(matrice) <- names
sigtab <- mergia.dataframe(matrice, as(tax_table_original, "matrix"))
write.table(sigtab,file=out3,sep="\t",quote=FALSE)