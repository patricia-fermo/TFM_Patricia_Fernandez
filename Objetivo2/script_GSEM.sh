### SCRIPT PARA GENOMIC SEM ###
##
## Marzo 2023
##
## Este script está pensado para obtener factores latentes de SCZ sin ISO (Isolation) y SCZ con ISO, mediante la metodología de GenomicSEM, y aplicando GWAS by substraction
##Tutorial aqui: https://rpubs.com/MichelNivard/565885

## Archivos necesarios:
## 1) Summary data que vayamos a usar: En este caso, SCZ-PGC3 y ISO UKBB 2018
## 2) rsID de las variantes de 1000G a frecuencia MAF > 0.05%: reference.1000G.maf.0.005.txt
## 3) Lista de SNPs que vayamos a emplear: w_hm3.snplist (SNPs de hapmap 3: 1173569 SNPs)
## 4) Archivos de LD score y LD weights para el cálculo de LDSC: eur_w_ld_chr


## 0) INSTALAMOS PROGRAMA Y DEPENDENCIAS
library(devtools)
install_github("GenomicSEM/GenomicSEM")
require(GenomicSEM)
library(GenomicSEM)

library("data.table")
library("dplyr")
library("ggplot2")

# Directorio de trabajo
setwd("/home/nodotea/PATRICIA/Nuevos_objetivos/Objetivo3")


################################################
## A) PREPARACIÓN DE SUMMARIES #################
################################################

## En primer lugar necesitamos pre-formatear los summary statistics que vayamos a usar para que tenga la info necesaria, y calcular N effective

### SCZ (PGC3 - Trubetskoy et al., 2022)
SCZ<-fread("PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv",data.table=FALSE)
SCZ$NCAS <-NULL
SCZ$NCON <-NULL
colnames(SCZ)<-c("CHR", "SNP", "BP", "A1", "A2", "FCAS", "FCON", "INFO", "BETA", "SE", "P", "Neff")  
write.table(SCZ, file = "SCZ_withNeff.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)



### Isolation (UKBB 2018)
ISO <-fread("1031.gwas.imputed_v3.both_sexes.tsv",data.table=FALSE)
head(ISO)

## Como no tiene rsID se lo pongo de 1000Genomas
ref<-fread("reference.1000G.maf.0.005.txt",data.table=FALSE)

ref$ID1<-paste(ref$CHR, ref$BP, ref$A1, ref$A2, sep=":")
ref$ID2<-paste(ref$CHR, ref$BP, ref$A2, ref$A1, sep=":")

ID1 = ref[,c("SNP", "ID1")]
ID2 = ref[,c("SNP", "ID2")]
colnames(ID1)<-c("SNP", "ID")
colnames(ID2)<-c("SNP", "ID")
ref_ID = rbind(ID1,ID2)

temp = merge(ISO,ref_ID,by.x = "variant", by.y = "ID")
head(temp)
ISO <- temp
rm(ID2,ID1,ref_ID,temp)

#separo la columna variant para obtener info de chr y bp
require(stringr)
ISO$new = ISO$variant
ISO$new2 = ISO$variant
ISO$new3 = ISO$variant
ISO$BP <-sapply(strsplit(ISO$variant, ":"), '[', 2)
ISO$CHR<-sapply(strsplit(ISO$new, ":"), '[', 1)
ISO$A1<-sapply(strsplit(ISO$new2, ":"), '[', 3)
ISO$A2<-sapply(strsplit(ISO$new3, ":"), '[', 4)

ISO<-ISO[,c("CHR", "BP", "A1", "A2", "SNP", "beta", "se", "pval", "n_complete_samples")]
colnames(ISO)<-c("CHR", "BP", "A2", "A1", "SNP", "BETA", "SE","P", "Neff")
write.table(ISO, file = "ISO_withNeff.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)



################################################
## B) Munge the summary statistics      ########
################################################
# En este paso usaremos los archivos creados para crear un data table conjunto de ambos, con mismas variantes

#create vector of the summary statistics files
files<-c("SCZ_withNeff.txt", "ISO_withNeff.txt")

#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3
hm3<-"eur_w_ld_chr/w_hm3.snplist"

#name the traits 
trait.names<-c("SCZ", "ISO")

#list the sample sizes. All but PTSD have SNP-specific sum of effective sample sizes so only its
#sample size is listed here
N=c(NA,NA)

#definte the imputation quality filter
info.filter=0.9

#define the MAF filter
maf.filter=0.01

#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)


################################################
## C) Run multivariable LDSC       #############
################################################
# Usamos esos archivos creados para aplicar LDSC (LD-score regression), que nos permite ver la estructura de covarianza de ambos rasgos teniendo en cuenta los LD blocks para evitar sesgos

#vector of munged summary statisitcs
traits<-c("SCZ.sumstats.gz","ISO.sumstats.gz")

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-c(.5,.5,.5,.5)

#vector of population prevalences
population.prev<-c(.01,.312)

#the folder of LD scores
ld<-"eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"eur_w_ld_chr/"

#name the traits
trait.names<-c("SCZ","ISO")

#run LDSC
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names,stand=TRUE)

#optional command to save the output as a .RData file for later use
save(LDSCoutput,file="LDSCoutput.RData")

## If you want to output the standard errors of the ld-score regression in the order they are listed in the genetic covariance (i.e., S) matrix, then you can run the three lines of code below. These SEs are also listed in the .log file produced by ldsc.
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))




##################################################
## D) Corremos el modelo  sin SNPs              ##
##################################################
## Empleando el archivos de LDSC generado podemos correr el modelo que vayamos a probar para ver las stats del mismo, antes de generar GWAS data para los factores latentes (SCZ_ISO y SCZ_only)
# Vemos si converge bien, y los parámetros son adecuados

load(file="LDSCoutput.RData")

model<-'SCZ_ISO=~NA*SCZ + start(0.1)*ISO
SCZ_only=~NA*SCZ
SCZ_only~~1*SCZ_only
SCZ_ISO~~1*SCZ_ISO
SCZ_ISO~~0*SCZ_only
ISO ~~ 0*SCZ
ISO~~0*ISO
SCZ~~0*SCZ'

output<-usermodel(LDSCoutput,estimation="DWLS",model=model)
output
save(output, file="SCZ_ISO_wo_SNPs.Rdata" )

#€ interpretación: https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects







################################################
###    E) MULTIVARIATE GWAS    #################
################################################
## escribimos el modelo en lenguaje Laavan

#user model
model<-'SCZ_ISO=~NA*SCZ + start(0.1)*SCZ + start(0.2)*ISO
SCZ_only=~NA*SCZ + start(0.7)*SCZ
SCZ_only~SNP
SCZ_ISO~SNP
SCZ_only~~1*SCZ_only
SCZ_ISO~~1*SCZ_ISO
SCZ_ISO~~0*SCZ_only
ISO ~~ 0*SCZ
ISO~~0*ISO
SCZ~~0*SCZ
SNP~~SNP'


################################################
## E1) Prepare the summary statistics for GWAS #
################################################

### necesitamos preparar los summaries para hacer el multiivariate GWAS. Para ello lo importante es corregir las betas que vengan de distintas medidas en el gwas de origen (logistic, lineal, etc...)

## The sumstats function takes 4 necessary arguments, along with 6 additional arguments whose default is NULL:
## files: The name of the summary statistics files.
## ref: The reference file used to calculate SNP variance across traits. (usamos 1000 genomas europea - provided file)
## trait.names: The names of the traits in the order that they are listed for the files.
## se.logit: Whether the SEs are on a logistic scale. Please note that for many traits, the SEs will be reported on a logistic scale even when the regression effects are reported as odds ratios (ORs). This information can usually be found in the README file for the summary statistics.
## OLS: Whether the phenotype was a continuous outcome analyzed using an observed least square (OLS; i.e., linear) estimator. If no information is provided the default is for the function to assume FALSE. 
## linprob: Whether the phenotype was a dichotomous outcome for which there are only Z-statistics in the summary statistics file -or- it was a dichotomous outcome analyzed using an OLS estimator
## N: A user provided N listed in the order the traits are listed for the files argument. If no information is provided, the default is for the function to assume NULL. When backing out a logistic beta using the linprob argument this requires the sum of effective sample sizes across the cohorts contributing to the GWAS
## betas: This argument can be used to directly supply the column name in the GWAS summary stats of the betas for continuous traits if the user knows the GWAS was run on an an already standardized phenotype.
## info.filter: The INFO filter to use. The package default is 0.6.
## maf.filter: The MAF filter to use. The package default is 0.01.
## keep.indel: Whether insertion deletions (indels) should be included in the output. Default = FALSE.
## parallel: Whether the function should be run in parallel. Default = FALSE.
## cores: When running in parallel, whether you want the computer to use a certain number of cores.

## Para saber cómo poner los parámetros según nuestro GWAS: https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-resources-and-key-information

## En nuestro caso:
## ISO se ha estimado con un modelo de least-squares linear model, así que asumimos que linprog es TRUE, mientras que en el caso de SCZ es un additive logistic regression model.
## SCZ : se.logit = T, OLS = F, linprob = F;  ISO : se.logit = F, OLS = F, linprob = T



######################################################################################
## E2) Combine the summary statistics and LDSC output and RUN USER-SPECIFIED MODEL ####
######################################################################################

## In order to run multivariate GWAS, the S and V matrices from LDSC and the newly prepared summary statistics files first need to be merged to create as many S and V matrices as there are SNPs.
## The userGWAS function does this by expanding the S and V matrices to include the SNP effects
## After constructing the S and V matrices, the function will run the user-specified model
## The userGWAS function takes the following arguments. Note that only covstruc, SNPs, and model are necessary arguments.

## MANDATORY
## covstruc: The output from LDSC.
## SNPs: The output from sumstats.
## model: The model that is being estimated (written in lavaan syntax)

## OPTIONAL
## estimation: Whether the models should be estimated using DWLS or ML estimation. The package default is DWLS.
## printwarn: An optional argument specifying whether the user wants individual lavaan warnings and errors to be included in the output. Default = TRUE
## sub: An optional argument specifying whether or not the user is requesting only specific components of the model output to be saved. 
## cores: When running in parallel, an optional argument specifying whether or not the user wants to use a certain number of cores for parallel processing. 
## toler: What tolerance level to use for matrix inversions. This is only something that needs to be of concern if warning/error messages are produced to the effect of "matrix could not be inverted".
## SNPSE: An optional argument specifying whether you want to set the SE of the SNP to a lower or higher value than the default of .0005. 
## parallel: An optional argument specifying whether you want the function to be run in parallel, or to be run serially. 
## GC: Level of Genomic Control (GC) you want the function to use. The default is 'standard' which adjusts the univariate GWAS standard errors by multiplying them by the square root of the univariate LDSC intercept.
## MPI: Whether the function should use multi-node processing (i.e., MPI). 
## smooth_check: Whether the function should save the largest Z-statistic difference between pre- and post-smoothed genetic covariance matrices.




## IMPORTANTE: corro tan solo los GWAS para los SNPs de hapmap3 (aprox 1,2 M SNPs). Son, por ejemplo, los que usa LDpred2 para el cálculo de PRS. Para mayor fluidez

## Creamos los archivos
SCZ<-fread("SCZ_withNeff.txt",data.table=FALSE)
#LNL<-fread("LNL_withNeff.txt",data.table=FALSE)
ISO<-fread("ISO_withNeff.txt",data.table=FALSE)

## SNPS EN HAPMAP 3
SNPlist<-fread("w_hm3.snplist",data.table=FALSE)

temp <- filter(SCZ, SNP  %in% SNPlist$SNP)
SCZ_hm3 <- temp

#temp <- filter(LNL, SNP  %in% SNPlist$SNP)
#LNL_hm3 <- temp

temp <- filter(ISO, SNP  %in% SNPlist$SNP)
ISO_hm3 <- temp

write.table(SCZ_hm3, file = "SCZ_hm3.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)
#write.table(LNL_hm3, file = "LNL_hm3.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)
write.table(ISO_hm3, file = "ISO_hm3.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)



## Creamos los archivos de ambos traits por cromosoma, para luego correr el multivariateGWAS de forma más rápida y eficiente computacionalmente
traitlist <- c("SCZ","ISO")
for(i in traitlist)
	{
	for (j in 1:22)
		{
		trait = fread(paste0(i,"_hm3.txt"), data.table=FALSE)
		temp = subset(trait, CHR == j)
		write.table(temp, file = paste0(i,".chrom.",j),sep="\t", quote=FALSE, row.names=FALSE)	
		}	
	}


## junto sumstats por chr y corro genomicSEM
for (j in 1:22)
	{
	#junto sumstats
	files<-c(paste0("SCZ.chrom.",j),paste0("ISO.chrom.",j))
	ref= "reference.1000G.maf.0.005.txt"
	trait.names<-c("SCZ","ISO")
	se.logit <-c(T,F)
	linprob<-c(F,T)
	info.filter=.6
	maf.filter=0.01
	OLS=NULL
	N=NULL
	betas=NULL

	SCZ_ISO_sumstats <- sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS = NULL,linprob=linprob, N=N, betas=NULL, info.filter=.6, maf.filter=0.01,keep.indel=FALSE,parallel=TRUE,cores=NULL)

	# corro genomicSEM
	outputGWAS<-userGWAS(covstruc = LDSCoutput, SNPs = SCZ_ISO_sumstats, model = model, sub =c("SCZ_only~SNP","SCZ_ISO~SNP"))

	write.table(outputGWAS[[1]], file = paste0("SCZ_only.chrom.",j), sep="\t", quote=FALSE, row.names=FALSE)	
	write.table(outputGWAS[[2]], file = paste0("SCZ_ISO.chrom.",j), sep="\t", quote=FALSE, row.names=FALSE)	
	}


### Junto los archivos exportados para unir los sumstats generados SCZ_only y SCZ_ISO
## SCZ_only
DATA <- c() 
for (j in 1:22)
	{
	#junto sumstats
	file<-fread(paste0("SCZ_only.chrom.",j), data.table=FALSE)
	DATA = rbind(DATA, file)
	}

temp<-subset(DATA, DATA$MAF <= .4 & DATA$MAF >= .1)
Neff<-mean(1/((2*temp$MAF*(1-temp$MAF))*temp$SE^2))
DATA$Neff = Neff
data <- DATA[,c("SNP","CHR","BP","MAF","A1","A2","est","SE","Z_Estimate","Pval_Estimate","chisq_pval","Neff")]
colnames(data) = c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","Z","P","chisq_P","Neff")
write.table(data, file = "SCZ_only_GSEM.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)

## pintamos manhattan
library(qqman)
manhattan(data, main = "Manhattan Plot", ylim = c(0, 38), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c(1:22))

# Guardamos el gráfico
ggsave("SCZ_only_Manhattan_Plot.png", height = 16, width = 12)


## SCZ_ISO
DATA <- c() 
for (j in 1:22)
	{
	#junto sumstats
	file<-fread(paste0("SCZ_ISO.chrom.",j), data.table=FALSE)
	DATA = rbind(DATA, file)
	}

temp<-subset(DATA, DATA$MAF <= .4 & DATA$MAF >= .1)
Neff<-mean(1/((2*temp$MAF*(1-temp$MAF))*temp$SE^2))
DATA$Neff = Neff
data <- DATA[,c("SNP","CHR","BP","MAF","A1","A2","est","SE","Z_Estimate","Pval_Estimate","chisq_pval","Neff")]
colnames(data) = c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","Z","P","chisq_P","Neff")
write.table(data, file = "SCZ_ISO_GSEM.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)

## pintamos manhattan
library(qqman)
manhattan(data, main = "Manhattan Plot", ylim = c(0, 38), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c(1:22))

# Guardamos el gráfico
ggsave("SCZ_ISO_Manhattan_Plot.png", height = 16, width = 12)

library(qqman)
manhattan(data, main = "Manhattan Plot", ylim = c(0, 38), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c(1:22))

# Guardamos el gráfico
ggsave("SCZ_ISO_Manhattan_Plot.png", height = 16, width = 12)





