# Script para cálculo de PRS empleando el método clásico de podado informativo de variantes y umbrales de asociación
# Enero 2022 - Patricia Fernández
# ==============================================
# Description -
#   El script utiliza los datos de genotipado de una muestra objetivo (cohorte de casos - controles EUGEI) que se quiere analizar (dataTarget) y los summary statistics de GWAS de meta-análisis de soledad/aislamiento  (dataDiscovery).
#
# ==============================================


library(base)

### A) LEEMOS LOS ARCHIVOS QUE VAMOS A UTILIZAR
# summary statistics (LNL_ISO)
setwd("/home/nodotea/PATRICIA/Objetivos/Datos/Discovery")
dataDiscovery <- read.table('MTAG_results.txt', header=TRUE, stringsAsFactors = FALSE)  ## 7583883 variantes

# Cohorte EUGEI WP6 (target sample)
setwd("/home/nodotea/PATRICIA/Objetivos/Datos/Target/WP6")
dataTarget <- read.table('EUGEI_WP6_Case_control_noMHC.bim', header = FALSE, stringsAsFactors = FALSE) ## 9340561 variantes



### B) FILTRAMOS DISCOVERY ###
head(dataDiscovery)

# Filtro de variantes imputadas a mala calidad
temp_discovery1<-dataDiscovery  ## 6938323 variantes

# Eliminamos cvariantes con strand ambiguity (las TA y CG). En este paso eliminamos tb indels. 
temp_discovery1$POS<-paste(temp_discovery1$CHR, temp_discovery1$BP_hg19, sep=":")
temp_discovery1$C<-paste(temp_discovery1$REF, temp_discovery1$ALT, sep="")
temp_discovery2<-subset(temp_discovery1, !(C %in% c("AT","TA", "CG", "GC"))) ## 7524527 variantes

# Creo un archivo con frecuencia de cada variante por si aparece repetida, y eliminamos repetidas
FRECUENCIAS_discovery<-as.data.frame(table(temp_discovery2$POS))
temp_discovery2_UNIQ<-temp_discovery2[temp_discovery2$POS %in% FRECUENCIAS_discovery$Var1[FRECUENCIAS_discovery$ALT_AF < 2], ]
temp_discovery2_UNIQ$POS<-sub("CHR", "", temp_discovery2_UNIQ$POS)
dataDiscovery <- temp_discovery2_UNIQ
head(dataDiscovery)              ## 7524527 variantes

dataDiscovery <- temp_discovery2
# Filtro columnas que no necesito. Me quedo con la informacion relevante, y renombro columnas
temp_discovery3<-dataDiscovery[,c("CHR", "BP_hg19", "REF", "ALT", "rsID", "BETA", "BETA_SE", "P_value", "POS")]
colnames(temp_discovery3)<-c("CHR", "BP", "A1", "A2", "SNP_ID", "BETA", "SE", "P", "SNP")
dataDiscovery <- temp_discovery3
rm(list=ls(pattern="^temp"))     ## Elimino objetis temporales

#borrar variante duplicada
dataDiscovery <- dataDiscovery[!duplicated(dataDiscovery$SNP), ]



### C) FILTRAMOS TARGET ##
head(dataTarget)

# Cambio nombre de columnas
colnames(dataTarget) <- c("V1", "varID", "V3", "pos", "A1", "A2")

# Creo un archivo con frecuencia de cada variante por si aparece repetida, y eliminamos repetidas
FRECUENCIAS_target<-as.data.frame(table(dataTarget$varID))
temp_target_UNIQ<-dataTarget[dataTarget$varID %in% FRECUENCIAS_target$Var1[FRECUENCIAS_target$Freq < 2], ]
dataTarget <- temp_target_UNIQ  ## 8220796 variantes
rm(list=ls(pattern="^temp"))         ## Elimino objetis temporales
rm(list=ls(pattern="^FRECUENCIAS"))  ## Elimino objetis temporales





### D) HACEMOS LA INTERSECCIÓN ENTRE LOS DOS ARCHIVOS PARA QUEDARNOS CON LAS VARIANTES PRESENTES EN AMBOS ARCHIVOS

# Me quedo con las variantes de un archivo presentes en el otro, y les cambio el nombre para distinguirlos
dataTarget_filtered <- dataTarget[which(dataTarget$varID %in% dataDiscovery$SNP),]                 ## 5541418       variantes
dataDiscovery_filtered <- dataDiscovery[which(dataDiscovery$SNP %in% dataTarget_filtered$varID),]  ## 5541418       variantes





### E) CLUNPING - Eliminar variantes en desequilibrio de ligamiento (LD; es decir, no independientes)

# Exportamos el archivo de discovery antes de hacer clumping
setwd("/home/nodotea/PATRICIA/Objetivos/Calculos_PRS")
write.table(dataDiscovery_filtered, "discovery_paraCLUMP.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Exportamos las variantes filtradas en la target
varIDs <- data.frame(dataTarget_filtered$varID)
# write.table(dataDiscovery_filtered, 'dataDiscovery_filtered.txt', row.names = FALSE, col.names = TRUE)
write.table(varIDs, 'variableIDs_forFiltering.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

# Creamos nuevos archivos de la target sample en formato PLINK .bed .bim and .fam solo con las variantes antes filtradas
#Eliminamos también las variantes dentro de la región del MHC (chr6:26Mb - 33Mb)
setwd("/home/nodotea/PATRICIA/Objetivos/Calculos_PRS")
system2("./plink", args=c("--bfile /home/nodotea/PATRICIA/Objetivos/Datos/Target/WP6/EUGEI_WP6_Case_control_noMHC", "--extract variableIDs_forFiltering.txt", "--exclude range MHC_filter.txt", "--make-bed", "--out EUGEI_WP6_Case_control_noMHC_new_filtered"))

# Hacemos CLUMPING usando  plink con el arechivo de la target ya filtrado, y usando como muestra discovery el exportado en los pasos anteriores
# Mediante el clumping vamos a quedarnos con las variantes independientes (no en LD - usando la propia target sample para ello) priorizando las más asociadas con SCZ, a partir de la información de la muestra discovery
# En un primer paso exportamos la lista de variantes resultante y posteriormente creamos un nuevo archivo de la target solo con esas variantes, que serán las únicas que se considerarán para el cálculo de PRS posterior
# Como parámetros de clumping usamos r2 < 0.1 y en una ventana de 500 Kb.
# IMPORTANTE: PLINK necesita reconocer el ID de la variante (columna SNP) y el valor de asociación de la muestra discovery (columna P). Aquí usamos como ID la columna chr:posición que creamos antes y es común a la discovery y a la target. Si no coinciden los IDs va a dar error
system2("./plink", args=c("--bfile EUGEI_WP6_Case_control_noMHC_new_filtered", "--clump discovery_paraCLUMP.txt", "--clump-p1 1", "--clump-p2 1", "--clump-kb 500", "--clump-r2 0.1", "--out EUGEI_WP6_Case_control_noMHC_filtered_CLUMPED"))
CLUMPED_PLINK<-read.table("EUGEI_WP6_Case_control_noMHC_filtered_CLUMPED.clumped", header=T)
write.table(CLUMPED_PLINK$SNP, "CLUMPED_PLINK.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
system2("./plink", args=c("--bfile EUGEI_WP6_Case_control_noMHC_new_filtered", "--extract CLUMPED_PLINK.txt", "--make-bed", "--out EUGEI_WP6_Case_control_noMHC_filtered_CLUMPED"))  ## 210691 variantes 






### F)  CALCULO DE PUNTUCIONES PRS  #####
head(dataDiscovery_filtered)
head(dataTarget_filtered)

### este paso no se hace en LNL_ISO
# IMPORTANTE! Si hay OR como variable hay que cambiara a Beta, usando logOR para convertirla a Beta
# dataDiscovery_filtered$BETA = log(dataDiscovery_filtered$OR)

setwd("/home/nodotea/PATRICIA/Objetivos/Calculos_PRS")
# colnames(dataDiscovery_filtered)<-c("SNP_ID", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "SNP")


# CReamos los archivos necesarios
# IMPORTANTE: Es necesario conocer cuál es el alelo sobre el que está calculado el efecto - si nos equivocamos veremos unos cálculos a la inversa. Consultar en el readme qué alelo es al que sde refiere en el OR
LNL_ISO_SCORE <- dataDiscovery_filtered[,c("SNP", "A1", "BETA")]



LNL_ISO_P <- dataDiscovery_filtered[,c("SNP", "P")]
write.table(LNL_ISO_SCORE, "LNL_ISO_SCORE.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(LNL_ISO_P, "LNL_ISO_P.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Calculamos PRS
# Empleamos el archivo SCORE, el P valor para seleccionar variantes y calcular puntuaciones
# Usamos un archivo de umbrales de P valor (Q_RANGES): En el definimos los umbrales que queremos emplear para distintos cálculos
setwd("/home/nodotea/PATRICIA/Objetivos/Calculos_PRS")
system2("./plink", args=c("--bfile EUGEI_WP6_Case_control_noMHC_filtered_CLUMPED", "--score LNL_ISO_SCORE.txt", "--q-score-file LNL_ISO_P.txt", "--q-score-range Q_RANGES.txt", "--out LNL_ISO_PRS_SCORE"))







### G) JUNTAR VCOVARIABLES DE ANÁLISIS Y PRS A DISTINTOS UMBRALES DE P VALOR.

# Calculamos mediante PCA las componentes de ancestralidad.
setwd("/home/nodotea/PATRICIA/Objetivos/Calculos_PRS")
system2("./plink", args=c("--bfile EUGEI_WP6_Case_control_noMHC_new_filtered", "--indep-pairwise 500 1 0.1", "--out PCA_PRUNED"))
system2("./plink", args=c("--bfile EUGEI_WP6_Case_control_noMHC_new_filtered", "--extract PCA_PRUNED.prune.in", "--make-bed", "--out EUGEI_WP6_Case_control_noMHC_filtered_PRUNED"))
system2("./plink", args=c("--bfile EUGEI_WP6_Case_control_noMHC_filtered_PRUNED", "--pca", "--out PCA"))

# Seleccionamos el archivo de eigenvectors (valores de PCA para las primeras 20 PC en los sujetos de la target sample) y exportamos las 10 primeras
PCA <- read.table('PCA.eigenvec', header=FALSE, stringsAsFactors = FALSE)  
colnames(PCA) = c("FID", "IID", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20")
PCA_final <-PCA[ ,c("FID", "IID", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10")]
write.table(PCA_final, "PCA.txt", sep="\t", quote=FALSE, row.names=FALSE)


# Usamos la información fenotípica en el archivo FAM. Y creamos un archivo con las covariables (Edad, sexo, componentes principales)
familyInfo <- read.table("EUGEI_WP6_Case_control_noMHC_new_filtered.fam", header=FALSE)
familyInfo <- data.frame(FID= familyInfo$V1, IID = familyInfo$V2, scode= familyInfo$V5, dcode = familyInfo$V6)

# scode = 1 (Males), scode = 2 (Females)
# dcode = 1 (Control),  dcode = 2 (Case), dcode = -9 (No Diagnosis Info)

head(familyInfo)

PCA       <- read.table("PCA.txt", header=T, stringsAsFactors = FALSE)
covars    <- merge(PCA, familyInfo, by = c("IID", "FID"))
age       <- read.table("age.txt", header=T, stringsAsFactors = FALSE)
covars    <- merge(covars, age, by = c("IID", "FID"))
country   <- read.table("country.txt", header=T, stringsAsFactors = TRUE)
covars    <- merge(covars, country, by = c("IID", "FID"))
write.table(covars, "covariables_LNL_ISO.txt", row.names = FALSE, col.names = TRUE, dec = ".")

covars <- read.table("covariables_LNL_ISO.txt", header = TRUE, stringsAsFactors = FALSE)
covars_filter <- covars[covars$dcode %in% c(1, 2),] # Los pheno -9 se quitan. 
variables_LNL_ISO <- covars_filter

for (j in 1:12)  # Pongo 12 porque he generado 12 PRS para los 12 umbrales del archivo Q_RANGES
	{
  	fin <- read.table(paste0("LNL_ISO_PRS_SCORE.S",j,".profile"), header = TRUE, stringsAsFactors =  FALSE)
  	score <- data.frame(FID= fin$FID, IID= fin$IID, fin$SCORE)
  	colnames(score)[3] <- paste0("SCORE_S",j)
  	variables_LNL_ISO <- merge(variables_LNL_ISO, score, by = c("IID", "FID"))
	}

write.table(variables_LNL_ISO, "LNL_ISO.txt", dec= ".", row.names = FALSE)  #### GUARDAR A RESULTADOS ####






### H) REGRESIÓN LOGÍSTICA PARA EVALULAR PRS 

library(rms)
library(pROC)
variables_LNL_ISO <- read.table('LNL_ISO.txt', header=TRUE, stringsAsFactors = FALSE)  

# Definimos variables de regresión logística. Metemos las 10 primeras PCs, sexo y edad, además de PRS.
pheno <- as.numeric(variables_LNL_ISO$dcode)
PC1  <- variables_LNL_ISO$C1
PC2  <- variables_LNL_ISO$C2
PC3  <- variables_LNL_ISO$C3
PC4  <- variables_LNL_ISO$C4
PC5  <- variables_LNL_ISO$C5
PC6  <- variables_LNL_ISO$C6
PC7  <- variables_LNL_ISO$C7
PC8  <- variables_LNL_ISO$C8
PC9  <- variables_LNL_ISO$C9
PC10 <- variables_LNL_ISO$C10
SEX   <- variables_LNL_ISO$scode
AGE   <- variables_LNL_ISO$AGE
COUNTRY   <- variables_LNL_ISO$COUNTRY


# Comparamos un modelo de regresión solo incluyendo covariables frente a otro que incluye además el PRS. De esa comparación extraemos la varianza explicada (pseudoR2) por la varible PRS
# El P valor de la varible PRS la extraemos del modelo logístico (glm)
# Hacemos un bucle para la comparación en cada uno de los umbrales de significación seleccionados. habrá un valor de R2 y P para cada umbral. Y además, en cada permutación
# Sacamos el valor de AUC tb, pero en ese caso usamos un modelo sin covaribles

H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + COUNTRY, scale = TRUE)
# H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX , scale = TRUE) #si no quisiera controlar por edad
print(H0)

dataRanges <- c() 
for (j in 1:12) 
	{
  	score <- variables_LNL_ISO[,grep(paste0("^SCORE_S",j,"$"), colnames(variables_LNL_ISO))]  # tenemos que ponerle esta sintaxis para que coja score1 y no score10,11,12 a la vez
  	H1 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score + COUNTRY, scale = TRUE)
 	DELTA_R2 <- H1$stats[10] - H0$stats[10]
 	R2 <- DELTA_R2 * 100  
  
  	phenoN <- pheno - 1 #  cmabiamos de 2/1 a 1/0 para el modelo de regresión
  	model <- glm(phenoN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score + COUNTRY, family = binomial())
 	P <- summary(model)

 	AUC = roc(pheno, score)
 	
  	# Guardamos los rsultados de la regresión en el dataframe
  	df <- data.frame(R2 = R2, P = P$coefficients[14,4], AUC = AUC$auc, Range = paste0( "Q_RANGE_S", j))
  	dataRanges <- rbind(dataRanges, df)
	}

data_LNL_ISO_lrm = dataRanges # Mejor umbral S9 P < 0.1

write.table(data_LNL_ISO_lrm, "data_LNL_ISO_lrm.txt", dec= ",", row.names = FALSE)  #### GUARDAR A RESULTADOS ####


# Le añadimos el valor de P de la comparación inicial mediante glm si bootstrap.
table2 <- read.table("data_LNL_ISO_lrm.txt", header=T, stringsAsFactors = FALSE, dec = ",")
threshold <- read.table("Q_RANGES.txt", header=FALSE, stringsAsFactors = FALSE, dec = ".")

TABLA_FINAL = cbind(threshold[,3],table2)
names(TABLA_FINAL)[names(TABLA_FINAL) == "threshold[, 3]"] <- 'Threshold' # renombreamos la columna que acabamos de añadir

write.table(TABLA_FINAL, "FINAL_data_LNL_ISO_lrm.txt", dec= ",", row.names = FALSE)





### K) RESULTADOS DE LA COMPARACIÓN POR DECILES.
## Para ello, asignamos un decil a cada sujeto dentro del umbral de significación de mayor R2 (en este caso P < v0.1).
## Después hacemos una comparación de la probabilidad de ser un caso con LNL_ISO dentro de cada decilr respecto al medio (5º) o respecto al decil más bajo (1º)

library(base)
library(rms)
library(dplyr)

## Cargamos la tabla con todos los scores y covariables, y creamos la varible decil asignando a cada sujeto un valor del 1..10 en función de su score en P<0.1
LNL_ISO <- read.table('LNL_ISO.txt', header=TRUE, stringsAsFactors = FALSE)  
LNL_ISO$decile <- ntile(LNL_ISO$SCORE_S3, 4) ### Asignar el umbral más significativo


## COMPARACIÓN FRRENTE AL DECIL 5 (PROMEDIO)
## Hacemos una subselección de sujetos en base a su decil, inlcuyendo siempre el 5º para hacer la comparación
y <- c(1,2,3,4,6,7,8,9,10)
for(i in y)
	{
  	variables_LNL_ISO <- LNL_ISO[ which(LNL_ISO$decile == 5 | LNL_ISO$decile == i), ]
 	write.table(variables_LNL_ISO, paste0("data",i), dec= ",", row.names = FALSE)
	}

## Cargamos cada archivo y calculamos el OR de la comparación de porbabilidad de ser caso en el decil seeleccionado en comparado con el decil 5
dataRanges <- c() 
for(i in y)
	{
  	variables_LNL_ISO <- read.table(paste0("data",i), header = TRUE, stringsAsFactors = FALSE, dec = ",")
  	variables_LNL_ISO$comparison <- ifelse(variables_LNL_ISO$decile == 5, 0, 1)

  	# Definimos variables de regresión logística. Metemos las 10 primeras PCs, sexo y edad. Metemos la variable PRS del siguiente modo:
  	# decil X == 1
  	# decil 5 == 0

	pheno <- as.numeric(variables_LNL_ISO$dcode)
	PC1  <- variables_LNL_ISO$C1
	PC2  <- variables_LNL_ISO$C2
	PC3  <- variables_LNL_ISO$C3
	PC4  <- variables_LNL_ISO$C4
	PC5  <- variables_LNL_ISO$C5
	PC6  <- variables_LNL_ISO$C6
	PC7  <- variables_LNL_ISO$C7
	PC8  <- variables_LNL_ISO$C8
	PC9  <- variables_LNL_ISO$C9
	PC10 <- variables_LNL_ISO$C10
	SEX   <- variables_LNL_ISO$scode
	AGE   <- variables_LNL_ISO$AGE
	COUNTRY <- variables_LNL_ISO$COUNTRY
 	SCORE <- variables_LNL_ISO$comparison

  	phenoN <- pheno - 1 #Coded for glm
  	model <- glm(phenoN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + SCORE + COUNTRY, family = binomial())
  	P <- summary(model)

  	coef <- summary(model)$coeff
  	OR = exp(coef[14,1])
  	CI_inf = exp(coef[14,1] - 1.96*coef[14,2])
  	CI_sup = exp(coef[14,1] + 1.96*coef[14,2])

  	# Guardamos en archivo y exportamos
  	df <- data.frame(OR = OR, CI_inf = CI_inf, CI_sup = CI_sup, P = P$coefficients[14,4], Decile = paste0( "Decile_", i))
  	dataRanges <- rbind(dataRanges, df)
  	}
data_decils_lgm = dataRanges
write.table(data_decils_lgm, "data_LNL_ISO_lrm_DECIL_VS_5th.txt", dec= ",", row.names = FALSE)  #### GUARDAR A RESULTADOS ####




## COMPARACIÓN FRRENTE AL DECIL UNO (MENOR RIESGO POSIBLE)
## Hacemos una subselección de sujetos en base a su decil, inlcuyendo siempre el 5º para hacer la comparación
y <- c(2,3,4,5,6,7,8,9,10)
for(i in y)
	{
  	variables_LNL_ISO <- LNL_ISO[ which(LNL_ISO$decile == 1 | LNL_ISO$decile == i), ]
 	write.table(variables_LNL_ISO, paste0("data",i), dec= ",", row.names = FALSE)
	}

## Cargamos cada archivo y calculamos el OR de la comparación de porbabilidad de ser caso en el decil seeleccionado en comparado con el decil 1
dataRanges <- c() 
for(i in y)
	{
  	variables_LNL_ISO <- read.table(paste0("data",i), header = TRUE, stringsAsFactors = FALSE, dec = ",")
  	variables_LNL_ISO$comparison <- ifelse(variables_LNL_ISO$decile == 1, 0, 1)

  	# Definimos variables de regresión logística. Metemos las 10 primeras PCs, sexo y edad. Metemos la variable PRS del siguiente modo:
  	# decil X == 1
  	# decil 1 == 0
  	
	pheno <- as.numeric(variables_LNL_ISO$dcode)
	PC1  <- variables_LNL_ISO$C1
	PC2  <- variables_LNL_ISO$C2
	PC3  <- variables_LNL_ISO$C3
	PC4  <- variables_LNL_ISO$C4
	PC5  <- variables_LNL_ISO$C5
	PC6  <- variables_LNL_ISO$C6
	PC7  <- variables_LNL_ISO$C7
	PC8  <- variables_LNL_ISO$C8
	PC9  <- variables_LNL_ISO$C9
	PC10 <- variables_LNL_ISO$C10
	SEX   <- variables_LNL_ISO$scode
	AGE   <- variables_LNL_ISO$AGE
	COUNTRY <- variables_LNL_ISO$COUNTRY
 	SCORE <- variables_LNL_ISO$comparison

  	phenoN <- pheno - 1 #Coded for glm
  	model <- glm(phenoN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + SCORE + COUNTRY, family = binomial())
  	P <- summary(model)

  	coef <- summary(model)$coeff
  	OR = exp(coef[14,1])
  	CI_inf = exp(coef[14,1] - 1.96*coef[14,2])
  	CI_sup = exp(coef[14,1] + 1.96*coef[14,2])

  	# Guardamos en archivo y exportamos
  	df <- data.frame(OR = OR, CI_inf = CI_inf, CI_sup = CI_sup, P = P$coefficients[14,4], Decile = paste0( "Decile_", i))
  	dataRanges <- rbind(dataRanges, df)
  	}
data_decils_lgm = dataRanges
write.table(data_decils_lgm, "data_LNL_ISO_lrm_DECIL_VS_1st.txt", dec= ",", row.names = FALSE)   #### GUARDAR A RESULTADOS ####









### L) GRÁFICOS DE LOS ANÁLISIS DE VARIANZA Y COMPARCIÓN POR DECILES

library(ggplot2)

####################################
## GRÁFICOS DE VARIANZA EXPLICADA ##
####################################

## GRÁFICO PARA PSEUDO R2 DE NAGGELKERKE.

# Cargamos el archivo final con todos los valores para exportar
prs.result <- read.table("FINAL_data_LNL_ISO_lrm.txt", header=T, stringsAsFactors = FALSE, dec = ",")

# Generamos un buen formato para lo P-valores del modeo, y ponerlos sobre el gráfico de R2. Redondeamos a 3 dígitos, menos cuando sea exponencial que sustituímos la e por x10
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) & prs.result$print.p == 0] <- format(prs.result$P[!is.na(prs.result$print.p) & prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)


## GRÁFICO PARA R2 EN LIABILITY SCALE
# cargamos el gráfico en ggplot. 
# Cambiamos el valor de threshold a FACTOR para que se pueda mantener el orden creciente de P valor.
ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +

    # Especificamos que queremos el valor de P-value sobre el límite superior de la errorbar, e inclinado
    geom_text(aes(label = paste(print.p), y = R2), vjust = -0.5, hjust = 0, angle = 45, cex = 4, parse = T) +

    # Especificamos el rango de la escala Y como 1.25 por el máximo de R2_liab, para dejar suficiente espacio para barras y p valores.
    scale_y_continuous(limits = c(-0.10, max(prs.result$R2) * 3)) +

    # Especificamos los labels de los ejes X e Y en itálica.
    xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
    ylab(expression(paste("PRS model fit: pseudo- ", R ^ 2))) +

    # configuramos la barplot que es el tipo de gráfico que queremos. Configuramos en fill para que el clor sea funció de su grado de asociación.
    geom_bar(aes(fill = -log10(P)), stat = "identity") +

    # Especificamos los colores de las barras en funciòn del grado de asociación (según pusimos arriba en fill een barplot). Ponemos un midpoint según el grado de asociación que tengamos. Y la leyenda relativa
    scale_fill_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-10,
        name = bquote(atop(-log[10] ~ model, italic(P) - value),)) +

    # Cambios estéticos del gráfico según queramos
    theme_classic() + theme(
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Guardamos el gráfico
ggsave("R2_NAG_PLOT.png", height = 16, width = 12)





#########################################
## GRÁFICOS DE COMPARACIÓN POR DECILES ##
#########################################

###  COMPARACIONES FRENTE AL 5 DECIL

table2=read.table("data_LNL_ISO_lrm_DECIL_VS_5th.txt", header = T, dec = ",", stringsAsFactors = FALSE)
table2$Decile <- sub("Decile_", "", table2$Decile)
table2$Decile <- factor(table2$Decile, levels = table2$Decile)


# Generamos un buen formato para lo P-valores del modeo, y ponerlos sobre el gráfico de R2. Redondeamos a 3 dígitos, menos cuando sea exponencial que sustituímos la e por x10
table2$print.p <- round(table2$P, digits = 3)
table2$print.p[!is.na(table2$print.p) & table2$print.p == 0] <- format(table2$P[!is.na(table2$print.p) & table2$print.p == 0], digits = 2)
table2$print.p <- sub("e", "*x*10^", table2$print.p)

# Cambiamos el valor de decil a FACTOR para que se pueda mantener el orden creciente de P valor.
ggplot(data = table2, aes(x = factor(Decile), y = OR, fill = -log10(P), color = -log10(P))) +

    # Especificamos que queremos el valor de P-value sobre el límite superior de la errorbar, e inclinado
    geom_text(aes(label = paste(print.p), y = CI_sup), vjust = -1.5, hjust = 0, angle = 45, cex = 4, parse = T, color = "black") +

	# Especificamos el rango de la escala Y como 1.6 por el máximo de R2, para dejar suficiente espacio para barras y p valores (esto puede tener que retocarse).
    scale_y_continuous(limits = c(0, max(table2$OR) * 1.6)) +

	# Especificamos los labels de los ejes X e Y.
    xlab(expression(paste("PRS decile"))) +
    ylab(expression(paste("OR (comparison against 5th decile)"))) +

	# configuramos la errorbarplot especificando su posiciómn en eje X e Y, así como formato.
	geom_errorbar( aes( x=factor(Decile), ymin=CI_inf, ymax=CI_sup, color = -log10(P)), width=0.4, alpha=0.9, size=1.3) + geom_point() +

	# Añadimos una linea horizontal indicando el valor neutro (OR = 1)
	geom_hline(yintercept=1, linetype="dashed", color = "grey", size=0.8) +

	# Especificamos los colores de las barras en funciòn del grado de asociación (según pusimos arriba en fill een barplot). Ponemos un midpoint según el grado de asociación que tengamos. Y la leyenda relativa
    scale_colour_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-3,
        name = bquote(atop(-log[10] ~ model, italic(P) - value),)) +

	# Le añado colores a los puntos tb del mismo modo que las barras, pero le quito la leyenda para que no aparezcan las dos (guide = "none")
        scale_fill_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-3,
        guide = "none") +

    # Cambios estéticos del gráfico según queramos
    theme_classic() + theme(
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Guardamos el gráfico
ggsave("DECILES_COMP_5th.png", height = 10, width = 12)




###  COMPARACIONES FRENTE AL PRIMER DECIL

table2=read.table("data_LNL_ISO_lrm_DECIL_VS_1st.txt", header = T, dec = ",", stringsAsFactors = FALSE)
table2$Decile <- sub("Decile_", "", table2$Decile)
table2$Decile <- factor(table2$Decile, levels = table2$Decile)


# Generamos un buen formato para lo P-valores del modeo, y ponerlos sobre el gráfico de R2. Redondeamos a 3 dígitos, menos cuando sea exponencial que sustituímos la e por x10
table2$print.p <- round(table2$P, digits = 3)
table2$print.p[!is.na(table2$print.p) & table2$print.p == 0] <- format(table2$P[!is.na(table2$print.p) & table2$print.p == 0], digits = 2)
table2$print.p <- sub("e", "*x*10^", table2$print.p)

# Cambiamos el valor de decil a FACTOR para que se pueda mantener el orden creciente de P valor.
ggplot(data = table2, aes(x = factor(Decile), y = OR, fill = -log10(P), color = -log10(P))) +

    # Especificamos que queremos el valor de P-value sobre el límite superior de la errorbar, e inclinado
    geom_text(aes(label = paste(print.p), y = CI_sup), vjust = -1.5, hjust = 0, angle = 45, cex = 4, parse = T, color = "black") +

	# Especificamos el rango de la escala Y como 1.7 por el máximo de R2, para dejar suficiente espacio para barras y p valores (esto puede tener que retocarse).
    scale_y_continuous(limits = c(0, max(table2$OR) * 1.7)) +

	# Especificamos los labels de los ejes X e Y.
    xlab(expression(paste("PRS decile"))) +
    ylab(expression(paste("OR (comparison against 1st decile)"))) +

	# configuramos la errorbarplot especificando su posiciómn en eje X e Y, así como formato.
	geom_errorbar( aes( x=factor(Decile), ymin=CI_inf, ymax=CI_sup, color = -log10(P)), width=0.4, alpha=0.9, size=1.3) + geom_point() +

	# Añadimos una linea horizontal indicando el valor neutro (OR = 1)
	geom_hline(yintercept=1, linetype="dashed", color = "grey", size=0.8) +

	# Especificamos los colores de las barras en funciòn del grado de asociación (según pusimos arriba en fill een barplot). Ponemos un midpoint según el grado de asociación que tengamos. Y la leyenda relativa
    scale_colour_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-3,
        name = bquote(atop(-log[10] ~ model, italic(P) - value),)) +

	# Le añado colores a los puntos tb del mismo modo que las barras, pero le quito la leyenda para que no aparezcan las dos (guide = "none")
        scale_fill_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-3,
        guide = "none") +

    # Cambios estéticos del gráfico según queramos
    theme_classic() + theme(
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Guardamos el gráfico
ggsave("DECILES_COMP_1st.png", height = 10, width = 12)








