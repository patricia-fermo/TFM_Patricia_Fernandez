#Script para análisis gen-ambiente de aislamiento (WP6) con variables ambientales.
#Marzo 2023 - Patricia Fernández

library(rms)
library(pROC)
library("dplyr")

setwd("/home/nodotea/PATRICIA/Nuevos_objetivos/Objetivo2")

## leemos la info de genética y ambiente, y juntamos, y guardamos
variables_ISO <- read.table('ISO.txt', header=TRUE, stringsAsFactors = FALSE) 
variables_ISO <- variables_ISO %>% select(-c(SCORE_S1, SCORE_S2, SCORE_S3, SCORE_S4, SCORE_S5, SCORE_S6, SCORE_S7, SCORE_S9, SCORE_S10, SCORE_S11, SCORE_S12))

variables_WP6 <- read.table('WP6_env_variables', header=TRUE, stringsAsFactors = FALSE)
variables_WP6 <- variables_WP6 %>% select(c(Beadchip_Location, emot_abu_bin1, phys_abu_bin1, sex_abu_bin1, emot_negl_bin1, phys_negl_bin1, bully_non_rest, winter_birth, hear_imp, canreg))

data = merge(variables_ISO, variables_WP6, by.x = "IID", by.y = "Beadchip_Location")

write.table(data, "ISO_GxE.txt", sep="\t", quote=FALSE, row.names=FALSE)  #### GUARDAR A RESULTADOS ####



## creamos la varaible genética de forma dicotomizada ( > 75%), tal cual hace Guloksuz (2019).
## En primer lugar, seleccionamos la población control, y miramos el corte del 75% de PRS, luego recodificamos en el total de individuos
data <- read.table('ISO_GxE.txt', header=TRUE, stringsAsFactors = FALSE) 
control_data = subset(data, dcode = 1)
cutoff = quantile(control_data$SCORE_S8, 0.75)
data$PGS = ifelse(data$SCORE_S8 < cutoff, 0, 1)

## cambio phenotipo pra logistic model
data$pheno <- ifelse(data$dcode == 1, 0, 1)



## Evaluamos para las variables ambientales diferentes:
library(lme4)
dataRanges <- c() 
for (j in c("emot_abu_bin1","phys_abu_bin1","sex_abu_bin1","emot_negl_bin1","phys_negl_bin1","bully_non_rest","winter_birth","hear_imp","canreg")) 
	{
	# Create three dummy variables based on "genetic" and "environmental" variables
	data$GE1 <- ifelse(data$PGS == 1 & data[,j] == 0, 1, 0)
	data$GE2 <- ifelse(data$PGS == 0 & data[,j] == 1, 1, 0)
	data$GE3 <- ifelse(data$PGS == 1 & data[,j] == 1, 1, 0)

	# Fit a multilevel logistic regression model with genetic and environmental variables
	model <- glmer(pheno ~ scale(C1) + scale(C2) + scale(C3) + scale(C4) + scale(C5) + scale(C6) + scale(C7) + scale(C8) + scale(C9) + scale(C10) + scode + scale(AGE) + GE1 + GE2 + GE3 + (1 | COUNTRY), data = data, family = binomial())
	P <- summary(model)
	OR_PGS = exp(P$coefficients[14,1])
	OR_Env = exp(P$coefficients[15,1])
	OR_PGS_x_Env = exp(P$coefficients[16,1])
	RERI = (exp(P$coefficients[16,1]) - exp(P$coefficients[14,1]) - exp(P$coefficients[15,1]) + 1)
  	
  	# Guardamos los rsultados de la regresión en el dataframe
  	df <- data.frame(environment = paste0( "variable ", j), RERI = RERI, OR_PGS = OR_PGS, OR_Env = OR_Env, OR_PGS_x_Env = OR_PGS_x_Env)
  	dataRanges <- rbind(dataRanges, df)
	}
write.table(dataRanges, "RERI_results.txt", sep="\t", quote=FALSE, row.names=FALSE)  #### GUARDAR A RESULTADOS ####






## CALCULAMOS RERIS POR PERMUTACIONES MEDIANTE BOOTSTRAP - HOUSEMADE

control_data = subset(data, dcode = 1)
cutoff = quantile(control_data$SCORE_S8, 0.75)
data$PGS = ifelse(data$SCORE_S8 < cutoff, 0, 1)

## cambio phenotipo pra logistic model
data$pheno <- ifelse(data$dcode == 1, 0, 1)

## Evaluamos para las variables ambientales diferentes:
library(lme4)

## Establecemos el número de simulaciones que queremos y hacemos el loop
n.simul <- 1000               

dataRanges <- c() 
# Seleccionamos los 1100 casos y 1212 controles pero al azar con remplazamiento (bootstrap). Por tanto, se repetirán algunos y faltarán otros en cada permutación. Sirve como análisis de sensibilidad
for(i in 1:n.simul) 
	{
	cases = subset(data, dcode == 2)
	controls = subset(data, dcode == 1)
	pseudocases = cases[sample(nrow(cases),size = 1100 , replace = TRUE),]
	pseudocontrols = controls[sample(nrow(controls),size = 1212, replace = TRUE),]
	table = rbind(pseudocases,pseudocontrols)


	for (j in c("emot_abu_bin1","phys_abu_bin1","sex_abu_bin1","emot_negl_bin1","phys_negl_bin1","bully_non_rest","winter_birth","hear_imp","canreg")) 
		{
		# Create three dummy variables based on "genetic" and "environmental" variables
		table$GE1 <- ifelse(table$PGS == 1 & table[,j] == 0, 1, 0)
		table$GE2 <- ifelse(table$PGS == 0 & table[,j] == 1, 1, 0)
		table$GE3 <- ifelse(table$PGS == 1 & table[,j] == 1, 1, 0)

		# Fit a multilevel logistic regression model with genetic and environmental variables
		model <- glmer(pheno ~ scale(C1) + scale(C2) + scale(C3) + scale(C4) + scale(C5) + scale(C6) + scale(C7) + scale(C8) + scale(C9) + scale(C10) + scode + scale(AGE) + GE1 + GE2 + GE3 + (1 | COUNTRY), data = table, family = binomial())
		P <- summary(model)
		OR_PGS = exp(P$coefficients[14,1])
		OR_Env = exp(P$coefficients[15,1])
		OR_PGS_x_Env = exp(P$coefficients[16,1])
		RERI = (exp(P$coefficients[16,1]) - exp(P$coefficients[14,1]) - exp(P$coefficients[15,1]) + 1)
  	
  		# Guardamos los rsultados de la regresión en el dataframe
  		df <- data.frame(environment = paste0( "variable ",j), RERI = RERI, OR_PGS = OR_PGS, OR_Env = OR_Env, OR_PGS_x_Env = OR_PGS_x_Env, perm = paste0("permutacion ",i))
  		dataRanges <- rbind(dataRanges, df)
		}
	}
write.table(dataRanges, "RERI_results_PERMS_BOOTSTRAP.txt", sep="\t", quote=FALSE, row.names=FALSE)  #### GUARDAR A RESULTADOS ####







## CALCULAMOS INTERVALOS DE CONFIANZA y P VALUES MEDIANTE BOOTSTRAP PERCENTILE METHOD

# leemos los datos generados anteriormente
table <- read.table("RERI_results_PERMS_BOOTSTRAP.txt", sep="\t", header=T, stringsAsFactors = FALSE)
table_original <- read.table("RERI_results.txt", sep="\t", header=T, stringsAsFactors = FALSE)

datos_nuevos <- c() 
for (j in c("emot_abu_bin1","phys_abu_bin1","sex_abu_bin1","emot_negl_bin1","phys_negl_bin1","bully_non_rest","winter_birth","hear_imp","canreg")) 
	{
 	VARIABLE <- table[table$environment == paste0("variable ",j),]
 	VARIABLE_original <- table_original[table_original$environment == paste0("variable ",j),]


 	# OR_PGS
 	OR_PGS = VARIABLE_original$OR_PGS
 	mean_OR = mean(VARIABLE$OR_PGS)

 	# intervalos de confinaza estimados por el bootstrap percentile method
 	IC_OR_PGS <- quantile(VARIABLE$OR_PGS, c(0.025, 0.975))
 	
	## Calcular el p-valor como la proporción de veces en las que la estadística de prueba (en este caso OR) es superior al valor deseado (en este caso 1) en las muestras generadas por remuestreo
	P_PGS <- sum(VARIABLE$OR_PGS < 1)/(1000 + 1)



 	# OR_ENV
 	OR_Env = VARIABLE_original$OR_Env
 	mean_OR = mean(VARIABLE$OR_Env)

 	# intervalos de confinaza estimados por el bootstrap percentile method
 	IC_OR_Env <- quantile(VARIABLE$OR_Env, c(0.025, 0.975))
 	
	## Calcular el p-valor como la proporción de veces en las que la estadística de prueba (en este caso OR) es superior al valor deseado (en este caso 1) en las muestras generadas por remuestreo
	P_ENV <- sum(VARIABLE$OR_Env < 1)/(1000 + 1)



 	# OR_PGS_x_Env
 	OR_PGS_x_Env = VARIABLE_original$OR_PGS_x_Env
 	mean_OR = mean(VARIABLE$OR_PGS_x_Env)

 	# intervalos de confinaza estimados por el bootstrap percentile method
 	IC_OR_PGS_x_Env <- quantile(VARIABLE$OR_PGS_x_Env, c(0.025, 0.975))
 	
 	## Calcular el p-valor como la proporción de veces en las que la estadística de prueba (en este caso OR) es superior al valor deseado (en este caso 1) en las muestras generadas por remuestreo
	P_PGS_x_Env <- sum(VARIABLE$OR_PGS_x_Env < 1)/(1000 + 1)



 	# RERI
 	RERI = VARIABLE_original$RERI
 	mean_RERI = mean(VARIABLE$RERI)

 	# intervalos de confinaza estimados por el bootstrap percentile method
 	IC_RERI <- quantile(VARIABLE$RERI, c(0.025, 0.975))
 	
 	## Calcular el p-valor como la proporción de veces en las que la estadística de prueba (en este caso RERI) es superior al valor deseado (en este caso 0) en las muestras generadas por remuestreo
	P_RERI <- sum(VARIABLE$RERI < 0)/(1000 + 1)



	## Juntamos todo en una tabla
  	df <- data.frame(environment = paste0( "variable ",j), OR_PGS = OR_PGS, OR_PGS_INF = IC_OR_PGS[1], OR_PGS_SUP = IC_OR_PGS[2], P_PGS = P_PGS, OR_Env = OR_Env, OR_Env_INF = IC_OR_Env[1], OR_Env_SUP = IC_OR_Env[2], P_ENV = P_ENV, OR_PGS_x_Env = OR_PGS_x_Env, OR_PGS_x_Env_INF = IC_OR_PGS_x_Env[1], OR_PGS_x_Env_SUP = IC_OR_PGS_x_Env[2], P_PGS_x_Env = P_PGS_x_Env, RERI = RERI, RERI_INF = IC_RERI[1], RERI_SUP = IC_RERI[2], P_RERI = P_RERI)
  	datos_nuevos = rbind(datos_nuevos, df)
	}

write.table(datos_nuevos, "RERI_results_CI_P_BOOTSTRAP.txt", sep="\t", quote=FALSE, row.names=FALSE, dec = ",")  #### GUARDAR A RESULTADOS ####










