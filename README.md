# Analysis of the genetic overlap of loneliness and schizophrenia

This repository includes all the scripts used to perform the analyses in the Master thesis, divided by the three main objectives. All code is written in R programming language and run in Ubuntu operative system. 

# Table of contents
1. [Abstract](#abstract)
2. [Objective 1](#objective1)
3. [Objective 2](#objective2)
4. [Objective 3](#objective3)
5. [References](#references)

## Abstract <a name="abstract"></a>

Schizophrenia is a polygenic and heritable severe mental disorder. It has been associated with social isolation and loneliness in previous research by a bidirectional causal relationship. Loneliness is a known cause for the onset and worsening of mental disorders as well as a risk factor for mortality. 

In this study, a multidisciplinary approach is taken to study the association between loneliness and psychosis with the use of genome-wide association studies (GWAS) data. Three main bioinformatics analyses are performed including polygenic risk score estimation, genomic structural equation modelling and gene-environment interaction analysis. The objective is to replicate the association found in previous literature in a European cohort. We further advance this objective by studying the potential overlap at different stages of psychosis, including schizophrenia and first episode psychosis data. We also perform a genomic dissection methodology and a study of the interaction of several environmental variables with schizophrenia subjected to isolation. 

We find a statistically significant genetic overlap between social isolation and schizophrenia, indicating an association with later stages of psychosis. We as well find a gene-environment interaction between cannabis use and schizophrenia subjected to isolation. These results illustrate the application of bioinformatics tools in translational medicine research and their possible use for innovation in the clinical practice in the form of adequate socialization tools in the treatment of psychosis patients.

## Objetive 1 <a name="objective1"></a>

* The script for PRS calculation is divided in three subfolders depending on the dataset analysed. WP6 includes the analysis of schizophrenia case/control samples, WP2 for first episode psychosis samples and GSEM for the mathematically-derived GWAS files from the Genomic SEM analysis (objective 2). Files `SCRIPT_PRS.R` include the code followed for all nine steps of polygenic risk score analysis. Four different phenotypes were analysed for WP6 and WP2 (SCZ, LNL, ISO and LNL-ISO). Two phenotypes were analysed for GSEM (SCZ_ISO and SCZ_only).

## Objetive 2 <a name="objective2"></a>

* `script_GSEM.R` contains all the code used for the Genomic SEM analysis following the methodology described by Grotzinger et al. (2019).

## Objetive 3 <a name="objective3"></a>

* Files `script_GxE.R` include the gene-environment analysis performed of three genetic variables: ISO, SCZ_ISO and SCZ_only. The methodology by Guloksuz et al. (2019) was followed for the analysis. 

## References <a name="references"></a>

* Grotzinger, A. D., Rhemtulla, M., de Vlaming, R., Ritchie, S. J., Mallard, T. T., Hill, W. D., Ip, H. F., Marioni, R. E., McIntosh, A. M., Deary, I. J., Koellinger, P. D., Harden, K. P., Nivard, M. G., & Tucker-Drob, E. M. (2019). Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits. Nature human behaviour, 3(5), 513–525. https://doi.org/10.1038/s41562-019-0566-x

* Guloksuz, S., Pries, L. K., Delespaul, P., Kenis, G., Luykx, J. J., Lin, B. D., Richards, A. L., Akdede, B., Binbay, T., Altınyazar, V., Yalınçetin, B., Gümüş-Akay, G., Cihan, B., Soygür, H., Ulaş, H., Cankurtaran, E., Kaymak, S. U., Mihaljevic, M. M., Petrovic, S. A., Mirjanic, T., … van Os, J. (2019). Examining the independent and joint effects of molecular genetic liability and environmental exposures in schizophrenia: results from the EUGEI study. World psychiatry : official journal of the World Psychiatric Association (WPA), 18(2), 173–182. https://doi.org/10.1002/wps.20629


