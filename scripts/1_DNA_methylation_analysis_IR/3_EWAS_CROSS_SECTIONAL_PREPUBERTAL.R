


####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

"Chunk 0: Installing OR loading required packages"

####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################



library(limma)
library(RColorBrewer)
library(plyr)
library(ggplot2)
library(reshape2)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #Specific for the 850k
library(IlluminaHumanMethylationEPICmanifest) #Specific for the 850k
library(DMRcate)
require(Gviz)
require(missMethyl)
require(GOplot)
library(RColorBrewer)
library(ggplot2)



####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

"Chunk 1: Loading Datasets and formating data"

####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################



#CARGA DE DATOS DE METILACION:
#mVals_QUANTILE <- read.csv2("/home/augusto/Escritorio/ANALISIS_EWAS/PRIMER_ANALISIS/BASES_DATOS_METILACION_GENERADAS_PREPROCESADAS/mVals_QUANTILE.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1)

#bVals_QUANTILE  <- read.csv2("/home/augusto/Escritorio/ANALISIS_EWAS/PRIMER_ANALISIS/BASES_DATOS_METILACION_GENERADAS_PREPROCESADAS/bVals_QUANTILE.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1)

mVals_BMIQ <- read.csv2("/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/BASES_DATOS_METILACION_GENERADAS_PREPROCESADAS/mVals_BMIQ_NEW.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1)

#bVals_BMIQ <- read.csv2("/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/BASES_DATOS_METILACION_GENERADAS_PREPROCESADAS/bVals_BMIQ.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1)

#para hacerlo solo con S100A4
#mVals_BMIQ <- read.csv2("/home/augusto/Escritorio/ANALISIS_EWAS/PRIMER_ANALISIS/S100A4/MVALS_S100A4.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1)
#bVals_BMIQ <- read.csv2("/home/augusto/Escritorio/ANALISIS_EWAS/PRIMER_ANALISIS/S100A4/BVALS_S100A4.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1)

########## CARGA DE COVARIABLES DEL DATASET:



dataDirectorys <- "/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/"
targets <- read.metharray.sheet(dataDirectorys, pattern="MethylationEPIC_Sample_Sheet_ConcepcionAguilera_ALL_SAMPLES_to_1_1_2020.csv")
targets

#targets$Sample_Group <- as.character(targets$Sample_Group) 

targets$Sample_Name <- gsub("_",".",targets$Sample_Name)

OLD <- c("R04767","R06071","R15536","R15604","R19729","R21247","R23281","R24959","R33384","R33551","R38002","R38183","R39751","R39962","R41181","R48442","R50197","R66312","R68913","R73191","R75314","R79081","R84145","R87905","R91724","R94897","R97793","S00681","S04885","S07203","S07989","S08206","S09513","S10605","S10751","S11300","S12854","S17163","S17875","S18220","S18312","S20789","S22181","S23321","S24114","S24971","S26140","S28200","S28967","S29760","S34730","S35116","S35901","S37056","S37925","S39280","S39767","S40403","S40754","S40755","S41894","S42693","S43108","S43120","S43130","S43518","S43541","S44237","S44251","S44281","S45131","S47677","S47760","S48108","S49990","S50050","S50144","S52236","S52889","S53656","S54613","S56185","S58315","S58671","S60493","S60709","S60829","S64632","S66605","S66870","S67718","S69608","S70326","S71367","S72403","S72647","S72682","S73741","S74848","S74908","S755143","S75679","S75961","S77036","S77980","S78405","S78549","S79522","S795613","S79880","S81492","S82063","S82098","S825037","S82734","S85555","S86176","S86393","S87205","S87406","S88575","S88810","S89592","S90413","S91084","S92224","S92904","S93897","S95885","S98047","S98733","Z004.2","Z005.2","Z006.2","Z010.2","Z011.2","Z013.2","Z014.2","Z015.2","Z016.2","Z021.1","Z021.2","Z022.2","Z027.2","Z030.2","Z032.1","Z033.2","Z034.2","Z035.1","Z041.2","Z043.2","Z044.2","Z048.2","Z051.2","Z055.2","Z070.1","Z070.2","Z074.1","Z074.2","Z077.1","Z079.2","Z080.1","Z081.1","Z091.1","Z093.1","Z093.2","Z094.2","Z095.1","Z099.2","Z101.1","Z101.2","Z102.1","Z102.2","Z104.2","Z105.1","Z105.2","Z106.1","Z112.1","Z114.1","Z114.2","Z115.2","Z118.2","Z120.1","Z120.2","Z124.2","Z127.1","Z129.2","Z130.1","Z131.1","Z133.2","Z136.1","Z136.2","Z137.1","Z137.2","Z138.2","Z139.2","Z140.1","Z140.2","Z141.2","Z142.2","Z143.2","Z144.2","Z147.2","Z148.2","Z149.2","Z153.2","Z156.2","Z160.2","Z166.2","Z167.2","Z176.2","Z180.2","Z42")

NEW <- c("MR021","MR017","MR020","MR007","MR001","MR004","MR025","MR012","MR015","MR014","MR009","MR027","MR003","MR011","MR002","MR024","MR023","MR019","MR010","MR022","MR018","MR008","MR016","MR013","MR005","MR006","MR026","MS265","MS215","MS243","MS282","MS273","MS272","MS224","MS277","MS280","MS226","MS248","MS211","MS276","MS297","MS299","MS236","MS278","MS255","MS216","MS214","MS250","MS202","MS285","MS261","MS252","MS293","MS225","MS269","MS229","MS271","MS246","MS222","MS223","MS233","MS258","MS207","MS205","MS208","MS203","MS230","MS206","MS200","MS201","MS305","MS218","MS260","MS267","MS303","MS228","MS295","MS247","MS249","MS204","MS279","MS227","MS244","MS296","MS212","MS289","MS287","MS270","MS283","MS291","MS251","MS300","MS240","MS241","MS235","MS275","MS210","MS259","MS268","MS254","MS217","MS231","MS221","MS298","MS257","MS237","MS262","MS288","MS238","MS274","MS239","MS302","MS234","MS281","MS264","MS232","MS256","MS301","MS266","MS253","MS220","MS284","MS242","MS213","MS286","MS290","MS245","MS219","MS304","MS292","MS263","Z478","Z437","Z436","Z476","Z446","Z408","Z475","Z465","Z469","Z472","Z440","Z471","Z473","Z477","Z459","Z442","Z443","Z441","Z474","Z454","Z455","Z420","Z480","Z448","Z466","Z410","Z401","Z470","Z416","Z483","Z402","Z406","Z461","Z407","Z434","Z435","Z419","Z444","Z457","Z424","Z400","Z423","Z452","Z427","Z462","Z415","Z458","Z404","Z445","Z414","Z439","Z453","Z484","Z421","Z468","Z438","Z433","Z403","Z405","Z417","Z431","Z432","Z422","Z426","Z425","Z451","Z467","Z464","Z463","Z486","Z487","Z413","Z412","Z430","Z479","Z447","Z482","Z449","Z450","Z460","Z411","Z456")

targets$Sample_Time
targets$Sample_Time[targets$Sample_Name %in% NEW] <- 2
targets$Sample_Time[targets$Sample_Name %in% OLD] <- 1
targets$Sample_Time

targets$Sample_Group <- as.factor(targets$Sample_Group)
targets$Sample_Group <- factor(targets$Sample_Group,levels(targets$Sample_Group)[c(1,6,4,3,5,2)])
table(targets$Sample_Group)

targets$BATCH <- rep(NA,nrow(targets))
targets$BATCH[which(is.na(targets$Sample_Center))] <- 2
targets$BATCH[which(!is.na(targets$Sample_Center))] <- 1
table(targets$BATCH)

targets <- targets[order(targets$Sample_Name),]
COVARS_ALL <- targets
str(COVARS_ALL)
names(COVARS_ALL)[1] <- "Code"


mVals_BMIQ <- mVals_BMIQ[,sort(colnames(mVals_BMIQ))]

colnames(mVals_BMIQ) == COVARS_ALL$Code


COVARS_ALL$ORDEN <- 1:nrow(COVARS_ALL)

########## PARA FENOTIPOS CONTINUOS



melted <- read.csv2("/home/augusto/Escritorio/BBDD_PUBMEP/BBDD_PubmeP_RESULTADO/BASES_METILACION/BASE_PUBMEP_LONGITUDINAL_LONGformat_NIÑAS_Y_NIÑOS_213inds_02_06_2020_EPIC.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"))
melted$Code <- gsub("-",".",melted$Code)
melted$Code_old <- melted$Code
melted$Code[214:426] <- melted$Code_new[214:426]

melted$Tanner_Time[1:213]
melted$Tanner_Time[214:426]
melted$Tanner_Time[214:426][ (melted$Tanner_Time[214:426] == 1 & !is.na(melted$Tanner_Time[214:426])) ] <- 0
melted$Tanner_Time[214:426]


melted$TIEMPO <- rep(NA,nrow(melted))
melted$TIEMPO[which(melted$Code %in% melted$Code_new)] <- 2
melted$TIEMPO[which(melted$Code %in% melted$Code_old)] <- 1
table(melted$TIEMPO)
melted$TIEMPO

dim(COVARS_ALL)
COVARS_FENOS_CONTINUOS <- merge(COVARS_ALL,melted, by= "Code",all.x=TRUE,all.y=FALSE,sort=FALSE)
dim(COVARS_FENOS_CONTINUOS)

COVARS_FENOS_CONTINUOS <- COVARS_FENOS_CONTINUOS[order(COVARS_FENOS_CONTINUOS$ORDEN),]

colnames(mVals_BMIQ) == COVARS_FENOS_CONTINUOS$Code

ESTIMATE_CELL_COUNTS <- read.csv2("/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/ESTIMACION_CELULAS_BLANCAS/ESTIMACION_CELULAS_BLANCAS.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"))
names(ESTIMATE_CELL_COUNTS)[1] <- "Code"
ESTIMATE_CELL_COUNTS$Code <- gsub("-",".",ESTIMATE_CELL_COUNTS$Code)
COVARS_FENOS_CONTINUOS$Code

COVARS_FENOS_CONTINUOS <- merge(ESTIMATE_CELL_COUNTS,COVARS_FENOS_CONTINUOS, by= "Code",all.x=TRUE,all.y=FALSE,sort=FALSE)

COVARS_FENOS_CONTINUOS <- COVARS_FENOS_CONTINUOS[order(COVARS_FENOS_CONTINUOS$ORDEN),]

colnames(mVals_BMIQ) == COVARS_FENOS_CONTINUOS$Code

#CPGs_MAPPING_S100A4 <- read.csv2("/home/augusto/Escritorio/PROVISIONAL_18_ENERO/ENERO_18/working_NOVIEMBRE_excepcion_/S100A4/CPGs_MAPPING_S100A4.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"))
#CPGs_MAPPING_S100A4[,4]
#SAVE_DATASETS_S100A4:
#write.csv2(mVals_BMIQ[which(rownames(mVals_BMIQ) %in% CPGs_MAPPING_S100A4[,4]),], file="/home/augusto/Escritorio/ANALISIS_EWAS/PRIMER_ANALISIS/S100A4/MVALS_S100A4.csv", row.names=TRUE)
#write.csv2(bVals_BMIQ[which(rownames(bVals_BMIQ) %in% CPGs_MAPPING_S100A4[,4]),], file="/home/augusto/Escritorio/ANALISIS_EWAS/PRIMER_ANALISIS/S100A4/BVALS_S100A4.csv", row.names=TRUE)



####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

"Chunk 1.1: Generacion de datasets que se usaran en cada analisiss"

####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################



##### ##### ##### NIÑxS CON 2 TIEMPOS ANALIZADOS



length(which(table(COVARS_FENOS_CONTINUOS$Code_old)==2))
[1] 102
CODES_CHILDREN_WITH_TWOTIMES <- COVARS_FENOS_CONTINUOS$Code_old[which(duplicated(COVARS_FENOS_CONTINUOS$Code_old))]
length(CODES_CHILDREN_WITH_TWOTIMES)
SELECCION_CHILDREN_WITH_TWO_TIMES <- COVARS_FENOS_CONTINUOS[which(COVARS_FENOS_CONTINUOS$Code_old %in% CODES_CHILDREN_WITH_TWOTIMES),]
SELECCION_CHILDREN_WITH_TWO_TIMES

dim(SELECCION_CHILDREN_WITH_TWO_TIMES)
[1]  204 3542

table(SELECCION_CHILDREN_WITH_TWO_TIMES$TIEMPO)
  1   2 
102 102 

table(SELECCION_CHILDREN_WITH_TWO_TIMES$GRUPOS_exp_interes_NIÑAS)/2

table(SELECCION_CHILDREN_WITH_TWO_TIMES$GRUPOS_exp_interes_NIÑOS)/2

table(SELECCION_CHILDREN_WITH_TWO_TIMES$Sample_Group)/2

#NIÑXS DE ESOS CON DOS TIEMPOS QUE NO TIENEN INFORMACION DE GRUPO EXPERIMENTAL LONGITUDINAL
SELECCION_CHILDREN_WITH_TWO_TIMES$Code_old[which(is.na(SELECCION_CHILDREN_WITH_TWO_TIMES$Sample_Group))]
"S755143"
"Z143-2"
"Z006-2"
"Z070-2"



##### ##### ##### NIÑxS CON SOLO 1 TIEMPO ANALIZADO (EL PÚBER)



length(which(table(COVARS_FENOS_CONTINUOS$Code_old)==1))
[1] 36
CODES_CHILDREN_WITH_ONETIME <- COVARS_FENOS_CONTINUOS$Code_old[which(!COVARS_FENOS_CONTINUOS$Code_old %in% CODES_CHILDREN_WITH_TWOTIMES)]

SELECCION_CHILDREN_WITH_ONE_TIME <- COVARS_FENOS_CONTINUOS[which(COVARS_FENOS_CONTINUOS$Code_old %in% CODES_CHILDREN_WITH_ONETIME),]
SELECCION_CHILDREN_WITH_ONE_TIME

dim(SELECCION_CHILDREN_WITH_ONE_TIME)
[1]   36 3542

table(SELECCION_CHILDREN_WITH_ONE_TIME$TIEMPO)
 2 
36 

table(SELECCION_CHILDREN_WITH_ONE_TIME$Sex_Time)
 Niña_T2 Varon_T2 
      17       19 



##### ##### ##### SELECCION POBLACION EFECTIVA ANALISIS LONGITUDINALES



codes102a <- SELECCION_CHILDREN_WITH_TWO_TIMES$Code_old[ which( SELECCION_CHILDREN_WITH_TWO_TIMES$Tanner_Time !=0 & SELECCION_CHILDREN_WITH_TWO_TIMES$time == 1 ) ]

codes102b <- SELECCION_CHILDREN_WITH_TWO_TIMES$Code_old[ which( SELECCION_CHILDREN_WITH_TWO_TIMES$Tanner_Time ==0 & SELECCION_CHILDREN_WITH_TWO_TIMES$time == 2 ) ]

codes_old_toExclude <- c(codes102a,codes102b)
codes_old_toExclude 

codes_all_toExclude <- SELECCION_CHILDREN_WITH_TWO_TIMES$Code[ which(SELECCION_CHILDREN_WITH_TWO_TIMES$Code_old %in% codes_old_toExclude ) ]


POBLACION_EFECTIVA_LONGI <- SELECCION_CHILDREN_WITH_TWO_TIMES[-which( SELECCION_CHILDREN_WITH_TWO_TIMES$Code %in% codes_all_toExclude ),]
dim(POBLACION_EFECTIVA_LONGI)


##### ##### ##### SELECCION POBLACION EFECTIVA ANALISIS TRANSVERSAL T0 PREPUBER



codes102 <- SELECCION_CHILDREN_WITH_TWO_TIMES$Code[ which( SELECCION_CHILDREN_WITH_TWO_TIMES$Tanner_Time == 0 & SELECCION_CHILDREN_WITH_TWO_TIMES$time == 1 ) ]

codes36 <- SELECCION_CHILDREN_WITH_ONE_TIME$Code[ which( SELECCION_CHILDREN_WITH_ONE_TIME$Tanner_Time == 0 ) ]

codes_all_toInclude <- c(codes102,codes36)
codes_all_toInclude 

POBLACION_EFECTIVA_PREPUBER_CROSS <- COVARS_FENOS_CONTINUOS[which( COVARS_FENOS_CONTINUOS$Code %in% codes_all_toInclude ),]
dim(POBLACION_EFECTIVA_PREPUBER_CROSS)



##### ##### ##### SELECCION POBLACION EFECTIVA ANALISIS TRANSVERSAL T1 PUBER



codes102 <- SELECCION_CHILDREN_WITH_TWO_TIMES$Code[ which( SELECCION_CHILDREN_WITH_TWO_TIMES$Tanner_Time != 0 & SELECCION_CHILDREN_WITH_TWO_TIMES$time == 2 ) ]

codes36 <- SELECCION_CHILDREN_WITH_ONE_TIME$Code[ which( SELECCION_CHILDREN_WITH_ONE_TIME$Tanner_Time != 0 ) ]

codes_all_toInclude <- c(codes102,codes36)
codes_all_toInclude 

POBLACION_EFECTIVA_PUBERtal_CROSS <- COVARS_FENOS_CONTINUOS[which( COVARS_FENOS_CONTINUOS$Code %in% codes_all_toInclude ),]
dim(POBLACION_EFECTIVA_PUBERtal_CROSS)





##### ##### ##### DATASETS RESULTANTES:



dim(POBLACION_EFECTIVA_LONGI)
[1]  196 3542

dim(POBLACION_EFECTIVA_PREPUBER_CROSS)
[1]  104 3542

dim(POBLACION_EFECTIVA_PUBERtal_CROSS)
[1]  132 3542


sort(POBLACION_EFECTIVA_LONGI$time)

sort(POBLACION_EFECTIVA_PREPUBER_CROSS$time)

sort(POBLACION_EFECTIVA_PUBERtal_CROSS$time)



sort(POBLACION_EFECTIVA_LONGI$Tanner_Time)

sort(POBLACION_EFECTIVA_PREPUBER_CROSS$Tanner_Time)

sort(POBLACION_EFECTIVA_PUBERtal_CROSS$Tanner_Time)


POBLACION_EFECTIVA_LONGI$Sex_Time <- gsub("_T1","",POBLACION_EFECTIVA_LONGI$Sex_Time)
POBLACION_EFECTIVA_LONGI$Sex_Time <- gsub("_T2","",POBLACION_EFECTIVA_LONGI$Sex_Time)
POBLACION_EFECTIVA_LONGI$Sex_Time <- as.factor(POBLACION_EFECTIVA_LONGI$Sex_Time)

POBLACION_EFECTIVA_PREPUBER_CROSS$Sex_Time <- gsub("_T1","",POBLACION_EFECTIVA_PREPUBER_CROSS$Sex_Time)
POBLACION_EFECTIVA_PREPUBER_CROSS$Sex_Time <- gsub("_T2","",POBLACION_EFECTIVA_PREPUBER_CROSS$Sex_Time)
POBLACION_EFECTIVA_PREPUBER_CROSS$Sex_Time <- as.factor(POBLACION_EFECTIVA_PREPUBER_CROSS$Sex_Time)

POBLACION_EFECTIVA_PUBERtal_CROSS$Sex_Time <- gsub("_T1","",POBLACION_EFECTIVA_PUBERtal_CROSS$Sex_Time)
POBLACION_EFECTIVA_PUBERtal_CROSS$Sex_Time <- gsub("_T2","",POBLACION_EFECTIVA_PUBERtal_CROSS$Sex_Time)
POBLACION_EFECTIVA_PUBERtal_CROSS$Sex_Time <- as.factor(POBLACION_EFECTIVA_PUBERtal_CROSS$Sex_Time)



COVARS_FENOS_CONTINUOS$Sex_Time <- gsub("_T1","",COVARS_FENOS_CONTINUOS$Sex_Time)
COVARS_FENOS_CONTINUOS$Sex_Time <- gsub("_T2","",COVARS_FENOS_CONTINUOS$Sex_Time)
COVARS_FENOS_CONTINUOS$Sex_Time <- as.factor(COVARS_FENOS_CONTINUOS$Sex_Time)






####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

"Chunk 2: QUALITY CONTROL PLOTS PRIOR ANALYSIS"

####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################



COVARS_ALL <- POBLACION_EFECTIVA_PREPUBER_CROSS

dim(COVARS_ALL)

# plot cell type proportions by TIME
#pdf("/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/ESTIMACION_CELULAS_BLANCAS/CELL_ESTIMATIO_BOXplot_PUBERTAL_STAGE.pdf", width=40/2.54, height=28/2.54)
par(mfrow=c(1,1))
age.pal <- brewer.pal(8,"Set1")
a = COVARS_ALL[COVARS_ALL$TIEMPO == 1,c(2:7)]
b = COVARS_ALL[COVARS_ALL$TIEMPO == 2,c(2:7)]
boxplot(a, at=0:5*3 + 1, xlim=c(0, 18), xaxt="n",
col=age.pal[1], main="", ylab="Cell type proportion")
boxplot(b, at=0:5*3 + 2, xaxt="n", add=TRUE, col=age.pal[2])
axis(1, at=0:5*3 + 1.5, labels=colnames(a), tick=TRUE)
legend("topleft", legend=c("Pre-pubertal","Pubertal"), fill=age.pal)
#dev.off()
#
#
## plot cell type proportions by EXPERIMENTAL GROUP
dataset <- na.omit(melt(COVARS_ALL[,c(2:7,10)]))
pal <- brewer.pal(6,"OrRd")
ggplot(dataset,aes(x=variable,y=value,fill=Sample_Group))+geom_boxplot()+ scale_fill_manual(values=pal)+ labs(x = "Cell Type", y = "Cell Counts")
#ggsave("/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/ESTIMACION_CELULAS_BLANCAS/CELL_ESTIMATION_BOXplot_GROUPS.pdf")
#
#
#
## plot cell type proportions by ANALYSIS BATCH
ji <- COVARS_ALL[,c(2:7,17)]
ji$BATCH<- as.factor(ji$BATCH)
dataset <- na.omit(melt(ji))
pal <- brewer.pal(5,"OrRd")
ggplot(dataset,aes(x=variable,y=value,fill=BATCH))+geom_boxplot()+ scale_fill_manual(values=pal)+ labs(x = "Cell Type", y = "Cell Counts")
#ggsave("/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/ESTIMACION_CELULAS_BLANCAS/CELL_ESTIMATION_BOXplot_BATCH.pdf")
#
#
#
#
## plot PCA by sex
#pdf("/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/plot_qc/PCA_BMIQ_NORMALIZED_SEX.pdf", width=35/2.54, height=18/2.54)
pal <- brewer.pal(11,"RdBu")[c(2,9)]
plotMDS(mVals_BMIQ[,which(colnames(mVals_BMIQ) %in% POBLACION_EFECTIVA_PREPUBER_CROSS$Code)], top=1000, gene.selection="common",col=pal[COVARS_ALL$Sex_Time])
legend("topleft", legend=levels(COVARS_ALL$Sex_Time), text.col=pal,bg="white",cex=0.7,title="Sex",title.col="black")
#dev.off()
#
#
#
#
#
## plot PCA by batch
# pdf("/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/plot_qc/PCA_BMIQ_NORMALIZED_batch.pdf", width=35/2.54, height=18/2.54)
pal <- brewer.pal(11,"RdBu")[c(2,9)]
plotMDS(mVals_BMIQ[,which(colnames(mVals_BMIQ) %in% POBLACION_EFECTIVA_PREPUBER_CROSS$Code)], top=1000, gene.selection="common",col=pal[factor(COVARS_ALL$BATCH)])
legend("topleft", legend=levels(factor(COVARS_ALL$BATCH)), text.col=pal,bg="white",cex=0.7,title="Batch",title.col="black")
#dev.off()
#



####################################################################################################################################################################################################################################################################################################################

#"FACTOR ANALYSIS SOBRE COVARIABLES PRINCIPALES"

####################################################################################################################################################################################################################################################################################################################

#install.packages("swamp")
#library(swamp)
#
#COVAR <- COVARS_ALL[,c("Code","Sample_Group","Sex_Time","Age_Time","Habitos_toxicos_en_el_embarazo_Time","Diabetes_gestacional_Time","Tanner_Time","Pubarquia_T2","Telarquia_T2","Axilarquia_T2","TIEMPO","BATCH","Presencia_de_Menarquia_Time","Si_es_si_edad_de_inicio_Time","counts.CD8T","counts.CD4T","counts.NK","counts.Bcell","counts.Mono","counts.Neu","Origen_T1")]
#
#COVAR$Code == colnames(mVals_BMIQ[,which(colnames(mVals_BMIQ) %in% POBLACION_EFECTIVA_LONGI$Code)])
#
#COVAR <- COVAR[,-1]
#
# calculate M-values IN BMIQ DATASET for statistical analysis
#
#mVals_BMIQ_QUANTILE <- swamp::prince(as.matrix(mVals_BMIQ[,which(colnames(mVals_BMIQ) %in% POBLACION_EFECTIVA_LONGI$Code)]), COVAR, top = 10, imputeknn = T, center = T, permute = F)
#
#pdf("/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/plot_qc/mVals_BMIQ_prince_plot.pdf", width=50/2.54, height=18/2.54)
#swamp::prince.plot(mVals_BMIQ_QUANTILE,note=T,Rsquared=T,cexRow=0.6)
#dev.off()



####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

"Chunk 3: ANALISIS LONGITUDINAL (within and between)"
#EXPLICADO EN: https://support.bioconductor.org/p/74247/

####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################



COVARS_ALL <- POBLACION_EFECTIVA_PREPUBER_CROSS

dim(COVARS_ALL)
# 104 3542

table(COVARS_ALL$TIEMPO)
#  1   2 
# 101   3 

table(COVARS_ALL$TIEMPO,COVARS_ALL$Tanner_Time)
#      0
#  1 101
#  2   3

table(COVARS_ALL$Sex_Time,COVARS_ALL$Sample_Group)
#        NW non-IR no change OB/OW non-IR to NW non-IR OB/OW non-IR no change
#  Niña                   12                         4                      9
#  Varon                  15                         0                     17
#       
#        OB/OW IR to non-IR OB/OW non-IR to IR OB/OW IR no change
#  Niña                   9                 11                  6
#  Varon                  5                  5                  4




table(COVARS_ALL$Sample_Group)/2
#      NW non-IR no change OB/OW non-IR to NW non-IR    OB/OW non-IR no change 
#                       27                         4                        26 
#       OB/OW IR to non-IR        OB/OW non-IR to IR        OB/OW IR no change 
#                       14                        16                        10        



COVARS_ALL$GRS_SI_O_NO <- rep(NA,nrow(COVARS_ALL))
COVARS_ALL$GRS_SI_O_NO[!is.na(COVARS_ALL$Genetic_Risk_Score_Locke44SNPs_T1)] <- 1
table(COVARS_ALL$GRS_SI_O_NO)
# 1 
#68 



table(COVARS_ALL$Sample_Group,COVARS_ALL$GRS_SI_O_NO)
#  NW non-IR no change       18
#  OB/OW non-IR to NW non-IR  1
#  OB/OW non-IR no change    16
#  OB/OW IR to non-IR         8
#  OB/OW non-IR to IR        11
#  OB/OW IR no change         9


dim(COVARS_ALL)
COVARS_ALL <- COVARS_ALL[-which((COVARS_ALL$Cole_Time=="Normopeso_T1" & COVARS_ALL$HOMA_0_AUG_Time==1) | (COVARS_ALL$Cole_Time=="Normopeso_T2" & COVARS_ALL$HOMA_0_AUG_Time==1)),]
dim(COVARS_ALL)





MICROARRAY <- mVals_BMIQ 

MICROARRAY <- MICROARRAY[,which(colnames(MICROARRAY) %in% COVARS_ALL$Code)] #SELECCIONAR EN ARRAY SOLO AQUELLOS NIÑOS SUSCEPTIBLES DE ANALISIS
dim(MICROARRAY)
#

COVARS_FOR_ANALYSIS <- COVARS_ALL

rownames(COVARS_FOR_ANALYSIS) <- COVARS_FOR_ANALYSIS$Code
colnames(MICROARRAY) == rownames(COVARS_FOR_ANALYSIS)
#TRUE




####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

"Chunk 3: PREPARANDO DATASET ANALISIS CROSS-SECTIONAL"
#EXPLICADO EN: https://support.bioconductor.org/p/74247/

####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################


####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

"Chunk 3.1: Seleccionar TIEMPO 0"

####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################


#SELECCION DE NIÑOS PREPUBERES
COVARS_ALL <- COVARS_FOR_ANALYSIS


COVARS_ALL$Cole_Time[ COVARS_ALL$Cole_Time == "Normopeso_T2"] <- "Normopeso_T1"
COVARS_ALL$Cole_Time[ COVARS_ALL$Cole_Time == "Sobrepeso_T2"] <- "Sobrepeso_T1"
COVARS_ALL$Cole_Time[ COVARS_ALL$Cole_Time == "Obeso_T2"] <- "Obeso_T1"
table(COVARS_ALL$Cole_Time)
#NOTA: NO HAY NINGUN INDIVIDUO SIN DATO
#Normopeso_T1     Obeso_T1 Sobrepeso_T1 
#          30           51           21 

table(COVARS_ALL$HOMA_0_AUG_Time)
#NOTA: HAY 2 INDIVIDUOS QUE NO TIENEN GRUPO EXPERIMENTAL HOMA
# 0  1 
#75  25


dim(COVARS_ALL)
#[1]  102 3543

COVARS_FOR_ANALYSIS <- COVARS_ALL

rownames(COVARS_FOR_ANALYSIS) <- COVARS_FOR_ANALYSIS$Code
colnames(MICROARRAY) == rownames(COVARS_FOR_ANALYSIS)
#TRUE

COVARS_FENOS_CONTINUOS_t0 <- COVARS_FOR_ANALYSIS
COVARS_FENOS_CONTINUOS_t0$HOMA_0_AUG_Time[ COVARS_FENOS_CONTINUOS_t0$HOMA_0_AUG_Time == 0 ] <- "No-IR"
COVARS_FENOS_CONTINUOS_t0$HOMA_0_AUG_Time[ COVARS_FENOS_CONTINUOS_t0$HOMA_0_AUG_Time == 1 ] <- "IR"


######################
######################
###################### UN SOLO FACTOR (GRUPOS LONGITUDINALES)
############################################
############################################

####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

97 INDIVIDUOS

####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

MICROARRAY_PREPUBER <- MICROARRAY[,-which(is.na(COVARS_FENOS_CONTINUOS_t0$Sample_Group))] 


COVARS_FENOS_CONTINUOS_PREPUBER <- COVARS_FENOS_CONTINUOS_t0[-which(is.na(COVARS_FENOS_CONTINUOS_t0$Sample_Group)),]


colnames(MICROARRAY_PREPUBER) == rownames(COVARS_FENOS_CONTINUOS_PREPUBER)

dim(COVARS_FENOS_CONTINUOS_PREPUBER)
#[1]  97 3543

GROUP_CONJUNTO <-  factor(COVARS_FENOS_CONTINUOS_PREPUBER$Sample_Group)
GROUP_CONJUNTO


OBESIDAD <- as.factor(COVARS_FENOS_CONTINUOS_PREPUBER$Cole_Time)
#OBESIDAD[ which(OBESIDAD == "Sobrepeso_T1") ] <- "Obeso_T1"
OBESIDAD <- factor(OBESIDAD,levels(OBESIDAD)[c(1,3,2)])
table(OBESIDAD)
OBESIDAD

HOMA <- as.factor(COVARS_FENOS_CONTINUOS_PREPUBER$HOMA_0_AUG_Time)
HOMA <- factor(HOMA,levels(HOMA)[c(2,1)])
HOMA


cor(as.numeric(COVARS_FENOS_CONTINUOS_PREPUBER$Sex_Time),as.numeric(COVARS_FENOS_CONTINUOS_PREPUBER$BATCH))
[1] 0.9019608 #CON ESTA CORRELACION ENTRE SEXO Y BATCH NO VAMOS A METER BATCH COMO UN FACTOR CONFUSOR


####################################################################

design = model.matrix(~ 0 + GROUP_CONJUNTO + COVARS_FENOS_CONTINUOS_PREPUBER$counts.CD8T + COVARS_FENOS_CONTINUOS_PREPUBER$counts.CD4T + COVARS_FENOS_CONTINUOS_PREPUBER$counts.NK + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Bcell + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Mono + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Neu + COVARS_FENOS_CONTINUOS_PREPUBER$Sex_Time + COVARS_FENOS_CONTINUOS_PREPUBER$Age_Time + as.factor(COVARS_FENOS_CONTINUOS_PREPUBER$Origen_T1))
colnames(design)[1:6] = levels(GROUP_CONJUNTO)
colnames(design)[7:16] = c("CD8T","CD4T","NK","BCELL","MONO","NEU","SEX","AGE","ORIGEN1","ORIGEN2")
colnames(design)

 colnames(design) <- gsub("/","", colnames(design))

 colnames(design) <- gsub("-","", colnames(design))

 colnames(design) <- gsub(" ","", colnames(design))

 colnames(design) 
 [1] "NWnonIRnochange"    "OBOWnonIRtoNWnonIR" "OBOWnonIRnochange" 
 [4] "OBOWIRtononIR"      "OBOWnonIRtoIR"      "OBOWIRnochange"    
 [7] "CD8T"               "CD4T"               "NK"                
[10] "BCELL"              "MONO"               "NEU"               
[13] "SEX"                "AGE"                "ORIGEN1"           
[16] "ORIGEN2"    


data.fit = lmFit(MICROARRAY_PREPUBER, design)
colnames(data.fit)
 [1] "NWnonIRnochange"    "OBOWnonIRtoNWnonIR" "OBOWnonIRnochange" 
 [4] "OBOWIRtononIR"      "OBOWnonIRtoIR"      "OBOWIRnochange"    
 [7] "CD8T"               "CD4T"               "NK"                
[10] "BCELL"              "MONO"               "NEU"               
[13] "SEX"                "AGE"  



CONTRASTES_BETWEEN <- rbind(c(1,0,-1,0,0,0,rep(0,10)),c(1,0,0,-1,0,0,rep(0,10)),c(1,0,0,0,-1,0,rep(0,10)),c(1,0,0,0,0,-1,rep(0,10)),c(0,0,1,-1,0,0,rep(0,10)),c(0,0,1,0,-1,0,rep(0,10)),c(0,0,1,0,0,-1,rep(0,10)),c(0,0,0,1,-1,0,rep(0,10)),c(0,0,0,1,0,-1,rep(0,10)),c(0,0,0,0,1,-1,rep(0,10)),c(0,1,0,0,0,-1,rep(0,10)),c(0,1,0,0,-1,0,rep(0,10)),c(0,1,0,-1,0,0,rep(0,10)),c(0,1,-1,0,0,0,rep(0,10)),c(-1,1,0,0,0,0,rep(0,10)))
colnames(CONTRASTES_BETWEEN) <- colnames(data.fit)
rownames(CONTRASTES_BETWEEN) <- as.character(c(1:15))
rownames(CONTRASTES_BETWEEN)[1] <- "NW non-IR no change VS OB/OW non-IR no change"
rownames(CONTRASTES_BETWEEN)[2] <- "NW non-IR no change VS OB/OW IR to non-IR"
rownames(CONTRASTES_BETWEEN)[3] <- "NW non-IR no change VS OB/OW non-IR to IR"
rownames(CONTRASTES_BETWEEN)[4] <- "NW non-IR no change VS  OB/OW IR no change"
rownames(CONTRASTES_BETWEEN)[5] <- "OB/OW non-IR no change VS OB/OW IR to non-IR"
rownames(CONTRASTES_BETWEEN)[6] <- "OB/OW non-IR no change VS OB/OW non-IR to IR"
rownames(CONTRASTES_BETWEEN)[7] <- "OB/OW non-IR no change VS OB/OW IR no change"
rownames(CONTRASTES_BETWEEN)[8] <- "OB/OW IR to non-IR VS OB/OW non-IR to IR"
rownames(CONTRASTES_BETWEEN)[9] <- "OB/OW IR to non-IR VS OB/OW IR no change"
rownames(CONTRASTES_BETWEEN)[10] <- "OB/OW non-IR to IR VS OB/OW IR no change"
rownames(CONTRASTES_BETWEEN)[11] <- "OB/OW non-IR to NW non-IR VS OB/OW IR no change"
rownames(CONTRASTES_BETWEEN)[12] <- "OB/OW non-IR to NW non-IR VS OB/OW non-IR to IR"
rownames(CONTRASTES_BETWEEN)[13] <- "OB/OW non-IR to NW non-IR VS OB/OW IR to non-IR"
rownames(CONTRASTES_BETWEEN)[14] <- "OB/OW non-IR to NW non-IR VS OB/OW non-IR no change"
rownames(CONTRASTES_BETWEEN)[15] <- "OB/OW non-IR to NW non-IR VS NW non-IR no change"


CONTRASTES_BETWEEN <- t(CONTRASTES_BETWEEN)

lmfit.cont <- contrasts.fit(data.fit, CONTRASTES_BETWEEN)
lmfit.cont.ebayes <- eBayes(lmfit.cont)
lmfit.cont.ebayes$p.value
summary(decideTests(lmfit.cont.ebayes))



############### PARA SACAR LISTA DE SIGNIFICATIVAS

colnames(summary(decideTests(lmfit.cont.ebayes)))
 [1] "NW non-IR no change VS OB/OW non-IR no change"      
 [2] "NW non-IR no change VS OB/OW IR to non-IR"          
 [3] "NW non-IR no change VS OB/OW non-IR to IR"          
 [4] "NW non-IR no change VS  OB/OW IR no change"         
 [5] "OB/OW non-IR no change VS OB/OW IR to non-IR"       
 [6] "OB/OW non-IR no change VS OB/OW non-IR to IR"       
 [7] "OB/OW non-IR no change VS OB/OW IR no change"       
 [8] "OB/OW IR to non-IR VS OB/OW non-IR to IR"           
 [9] "OB/OW IR to non-IR VS OB/OW IR no change"           
[10] "OB/OW non-IR to IR VS OB/OW IR no change"           
[11] "OB/OW non-IR to NW non-IR VS OB/OW IR no change"    
[12] "OB/OW non-IR to NW non-IR VS OB/OW non-IR to IR"    
[13] "OB/OW non-IR to NW non-IR VS OB/OW IR to non-IR"    
[14] "OB/OW non-IR to NW non-IR VS OB/OW non-IR no change"
[15] "OB/OW non-IR to NW non-IR VS NW non-IR no change"  

colnames(data.fit)
 [1] "NWnonIRnochange"    "OBOWnonIRtoNWnonIR" "OBOWnonIRnochange" 
 [4] "OBOWIRtononIR"      "OBOWnonIRtoIR"      "OBOWIRnochange"    
 [7] "CD8T"               "CD4T"               "NK"                
[10] "BCELL"              "MONO"               "NEU"               
[13] "SEX"                "AGE"                "ORIGEN1"           
[16] "ORIGEN2"          


#A CONTINUACION REPITO LO ANTERIOR DE OTRA FORMA PARA ASEGURARME QUE TODO OK.

#Now we can make any comparisons between the experimental conditions in the usual way, example:
cm <- makeContrasts(NWvsOBNW=NWnonIRnochange-OBOWnonIRtoNWnonIR,NWvsOBsano=NWnonIRnochange-OBOWnonIRnochange,NWvsOBIRnoir=NWnonIRnochange-OBOWIRtononIR,NWvsOBnoirIR=NWnonIRnochange-OBOWnonIRtoIR,NWvsOBIRnochange=NWnonIRnochange-OBOWIRnochange,OBNWvsOBsano=OBOWnonIRtoNWnonIR-OBOWnonIRnochange,OBNWvsOBIRnoir=OBOWnonIRtoNWnonIR-OBOWIRtononIR,OBNWvsOBnoirIR=OBOWnonIRtoNWnonIR-OBOWnonIRtoIR,OBNWvsOBIRnochange=OBOWnonIRtoNWnonIR-OBOWIRnochange,OBsanovsOBIRnoir=OBOWnonIRnochange-OBOWIRtononIR,OBsanovsOBnoirIR=OBOWnonIRnochange-OBOWnonIRtoIR,OBsanovsOBIRnochange=OBOWnonIRnochange-OBOWIRnochange,OBIRnoirVSOBnoirIR=OBOWIRtononIR-OBOWnonIRtoIR,OBIRnoirVSOBIRnochange=OBOWIRtononIR-OBOWIRnochange,OBnoirIRVSOBIRnochange=OBOWnonIRtoIR-OBOWIRnochange, levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(data.fit, cm)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))
       NWvsOBNW NWvsOBsano NWvsOBIRnoir NWvsOBnoirIR NWvsOBIRnochange
Down       6475          0            0            0                0
NotSig   799100     809380       809381       809381           809381
Up         3806          1            0            0                0
       OBNWvsOBsano OBNWvsOBIRnoir OBNWvsOBnoirIR OBNWvsOBIRnochange
Down           9869          16517           5971               3104
NotSig       785582         780258         789644             801333
Up            13930          12606          13766               4944
       OBsanovsOBIRnoir OBsanovsOBnoirIR OBsanovsOBIRnochange
Down                  0                0                    0
NotSig           809381           809381               809381
Up                    0                0                    0
       OBIRnoirVSOBnoirIR OBIRnoirVSOBIRnochange OBnoirIRVSOBIRnochange
Down                    0                      0                      0
NotSig             809381                 809381                 809381
Up                      0                      0                      0



############### PARA SACAR LISTA DE SIGNIFICATIVAS


colnames(summary(decideTests(data.fit.eb)))
 [1] "NWvsOBNW"               "NWvsOBsano"             "NWvsOBIRnoir"          
 [4] "NWvsOBnoirIR"           "NWvsOBIRnochange"       "OBNWvsOBsano"          
 [7] "OBNWvsOBIRnoir"         "OBNWvsOBnoirIR"         "OBNWvsOBIRnochange"    
[10] "OBsanovsOBIRnoir"       "OBsanovsOBnoirIR"       "OBsanovsOBIRnochange"  
[13] "OBIRnoirVSOBnoirIR"     "OBIRnoirVSOBIRnochange" "OBnoirIRVSOBIRnochange"


# seleccion de significativas
ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)

ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)


#DMPs_significants <- subset(DMPs, P.Value < 0.000001)
DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/NWvsOBNW_PREPUBER_adjusted_LEUCOS_AGE_SEX_FDR.csv", row.names=FALSE)


# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=2, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/NWvsOBsano_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001FDR05.csv", row.names=FALSE)




# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=3, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/NWvsOBIRnoir_PREPUBER_adjusted_LEUCOS_AGE_SEXandrawP0001.csv", row.names=FALSE)



# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=4, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/NWvsOBnoirIR_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)





# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=5, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/NWvsOBIRnochange_PREPUBER_adjusted_LEUCOS_AGE_SEXrawP0001.csv", row.names=FALSE)





# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=6, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

#DMPs_significants <- subset(DMPs, P.Value < 0.000001)
DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/OBNWvsOBsano_PREPUBER_adjusted_LEUCOS_AGE_SEX_FDR.csv", row.names=FALSE)




# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=7, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

#DMPs_significants <- subset(DMPs, P.Value < 0.000001)
DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/OBNWvsOBIRnoir_PREPUBER_adjusted_LEUCOS_AGE_SEX_FDR.csv", row.names=FALSE)



# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=8, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

#DMPs_significants <- subset(DMPs, P.Value < 0.000001)
DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/OBNWvsOBnoirIR_PREPUBER_adjusted_LEUCOS_AGE_SEX_FDR.csv", row.names=FALSE)




# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=9, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

#DMPs_significants <- subset(DMPs, P.Value < 0.000001)
DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/OBNWvsOBIRnochange_PREPUBER_adjusted_LEUCOS_AGE_SEX_FDR.csv", row.names=FALSE)





# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=10, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/OBsanovsOBIRnoir_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)






# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=11, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/OBsanovsOBnoirIR_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)



# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=12, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/OBsanovsOBIRnochange_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)




# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=13, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/OBIRnoirVSOBnoirIR_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)




# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=14, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/OBIRnoirVSOBIRnochange_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)




# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=15, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_LONGITUDINALES/OBnoirIRVSOBIRnochange_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)






















######################
######################
###################### TRES FACTORES : NORMO vs OW vs OB
############################################
############################################


####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

"SELECCION DE POBLACION"

####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

MICROARRAY_PREPUBER <- MICROARRAY[,-which(is.na(COVARS_FENOS_CONTINUOS_t0$HOMA_0_AUG_Time))]

COVARS_FENOS_CONTINUOS_PREPUBER <- COVARS_FENOS_CONTINUOS_t0[-which(is.na(COVARS_FENOS_CONTINUOS_t0$HOMA_0_AUG_Time)),]

dim(COVARS_FENOS_CONTINUOS_PREPUBER)
#[1]  100 3543

colnames(MICROARRAY_PREPUBER) == rownames(COVARS_FENOS_CONTINUOS_PREPUBER)


OBESIDAD <- as.factor(COVARS_FENOS_CONTINUOS_PREPUBER$Cole_Time)
#OBESIDAD[ which(OBESIDAD == "Sobrepeso_T1") ] <- "Obeso_T1"
OBESIDAD <- factor(OBESIDAD,levels(OBESIDAD)[c(1,3,2)])
table(OBESIDAD)
OBESIDAD

HOMA <- as.factor(COVARS_FENOS_CONTINUOS_PREPUBER$HOMA_0_AUG_Time)
HOMA <- factor(HOMA,levels(HOMA)[c(2,1)])
HOMA

####################################################################

factored.design = model.matrix(~ 0 + OBESIDAD + COVARS_FENOS_CONTINUOS_PREPUBER$counts.CD8T + COVARS_FENOS_CONTINUOS_PREPUBER$counts.CD4T + COVARS_FENOS_CONTINUOS_PREPUBER$counts.NK + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Bcell + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Mono + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Neu + COVARS_FENOS_CONTINUOS_PREPUBER$Sex_Time + COVARS_FENOS_CONTINUOS_PREPUBER$Age_Time + factor(COVARS_FENOS_CONTINUOS_PREPUBER$Origen_T1) + HOMA)



 colnames(factored.design) <- gsub("OBESIDAD","", colnames(factored.design))

 colnames(factored.design) <- gsub("_T1","", colnames(factored.design))

 colnames(factored.design)[4:14]  <- c("CD8T","CD4T","NK","BCELL","MONO","NEU","SEX","AGE","ORIGEN1","ORIGEN2","HOMA")

 colnames(factored.design)

data.fit = lmFit(MICROARRAY_PREPUBER,factored.design)



#Now we can make any comparisons between the experimental conditions in the usual way, example:
cm <- makeContrasts(NWvsOW=Normopeso-Sobrepeso,NWvsOB=Normopeso-Obeso,OWvsOB=Sobrepeso-Obeso, levels=factored.design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(data.fit, cm)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))


ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBESIDAD_NO_solapando_OBOW/NormopesoVSSobrepeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)




# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=2, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBESIDAD_NO_solapando_OBOW/NormopesoVSObeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)





# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=3, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBESIDAD_NO_solapando_OBOW/SobrepesoVSObeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)






















######################
######################
###################### DOS FACTORES : NORMO vs OWOB
############################################
############################################

dim(COVARS_FENOS_CONTINUOS_PREPUBER)
#[1]  100 3543

OBESIDAD__ <- as.factor(COVARS_FENOS_CONTINUOS_PREPUBER$Cole_Time)
OBESIDAD__[ which(OBESIDAD__ == "Sobrepeso_T1") ] <- "Obeso_T1"
table(OBESIDAD__)
OBESIDAD__ <- factor(OBESIDAD__)
OBESIDAD__


factored.design = model.matrix(~ OBESIDAD__ + COVARS_FENOS_CONTINUOS_PREPUBER$counts.CD8T + COVARS_FENOS_CONTINUOS_PREPUBER$counts.CD4T + COVARS_FENOS_CONTINUOS_PREPUBER$counts.NK + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Bcell + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Mono + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Neu + COVARS_FENOS_CONTINUOS_PREPUBER$Sex_Time + COVARS_FENOS_CONTINUOS_PREPUBER$Age_Time + factor(COVARS_FENOS_CONTINUOS_PREPUBER$Origen_T1) + HOMA)
data.fit = lmFit(MICROARRAY_PREPUBER,factored.design)
data.fit.eb = eBayes(data.fit)
summary(decideTests(data.fit.eb))


colnames(summary(decideTests(data.fit.eb)))

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=2, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBESIDAD_solapando_OBOW/NormopesoObeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)



# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=13, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBESIDAD_solapando_OBOW/IR_VS_noIR_PREPUBER_adjusted_LEUCOS_AGE_SEX_FDR.csv", row.names=FALSE)




















######################
######################
###################### CINCO FACTORES : NORMO vs OWnoIR vs OWIR vs OBnoIR vs OBIR
############################################
############################################



####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

"SELECCION DE POBLACION"

####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################

MICROARRAY_PREPUBER <- MICROARRAY[,-which(is.na(COVARS_FENOS_CONTINUOS_t0$HOMA_0_AUG_Time
))] 


COVARS_FENOS_CONTINUOS_PREPUBER <- COVARS_FENOS_CONTINUOS_t0[-which(is.na(COVARS_FENOS_CONTINUOS_t0$HOMA_0_AUG_Time
)),]

dim(COVARS_FENOS_CONTINUOS_PREPUBER)
#[1]  100 3543

colnames(MICROARRAY_PREPUBER) == rownames(COVARS_FENOS_CONTINUOS_PREPUBER)




OBESIDAD <- as.factor(COVARS_FENOS_CONTINUOS_PREPUBER$Cole_Time)
#OBESIDAD[ which(OBESIDAD == "Sobrepeso_T1") ] <- "Obeso_T1"
OBESIDAD <- factor(OBESIDAD,levels(OBESIDAD)[c(1,3,2)])
table(OBESIDAD)
OBESIDAD

HOMA <- as.character(COVARS_FENOS_CONTINUOS_PREPUBER$HOMA_0_AUG_Time)
HOMA[ HOMA == "No-IR" ] <- "NoIR"
HOMA <- as.factor(HOMA)
HOMA <- factor(HOMA,levels(HOMA)[c(2,1)])
HOMA


HOMAOB <- paste(HOMA,OBESIDAD,sep="_")

HOMAOB <- factor(HOMAOB)

HOMAOB <- factor(HOMAOB,levels(HOMAOB)[c(3,5,2,4,1)])
HOMAOB
#NoIR_Normopeso_T1 NoIR_Sobrepeso_T1   IR_Sobrepeso_T1     NoIR_Obeso_T1 
#               29                15                 5                31 
#      IR_Obeso_T1 
#               20 

####################################################################

factored.design = model.matrix(~ 0 + HOMAOB + COVARS_FENOS_CONTINUOS_PREPUBER$counts.CD8T + COVARS_FENOS_CONTINUOS_PREPUBER$counts.CD4T + COVARS_FENOS_CONTINUOS_PREPUBER$counts.NK + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Bcell + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Mono + COVARS_FENOS_CONTINUOS_PREPUBER$counts.Neu + COVARS_FENOS_CONTINUOS_PREPUBER$Sex_Time + COVARS_FENOS_CONTINUOS_PREPUBER$Age_Time + factor(COVARS_FENOS_CONTINUOS_PREPUBER$Origen_T1))
colnames(factored.design)[1:5] <- c("NoIRNormopeso","NoIRSobrepeso","IRSobrepeso","NoIRObeso","IRObeso")
 colnames(factored.design)[6:15]  <- c("CD8T","CD4T","NK","BCELL","MONO","NEU","SEX","AGE","ORIGEN1","ORIGEN2")

data.fit = lmFit(MICROARRAY_PREPUBER,factored.design)
colnames(data.fit)



#Now we can make any comparisons between the experimental conditions in the usual way, example:
cm <- makeContrasts(NoIRNormopesoVSNoIRSobrepeso=NoIRNormopeso-NoIRSobrepeso,NoIRNormopesoVSIRSobrepeso=NoIRNormopeso-IRSobrepeso,NoIRNormopesoVSNoIRObeso=NoIRNormopeso-NoIRObeso,NoIRNormopesoVSIRObeso=NoIRNormopeso-IRObeso,NoIRSobrepesoVSIRSobrepeso=NoIRSobrepeso-IRSobrepeso,NoIRSobrepesoVSNoIRObeso=NoIRSobrepeso-NoIRObeso,NoIRSobrepesoVSIRObeso=NoIRSobrepeso-IRObeso,IRSobrepesoVSNoIRObeso=IRSobrepeso-NoIRObeso,IRSobrepesoVSIRObeso=IRSobrepeso-IRObeso,NoIRObesoVSIRObeso=NoIRObeso-IRObeso, levels=factored.design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(data.fit, cm)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))


ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)




# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBIR_UNIDOS/NoIRNormopesoVSNoIRSobrepeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)



# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=2, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBIR_UNIDOS/NoIRNormopesoVSIRSobrepeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_FDR05rawP0001.csv", row.names=FALSE)




# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=3, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBIR_UNIDOS/NoIRNormopesoVSNoIRObeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)



# seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=4, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBIR_UNIDOS/NoIRNormopesoVSIRObeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)




  # seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=5, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBIR_UNIDOS/NoIRSobrepesoVSIRSobrepeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_FDR05andrawP0001.csv", row.names=FALSE)




  # seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=6, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBIR_UNIDOS/NoIRSobrepesoVSNoIRObeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)

 
  # seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=7, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBIR_UNIDOS/NoIRSobrepesoVSIRObeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_FDR05rawP0001.csv", row.names=FALSE)



  # seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=8, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBIR_UNIDOS/IRSobrepesoVSNoIRObeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_FDR05rawP0001.csv", row.names=FALSE)


  # seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=9, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)
DMPs_significants$UCSC_RefGene_Name

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBIR_UNIDOS/IRSobrepesoVSIRObeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_FDR05rawP0001.csv", row.names=FALSE)



  # seleccion de significativas
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=10, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT3_PREPUBERTAL/GRUPOS_OBIR_UNIDOS/NoIRObesoVSIRObeso_PREPUBER_adjusted_LEUCOS_AGE_SEX_rawP0001.csv", row.names=FALSE)



