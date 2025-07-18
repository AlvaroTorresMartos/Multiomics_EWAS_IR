





AJUSTANDO por edad:





######################################################################################################
	WITHIN GROUPS COMPARISONS (T1 vs T0), de acuerdo a manual LIMMA
######################################################################################################

#continuar con este script tras lanzar el principal hasta la linea 550 de codigo DEL SCRIPT #1.



Treat <- factor(paste(Time,Group,sep="."))
counts.CD8T <- COVARS_FOR_ANALYSIS$counts.CD8T
counts.CD4T <- COVARS_FOR_ANALYSIS$counts.CD4T
counts.NK <- COVARS_FOR_ANALYSIS$counts.NK
counts.Bcell <- COVARS_FOR_ANALYSIS$counts.Bcell
counts.Mono <- COVARS_FOR_ANALYSIS$counts.Mono
counts.Neu <- COVARS_FOR_ANALYSIS$counts.Neu
Sex <- COVARS_FOR_ANALYSIS$Sex_Time
Sex
Age <- COVARS_FOR_ANALYSIS$Age_Time
origen <- as.factor(COVARS_FOR_ANALYSIS$Origen_T1)


design <- model.matrix(~0+Treat+counts.CD8T+counts.CD4T+counts.NK+counts.Bcell+counts.Mono+counts.Neu+Sex+Age+origen) 
colnames(design)[1:12] <- levels(Treat)


 colnames(design)
 [1] "Baseline.NW non-IR no change"       "Baseline.OB/OW IR no change"       
 [3] "Baseline.OB/OW IR to non-IR"        "Baseline.OB/OW non-IR no change"   
 [5] "Baseline.OB/OW non-IR to IR"        "Baseline.OB/OW non-IR to NW non-IR"
 [7] "T1.NW non-IR no change"             "T1.OB/OW IR no change"             
 [9] "T1.OB/OW IR to non-IR"              "T1.OB/OW non-IR no change"         
[11] "T1.OB/OW non-IR to IR"              "T1.OB/OW non-IR to NW non-IR"      
[13] "counts.CD8T"                        "counts.CD4T"                       
[15] "counts.NK"                          "counts.Bcell"                      
[17] "counts.Mono"                        "counts.Neu"                        
[19] "SexVaron"                           "Age"                               
[21] "origen1"                            "origen2"              

 colnames(design) <- gsub("/","", colnames(design))

 colnames(design) <- gsub("-","", colnames(design))

 colnames(design) <- gsub(" ","", colnames(design))

 colnames(design) <- gsub("Treat","", colnames(design))

 colnames(design)
 [1] "Baseline.NWnonIRnochange"    "Baseline.OBOWIRnochange"    
 [3] "Baseline.OBOWIRtononIR"      "Baseline.OBOWnonIRnochange" 
 [5] "Baseline.OBOWnonIRtoIR"      "Baseline.OBOWnonIRtoNWnonIR"
 [7] "T1.NWnonIRnochange"          "T1.OBOWIRnochange"          
 [9] "T1.OBOWIRtononIR"            "T1.OBOWnonIRnochange"       
[11] "T1.OBOWnonIRtoIR"            "T1.OBOWnonIRtoNWnonIR"      
[13] "counts.CD8T"                 "counts.CD4T"                
[15] "counts.NK"                   "counts.Bcell"               
[17] "counts.Mono"                 "counts.Neu"                 
[19] "SexVaron"                    "Age"                        
[21] "origen1"                     "origen2"      

#Then we estimate the correlation between measurements made on the same subject:
corfit <- duplicateCorrelation(MICROARRAY,design,block=Patient)
corfit$consensus
#Then this inter-subject correlation is input into the linear model fit:
fit <- lmFit(MICROARRAY,design,block=Patient,correlation=corfit$consensus)

#Now we can make any comparisons between the experimental conditions in the usual way, example:
cm <- makeContrasts(T1_VS_BASELINE_NW = T1.NWnonIRnochange-Baseline.NWnonIRnochange, T1_VS_BASELINE_obIRnochange = T1.OBOWIRnochange-Baseline.OBOWIRnochange, T1_VS_BASELINE_obIRtoNOIR = T1.OBOWIRtononIR-Baseline.OBOWIRtononIR, T1_VS_BASELINE_obNOIRnochange = T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange, T1_VS_BASELINE_obNOIRtoIR = T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR, T1_VS_BASELINE_OBOWnonIRtoNWnonIR= T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR, levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))



############### PARA SACAR LISTA DE SIGNIFICATIVAS


colnames(summary(decideTests(data.fit.eb)))
[1] "T1_VS_BASELINE_NW"             "T1_VS_BASELINE_obIRnochange"  
[3] "T1_VS_BASELINE_obIRtoNOIR"     "T1_VS_BASELINE_obNOIRnochange"
[5] "T1_VS_BASELINE_obNOIRtoIR"     "T1_VS_BASELINE_OBOWnonIRtoNWnonIR"


############PARA VER QUE SONDAS HAN SALIDO EN CADA COMPARACION (CAMBIAR EL COEF)


# get the 850k annotation data
ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1,genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_NW_non_IR_no_change_WITHIN_RAW0001_age.csv", row.names=FALSE)



 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=2, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_OB_OW_IR_no_change_WITHIN_FDR05RAW0001_age.csv", row.names=FALSE)





 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=3, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_OB_OW_IR_to_non_IR_WITHIN_raw0001_age.csv", row.names=FALSE)




 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=4, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_OB_OW_non_IR_no_change_WITHIN_raw0001_age.csv", row.names=FALSE)



 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=5, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_OB_OW_non_IR_to_IR_WITHIN_raw0001_age.csv", row.names=FALSE)




 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=6, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_OBOWnonIRtoNWnonIR_WITHIN_FDR05_and_raw0001_age.csv", row.names=FALSE)














 colnames(design)
 [1] "Baseline.NWnonIRnochange"    "Baseline.OBOWIRnochange"    
 [3] "Baseline.OBOWIRtononIR"      "Baseline.OBOWnonIRnochange" 
 [5] "Baseline.OBOWnonIRtoIR"      "Baseline.OBOWnonIRtoNWnonIR"
 [7] "T1.NWnonIRnochange"          "T1.OBOWIRnochange"          
 [9] "T1.OBOWIRtononIR"            "T1.OBOWnonIRnochange"       
[11] "T1.OBOWnonIRtoIR"            "T1.OBOWnonIRtoNWnonIR"      
[13] "counts.CD8T"                 "counts.CD4T"                
[15] "counts.NK"                   "counts.Bcell"               
[17] "counts.Mono"                 "counts.Neu"                 
[19] "SexVaron"                    "Age"                        
[21] "origen1"                     "origen2" 









############ NORMO VS OBIRNOCHANGE 

con.t1vt2<-makeContrasts(((T1.NWnonIRnochange-Baseline.NWnonIRnochange)-(T1.OBOWIRnochange-Baseline.OBOWIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_NORMO VS OBIRNOCHANGE_FDR05RAW0001_AGE.csv", row.names=FALSE)




############ NORMO VS obIRtoNOIR 

con.t1vt2<-makeContrasts(((T1.NWnonIRnochange-Baseline.NWnonIRnochange)-(T1.OBOWIRtononIR-Baseline.OBOWIRtononIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_NORMO VS obIRtoNOIR_RAW0001_AGE.csv", row.names=FALSE)





############ NORMO VS obnonIRnochange 

con.t1vt2<-makeContrasts(((T1.NWnonIRnochange-Baseline.NWnonIRnochange)-(T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_NORMO VS obnonIRnochange_RAW0001_AGE.csv", row.names=FALSE)






############ NORMO VS OBOWnonIRtoIR 

con.t1vt2<-makeContrasts(((T1.NWnonIRnochange-Baseline.NWnonIRnochange)-(T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_NORMO VS OBOWnonIRtoIR_RAW0001_AGE.csv", row.names=FALSE)













*************************************************************
*************************************************************
*************************************************************

############ OBIRNOCHANGE VS OBOWIRtononIR 

con.t1vt2<-makeContrasts(((T1.OBOWIRnochange-Baseline.OBOWIRnochange)-(T1.OBOWIRtononIR-Baseline.OBOWIRtononIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

#DMPs_significants <- subset(DMPs, P.Value < 0.0001)
DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBIRNOCHANGE VS OBIRtononIR_RAW0001andFDR05_AGE.csv", row.names=FALSE)




*************************************************************
*************************************************************
*************************************************************


############ OBIRNOCHANGE VS OBOWnonIRnochange

con.t1vt2<-makeContrasts(((T1.OBOWIRnochange-Baseline.OBOWIRnochange)-(T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBIRNOCHANGE VS OBnoIRNOCHANGE_RAW0001andFDR05_AGE.csv", row.names=FALSE)









############ OBIRnoCHANGE VS OBOWnonIRtoIR 

con.t1vt2<-makeContrasts(((T1.OBOWIRnochange-Baseline.OBOWIRnochange)-(T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBIRNOCHANGE VS OBOWnonIRtoIR_RAW0001_AGE.csv", row.names=FALSE)
















############ OBOWIRtononIR VS OBOWnonIRnochange 

con.t1vt2<-makeContrasts(((T1.OBOWIRtononIR-Baseline.OBOWIRtononIR)-(T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWIRtononIR VS OBOWnonIRnochange_fdrRAW0001_AGE.csv", row.names=FALSE)






############ OBOWIRtononIR VS OBOWnonIRtoIR 

con.t1vt2<-makeContrasts(((T1.OBOWIRtononIR-Baseline.OBOWIRtononIR)-(T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))


ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWIRtononIR VS OBOWnonIRtoIR_RAW0001_AGE.csv", row.names=FALSE)













############ OBOWnonIRnochange VS OBOWnonIRtoIR 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange)-(T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))


ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRnochange VS OBOWnonIRtoIR_RAW0001_AGE.csv", row.names=FALSE)


















############ OBOWnonIRtoNWnonIR VS NWnonIRnochange 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR)-(T1.NWnonIRnochange-Baseline.NWnonIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))


ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRtoNWnonIR VS NWnonIRnochange_RAW0001_FDR05_AGE.csv", row.names=FALSE)












############ OBOWnonIRtoNWnonIR VS OBOWIRnochange 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR)-(T1.OBOWIRnochange-Baseline.OBOWIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))


ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

#DMPs_significants <- subset(DMPs, P.Value < 0.0001)
DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRtoNWnonIR VS OBOWIRnochange_FDR05_AGE.csv", row.names=FALSE)











############ OBOWnonIRtoNWnonIR VS OBOWIRtononIR 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR)-(T1.OBOWIRtononIR-Baseline.OBOWIRtononIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))


ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRtoNWnonIR VS OBOWIRtononIR_RAW0001_FDR05_AGE.csv", row.names=FALSE)












############ OBOWnonIRtoNWnonIR VS OBOWnonIRnochange 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR)-(T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))


ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)


 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRtoNWnonIR VS OBOWnonIRnochange_RAW0001_FDR05_AGE.csv", row.names=FALSE)






############ OBOWnonIRtoNWnonIR VS OBOWnonIRtoIR 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR)-(T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))



ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)



 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRtoNWnonIR VS OBOWnonIRtoIR_RAW0001_FDR05_AGE.csv", row.names=FALSE)


































































sin ajustar por edad:





######################################################################################################
	WITHIN GROUPS COMPARISONS (T1 vs T0), de acuerdo a manual LIMMA
######################################################################################################

#continuar con este script tras lanzar el principal hasta la linea 524 de codigo DEL SCRIPT #1.



Treat <- factor(paste(Time,Group,sep="."))
counts.CD8T <- COVARS_FOR_ANALYSIS$counts.CD8T
counts.CD4T <- COVARS_FOR_ANALYSIS$counts.CD4T
counts.NK <- COVARS_FOR_ANALYSIS$counts.NK
counts.Bcell <- COVARS_FOR_ANALYSIS$counts.Bcell
counts.Mono <- COVARS_FOR_ANALYSIS$counts.Mono
counts.Neu <- COVARS_FOR_ANALYSIS$counts.Neu
Sex <- COVARS_FOR_ANALYSIS$Sex_Time

origen <- as.factor(COVARS_FOR_ANALYSIS$Origen_T1)





design <- model.matrix(~0+Treat+counts.CD8T+counts.CD4T+counts.NK+counts.Bcell+counts.Mono+counts.Neu+Sex+origen) 
colnames(design)[1:12] <- levels(Treat)


 colnames(design)
 [1] "Baseline.NW non-IR no change"       "Baseline.OB/OW IR no change"       
 [3] "Baseline.OB/OW IR to non-IR"        "Baseline.OB/OW non-IR no change"   
 [5] "Baseline.OB/OW non-IR to IR"        "Baseline.OB/OW non-IR to NW non-IR"
 [7] "T1.NW non-IR no change"             "T1.OB/OW IR no change"             
 [9] "T1.OB/OW IR to non-IR"              "T1.OB/OW non-IR no change"         
[11] "T1.OB/OW non-IR to IR"              "T1.OB/OW non-IR to NW non-IR"      
[13] "counts.CD8T"                        "counts.CD4T"                       
[15] "counts.NK"                          "counts.Bcell"                      
[17] "counts.Mono"                        "counts.Neu"                        
[19] "SexVaron"                           "origen1"                           
[21] "origen2"                                  

 colnames(design) <- gsub("/","", colnames(design))

 colnames(design) <- gsub("-","", colnames(design))

 colnames(design) <- gsub(" ","", colnames(design))

 colnames(design) <- gsub("Treat","", colnames(design))

 colnames(design)
[1] "Baseline.NWnonIRnochange"    "Baseline.OBOWIRnochange"    
 [3] "Baseline.OBOWIRtononIR"      "Baseline.OBOWnonIRnochange" 
 [5] "Baseline.OBOWnonIRtoIR"      "Baseline.OBOWnonIRtoNWnonIR"
 [7] "T1.NWnonIRnochange"          "T1.OBOWIRnochange"          
 [9] "T1.OBOWIRtononIR"            "T1.OBOWnonIRnochange"       
[11] "T1.OBOWnonIRtoIR"            "T1.OBOWnonIRtoNWnonIR"      
[13] "counts.CD8T"                 "counts.CD4T"                
[15] "counts.NK"                   "counts.Bcell"               
[17] "counts.Mono"                 "counts.Neu"                 
[19] "SexVaron"                    "origen1"                    
[21] "origen2"                  

#Then we estimate the correlation between measurements made on the same subject:
corfit <- duplicateCorrelation(MICROARRAY,design,block=Patient)
corfit$consensus
#Then this inter-subject correlation is input into the linear model fit:
fit <- lmFit(MICROARRAY,design,block=Patient,correlation=corfit$consensus)

#Now we can make any comparisons between the experimental conditions in the usual way, example:
cm <- makeContrasts(T1_VS_BASELINE_NW = T1.NWnonIRnochange-Baseline.NWnonIRnochange, T1_VS_BASELINE_obIRnochange = T1.OBOWIRnochange-Baseline.OBOWIRnochange, T1_VS_BASELINE_obIRtoNOIR = T1.OBOWIRtononIR-Baseline.OBOWIRtononIR, T1_VS_BASELINE_obNOIRnochange = T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange, T1_VS_BASELINE_obNOIRtoIR = T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR, T1_VS_BASELINE_OBOWnonIRtoNWnonIR= T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR, levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))



############### PARA SACAR LISTA DE SIGNIFICATIVAS


colnames(summary(decideTests(data.fit.eb)))
[1] "T1_VS_BASELINE_NW"             "T1_VS_BASELINE_obIRnochange"  
[3] "T1_VS_BASELINE_obIRtoNOIR"     "T1_VS_BASELINE_obNOIRnochange"
[5] "T1_VS_BASELINE_obNOIRtoIR"     "T1_VS_BASELINE_OBOWnonIRtoNWnonIR"


############PARA VER QUE SONDAS HAN SALIDO EN CADA COMPARACION (CAMBIAR EL COEF)

 # get the 850k annotation data
ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(data.fit.eb$coefficients),]

ann850k[,4] == rownames(data.fit.eb$coefficients)




 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1,genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

#DMPs_significants <- subset(DMPs, P.Value < 0.0001)
DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_NW_non_IR_no_change_WITHIN_FDR05.csv", row.names=FALSE)



 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=2, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_OB_OW_IR_no_change_WITHIN_FDR05RAW0001.csv", row.names=FALSE)





 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=3, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_OB_OW_IR_to_non_IR_WITHIN_andFDR05.csv", row.names=FALSE)




 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=4, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)

write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_OB_OW_non_IR_no_change_WITHIN_FDR05.csv", row.names=FALSE)



 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=5, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_OB_OW_non_IR_to_IR_WITHIN_FDR05raw0001.csv", row.names=FALSE)




 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=6, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/WITHIN_SCRIPT2/AJUSTE_EDAD_SEXO_LEUCOS/DMPs_significants_adjusted_OBOWnonIRtoNWnonIR_WITHIN_raw0001andFDR05.csv", row.names=FALSE)

























############ NORMO VS OBIRNOCHANGE 

con.t1vt2<-makeContrasts(((T1.NWnonIRnochange-Baseline.NWnonIRnochange)-(T1.OBOWIRnochange-Baseline.OBOWIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_NORMO VS OBIRNOCHANGE_FDR05RAW0001.csv", row.names=FALSE)




############ NORMO VS obIRtoNOIR 

con.t1vt2<-makeContrasts(((T1.NWnonIRnochange-Baseline.NWnonIRnochange)-(T1.OBOWIRtononIR-Baseline.OBOWIRtononIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_NORMO VS obIRtoNOIR_RAW0001.csv", row.names=FALSE)





############ NORMO VS obnonIRnochange 

con.t1vt2<-makeContrasts(((T1.NWnonIRnochange-Baseline.NWnonIRnochange)-(T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_NORMO VS obnonIRnochange_RAW0001.csv", row.names=FALSE)






############ NORMO VS OBOWnonIRtoIR 

con.t1vt2<-makeContrasts(((T1.NWnonIRnochange-Baseline.NWnonIRnochange)-(T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_NORMO VS OBOWnonIRtoIR_RAW0001.csv", row.names=FALSE)

















############ OBIRNOCHANGE VS OBOWIRtononIR 

con.t1vt2<-makeContrasts(((T1.OBOWIRnochange-Baseline.OBOWIRnochange)-(T1.OBOWIRtononIR-Baseline.OBOWIRtononIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

#DMPs_significants <- subset(DMPs, P.Value < 0.0001)
DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBIRNOCHANGE VS OBOWIRtononIR_FDR05.csv", row.names=FALSE)













############ OBIRNOCHANGE VS OBOWnonIRnochange 

con.t1vt2<-makeContrasts(((T1.OBOWIRnochange-Baseline.OBOWIRnochange)-(T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBIRNOCHANGE VS OBOWnonIRnochange_RAW0001andFDR05.csv", row.names=FALSE)









############ OBIRnoCHANGE VS OBOWnonIRtoIR 

con.t1vt2<-makeContrasts(((T1.OBOWIRnochange-Baseline.OBOWIRnochange)-(T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBIRNOCHANGE VS OBOWnonIRtoIR_0001RAW.csv", row.names=FALSE)

















############ OBOWIRtononIR VS OBOWnonIRnochange 

con.t1vt2<-makeContrasts(((T1.OBOWIRtononIR-Baseline.OBOWIRtononIR)-(T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWIRtononIR VS OBOWnonIRnochange_RAW0001.csv", row.names=FALSE)






############ OBOWIRtononIR VS OBOWnonIRtoIR 

con.t1vt2<-makeContrasts(((T1.OBOWIRtononIR-Baseline.OBOWIRtononIR)-(T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWIRtononIR VS OBOWnonIRtoIR_RAW0001.csv", row.names=FALSE)













############ OBOWnonIRnochange VS OBOWnonIRtoIR 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange)-(T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRnochange VS OBOWnonIRtoIR_RAW0001.csv", row.names=FALSE)


















############ OBOWnonIRtoNWnonIR VS NWnonIRnochange 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR)-(T1.NWnonIRnochange-Baseline.NWnonIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRtoNWnonIR VS NWnonIRnochange_RAW0001_FDR05.csv", row.names=FALSE)












############ OBOWnonIRtoNWnonIR VS OBOWIRnochange 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR)-(T1.OBOWIRnochange-Baseline.OBOWIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

#DMPs_significants <- subset(DMPs, P.Value < 0.0001)
DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRtoNWnonIR VS OBOWIRnochange_FDR05.csv", row.names=FALSE)











############ OBOWnonIRtoNWnonIR VS OBOWIRtononIR 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR)-(T1.OBOWIRtononIR-Baseline.OBOWIRtononIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRtoNWnonIR VS OBOWIRtononIR_RAW0001_FDR05.csv", row.names=FALSE)












############ OBOWnonIRtoNWnonIR VS OBOWnonIRnochange 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR)-(T1.OBOWnonIRnochange-Baseline.OBOWnonIRnochange)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRtoNWnonIR VS OBOWnonIRnochange_RAW0001_FDR05.csv", row.names=FALSE)








############ OBOWnonIRtoNWnonIR VS OBOWnonIRtoIR 

con.t1vt2<-makeContrasts(((T1.OBOWnonIRtoNWnonIR-Baseline.OBOWnonIRtoNWnonIR)-(T1.OBOWnonIRtoIR-Baseline.OBOWnonIRtoIR)), levels=design)
#Then compute these contrasts and moderated t-tests:
fit2 <- contrasts.fit(fit, con.t1vt2)
fit2 <- eBayes(fit2)
data.fit.eb <- fit2
summary(decideTests(data.fit.eb))

 # ANNOT PROBES AND GET SIGNIFICANT cpg RESULTS ORDERED
ann850kSub <- ann850k[,
c(1:4,12:19,22,23,25:ncol(ann850k))]
DMPs <- topTable(data.fit.eb, num=Inf, coef=1, genelist=ann850kSub,adjust.method='fdr')
dim(DMPs)
head(DMPs)

DMPs_significants <- subset(DMPs, P.Value < 0.0001)
#DMPs_significants <- subset(DMPs, adj.P.Val < 0.05)
dim(DMPs_significants)


write.csv2(DMPs_significants, file="/home/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/RESULTADOS_ANALISIS/SCRIPT2/BETWEEN_SCRIPT2/between_OBOWnonIRtoNWnonIR VS OBOWnonIRtoIR_RAW0001_FDR05.csv", row.names=FALSE)






























#Treat <- factor(paste(Time,Group,sep="."))
#counts.CD8T <- COVARS_FOR_ANALYSIS$counts.CD8T
#counts.CD4T <- COVARS_FOR_ANALYSIS$counts.CD4T
#counts.NK <- COVARS_FOR_ANALYSIS$counts.NK
#counts.Bcell <- COVARS_FOR_ANALYSIS$counts.Bcell
#counts.Mono <- COVARS_FOR_ANALYSIS$counts.Mono
#counts.Neu <- COVARS_FOR_ANALYSIS$counts.Neu
#Sex <- COVARS_FOR_ANALYSIS$Sex_Time
#Sex
#Age <- COVARS_FOR_ANALYSIS$Age_Time
#origen <- as.factor(COVARS_FOR_ANALYSIS$Origen_T1)


#design <- model.matrix(~0+Treat+counts.CD8T+counts.CD4T+counts.NK+counts.Bcell+counts.Mono+counts.Neu+Sex+Age+origen) 
#colnames(design)[1:12] <- levels(Treat)


# colnames(design)

#design[1:5,1:9]




#designS <- model.matrix(~0+Time:Group+counts.CD8T+counts.CD4T+counts.NK+counts.Bcell+counts.Mono+counts.Neu+Sex+Age+origen) 
#colnames(designS)

#designS <- designS[,c(12:23,1:11)]
#colnames(designS)
#designS <- designS[,c(1,11,7,5,9,3,2,12,8,6,10,4,13:23)]
#colnames(designS)
#design[1:5,1:12] == designS[1:5,1:12]









