
# Packages -----
library(dplyr)
library(janitor)
library(Hmisc)
library(CMplot)
library(ggrepel)

# Import data -----
lista_def = readxl::read_xlsx("./2024_03_13_lista_def_CpGs_augusto_ampliada_288filas_193CpGs_89genes.xlsx")

lista_def2 = lista_def %>% 
  mutate(
    gene = paste0(Name, "_", UCSC_RefGene_Name),
    gene = gsub("sep-01", "SEPT1", gene),
    gene = gsub("_NA", "", gene)
  ) %>% 
  clean_names() %>% 
  select(chr, gene, pos, p_value, approach_group_comparison)


g5vsg3_between = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Longitudinal G5 vs. G3") %>% 
  dplyr::select(chr, gene, pos, p_value) 

g3vsg4_between = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Longitudinal G3 vs. G4") %>% 
  dplyr::select(chr, gene, pos, p_value) 

g2vsg4_between = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Longitudinal G2 vs. G4") %>% 
  dplyr::select(chr, gene, pos, p_value) 

g3_within = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Longitudinal G3 (within)") %>% 
  dplyr::select(chr, gene, pos, p_value) 

g4_within = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Longitudinal G4 (within)") %>% 
  dplyr::select(chr, gene, pos, p_value) 

nw_owob_cross = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Pubertal non-IR NW vs. IR Obese & Overweight") %>% 
  dplyr::select(chr, gene, pos, p_value) 

owob_owb_cross = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Pubertal non-IR Obese & Overweight vs. IR Obese & Overweight") %>% 
  dplyr::select(chr, gene, pos, p_value) 

noir_ir_cross = lista_def2 %>% 
  dplyr::filter(approach_group_comparison == "Pubertal IR vs. non-IR") %>% 
  dplyr::select(chr, gene, pos, p_value)

allgenes2 = read.csv("./allgenes_manhattanplot.csv", row.names = 1)


# Longitudinal within G3 ------

df_original = rbind(g3_within, allgenes2)

df_original = df_original %>% 
  dplyr::select(gene, chr, pos, p_value) %>%
  dplyr::mutate(chr = gsub("chr", "", chr), 
                chr = as.numeric(chr)) %>% 
  dplyr::arrange(chr)

# Definimos los umbrales de significancia
genomewide_threshold = 5e-8
suggestive_threshold = 1e-5      

highlight_genes = df_original %>%
  filter(p_value < suggestive_threshold) 


CMplot(df_original, plot.type = "m", col=c("grey30","grey60"),
       LOG10 = TRUE, band = 0, 
       threshold = suggestive_threshold, # c(suggestive_threshold, genomewide_threshold),
       threshold.col = "blue", # c("blue", "red"),  
       threshold.lty = 1, #c(1, 1),           
       highlight = highlight_genes$gene,
       highlight.text = highlight_genes$gene,
       highlight.text.col = "black",
       highlight.text.cex = 1,
       highlight.col = "red",    
       amplify = TRUE,                    
       ylim = c(0, 7),                   
       signal.col = "red", 
       signal.cex = 1,
       cex = 0.3, 
       file = "png", 
       file.name = "_g3_within", 
       main = "Longitudinal Changes within G3"
)


# Longitudinal within G4 ------

df_original = rbind(g4_within, allgenes2)

df_original = df_original %>% 
  dplyr::select(gene, chr, pos, p_value) %>%
  dplyr::mutate(chr = gsub("chr", "", chr), 
                chr = as.numeric(chr)) %>% 
  dplyr::arrange(chr)

# Definimos los umbrales de significancia
genomewide_threshold = 5e-8
suggestive_threshold = 1e-5      

highlight_genes = df_original %>%
  filter(p_value < suggestive_threshold) 


CMplot(df_original, plot.type = "m", col=c("grey30","grey60"),
       LOG10 = TRUE,  band = 0, 
       threshold = suggestive_threshold, # c(suggestive_threshold, genomewide_threshold),
       threshold.col = "blue", # c("blue", "red"),  
       threshold.lty = 1, #c(1, 1),           
       highlight = highlight_genes$gene,
       highlight.text = highlight_genes$gene,
       highlight.text.col = "black",
       highlight.text.cex = 1,
       highlight.col = "red",    
       amplify = TRUE,                    
       ylim = c(0, 7),                   
       signal.col = "red", 
       signal.cex = 1,
       cex = 0.3, 
       file = "png", 
       file.name = "_g4_within", 
       main = "Longitudinal Changes within G4"
)


# Longitudinal G5 vs G3 between------

df_original = rbind(g5vsg3_between, allgenes2)

df_original = df_original %>% 
  dplyr::select(gene, chr, pos, p_value) %>%
  dplyr::mutate(chr = gsub("chr", "", chr), 
                chr = as.numeric(chr)) %>% 
  dplyr::arrange(chr)

# Definimos los umbrales de significancia
genomewide_threshold = 5e-8
suggestive_threshold = 1e-5      

highlight_genes = df_original %>%
  filter(p_value < suggestive_threshold) 


CMplot(df_original, plot.type = "m", col=c("grey30","grey60"),
       LOG10 = TRUE, band = 0, 
       threshold = suggestive_threshold, # c(suggestive_threshold, genomewide_threshold),
       threshold.col = "blue", # c("blue", "red"),  
       threshold.lty = 1, #c(1, 1),           
       highlight = highlight_genes$gene,
       highlight.text = highlight_genes$gene,
       highlight.text.col = "black",
       highlight.text.cex = 1,
       highlight.col = "red",    
       amplify = TRUE,                    
       ylim = c(0, 7),                   
       signal.col = "red", 
       signal.cex = 1,
       cex = 0.3, 
       file = "png", 
       file.name = "_g5vsg3_between", 
       main = "Longitudinal Changes between G3 vs G5"
)


# Longitudinal G3 vs G4 between ------

df_original = rbind(g3vsg4_between, allgenes2)

df_original = df_original %>% 
  dplyr::select(gene, chr, pos, p_value) %>%
  dplyr::mutate(chr = gsub("chr", "", chr), 
                chr = as.numeric(chr)) %>% 
  dplyr::arrange(chr)

# Definimos los umbrales de significancia
genomewide_threshold = 5e-8
suggestive_threshold = 1e-5      

highlight_genes = df_original %>%
  filter(p_value < suggestive_threshold) 


CMplot(df_original, plot.type = "m", col=c("grey30","grey60"),
       LOG10 = TRUE, band = 0, 
       threshold = suggestive_threshold, # c(suggestive_threshold, genomewide_threshold),
       threshold.col = "blue", # c("blue", "red"),  
       threshold.lty = 1, #c(1, 1),           
       highlight = highlight_genes$gene,
       highlight.text = highlight_genes$gene,
       highlight.text.col = "black",
       highlight.text.cex = 1,
       highlight.col = "red",    
       amplify = TRUE,                    
       ylim = c(0, 7),                   
       signal.col = "red", 
       signal.cex = 1,
       cex = 0.3, 
       file = "png", 
       file.name = "_g3vsg4_between", 
       main = "Longitudinal Changes between G3 vs G4"
)


# Longitudinal G2 vs G4 between ------

df_original = rbind(g2vsg4_between, allgenes2)

df_original = df_original %>% 
  dplyr::select(gene, chr, pos, p_value) %>%
  dplyr::mutate(chr = gsub("chr", "", chr), 
                chr = as.numeric(chr)) %>% 
  dplyr::arrange(chr)

# Definimos los umbrales de significancia
genomewide_threshold = 5e-8
suggestive_threshold = 1e-5      

highlight_genes = df_original %>%
  filter(p_value < suggestive_threshold) 


CMplot(df_original, plot.type = "m", col=c("grey30","grey60"),
       LOG10 = TRUE, band = 0, 
       threshold = suggestive_threshold, # c(suggestive_threshold, genomewide_threshold),
       threshold.col = "blue", # c("blue", "red"),  
       threshold.lty = 1, #c(1, 1),           
       highlight = highlight_genes$gene,
       highlight.text = highlight_genes$gene,
       highlight.text.col = "black",
       highlight.text.cex = 1,
       highlight.col = "red",    
       amplify = TRUE,                    
       ylim = c(0, 7),                   
       signal.col = "red", 
       signal.cex = 1,
       cex = 0.3, 
       file = "png", 
       file.name = "_g2vsg4_between", 
       main = "Longitudinal Changes between G2 vs G4"
)


# Cross-sectional Pubertal noIR NW vs IR obesity and overweight ------

df_original = rbind(nw_owob_cross, allgenes2)

df_original = df_original %>% 
  dplyr::select(gene, chr, pos, p_value) %>%
  dplyr::mutate(chr = gsub("chr", "", chr), 
                chr = as.numeric(chr)) %>% 
  dplyr::arrange(chr)

# Definimos los umbrales de significancia
genomewide_threshold = 5e-8
suggestive_threshold = 1e-5      

highlight_genes = df_original %>%
  filter(p_value < suggestive_threshold) 


CMplot(df_original, plot.type = "m", col=c("grey30","grey60"),
       LOG10 = TRUE, band = 0, 
       threshold = suggestive_threshold, # c(suggestive_threshold, genomewide_threshold),
       threshold.col = "blue", # c("blue", "red"),  
       threshold.lty = 1, #c(1, 1),           
       highlight = highlight_genes$gene,
       highlight.text = highlight_genes$gene,
       highlight.text.col = "black",
       highlight.text.cex = 1,
       highlight.col = "red",    
       amplify = TRUE,                    
       ylim = c(0, 7),                   
       signal.col = "red", 
       signal.cex = 1,
       cex = 0.3, 
       file = "png", 
       file.name = "_nw_owob_cross", 
       main = "Pubertal Cross-sectional Changes between noIR (Normoweight) vs IR (Overweight/Obesity)"
)


# Cross-sectional Pubertal noIR Obesity and overweight vs IR Obesity and overweight------

df_original = rbind(owob_owb_cross, allgenes2)

df_original = df_original %>% 
  dplyr::select(gene, chr, pos, p_value) %>%
  dplyr::mutate(chr = gsub("chr", "", chr), 
                chr = as.numeric(chr)) %>% 
  dplyr::arrange(chr)

# Definimos los umbrales de significancia
genomewide_threshold = 5e-8
suggestive_threshold = 1e-5      

highlight_genes = df_original %>%
  filter(p_value < suggestive_threshold) 


CMplot(df_original, plot.type = "m", col=c("grey30","grey60"),
       LOG10 = TRUE, band = 0, 
       threshold = suggestive_threshold, # c(suggestive_threshold, genomewide_threshold),
       threshold.col = "blue", # c("blue", "red"),  
       threshold.lty = 1, #c(1, 1),           
       highlight = highlight_genes$gene,
       highlight.text = highlight_genes$gene,
       highlight.text.col = "black",
       highlight.text.cex = 1,
       highlight.col = "red",    
       amplify = TRUE,                    
       ylim = c(0, 7),                   
       signal.col = "red", 
       signal.cex = 1,
       cex = 0.3, 
       file = "png", 
       file.name = "_owob_owb_cross", 
       main = "Pubertal Cross-sectional Changes between noIR (Overweight/Obesity) vs IR (Overweight/Obesity)"
)


# Cross-sectional Pubertal noIR vs IR  ------

df_original = rbind(noir_ir_cross, allgenes2)

df_original = df_original %>% 
  dplyr::select(gene, chr, pos, p_value) %>%
  dplyr::mutate(chr = gsub("chr", "", chr), 
                chr = as.numeric(chr)) %>% 
  dplyr::arrange(chr)

# Definimos los umbrales de significancia
genomewide_threshold = 5e-8
suggestive_threshold = 1e-5      

highlight_genes = df_original %>%
  filter(p_value < suggestive_threshold) 


CMplot(df_original, plot.type = "m", col=c("grey30","grey60"),
       LOG10 = TRUE, band = 0, 
       threshold = suggestive_threshold, # c(suggestive_threshold, genomewide_threshold),
       threshold.col = "blue", # c("blue", "red"),  
       threshold.lty = 1, #c(1, 1),           
       highlight = highlight_genes$gene,
       highlight.text = highlight_genes$gene,
       highlight.text.col = "black",
       highlight.text.cex = 1,
       highlight.col = "red",    
       amplify = TRUE,                    
       ylim = c(0, 7),                   
       signal.col = "red", 
       signal.cex = 1,
       cex = 0.3, 
       file = "png", 
       file.name = "_noir_ir_cross", 
       main = "Pubertal Cross-sectional Changes between noIR vs IR"
       
)
