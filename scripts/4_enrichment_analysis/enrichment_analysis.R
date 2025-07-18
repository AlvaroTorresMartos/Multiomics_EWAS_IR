

if (!require("pacman")) install.packages("pacman")

pacman::p_load(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, 
               IlluminaHumanMethylationEPICmanifest, 
               org.Hs.eg.db, GO.db, KEGGREST, limma, missMethyl, tibble, dplyr
               )

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
all = ann850k@listData$Name

cpgs = readxl::read_xlsx("./data/2024_03_13_lista_def_CpGs_augusto_ampliada_288filas_193CpGs_89genes.xlsx")


# GO
go =  gometh(sig.cpg=cpgs$Name, all.cpg=all, 
               array.type="EPIC", collection="GO") 

# head(go[order(go$P.DE), ], n=10)
# go_sig = subset(go, P.DE < 0.01)
# dim(go_sig)
#[1] 1   6



# KEGG
kegg =  gometh(sig.cpg=cpgs$Name, all.cpg=all, 
             array.type="EPIC", collection="KEGG") 

# kegg_sig = subset(kegg, P.DE < 0.01)
# dim(kegg_sig)
#[1] 1   6

# Script go_kegg.R
go_genes = readRDS("./data/processed/GO_genes.RDS")
kegg_genes = readRDS("./data/processed/KEGG_genes.RDS")


go = go %>% 
  tibble::rownames_to_column(var = "Term") %>% 
  dplyr::left_join(go_genes, by = "Term")

kegg = kegg %>% 
  tibble::rownames_to_column(var = "Term") %>% 
  dplyr::left_join(kegg_genes, by = "Term")

rm(go_genes, kegg_genes, all, ann850k)


genes = cpgs %>% 
  dplyr::distinct(UCSC_RefGene_Name) %>% 
  dplyr::select(UCSC_RefGene_Name) %>% 
  na.omit()


find_differential_genes = function(gene_list, differential_genes = genes$UCSC_RefGene_Name) {
  if (is.na(gene_list) || gene_list == "") {
    return("")  
  }
  gene_list <- unlist(strsplit(as.character(gene_list), ",\\s*"))  
  matched_genes = intersect(gene_list, differential_genes)
                            
  if (length(matched_genes) == 0) {
    return("")
  }

  return(paste(matched_genes, collapse = ", "))
}


go_def = go %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Differential_genes = find_differential_genes(Genes)) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(P.DE < 0.1)

kegg_def = kegg %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Differential_genes = find_differential_genes(Gene_Symbols)) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(P.DE < 0.1)

go_def = go_def %>% 
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::select(Term:FDR, Differential_genes) %>% 
  dplyr::rename(GENES_DE = Differential_genes) %>% 
  dplyr::arrange(P.DE)

kegg_def = kegg_def %>% 
  dplyr::mutate(Term = paste0("path:", Term, sep = "")) %>% 
  dplyr::select(Term:FDR, Differential_genes) %>% 
  dplyr::rename(GENES_DE = Differential_genes) %>% 
  dplyr::arrange(P.DE)

# writexl::write_xlsx(kegg_def, "./results/Enrichment_analysis/Additional file 8.xlsx")
# writexl::write_xlsx(go_def, "./results/Enrichment_analysis/Additional file 9.xlsx")

