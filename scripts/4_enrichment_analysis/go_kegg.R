



# Instalar y cargar paquetes necesarios
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pacman::p_load(clusterProfiler, org.Hs.eg.db, KEGGREST, tidyr, dplyr)



# 1️ Obtener genes asociados a términos GO
go_terms = keys(org.Hs.eg.db, keytype = "GO")
go_mapping = AnnotationDbi::select(org.Hs.eg.db, keys = go_terms, 
                                    keytype = "GO", columns = c("SYMBOL", "GO"))

# Agrupar genes por término GO
go_genes = go_mapping %>%
  group_by(GO) %>%
  summarise(Genes = list(unique(SYMBOL))) %>%
  ungroup()

go_genes = go_genes %>%
  mutate(Genes = sapply(Genes, function(g) paste(g, collapse = ", ")))

# Renombrar columna GO -> Term
colnames(go_genes) = c("Term", "Gene_Symbols")

# 2️ Obtener genes asociados a vías KEGG
kegg_paths = keggList("pathway", "hsa")  # Obtener todas las rutas KEGG para humanos
kegg_ids = names(kegg_paths)

# Obtener genes asociados a cada vía KEGG
kegg_mapping = lapply(kegg_ids, function(pid) {
  genes <- keggGet(pid)[[1]]$GENE  # Extrae genes de la vía KEGG
  if (!is.null(genes)) {
    genes <- genes[seq(1, length(genes), 2)]  # Extrae solo los nombres de genes
  }
  data.frame(Term = pid, Genes = I(list(genes)))
})

# Convertir lista en dataframe
kegg_genes = bind_rows(kegg_mapping)

kegg_genes = kegg_genes %>%
  mutate(Genes = sapply(Genes, function(g) paste(g, collapse = ", ")))

# Función para convertir IDs de Entrez a nombres de genes
convert_entrez_to_symbol <- function(entrez_list) {
  if (is.na(entrez_list) || entrez_list == "") {
    return("")  # Si la fila está vacía, devolver cadena vacía
  }
  
  entrez_list <- unlist(strsplit(entrez_list, ", "))  # Separar IDs en una lista
  # Verificar que la lista no esté vacía después de separar
  if (length(entrez_list) == 0) {
    return("")
  }
  
  symbols <- mapIds(org.Hs.eg.db, keys = entrez_list, 
                    column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  symbols <- symbols[!is.na(symbols)]  # Eliminar valores NA

  # Si después de la conversión no hay genes, devolver cadena vacía
  if (length(symbols) == 0) {
    return("")
  }
  
  return(paste(symbols, collapse = ", "))  # Unir los nombres en un solo string
}

# Aplicar la conversión a cada fila de 'kegg_genes'
kegg_genes <- kegg_genes %>%
  rowwise() %>%
  mutate(Gene_Symbols = ifelse(Genes != "", convert_entrez_to_symbol(Genes), "")) %>%
  ungroup()

# Ver los primeros resultados
head(kegg_genes)

# saveRDS(go_genes, "./data/processed/GO_genes.RDS")
# saveRDS(kegg_genes, "./data/processed/KEGG_genes.RDS")
