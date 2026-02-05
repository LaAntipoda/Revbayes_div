#metadata en tabla para DEC 
library(ape)
library(ggtree)
library(treeio)
library(tidyverse)
library(phytools)

# Cargar datos
datos <- read.csv("/home/daniela/GSM/Alineamientos/Alineamientos/nuevos/clusters/metadata_clusters/label_metadatos_tree.csv")
subarbol <- read.nexus("/home/daniela/GSM/Revbayes/Rev_in/98_2601_ingroup.nex")

#filtrar a solo los OTUs del arbol
# Obtener los nombres de las terminales
nombres_tips <- subarbol$tip.label

datos <- datos %>%
  filter(OTU %in% nombres_tips)

faltantes <- nombres_tips[!nombres_tips %in% datos$OTU]

######
recodificacion <- c(
  "NAA" = "U",
  "NT" = "T",
  "IM" = "A",
  "AA" = "A",
  "AT" = "A",
  "PA" = "E"
  # Agrega más si necesitas: "viejo" = "nuevo"
)

for (viejo in names(recodificacion)) {
  nuevo <- recodificacion[viejo]
  datos$REALM <- gsub(paste0("\\b", viejo, "\\b"), nuevo, datos$REALM)
}


# PASO 2: Obtener realms únicos por OTU (sin strsplit)
datos_colapsados <- datos %>%
  group_by(OTU) %>%
  summarise(
    realms = paste(unique(REALM), collapse = ", "),
    n_realms = length(unique(REALM)),
    n_sitios = n()
  ) %>%
  ungroup()

# PASO 3: Identificar todos los realms únicos en el dataset
todos_realms <- unique(recodificacion)
todos_realms <- sort(todos_realms)  # Ordenar alfabéticamente

# PASO 4: Crear matriz binaria
n_taxa <- nrow(datos_colapsados)
n_realms <- length(todos_realms)

# Inicializar matriz con ceros
matriz_binaria <- matrix(0, nrow = n_taxa, ncol = n_realms)
rownames(matriz_binaria) <- datos_colapsados$OTU

# Llenar matriz - marcar presencias
for (i in 1:nrow(datos)) {
  otu <- datos$OTU[i]
  realm <- datos$REALM[i]
  
  # Encontrar índices
  fila_index <- which(datos_colapsados$OTU == otu)
  col_index <- which(todos_realms == realm)
  
  if (length(fila_index) > 0 && length(col_index) > 0) {
    matriz_binaria[fila_index, col_index] <- 1
  }
}

# Convertir a strings
matriz_strings <- apply(matriz_binaria, 1, paste, collapse = "")

# PASO 5: Crear archivo NEXUS
nexus_file <- "/home/daniela/GSM/Revbayes/DEC/data/matriz_dec_realms_98_2601_ingroup.nex"

cat("#NEXUS\n", file = nexus_file)
cat("Begin data;\n", file = nexus_file, append = TRUE)
cat(sprintf("Dimensions ntax=%d nchar=%d;\n", n_taxa, n_realms), 
    file = nexus_file, append = TRUE)
cat('Format datatype=Standard missing=? gap=- labels="01";\n', 
    file = nexus_file, append = TRUE)
cat("Matrix\n", file = nexus_file, append = TRUE)

# Escribir encabezado con códigos de realms
cat(sprintf("    [%s]\n", paste(todos_realms, collapse = "")), 
    file = nexus_file, append = TRUE)

# Escribir cada taxón con su codificación
for (i in 1:n_taxa) {
  taxon_name <- datos_colapsados$OTU[i]
  # Ajustar espaciado para alineación
  espacios <- paste(rep(" ", max(50 - nchar(taxon_name), 1)), collapse = "")
  cat(sprintf("    %s%s%s\n", taxon_name, espacios, matriz_strings[i]), 
      file = nexus_file, append = TRUE)
}

cat("    ;\n", file = nexus_file, append = TRUE)
cat("End;\n", file = nexus_file, append = TRUE)

print(paste("Archivo NEXUS creado:", nexus_file))
print(paste("Número de taxa:", n_taxa))
print(paste("Número de realms:", n_realms))

