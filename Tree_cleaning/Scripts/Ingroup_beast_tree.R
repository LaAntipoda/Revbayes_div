#preparar un arbol de beast para revbayes

# Cargar los paquetes necesarios
library(ape)
library(ggtree)
library(treeio)
library(tidyverse)
library(phytools)

# Definir la ruta del archivo NEXUS
f_beast <- "/home/daniela/GSM/Alineamientos/Alineamientos/nuevos/clusters/alineamientos_cluster/Beast/Calibrar_98_beast_2701_strict.tree"

# Cargar el árbol
arbol <- read.nexus(f_beast) #debe ser en nexus o no lo lee

# Calcular el Límite Máximo del Eje X
max_edge_length <- max(arbol$edge.length, na.rm = TRUE)

# Graficar el árbol
ggtree(arbol) +
  geom_tiplab(size = 3) +
  geom_text2(aes(subset = !isTip, label = node), hjust = -0.3) +  # Etiquetar nodos internos
  xlim(-0.5, max_edge_length * 1.4) +
  theme_tree2()



#Separar el ingroup

nodo_interes <- 693 #562 #693  # Reemplaza este valor por el nodo de tu interés

# Obtener los índices de los nodos descendientes
nodos_internos <- getDescendants(arbol, nodo_interes)

# Filtrar solo los índices que corresponden a las terminales (tips)
indices_tips <- nodos_internos[nodos_internos <= length(arbol$tip.label)]

# Obtener los nombres de las terminales
nombres_tips <- arbol$tip.label[indices_tips]

# Mostrar los nombres de las terminales
print(nombres_tips)

# Obtener el nodo más reciente en común del ingroup
mrca_ingroup <- getMRCA(arbol, nombres_tips)

# Extraer el subárbol del ingroup
subarbol <- extract.clade(arbol, mrca_ingroup)

# Calcular la longitud máxima de las ramas del subárbol
max_edge_length <- max(subarbol$edge.length, na.rm = TRUE)

# Graficar el subárbol
ggtree(subarbol) +
  geom_tiplab(size = 3) +
  xlim(-0.5, max_edge_length * 2) +
  theme_tree2()


#quitar especies duplicadas
# Separar los nombres de las especies (asumiendo el formato "ID_Especie")
species_info <- data.frame(
  tip_label = nombres_tips,
  species = str_extract(nombres_tips, "[A-Za-z]+_[A-Za-z]+(?:_[A-Za-z]+)?$")
)

# Eliminar las ssp _h y _p
# Seleccionar solo una muestra por especie
unique_species <- species_info %>%
  filter(!str_detect(species, "_[hp]$")) %>%
  group_by(species) %>% 
  slice(1) %>%  # Selecciona solo la primera aparición de cada especie
  ungroup()

# Extraer los nombres de los tips que queremos conservar
selected_tips <- unique_species$tip_label

# Remover duplicados
subarbol_final <- drop.tip(subarbol, setdiff(subarbol$tip.label, selected_tips))

ggtree(subarbol) + 
  geom_tiplab(size = 3) + 
  xlim(-0.5, max(subarbol$edge.length) * 2)  # Expande el espacio a la izquierda

#guardar el abrol para revbayes

write.nexus(subarbol, file="/home/daniela/GSM/Revbayes/Rev_in/98_2701_ingroup.nex")

# Guardar en formato Newick (.tre)
write.tree(subarbol, file = "/home/daniela/GSM/Revbayes/Rev_in/98_2701_ingroup.tre")

#Calibrar_97_beast_2101_strict.tree
