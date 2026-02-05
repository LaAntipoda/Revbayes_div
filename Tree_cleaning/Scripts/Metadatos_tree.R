
library(sf)
library(dplyr)
library(readr)
library(lwgeom)
library(tidyr)

####### Plots de GSM
#mis datos coords
tabla <- read_table("/home/daniela/GSM/originales/OG/Fungi_GSMc_sample_metadata.txt")

gsm_df <- tabla %>%
  select(plot,latitude, longitude)

###### Plots de EUK 

df_euk_og <- read.csv("/home/daniela/GSM/Alineamientos/Alineamientos/nuevos/clusters/metadata_clusters/Euk_metadatos.csv") %>%
  rename(plot = sample_ID)# renombrar las columna Sample por label

euk_df <- df_euk_og %>%
  select(plot,latitude, longitude)

#merge

# Convertir latitude a character en ambos dataframes
gsm_df$latitude <- as.numeric(gsm_df$latitude)
gsm_df$longitude <- as.numeric(gsm_df$longitude)
euk_df$latitude <- as.numeric(euk_df$latitude)
euk_df$longitude <- as.numeric(euk_df$longitude)

# Ahora sí combinar
df_combined <- bind_rows(gsm_df, euk_df)

#quitar muestras raras
df_combined <- df_combined %>%
  filter(!is.na(latitude))

#transformar sf
puntos_sf <- st_as_sf(df_combined, coords = c("longitude", "latitude"), crs = 4326)

#el sapefile
shapefile <- st_read("Ecoregion/official/wwf_terr_ecos.shp")
shapefile <- st_make_valid(shapefile)
st_crs(puntos_sf) <- st_crs(shapefile)

resultados <- st_join(puntos_sf, shapefile)
write_csv(st_drop_geometry(resultados), "/home/daniela/GSM/Alineamientos/Alineamientos/nuevos/clusters/metadata_clusters/metadatos_tree.csv")  # o solo tabla



##################################
#otu table
otu_gsm = read.delim("Stephanosporaceae/Stephanosporaceae_table.txt", check.names = FALSE)

# Crear tabla de presencia/ausencia (1/0) para euk
otu_euk  <- read.csv("/home/daniela/GSM/Alineamientos/Alineamientos/nuevos/clusters/metadata_clusters/Euk_metadatos.csv") %>%
  rename(plot = sample_ID) %>%
  rename(OTU = accession)

otu_euk <- otu_euk %>%
  select(OTU, plot) %>%           # Solo las columnas que necesitas
  distinct() %>%                   # Eliminar duplicados
  mutate(presencia = 1) %>%        # Crear columna de presencia
  pivot_wider(
    names_from = plot,
    values_from = presencia,
    values_fill = 0
  )

#otu table para todos los label 
otu_combined <- bind_rows(otu_euk, otu_gsm)%>%
  mutate(across(everything(), ~replace_na(., 0)))

# Metadatos por plot
plots <- read.csv("/home/daniela/GSM/Alineamientos/Alineamientos/nuevos/clusters/metadata_clusters/metadatos_tree.csv")

#asociar OTU-metadatos

library(tidyr)

# Primero convertir de ancho a largo
otu_plot <- otu_combined %>%
  pivot_longer(
    cols = -OTU,              # Todas las columnas excepto OTU
    names_to = "plot",
    values_to = "presencia"
  ) %>%
  filter(presencia == 1) %>%  # Solo donde está presente
  select(-presencia)

# Luego unir con metadatos
label_metadatos <- otu_plot %>%
  left_join(plots, by = "plot")

write_csv(label_metadatos, "/home/daniela/GSM/Alineamientos/Alineamientos/nuevos/clusters/metadata_clusters/label_metadatos_tree.csv")



#############

############################ cluster a partir de aqui 

############


# quitar labels repetidos 

label_metadatos <- read_csv("/home/daniela/GSM/Alineamientos/Alineamientos/nuevos/clusters/metadata_clusters/label_metadatos_tree.csv")


#crear categorias de distribucion 


label_metadatos_collapsed <- label_metadatos %>%
  filter(!is.na(REALM) & REALM != "") %>%
  group_by(OTU) %>%
  summarise(
    n_sitios = n(),
    realms = paste(unique(REALM), collapse = ", "),
    n_realms = n_distinct(REALM),
    # Crear categoría de distribución biogeográfica
    categoria_realm = case_when(
      n_distinct(REALM) == 1 ~ paste0("Endémico_", first(REALM)),
      n_distinct(REALM) == 2 ~ "Dos_realms",
      n_distinct(REALM) >= 3 ~ "Cosmopolita",
      TRUE ~ "Desconocido"
    ),
    # O categorías más específicas
    tipo_distribucion = case_when(
      all(REALM == "NT") ~ "Neotropical",
      all(REALM == "PA") ~ "Paleartico",
      all(REALM == "IM") ~ "Indiomalayo",
      all(REALM == "AA") ~ "Australasia",
      all(REALM == "OC") ~ "Oceania",
      all(REALM == "AN") ~ "Antartica",
      all(REALM == "AT") ~ "Afrotropico",
      all(REALM == "NAA") ~ "Neartico",
      all(REALM %in% c("NT", "NAA")) ~ "América",
      all(REALM %in% c("NT", "AT")) ~ "Trópico",
      all(REALM %in% c("PA", "IM")) ~ "Asia",
      all(REALM %in% c("PA", "AT")) ~ "Euroafricano",
      all(REALM %in% c("PA", "NAA")) ~ "Norteglobal",
      n_distinct(REALM) >= 3 ~ "Amplia_distribución",
      TRUE ~ "Multi-realm"
    ),
    .groups = 'drop'
  )

write.csv(label_metadatos_collapsed, "/home/daniela/GSM/Alineamientos/Alineamientos/nuevos/clusters/metadata_clusters/realm_metadatos_tree.csv")


