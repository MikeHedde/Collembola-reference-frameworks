#Ce script correspond à la mise en forme des différents jeux de données que j' ai utilisé pour mon stage M2 en 2025. Les jeux de données orginaux étant peu ou mal structurés, j'ai du faire quelques modifications sous Excel afin de gagner du temps. Ce script ne reprend que les étapes principales de la mise en forme, à savoir le calcul des densités de collemboles au m² (moyenne des échantillan x 0.002827) pour chaque site à partir des effectifs bruts de chaque carottier, l'extraction, la conversion et l'inclusion des coordonnées des sites dans les datasets, la vérification taxref, et la combinaison en un unique jeu de donnée pour les gouverner tous. 

  
  
  

##################Calcul des densités par site au m²##########

setwd("C:/UP/M2/Stage/Stats/Datasets")

#packages
library(readxl)
library(dplyr)
library(openxlsx)
library(stringr)

###################RMQS 2024##############


#Importation du dataset
RMQS_collembola <- read_excel("RMQS_2024_COLLEMBOLA.xlsx",sheet="Feuille stats")
RMQS_collembola <- RMQS_collembola %>% mutate(Site = sub("_lame[0-9]+$", "", Etiquettes))

# Vérifier si l'extraction a bien fonctionné
unique(RMQS_collembola$Site)

# calculer les moyennes au m²
RMQS_collembola <- RMQS_collembola %>%
  group_by(Site) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  mutate(across(-Site, ~ . / 0.002827))


# Charger le workbook
if (file.exists("RMQS_2024_COLLEMBOLA.xlsx")) {
  wb <- loadWorkbook("RMQS_2024_COLLEMBOLA.xlsx")
} else {
  wb <- createWorkbook()
}

# Ajouter une nouvelle feuille
addWorksheet(wb, "Moyennes collemboles par m²")

# Écrire les données
writeData(wb, sheet = "Moyennes collemboles par m²", RMQS_collembola)

# Sauvegarder
saveWorkbook(wb, "RMQS_2024_COLLEMBOLA.xlsx", overwrite = TRUE)


##################ANDRA#############################

ANDRA <- read_excel ("ANDRA_ARTHROPODA.xlsx", sheet="temporary")

# calculer les moyennes au m²
ANDRA <- ANDRA %>%
  group_by(Site, Annee) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(!c(Site, Annee), ~ . / 0.002827))

# Charger le workbook
if (file.exists("ANDRA_ARTHROPODA.xlsx")) {
  wb <- loadWorkbook("ANDRA_ARTHROPODA.xlsx")
} else {
  wb <- createWorkbook()
}

# Ajouter une nouvelle feuille
addWorksheet(wb, "Feuille stats")

# Écrire les données
writeData(wb, sheet = "Feuille stats", ANDRA)

# Sauvegarder
saveWorkbook(wb, "ANDRA_ARTHROPODA.xlsx", overwrite = TRUE)


###########TIGA rural#####################
TIGA_rural <- read_excel ("TIGA_rural_MESOFAUNA.xlsx", sheet="temporary")

# calculer les moyennes au m²
TIGA_rural <- TIGA_rural %>%
  mutate(Site = str_remove(Site, "_[1-3]$")) %>%  # Supprime _1, _2 ou _3 à la fin du nom
  group_by(Site) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(!Site, ~ . / 0.002827))

# Charger le workbook
if (file.exists("TIGA_rural_MESOFAUNA.xlsx")) {
  wb <- loadWorkbook("TIGA_rural_MESOFAUNA.xlsx")
} else {
  wb <- createWorkbook()
}

# Ajouter une nouvelle feuille
addWorksheet(wb, "Feuille stats")

# Écrire les données
writeData(wb, sheet = "Feuille stats", TIGA_rural)

# Sauvegarder
saveWorkbook(wb, "TIGA_rural_MESOFAUNA.xlsx", overwrite = TRUE)


###########BISES#########"

BISES <- read_excel ("BISES_MESOFAUNA.xlsx", sheet="temporary")

# calculer les moyennes au m²
BISES <- BISES %>%
  group_by(Site) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  mutate(across(-Site, ~ . / 0.002827))

# Charger le workbook
if (file.exists("BISES_MESOFAUNA.xlsx")) {
  wb <- loadWorkbook("BISES_MESOFAUNA.xlsx")
} else {
  wb <- createWorkbook()
}

# Ajouter une nouvelle feuille
addWorksheet(wb, "Feuille stats")

# Écrire les données
writeData(wb, sheet = "Feuille stats", BISES)

# Sauvegarder
saveWorkbook(wb, "BISES_MESOFAUNA.xlsx", overwrite = TRUE)


############RMQS_Bretagne###########
RMQS_Bretagne <- read_excel ("RMQS_Bretagne.xls", sheet="temporary") %>%  mutate(across(where(is.numeric), ~ replace_na(., 0)))

# calculer les moyennes au m²
RMQS_Bretagne <- RMQS_Bretagne %>%
  group_by(Site) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  mutate(across(-Site, ~ . / 0.002827))

# Charger le workbook
if (file.exists("RMQS_Bretagne.xlsx")) {
  wb <- loadWorkbook("RMQS_Bretagne.xlsx")
} else {
  wb <- createWorkbook()
}

# Ajouter une nouvelle feuille
addWorksheet(wb, "Feuille stats")

# Écrire les données
writeData(wb, sheet = "Feuille stats", RMQS_Bretagne)

# Sauvegarder
saveWorkbook(wb, "RMQS_Bretagne.xlsx", overwrite = TRUE)


#####################Bioindicateur2###########################
Bioindicateur2 <- read_excel ("Bioindicateur2_ARTHROPODA.xls", sheet="temporary")

#enlever les numéros échantillons qui perturbent le regroupement
Bioindicateur2$Site <- substr(Bioindicateur2$Site, 1, nchar(Bioindicateur2$Site) - 1)


# calculer les moyennes au m²
Bioindicateur2 <- Bioindicateur2 %>%
  group_by(Site, Annee, Nom_site, Usage, Traitement, Contamination) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(!c(Site, Annee, Nom_site, Usage, Traitement, Contamination), ~ . / 0.002827))

# Charger le workbook
if (file.exists("Bioindicateur2_ARTHROPODA.xlsx")) {
  wb <- loadWorkbook("Bioindicateur2_ARTHROPODA.xlsx")
} else {
  wb <- createWorkbook()
}

# Ajouter une nouvelle feuille
addWorksheet(wb, "Feuille stats")

# Écrire les données
writeData(wb, sheet = "Feuille stats", Bioindicateur2)

# Sauvegarder
saveWorkbook(wb, "Bioindicateur2_ARTHROPODA.xlsx", overwrite = TRUE)









#############ajouter les coordonnées et les convertir en WGS84 et Lambert 93####################

setwd("C:/UP/M2/Stage/Stats/Datasets")

#packages
library(readxl)
library(openxlsx)
library(sf)
library(dplyr)

#################RMQS 2024######################

RMQS_2024 <- read_excel("RMQS_2024_COLLEMBOLA.xlsx",sheet="Feuille stats")
RMQS_stations <- read_excel("RMQS_2024_COLLEMBOLA.xlsx",sheet="Stations RMQS")


# Sélectionner uniquement les colonnes utiles de RMQS_stations
coords <- RMQS_stations %>%
  select(NOM_STATION, `X (Latitude)`, `Y (Longitude)`) %>%
  rename(Site = NOM_STATION, Y_WGS84 = `X (Latitude)`, X_WGS84 = `Y (Longitude)`)

# Fusionner avec RMQS_collembola en gardant l'ordre correct des colonnes
RMQS_2024 <- RMQS_2024 %>%
  left_join(coords, by = "Site") %>%
  relocate(X_WGS84, Y_WGS84, .after = Annee)


# Convertir en objet sf avec les coordonnées actuelles (WGS84 EPSG:4326)
RMQS_2024_sf <- st_as_sf(RMQS_2024, coords = c("X_WGS84", "Y_WGS84"), crs = 4326)

# Transformer en Lambert 93 (EPSG:2154)
RMQS_2024_sf <- st_transform(RMQS_2024_sf, crs = 2154)

# Extraire les nouvelles coordonnées
coords_L93 <- st_coordinates(RMQS_2024_sf)
RMQS_2024 <- RMQS_2024 %>%
  mutate(X_L93 = coords_L93[,1], 
         Y_L93 = coords_L93[,2]) %>%
  relocate(X_L93, Y_L93, .after = Y_WGS84)


# Charger le workbook
if (file.exists("RMQS_2024_COLLEMBOLA.xlsx")) {
  wb <- loadWorkbook("RMQS_2024_COLLEMBOLA.xlsx")
} else {
  wb <- createWorkbook()
}

# Écrire les données
writeData(wb, sheet = "Feuille stats", RMQS_2024)

# Sauvegarder
saveWorkbook(wb, "RMQS_2024_COLLEMBOLA.xlsx", overwrite = TRUE)


###################RMQS 2021####################

RMQS_2021 <- read_excel("RMQS_2021_COLLEMBOLA.xlsx",sheet="Feuille stats")
RMQS_stations <- read_excel("RMQS_2021_COLLEMBOLA.xlsx",sheet="RMQS1_analyses_composites_18_11")


# Sélectionner uniquement les colonnes utiles de RMQS_stations
coords <- RMQS_stations %>%
  select(id_site, `x_theo`, `y_theo`) %>%
  rename(Site = id_site, X_L93 = `x_theo`, Y_L93 = `y_theo`) %>%
  distinct(Site, .keep_all = TRUE)  # Supprime les doublons


# Fusionner avec RMQS_collembola en gardant l'ordre correct des colonnes
RMQS_2021 <- RMQS_2021 %>%
  left_join(coords, by = "Site") %>%
  relocate(X_L93, Y_L93, .after = Annee)


# Convertir en objet sf avec les coordonnées actuelles (WGS84 EPSG:4326)
RMQS_2021_sf <- st_as_sf(RMQS_2021, coords = c("X_L93", "Y_L93"), crs = 2154)

# Transformer 
RMQS_2021_sf <- st_transform(RMQS_2021_sf, crs = 4326)

# Extraire les nouvelles coordonnées
coords_WGS84 <- st_coordinates(RMQS_2021_sf)
RMQS_2021 <- RMQS_2021 %>%
  mutate(X_WGS84 = coords_WGS84[,1], 
         Y_WGS84 = coords_WGS84[,2]) %>%
  relocate(X_WGS84, Y_WGS84, .after = Annee)


# Charger le workbook
if (file.exists("RMQS_2021_COLLEMBOLA.xlsx")) {
  wb <- loadWorkbook("RMQS_2021_COLLEMBOLA.xlsx")
} else {
  wb <- createWorkbook()
}

# Écrire les données
writeData(wb, sheet = "Feuille stats", RMQS_2021)

# Sauvegarder
saveWorkbook(wb, "RMQS_2021_COLLEMBOLA.xlsx", overwrite = TRUE)


#####################ANDRA##################

ANDRA <- read_excel("ANDRA_ARTHROPODA.xlsx",sheet="Feuille stats")
ANDRA_stations <- read_excel("ANDRA_ARTHROPODA.xlsx",sheet="Coord")

# Sélectionner uniquement les colonnes utiles de RMQS_stations
coords <- ANDRA_stations %>%
  select(Reference_Site, `X_L1`, `Y_L1`) %>%
  rename(Site = Reference_Site) %>%
  distinct(Site, .keep_all = TRUE)  # Supprime les doublons
ANDRA <- ANDRA %>%
  left_join(coords, by = "Site") %>%
  relocate(X_L1, Y_L1, .after = Annee)



# Convertir en objet sf avec les coordonnées actuelles (WGS84 EPSG:4326)
ANDRA_sf <- st_as_sf(ANDRA, coords = c("X_L1", "Y_L1"), crs = 27571)

# Transformer 
ANDRA_sf <- st_transform(ANDRA_sf, crs = 2154)

# Extraire les nouvelles coordonnées
coords_L93 <- st_coordinates(ANDRA_sf)
ANDRA <- ANDRA %>%
  mutate(X_L93 = coords_L93[,1], 
         Y_L93 = coords_L93[,2]) %>%
  relocate(X_L93, Y_L93, .after = Annee)

# Transformer 
ANDRA_sf <- st_transform(ANDRA_sf, crs = 4326)

# Extraire les nouvelles coordonnées
coords_WGS84 <- st_coordinates(ANDRA_sf)
ANDRA <- ANDRA %>%
  mutate(X_WGS84 = coords_WGS84[,1], 
         Y_WGS84 = coords_WGS84[,2]) %>%
  relocate(X_WGS84, Y_WGS84, .after = Annee)

#supprimer les colonnes L1
ANDRA <- ANDRA %>% select(-X_L1, -Y_L1)


# Charger le workbook
if (file.exists("ANDRA_ARTHROPODA.xlsx")) {
  wb <- loadWorkbook("ANDRA_ARTHROPODA.xlsx")
} else {
  wb <- createWorkbook()
}

# Écrire les données
writeData(wb, sheet = "Feuille stats", ANDRA)

# Sauvegarder
saveWorkbook(wb, "ANDRA_ARTHROPODA.xlsx", overwrite = TRUE)


############TIGA_rural####################

TIGA_rural <- read_excel("TIGA_rural_MESOFAUNA.xlsx",sheet="Feuille stats")
TIGA_stations <- read_excel("TIGA_rural_MESOFAUNA.xlsx",sheet="Info sites")

library(lubridate)

# Convertir la date en format Date et extraire l'année
TIGA_stations <- TIGA_stations %>%
  mutate(Annee_extraite = year(dmy(date_du_prelevement)))  # Si format jj/mm/aaaa

# Fusionner les jeux de données pour ajouter l'année
TIGA_rural <- TIGA_rural %>%
  left_join(TIGA_stations %>% select(code_ech, Annee_extraite), by = c("Site" = "code_ech")) %>%
  mutate(Annee = coalesce(Annee, Annee_extraite)) %>%  # Remplacer les valeurs manquantes
  select(-Annee_extraite)  # Supprimer la colonne temporaire


# Sélectionner uniquement les colonnes utiles de RMQS_stations
coords <- TIGA_stations %>%
  select(code_ech, `X_L93`, `Y_L93`) %>%
  rename(Site = code_ech) %>%
  distinct(Site, .keep_all = TRUE)  # Supprime les doublons


# Fusionner avec RMQS_collembola en gardant l'ordre correct des colonnes
TIGA_rural <- TIGA_rural %>%
  left_join(coords, by = "Site") %>%
  relocate(X_L93, Y_L93, .after = Annee)


# Convertir en objet sf avec les coordonnées actuelles (WGS84 EPSG:4326)
TIGA_rural_sf <- st_as_sf(TIGA_rural, coords = c("X_L93", "Y_L93"), crs = 2154)

# Transformer 
TIGA_rural_sf <- st_transform(TIGA_rural_sf, crs = 4326)

# Extraire les nouvelles coordonnées
coords_WGS84 <- st_coordinates(TIGA_rural_sf)
TIGA_rural <- TIGA_rural %>%
  mutate(X_WGS84 = coords_WGS84[,1], 
         Y_WGS84 = coords_WGS84[,2]) %>%
  relocate(X_WGS84, Y_WGS84, .after = Annee)


# Charger le workbook
if (file.exists("TIGA_rural_MESOFAUNA.xlsx")) {
  wb <- loadWorkbook("TIGA_rural_MESOFAUNA.xlsx")
} else {
  wb <- createWorkbook()
}

# Écrire les données
writeData(wb, sheet = "Feuille stats", TIGA_rural)

# Sauvegarder
saveWorkbook(wb, "TIGA_rural_MESOFAUNA.xlsx", overwrite = TRUE)




##################RMQS_Bretagne######################
RMQS_Bretagne <- read_excel("RMQS_Bretagne.xlsx",sheet="Feuille stats")
RMQS_stations <- read_excel("RMQS_2021_COLLEMBOLA.xlsx",sheet="RMQS1_analyses_composites_18_11")


# Sélectionner uniquement les colonnes utiles de RMQS_stations
coords <- RMQS_stations %>%
  select(id_site, `x_theo`, `y_theo`) %>%
  rename(Site = id_site, X_L93 = `x_theo`, Y_L93 = `y_theo`) %>%
  distinct(Site, .keep_all = TRUE)  # Supprime les doublons


# Fusionner avec RMQS_collembola en gardant l'ordre correct des colonnes
RMQS_Bretagne <- RMQS_Bretagne %>%
  left_join(coords, by = "Site") %>%
  relocate(X_L93, Y_L93, .after = Annee)


# Convertir en objet sf avec les coordonnées actuelles (WGS84 EPSG:4326)
RMQS_Bretagne_sf <- st_as_sf(RMQS_Bretagne, coords = c("X_L93", "Y_L93"), crs = 2154)

# Transformer 
RMQS_Bretagne_sf <- st_transform(RMQS_Bretagne_sf, crs = 4326)

# Extraire les nouvelles coordonnées
coords_WGS84 <- st_coordinates(RMQS_Bretagne_sf)
RMQS_Bretagne <- RMQS_Bretagne %>%
  mutate(X_WGS84 = coords_WGS84[,1], 
         Y_WGS84 = coords_WGS84[,2]) %>%
  relocate(X_WGS84, Y_WGS84, .after = Annee)


# Charger le workbook
if (file.exists("RMQS_Bretagne.xlsx")) {
  wb <- loadWorkbook("RMQS_Bretagne.xlsx")
} else {
  wb <- createWorkbook()
}

# Écrire les données
writeData(wb, sheet = "Feuille stats", RMQS_Bretagne)

# Sauvegarder
saveWorkbook(wb, "RMQS_Bretagne.xlsx", overwrite = TRUE)


##################Bioindicateur2######################
Bioindicateur2 <- read_excel("Bioindicateur2_ARTHROPODA.xls",sheet="Feuille stats")
Bioindicateur2_stations <- read_excel("Bioindicateur2_ARTHROPODA.xls",sheet="Coord")


# Sélectionner uniquement les colonnes utiles de RMQS_stations
coords <- Bioindicateur2_stations %>%
  select(NOM_SITE, `X`, `Y`) %>%
  rename(Nom_site = NOM_SITE, X_WGS84 = `Y`, Y_WGS84 = `X`) %>%
  distinct(Nom_site, .keep_all = TRUE)  # Supprime les doublons


# Fusionner en gardant l'ordre correct des colonnes
Bioindicateur2 <- Bioindicateur2 %>%
  left_join(coords, by = "Nom_site") %>%
  relocate(X_WGS84, Y_WGS84, .after = Traitement)


# Convertir en objet sf avec les coordonnées actuelles (WGS84 EPSG:4326)
Bioindicateur2_sf <- st_as_sf(Bioindicateur2, coords = c("X_WGS84", "Y_WGS84"), crs = 4326)

# Transformer 
Bioindicateur2_sf <- st_transform(Bioindicateur2_sf, crs = 2154)

# Extraire les nouvelles coordonnées
coords_L93 <- st_coordinates(Bioindicateur2_sf)
Bioindicateur2 <- Bioindicateur2 %>%
 mutate(X_L93 = coords_L93[,1], 
         Y_L93 = coords_L93[,2]) %>%
  relocate(X_L93, Y_L93, .after = Y_WGS84)


# Charger le workbook
if (file.exists("Bioindicateur2_ARTHROPODA.xlsx")) {
  wb <- loadWorkbook("Bioindicateur2_ARTHROPODA.xlsx")
} else {
  wb <- createWorkbook()
}

# Ajouter une nouvelle feuille
addWorksheet(wb, "Feuille stats")

# Écrire les données
writeData(wb, sheet = "Feuille stats", Bioindicateur2)

# Sauvegarder
saveWorkbook(wb, "Bioindicateur2_ARTHROPODA.xlsx", overwrite = TRUE)









######################## Vérification Taxref###########

setwd("C:/UP/M2/Stage/Stats/Datasets")

library(readxl)
library(openxlsx)

TaxRef<-read.csv("TaxRef18_Collembola.csv", sep=";")

#créer une colonne GENRE
TaxRef$LB_NOM <- iconv(TaxRef$LB_NOM, from = "", to = "UTF-8", sub = "byte")
TaxRef$GENRE <- sub(" .*", "", TaxRef$LB_NOM)

#créer une fonction 
check_in_taxref <- function(dataset, taxref, start_col = 13) {
  # Convertir les noms en UTF-8 pour éviter les problèmes d'encodage
  taxref$LB_NOM <- iconv(taxref$LB_NOM, from = "", to = "UTF-8", sub = "byte")
  
  # Extraire les noms d'espèces (colonnes après la 12e)
  species_names <- colnames(dataset)[start_col:ncol(dataset)]
  
  # Vérifier la correspondance directe avec LB_NOM de TaxRef
  direct_match <- species_names %in% taxref$LB_NOM
  
  # Identifier les noms au format "Genre sp."
  genre_sp <- grepl("^[A-Z][a-z]+ sp\\.$", species_names)
  
  # Extraire les genres de ces noms
  genres <- sub(" sp\\.$", "", species_names[genre_sp])
  
  # Vérifier si ces genres existent dans TaxRef$GENRE
  valid_genres <- genres %in% taxref$GENRE
  
  # Créer une validation globale
  valid_species <- direct_match | (genre_sp & (species_names %in% paste0(taxref$GENRE, " sp.")))
  
  # Filtrer les noms qui n'ont pas de correspondance
  species_not_in_taxref <- species_names[!valid_species]
  
  # Retourner les espèces non valides
  return(species_not_in_taxref)
}


#importation des datasets
RMQS_2024 <- read_excel("RMQS_2024_COLLEMBOLA.xlsx",sheet="Feuille stats")
RMQS_2021 <- read_excel("RMQS_2021_COLLEMBOLA.xlsx",sheet="Feuille stats")
RMQS_Bretagne <- read_excel("RMQS_Bretagne.xls",sheet="Feuille stats")
ANDRA <- read_excel("ANDRA_ARTHROPODA.xlsx",sheet="Feuille stats")
TIGA_rural <- read_excel("TIGA_rural_MESOFAUNA.xlsx",sheet="Feuille stats")
Bioindicateur2 <- read_excel("Bioindicateur2_ARTHROPODA.xls",sheet="Feuille stats")

#vérif TaxRef
check_in_taxref(RMQS_2024, taxref = TaxRef, start_col = 13) 
check_in_taxref(RMQS_2021, taxref = TaxRef, start_col = 13)
check_in_taxref(RMQS_Bretagne, taxref = TaxRef, start_col = 13)
check_in_taxref(ANDRA, taxref = TaxRef, start_col = 13)
check_in_taxref(TIGA_rural, taxref = TaxRef, start_col = 13)
check_in_taxref(Bioindicateur2, taxref = TaxRef, start_col = 11)

#on exclura de toute façon les niveaux taxonomiques supérieurs au genre



####################compiler tous les jeux de données################

setwd("C:/UP/M2/Stage/Stats/Datasets")

#packages
library(readxl)
library(dplyr)
library(openxlsx)
library(purrr)
library(stringr)

library(tidyverse)
library(readxl)

# Importation des jeux de données
RMQS_2024 <- read_excel("RMQS_2024_COLLEMBOLA.xlsx", sheet = "Feuille stats")
RMQS_2021 <- read_excel("RMQS_2021_COLLEMBOLA.xlsx", sheet = "Feuille stats")
RMQS_Bretagne <- read_excel("RMQS_Bretagne.xls", sheet = "Feuille stats")
ANDRA <- read_excel("ANDRA_ARTHROPODA.xlsx", sheet = "Feuille stats")
TIGA_rural <- read_excel("TIGA_rural_MESOFAUNA.xlsx", sheet = "Feuille stats")
#Bioindicateur2 <- read_excel("Bioindicateur2_ARTHROPODA.xls", sheet = "Feuille stats")

# Liste des datasets
datasets <- list(RMQS_2024, RMQS_2021, RMQS_Bretagne, ANDRA, TIGA_rural)

# Identifier les colonnes fixes communes à tous les datasets
colonnes_fixes <- c("Site", "Projet", "Strategie", "Echelle", "Annee", 
                    "X_WGS84", "Y_WGS84", "X_L93", "Y_L93", 
                    "CLC_niveau1", "CLC_niveau2", "CLC_niveau3")

# Harmoniser les colonnes fixes (éviter les problèmes de type)
datasets <- map(datasets, function(df) {
  df %>% mutate(across(all_of(colonnes_fixes), as.character)) # Tout convertir en texte
})

# Transformer chaque dataset en format long pour gérer les espèces correctement
datasets_long <- map(datasets, function(df) {
  df %>%
    pivot_longer(cols = -all_of(colonnes_fixes), names_to = "Espece", values_to = "Abondance")
})

# Fusionner tous les datasets en un seul grand tableau long
data_long <- bind_rows(datasets_long)

# Repasser en format large avec une colonne par espèce
Dataset_ultime <- data_long %>%
  pivot_wider(names_from = "Espece", values_from = "Abondance", values_fill = NA)

# Vérifier les dimensions et un aperçu du dataset final
print(dim(Dataset_ultime))
View(Dataset_ultime)


# Charger le workbook
if (file.exists("Stage_M2_Helio_Suarez.xlsx")) {
  wb <- loadWorkbook("Stage_M2_Helio_Suarez.xlsx")
} else {
  wb <- createWorkbook()
}

# Ajouter une nouvelle feuille
addWorksheet(wb, "Data")


# Écrire les données
writeData(wb, sheet = "Data", Dataset_ultime)

# Sauvegarder
saveWorkbook(wb, "Stage_M2_Helio_Suarez.xlsx", overwrite = TRUE)
