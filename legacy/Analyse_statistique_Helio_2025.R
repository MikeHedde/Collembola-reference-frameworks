setwd("C:/UP/M2/Stage/Stats")

#packages
library(readxl)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(e1071) #skewness et kurtosis
library(car)
library(cowplot)
library(cramer) #test de Cramer-von Mises
library(PMCMRplus) #test post-hoc de Nemenyi

#dataset
data <- read_excel("Stage_M2_Helio_Suarez.xlsx")
data[, 13:ncol(data)][is.na(data[, 13:ncol(data)])] <- 0 #remplacer les NA par des 0
data[,13:ncol(data)]<-round(data[,13:ncol(data)]) #arrondir les effectifs à l'unité

############1)Calcul des indices######

communities<-data[,13:205] #ici je crée un subset qui ne contient que les colonnes des espèces d'arthropodes
# Sélection des colonnes dont le nom contient un espace
communities_strict <- communities[, grepl(" ", names(communities))] #communities_strict ne contient que les taxons dont la résolution est au genre ou à l'espèce
length(which(colSums(communities)<119))#28 singleton
length(which(rowSums(communities)<1))#5 sites vides

data$RSb <- specnumber(communities_strict) #Richesse specifique brute
#rarecurve(communities, step = 10, col = "blue", label = FALSE, xlim = c(0, 3000), ylim = c(0, 30)) # Tracer la courbe avec de raréfaction
data$RSr <- rarefy(communities_strict, sample = 500) #Richesse specifique rarefie
data$ab <- rowSums(communities) #Abondance
data$Shannon <- diversity(communities_strict, "shannon") #Shannon (permet de calculer Pielou)
data$Pielou <- with(data, ifelse(RSr <= 1, NA, Shannon / log(RSr)))

#charger statsvar et multiCvM


#################Analyse #############

data <- data %>% 
  mutate(Projet = factor(Projet, levels = rev(c("RMQS_Biodiversite", "RMQS_BioDiv_Bretagne", "ANDRA","TIGA_rural"))))


#############Richessse specifique##########

#violin plot
plotRSr<- ggplot(data, aes(x = RSr, y = Projet, fill = Projet)) +
  geom_violin(trim = TRUE, alpha = 0.7) +  # Violin plot horizontal
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Boxplot centré
  geom_jitter(width = 0, height = 0.2, alpha = 0.5) +  # Points individuels (jitter vertical)
  scale_fill_manual(values = c("RMQS_Biodiversite" = "#56B4E9", "RMQS_BioDiv_Bretagne" = "#009E73","ANDRA" = "#F0E442", "TIGA_rural" = "#D55E00"),
                    labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural"), guide = guide_legend(reverse = TRUE))+
  labs(x = "Richesse spécifique", y ="") + 
  scale_y_discrete(labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural")) +
  theme_minimal()

#extraction de la légende pour le plot grid et retrait de RSr
legend <- get_legend(plotRSr + theme(legend.position = "right"))
plotRSr <- plotRSr + theme(legend.position = "none")
plotRSr

# Résumé statistique
stats_varall(data$RSr, data$Projet)

#Cramer von Mises
multiCvM(data$RSr, data$Projet, seed = 123)



##################Equitabilite de Pielou################

#violin plot
plotPielou<-ggplot(data, aes(x = Pielou, y = Projet, fill = Projet)) +
  geom_violin(trim = TRUE, alpha = 0.7) +  # Violin plot horizontal
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Boxplot centré
  geom_jitter(width = 0, height = 0.2, alpha = 0.5) +  # Points individuels (jitter vertical)
  scale_fill_manual(values = c("RMQS_Biodiversite" = "#56B4E9", "RMQS_BioDiv_Bretagne" = "#009E73","ANDRA" = "#F0E442", "TIGA_rural" = "#D55E00"),
                    labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural"), guide = guide_legend(reverse = TRUE))+
  labs(x = "Indice de Pielou", y ="") + 
  scale_y_discrete(labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural")) +
  theme_minimal() + theme(legend.position = "none")
plotPielou

# Résumé statistique
stats_varall(data$Pielou, data$Projet)

#Cramer von Mises
multiCvM(data$Pielou, data$Projet, seed = 123)




##################Diversité de Shannon################

#violin plot
plotShannon<-ggplot(data, aes(x = Shannon, y = Projet, fill = Projet)) +
  geom_violin(trim = TRUE, alpha = 0.7) +  # Violin plot horizontal
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Boxplot centré
  geom_jitter(width = 0, height = 0.2, alpha = 0.5) +  # Points individuels (jitter vertical)
  scale_fill_manual(values = c("RMQS_Biodiversite" = "#56B4E9", "RMQS_BioDiv_Bretagne" = "#009E73","ANDRA" = "#F0E442", "TIGA_rural" = "#D55E00"),
                    labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural"), guide = guide_legend(reverse = TRUE))+
  labs(x = "Indice de Shannon", y ="") + 
  scale_y_discrete(labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural")) +
  theme_minimal() + theme(legend.position = "none")
plotShannon

# Résumé statistique
stats_varall(data$Shannon, data$Projet)

#Cramer von Mises
multiCvM(data$Shannon, data$Projet, seed = 123)




################Densité au m²##########

#violin plot
plotdensity<- ggplot(data, aes(x = ab, y = Projet, fill = Projet)) +
  geom_violin(trim = TRUE, alpha = 0.7) +  # Violin plot horizontal
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Boxplot centré
  geom_jitter(width = 0, height = 0.2, alpha = 0.5) +  # Points individuels (jitter vertical)
  scale_fill_manual(values = c("RMQS_Biodiversite" = "#56B4E9", "RMQS_BioDiv_Bretagne" = "#009E73","ANDRA" = "#F0E442", "TIGA_rural" = "#D55E00"),
                    labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural"), guide = guide_legend(reverse = TRUE))+
  labs(x = "Densité au m²", y ="") + 
  scale_y_discrete(labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural")) +
  theme_minimal() + theme(legend.position = "none")
plotdensity

# Résumé statistique
stats_varall(data$ab, data$Projet)

#Cramer von Mises
multiCvM(data$ab, data$Projet, seed = 123)

#########################plot grid################
grid<-plot_grid(plotRSr, plotShannon, plotPielou, plotdensity, labels=c("A","B","C","D"))
add_margin <- function(p) {
  p + theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5))  # marge droite = 15 pts
}
grid<-add_margin(grid)
grid

#ggsave("all_sites_Shanonfirst.pdf", plot = grid, width = 12, height = 10)
#ggsave("agri_sites_Shanonfirst.pdf", plot = grid, width = 12, height = 10)

#final_plot <- plot_grid(grid, legend, rel_widths = c(1, 0.2), ncol = 2)
#print(final_plot)


######################filtrer par habitat################
table(data$CLC_niveau1)
table(data$CLC_niveau2)
table(Bioindicateur2$CLC_niveau1)
table(Bioindicateur2$CLC_niveau2)
#au vu de la diminution de puissance statistique engendrée par l'augmentation du niveau de CLC, on restera au niveau 1
table(data$CLC_niveau1,data$Projet)
data <- data[!is.na(data$CLC_niveau1) & data$CLC_niveau1 == 2, ]

table(data$Projet)






#######################simulation des scores###################

#mise en place de la base de données à tester
Bioindicateur2 <- read_excel("Bioindicateur2_ARTHROPODA.xls", sheet="Feuille stats")
Bioindicateur2[,15:ncol(Bioindicateur2)]<-round(Bioindicateur2[,15:ncol(Bioindicateur2)])
Bioindicateur2 <- Bioindicateur2[!is.na(Bioindicateur2$CLC_niveau1) & Bioindicateur2$CLC_niveau1 == 2, ]

communities2<-Bioindicateur2[,15:61] 
communities2_strict <- communities2[, grepl(" ", names(communities2))] 
length(which(rowSums(communities2)<1))#0 sites vides

Bioindicateur2$RSr <- rarefy(communities2_strict, sample = 500) #Richesse specifique rarefie
Bioindicateur2$ab <- rowSums(communities2) #Abondance
Bioindicateur2$Shannon <- diversity(communities2_strict, "shannon") #Shannon (permet de calculer Pielou)
Bioindicateur2$Pielou <- with(Bioindicateur2, ifelse(RSr <= 1, NA, Shannon / log(RSr)))



#mise en place des référentiels
gamme_RMQS<-(data[data$Projet=="RMQS_Biodiversite",])
gamme_RMQS_Bretagne<-(data[data$Projet=="RMQS_BioDiv_Bretagne",])
gamme_ANDRA<-(data[data$Projet=="ANDRA",])
gamme_TIGA_rural<-(data[data$Projet=="TIGA_rural",])

#appeler fonction GiveScores

#calcul des scores

Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_RMQS, var = "RSr", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_RMQS_Bretagne, var = "RSr", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_ANDRA, var = "RSr", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_TIGA_rural, var = "RSr", nb_class = 7, inv=F)

Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_RMQS, var = "Pielou", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_RMQS_Bretagne, var = "Pielou", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_ANDRA, var = "Pielou", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_TIGA_rural, var = "Pielou", nb_class = 7, inv=F)

Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_RMQS, var = "Pielou", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_RMQS_Bretagne, var = "Pielou", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_ANDRA, var = "Pielou", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_TIGA_rural, var = "Pielou", nb_class = 7, inv=F)

Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_RMQS, var = "Shannon", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_RMQS_Bretagne, var = "Shannon", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_ANDRA, var = "Shannon", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_TIGA_rural, var = "Shannon", nb_class = 7, inv=F)

Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_RMQS, var = "ab", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_RMQS_Bretagne, var = "ab", nb_class = 7, inv=F)
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_ANDRA, var = "ab", nb_class = 7, inv=F)
gamme_TIGA_rural<- gamme_TIGA_rural[gamme_TIGA_rural$ab < 75000, ]#gamme sans outlier d'abondance
Bioindicateur2 <- GiveScores(Bioindicateur2, ref = gamme_TIGA_rural, var = "ab", nb_class = 7, inv=F)

data.frame(
  Colonne = names(scores <- colSums(Bioindicateur2[, grep("^score", names(Bioindicateur2))], na.rm = TRUE)),
  Somme = as.vector(scores)
)


####################Richesse spécifique#####################
rsr_scores <- Bioindicateur2 %>%
  select(site = 1,  # Remplace "1" par le nom ou numéro de la colonne identifiant les sites
         RMQS_Biodiversite = score_RSr_gamme_RMQS,
         RMQS_BioDiv_Bretagne = score_RSr_gamme_RMQS_Bretagne,
         ANDRA = score_RSr_gamme_ANDRA,
         TIGA_rural = score_RSr_gamme_TIGA_rural) %>%
  drop_na()

# Passer au format long
rsr_long <- rsr_scores %>%
  pivot_longer(cols = -site, names_to = "Gamme", values_to = "Score")

# S'assurer que les colonnes sont bien des facteurs
rsr_long$site <- as.factor(rsr_long$site)
rsr_long$Gamme <- as.factor(rsr_long$Gamme)

# Forcer l’ordre des niveaux de Gamme
rsr_long$Gamme <- factor(rsr_long$Gamme,
                         levels = c("RMQS_Biodiversite", "RMQS_BioDiv_Bretagne", "ANDRA", "TIGA_rural"))


# Test de Friedman
friedman.test(Score ~ Gamme | site, data = rsr_long)

# Post-hoc de Nemenyi
frdAllPairsNemenyiTest(Score ~ Gamme | site, data = rsr_long)

# 1. Créer la table de fréquences
rsr_bar <- rsr_long %>%
  mutate(Score = as.factor(Score)) %>%
  count(Gamme, Score)

# 2. Forcer l’ordre croissant des scores (0 en bas, 6 en haut)
rsr_bar$Score <- factor(rsr_bar$Score, levels = 6:0)

# 3. Palette daltoniens-friendly
score_colors <- c(
  "0" = "#8B4513",  # marron
  "1" = "#D73027",  # rouge
  "2" = "#FC8D59",  # orange clair
  "3" = "#FEE08B",  # jaune
  "4" = "#D9EF8B",  # vert clair
  "5" = "#91BFDB",  # bleu clair
  "6" = "#4575B4"   # bleu foncé
)

# 4. Barplot empilé avec contours
plot_score_RSr<-ggplot(rsr_bar, aes(x = Gamme, y = n, fill = Score)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # contours noirs
  scale_fill_manual(values = score_colors, name = "Score") +
  labs(x = "Barème utilisé", y = "Nombre de sites", title = "Richesse spécifique") +
  theme_minimal() +
  scale_x_discrete(labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural"))+
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
plot_score_RSr <- plot_score_RSr + theme(legend.position = "none")
plot_score_RSr



####################Pielou#####################
Pielou_scores <- Bioindicateur2 %>%
  select(site = 1,  # Remplace "1" par le nom ou numéro de la colonne identifiant les sites
         RMQS_Biodiversite = score_Pielou_gamme_RMQS,
         RMQS_BioDiv_Bretagne = score_Pielou_gamme_RMQS_Bretagne,
         ANDRA = score_Pielou_gamme_ANDRA,
         TIGA_rural = score_Pielou_gamme_TIGA_rural) %>%
  drop_na()

# Passer au format long
Pielou_long <- Pielou_scores %>%
  pivot_longer(cols = -site, names_to = "Gamme", values_to = "Score")

# S'assurer que les colonnes sont bien des facteurs
Pielou_long$site <- as.factor(Pielou_long$site)
Pielou_long$Gamme <- as.factor(Pielou_long$Gamme)

# Forcer l’ordre des niveaux de Gamme
Pielou_long$Gamme <- factor(Pielou_long$Gamme,
                         levels = c("RMQS_Biodiversite", "RMQS_BioDiv_Bretagne", "ANDRA", "TIGA_rural"))


# Test de Friedman
friedman.test(Score ~ Gamme | site, data = Pielou_long)

# Post-hoc de Nemenyi
frdAllPairsNemenyiTest(Score ~ Gamme | site, data = Pielou_long)


# 1. Créer la table de fréquences
Pielou_bar <- Pielou_long %>%
  mutate(Score = as.factor(Score)) %>%
  count(Gamme, Score)

# 2. Forcer l’ordre croissant des scores (0 en bas, 6 en haut)
Pielou_bar$Score <- factor(Pielou_bar$Score, levels = 6:0)

# 4. Barplot empilé avec contours
plot_score_Pielou<-ggplot(Pielou_bar, aes(x = Gamme, y = n, fill = Score)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # contours noirs
  scale_fill_manual(values = score_colors, name = "Score") +
  labs(x = "Barème utilisé", y = "Nombre de sites", title = "Indice de Pielou") +
  theme_minimal() +
  scale_x_discrete(labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural"))+
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
plot_score_Pielou <- plot_score_Pielou + theme(legend.position = "none")
plot_score_Pielou


####################Shannon#####################
Shannon_scores <- Bioindicateur2 %>%
  select(site = 1,  # Remplace "1" par le nom ou numéro de la colonne identifiant les sites
         RMQS_Biodiversite = score_Shannon_gamme_RMQS,
         RMQS_BioDiv_Bretagne = score_Shannon_gamme_RMQS_Bretagne,
         ANDRA = score_Shannon_gamme_ANDRA,
         TIGA_rural = score_Shannon_gamme_TIGA_rural) %>%
  drop_na()

# Passer au format long
Shannon_long <- Shannon_scores %>%
  pivot_longer(cols = -site, names_to = "Gamme", values_to = "Score")

# S'assurer que les colonnes sont bien des facteurs
Shannon_long$site <- as.factor(Shannon_long$site)
Shannon_long$Gamme <- as.factor(Shannon_long$Gamme)

# Forcer l’ordre des niveaux de Gamme
Shannon_long$Gamme <- factor(Shannon_long$Gamme,
                         levels = c("RMQS_Biodiversite", "RMQS_BioDiv_Bretagne", "ANDRA", "TIGA_rural"))

# Test de Friedman
friedman.test(Score ~ Gamme | site, data = Shannon_long)

# Post-hoc de Nemenyi
frdAllPairsNemenyiTest(Score ~ Gamme | site, data = Shannon_long)

# Test post-hoc de Conover
frdAllPairsConoverTest(
  y = Shannon_long$Score,
  groups = Shannon_long$Gamme,
  blocks = Shannon_long$site,
  p.adjust.method = "holm"
)
# 1. Créer la table de fréquences
Shannon_bar <- Shannon_long %>%
  mutate(Score = as.factor(Score)) %>%
  count(Gamme, Score)

# 2. Forcer l’ordre croissant des scores (0 en bas, 6 en haut)
Shannon_bar$Score <- factor(Shannon_bar$Score, levels = 6:0)

# 4. Barplot empilé avec contours
plot_score_Shannon<-ggplot(Shannon_bar, aes(x = Gamme, y = n, fill = Score)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # contours noirs
  scale_fill_manual(values = score_colors, name = "Score") +
  labs(x = "Barème utilisé", y = "Nombre de sites", title = "Indice de Shannon") +
  theme_minimal() +
  scale_x_discrete(labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural"))+
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
plot_score_Shannon <- plot_score_Shannon + theme(legend.position = "none")
plot_score_Shannon


####################ab#####################
ab_scores <- Bioindicateur2 %>%
  select(site = 1,  # Remplace "1" par le nom ou numéro de la colonne identifiant les sites
         RMQS_Biodiversite = score_ab_gamme_RMQS,
         RMQS_BioDiv_Bretagne = score_ab_gamme_RMQS_Bretagne,
         ANDRA = score_ab_gamme_ANDRA,
         TIGA_rural = score_ab_gamme_TIGA_rural) %>%
  drop_na()

# Passer au format long
ab_long <- ab_scores %>%
  pivot_longer(cols = -site, names_to = "Gamme", values_to = "Score")

# S'assurer que les colonnes sont bien des facteurs
ab_long$site <- as.factor(ab_long$site)
ab_long$Gamme <- as.factor(ab_long$Gamme)

# Forcer l’ordre des niveaux de Gamme
ab_long$Gamme <- factor(ab_long$Gamme,
                         levels = c("RMQS_Biodiversite", "RMQS_BioDiv_Bretagne", "ANDRA", "TIGA_rural"))


# Test de Friedman
friedman.test(Score ~ Gamme | site, data = ab_long)

# Post-hoc de Nemenyi
frdAllPairsNemenyiTest(Score ~ Gamme | site, data = ab_long)

# 1. Créer la table de fréquences
ab_bar <- ab_long %>%
  mutate(Score = as.factor(Score)) %>%
  count(Gamme, Score)

# 2. Forcer l’ordre croissant des scores (0 en bas, 6 en haut)
ab_bar$Score <- factor(ab_bar$Score, levels = 6:0)

# 4. Barplot empilé avec contours
plot_score_ab<-ggplot(ab_bar, aes(x = Gamme, y = n, fill = Score)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # contours noirs
  scale_fill_manual(values = score_colors, name = "Score") +
  labs(x = "Barème utilisé", y = "Nombre de sites", title = "Densité au m²") +
  theme_minimal() +
  scale_x_discrete(labels = c("RMQS_Biodiversite" = "RMQS Biodiversité", "RMQS_BioDiv_Bretagne" = "RMQS BioDiv Bretagne", "TIGA_rural" = "TIGA rural"))+
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
plot_score_ab <- plot_score_ab + theme(legend.position = "none")
plot_score_ab


#################extraction légende#############

# Ajout manuel de lignes factices avec les scores manquants
missing_scores <- setdiff(as.character(6:0), unique(as.character(rsr_bar$Score)))
dummy_rows <- expand.grid(
  Gamme = unique(rsr_bar$Gamme),
  Score = missing_scores
) %>%
  mutate(n = 0)

# Fusionner avec les vraies données
rsr_bar_complete <- bind_rows(rsr_bar, dummy_rows)
rsr_bar_complete$Score <- factor(rsr_bar_complete$Score, levels = 6:0)  # refaire au cas où

plot_score_RSr_legend <- ggplot(rsr_bar_complete, aes(x = Gamme, y = n, fill = Score)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = score_colors, name = "Score") +
  theme_minimal() +
  theme(legend.position = "right")

legend_score <- get_legend(plot_score_RSr_legend)



#################plot grid final#################
plots_combined <- plot_grid(plot_score_RSr, plot_score_Shannon,
                            plot_score_Pielou, plot_score_ab,
                            labels = c("A", "B", "C", "D"),
                            ncol = 2)

final_plot <- plot_grid(plots_combined, legend_score,
                        ncol = 2, rel_widths = c(1, 0.1))  # ajuste selon la taille voulue

final_plot

#ggsave("notes_baremes_tous_sites.pdf", plot = final_plot, width = 14, height = 10)



######################obtenir et représenter les baremes#############

# Récupération des objets "bareme_*"
baremes_list <- mget(ls(pattern = "^bareme_"))

# Extraction des noms pour créer les colonnes "Referentiel" et "Indice"
combined_baremes <- bind_rows(
  lapply(names(baremes_list), function(name) {
    # Séparation du nom pour extraire le référentiel et l'indice
    parts <- strsplit(name, "_")[[1]]
    referentiel <- paste(parts[2:(length(parts)-1)], collapse = "_")
    indice <- parts[length(parts)]
    
    # Ajout des colonnes "Referentiel" et "Indice"
    baremes_list[[name]] %>%
      mutate(Referentiel = referentiel, Indice = indice)
  })
)

# Affichage du résultat
print(combined_baremes)


# Définir les limites supérieures pour chaque indice
limites_sup <- c(RSr = 25, Shannon = 3, Pielou = 1, ab = 60000)

# Remplacer les bornes infinies dans combined_baremes
combined_baremes <- combined_baremes %>%
  mutate(
    Minimum = ifelse(is.infinite(Minimum), 0, Minimum),
    Maximum = case_when(
      Indice == "RSr" ~ pmin(Maximum, limites_sup["RSr"]),
      Indice == "Shannon" ~ pmin(Maximum, limites_sup["Shannon"]),
      Indice == "Pielou" ~ pmin(Maximum, limites_sup["Pielou"]),
      Indice == "ab" ~ pmin(Maximum, limites_sup["ab"]),
      TRUE ~ Maximum
    )
  ) %>%
  # Supprimer les classes avec une étendue nulle
  filter(Maximum > Minimum)

combined_baremes_RSr <- combined_baremes %>% filter(Indice == "RSr")
combined_baremes_Shannon <- combined_baremes %>% filter(Indice == "Shannon")
combined_baremes_Pielou <- combined_baremes %>% filter(Indice == "Pielou")
combined_baremes_ab <- combined_baremes %>% filter(Indice == "ab")

####################Richesse spécifique####################


# Créer un vecteur de correspondance entre les niveaux numériques et les labels souhaités
y_labels <- c(
  "gamme_RMQS" = "RMQS Biodiversité",
  "gamme_RMQS_Bretagne" = "RMQS BioDiv Bretagne",
  "gamme_ANDRA" = "ANDRA",
  "gamme_TIGA_rural" = "TIGA rural"
)

# Extraire les niveaux (référentiels) dans l’ordre voulu
referentiel_levels <- rev(names(y_labels))

# Créer le graphique
p_RSr <- ggplot(combined_baremes_RSr) +
  geom_rect(aes(
    xmin = Minimum, xmax = Maximum,
    ymin = as.numeric(factor(Referentiel, levels = referentiel_levels)) - 0.4,
    ymax = as.numeric(factor(Referentiel, levels = referentiel_levels)) + 0.4,
    fill = as.factor(Score)
  ), color = "black", size = 0.2) +
  scale_y_continuous(
    breaks = 1:4,
    labels = rev(unname(y_labels))
  ) +
  scale_fill_manual(values = score_colors, name = "Score attribué", guide = guide_legend(title.position = "top")) +
  labs(x = "Valeur de richesse spécifique") +
  theme_minimal() +
  theme(legend.position = "right")

# Afficher
p_RSr <- p_RSr + theme(legend.position = "none")
print(p_RSr)


########################Shannon##############
# Créer le graphique
p_Shannon <- ggplot(combined_baremes_Shannon) +
  geom_rect(aes(
    xmin = Minimum, xmax = Maximum,
    ymin = as.numeric(factor(Referentiel, levels = referentiel_levels)) - 0.4,
    ymax = as.numeric(factor(Referentiel, levels = referentiel_levels)) + 0.4,
    fill = as.factor(Score)
  ), color = "black", size = 0.2) +
  scale_y_continuous(
    breaks = 1:4,
    labels = rev(unname(y_labels))
  ) +
  scale_fill_manual(values = score_colors, name = "Score attribué", guide = guide_legend(title.position = "top")) +
  labs(x = "Valeur de richesse spécifique") +
  theme_minimal() +
  theme(legend.position = "right")

# Afficher
p_Shannon <- p_Shannon + theme(legend.position = "none")
print(p_Shannon)


###########################Pielou######################
# Créer le graphique
p_Pielou <- ggplot(combined_baremes_Pielou) +
  geom_rect(aes(
    xmin = Minimum, xmax = Maximum,
    ymin = as.numeric(factor(Referentiel, levels = referentiel_levels)) - 0.4,
    ymax = as.numeric(factor(Referentiel, levels = referentiel_levels)) + 0.4,
    fill = as.factor(Score)
  ), color = "black", size = 0.2) +
  scale_y_continuous(
    breaks = 1:4,
    labels = rev(unname(y_labels))
  ) +
  scale_fill_manual(values = score_colors, name = "Score attribué", guide = guide_legend(title.position = "top")) +
  labs(x = "Valeur de richesse spécifique") +
  theme_minimal() +
  theme(legend.position = "right")

# Afficher
p_Pielou <- p_Pielou + theme(legend.position = "none")
print(p_Pielou)

########################densité#########################
# Créer le graphique
p_ab <- ggplot(combined_baremes_ab) +
  geom_rect(aes(
    xmin = Minimum, xmax = Maximum,
    ymin = as.numeric(factor(Referentiel, levels = referentiel_levels)) - 0.4,
    ymax = as.numeric(factor(Referentiel, levels = referentiel_levels)) + 0.4,
    fill = as.factor(Score)
  ), color = "black", size = 0.2) +
  scale_y_continuous(
    breaks = 1:4,
    labels = rev(unname(y_labels))
  ) +
  scale_fill_manual(values = score_colors, name = "Score attribué", guide = guide_legend(title.position = "top")) +
  labs(x = "Valeur de richesse spécifique") +
  theme_minimal() +
  theme(legend.position = "right")

# Afficher
p_ab <- p_ab + theme(legend.position = "none")
print(p_ab)



##############graphique final######################
# Combiner les graphiques en colonne
final_plot <- plot_grid(
  p_RSr, p_Shannon, p_Pielou, p_ab,
  labels = c("A", "B", "C", "D"),
  ncol = 1, 
  align = "v"
)

# Extraire la légende en tant que grob
legend_score_grob <- get_legend(
  ggplot(combined_baremes) +
    geom_rect(aes(
      xmin = Minimum, xmax = Maximum,
      ymin = 0, ymax = 1, 
      fill = as.factor(Score)
    )) +
    scale_fill_manual(values = score_colors, name = "Score attribué") +
    theme_minimal()
)

# Ajouter la légende sur la droite
final_with_legend <- plot_grid(
  final_plot,
  legend_score_grob,
  ncol = 2,
  rel_widths = c(4, 0.5)  # Ajuste la largeur relative des colonnes
)

# Afficher le graphique final
print(final_with_legend)

#ggsave("baremes_tous_sites.pdf", plot = final_with_legend, width = 14, height = 10)






###############méthodes pour barèmes#################
library(ggplot2)
library(dplyr)

# Définir les couleurs des scores
score_colors <- c(
  "0" = "#8B4513",  # marron
  "1" = "#D73027",  # rouge
  "2" = "#FC8D59",  # orange clair
  "3" = "#FEE08B",  # jaune
  "4" = "#D9EF8B",  # vert clair
  "5" = "#91BFDB",  # bleu clair
  "6" = "#4575B4"   # bleu foncé
)

# Étendre les limites d'affichage de l'axe des x
x_extension <- 2
x_vals <- seq(-4 - x_extension, 4 + x_extension, length.out = 1000)
dens_vals <- dnorm(x_vals)

# Définir les bornes des classes
quants <- qnorm(c(0.00001, 0.2, 0.4, 0.6, 0.8, 0.99999))
bounds <- c(min(x_vals), quants, max(x_vals))

# Créer les segments pour la barre
bar_df <- data.frame(
  xmin = bounds[-length(bounds)],
  xmax = bounds[-1],
  score = factor(0:6)
)

# Créer le graphique
method <- ggplot() +
  # Remplissage sous la courbe segmenté selon les classes
  lapply(1:nrow(bar_df), function(i) {
    x_seg <- seq(bar_df$xmin[i], bar_df$xmax[i], length.out = 300)
    y_seg <- dnorm(x_seg)
    geom_area(data = data.frame(x = x_seg, y = y_seg),
              aes(x = x, y = y),
              fill = score_colors[as.character(bar_df$score[i])],
              alpha = 0.2)
  }) +
  
  # Courbe de densité normale
  geom_line(aes(x = x_vals, y = dens_vals), color = "black", size = 1) +
  
  # Barre colorée des scores
  geom_rect(data = bar_df,
            aes(xmin = xmin, xmax = xmax, ymin = -0.02, ymax = -0.005, fill = score),
            color = "black") +
  
  # Scores dans la barre
  geom_text(data = bar_df,
            aes(x = (xmin + xmax) / 2, y = -0.0125, label = score),
            size = 4) +
  
  # Palette des couleurs
  scale_fill_manual(values = score_colors, name = "Score attribué") +
  
  # Affichage
  coord_cartesian(xlim = c(min(x_vals), max(x_vals)), ylim = c(-0.05, 0.45)) +
  labs(x = "Valeur de l’indice", y = "Distribution dans le référentiel") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none"
  )

method

ggsave("methode_baremes.pdf", plot = method, width = 14, height = 10)
