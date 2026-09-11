#multiCvM est une fonction qui permet de réaliser des tests de Cramer von Mises entre plus de deux distributions. Elle s'utilise ainsi : 
multiCvM(x, y, seed=, na.omit=)
#avec x la variable réponse et y la variable explicative
#seed= fixe la seed (fonction set.seed()) assurant la reproductibilité
#na.omit= permet d'exclure les NA (na.omit=T)


#Fonction : 
multiCvM <- function(x, y, seed = 1234, na.omit = TRUE) {
  set.seed(seed)  # Fixer la seed
  
  # Option pour supprimer les NA
  if (na.omit) {
    na_idx <- !is.na(x) & !is.na(y)
    x <- x[na_idx]
    y <- y[na_idx]
  }
  
  # Vérifier que y est un facteur
  y <- as.factor(y)
  
  # Obtenir toutes les paires de comparaisons
  combinaisons <- combn(levels(y), 2, simplify = FALSE)
  
  # Calcul des statistiques pour chaque paire
  results <- lapply(combinaisons, function(pair) {
    data1 <- x[y == pair[1]]
    data2 <- x[y == pair[2]]
    
    test_result <- cramer.test(data1, data2)
    p_val <- test_result$p.value
    W_val <- test_result$statistic
    
    # Format spécial pour les p-values très faibles
    p_label <- ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val))
    
    # Définir la catégorie de signification
    significance <- dplyr::case_when(
      p_val < 0.001 ~ "Très significatif",
      p_val < 0.05  ~ "Significatif",
      TRUE          ~ "Non significatif"
    )
    
    tibble::tibble(
      Groupe1 = pair[1], 
      Groupe2 = pair[2], 
      W_value = W_val,
      p_value = p_val, 
      p_label = p_label, 
      Significatif = significance
    )
  })
  
  # Fusionner les résultats en un dataframe
  results_df <- dplyr::bind_rows(results) 
  
  # Afficher le tableau complet dans la console
  print(results_df)
  
  # Définition des couleurs pour chaque niveau de significativité
  colors <- c("Non significatif" = "gray", 
              "Significatif" = "plum2", 
              "Très significatif" = "purple4")
  
  # Création du tile plot
  ggplot2::ggplot(results_df, ggplot2::aes(x = Groupe1, y = Groupe2, fill = Significatif)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = p_label), color = "white", size = 4) +
    ggplot2::scale_fill_manual(values = colors, name = "Significativité") +
    ggplot2::labs(title = "Test de Cramer-von Mises entre projets", x = "Projet 1", y = "Projet 2") +
    ggplot2::theme_minimal()
}

