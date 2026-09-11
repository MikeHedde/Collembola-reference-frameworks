# Fonction pour extraire les stats souhaitées
stats_var <- function(x) {
  c(Moyenne = mean(x, na.rm = TRUE),
    Ecart_type = sd(x, na.rm = TRUE),
    Skewness = skewness(x, na.rm = TRUE),
    Kurtosis = kurtosis(x, na.rm = TRUE)-3
  )
}

# Fonction pour automatiser stats_var
stats_varall <- function(x, y) {
  # Transformer y en facteur pour s'assurer de l'ordre des niveaux
  y <- as.factor(y)
  niveaux <- levels(y)
  
  # Appliquer stats_var à chaque modalité de y
  resultats <- lapply(niveaux, function(niv) {
    stats_var(x[y == niv])
  })
  
  # Nommer les résultats avec les niveaux
  noms <- niveaux
  names(resultats) <- noms
  
  # Combiner les résultats dans un data.frame
  # Si stats_var retourne un vecteur nommé
  df <- do.call(cbind, resultats)
  colnames(df) <- noms
  return(as.data.frame(df))
}