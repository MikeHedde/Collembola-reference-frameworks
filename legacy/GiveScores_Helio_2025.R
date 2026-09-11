##################GiveScores##############
#la fonction GiveScores permet d'attribuer des scores à un jeu de données en se basant sur ses indices par rapport à un jeu de données référentiel
#elle s'exécute de la sorte 
GiveScores(x, ref =, var =, nb_class = 5, inv = FALSE)
#avec x le jeu de données à noter
#ref = précise le jeu de données référentiel
#var = précise la variable (indice) à prendre en compte pour la notation
#nb_class = précise le nombre de classes de scores. La première et la dernière classes de scores correspondent à des valeurs d'indices en dehors de l'étendue du référentiel. Les classes intermédiaires se partagent l'étendue du référentiel de manière égale.
#inv= précise si l'attribution des scores doit se faire de façon inverse. inv = F induira des scores croissants lorsque la valeur d'indice est croissante, alors que inv = T induira des scores décroissants lorsque la valeur d'indice est croissante
  
#GiveScores crée également un data.frame nommé dynamiquement bareme_referentiel_indice qui s'affiche lorsque la fonction est exécutée. Ce data.frame spécifie les limites des classes de valeurs d'indices pour chaque note. 

#A noter que si une valeur de l'indiviu statistique noté tombe exactement sur une limite, la note de la classe supérieure lui est attribué.


#fonction
  GiveScores <- function(x, ref, var, nb_class = 5, inv = FALSE) {
    if (nb_class < 3) stop("nb_class doit être >= 3 pour inclure les bornes extrêmes + au moins une classe centrale.")
    
    # Vérification des variables
    if (!(var %in% colnames(ref))) stop(paste0("La variable '", var, "' n'existe pas dans le jeu de données référentiel."))
    if (!(var %in% colnames(x))) stop(paste0("La variable '", var, "' n'existe pas dans le jeu de données cible."))
    
    # Nettoyage NA dans le référentiel
    ref_noNA <- ref[!is.na(ref[[var]]), ]
    
    # Nombre de classes internes
    nb_internal <- nb_class - 2
    
    # Calcul des bornes internes selon l'ordre souhaité
    internal_breaks <- quantile(ref_noNA[[var]], 
                                probs = seq(0, 1, length.out = nb_internal + 1), 
                                na.rm = TRUE, names = FALSE)
    
    # Création des bornes finales
    breaks <- c(-Inf, internal_breaks, Inf)
    
    # Attribution des scores
    scores <- 0:(nb_class - 1)
    if (inv) {
      scores <- rev(scores)
    }
    
    # Attribution des scores aux valeurs de x
    var_values <- x[[var]]
    score <- cut(var_values,
                 breaks = breaks,
                 include.lowest = TRUE,
                 labels = scores,
                 right = FALSE)
    
    # Création du nom de la colonne
    ref_name <- deparse(substitute(ref))
    new_col_name <- paste0("score_", var, "_", ref_name)
    
    # Ajout au dataset
    x[[new_col_name]] <- as.integer(as.character(score))
    
    # Création du tableau des classes
    class_table <- data.frame(
      Minimum = head(breaks, -1),
      Maximum = head(tail(breaks, -1), nb_class),
      Score = scores
    )
    rownames(class_table) <- paste("Classe", scores)
    
    # Nom dynamique pour le tableau des classes
    bareme_name <- paste0("bareme_", ref_name, "_", var)
    assign(bareme_name, class_table, envir = .GlobalEnv)
    
    # Affichage du tableau
    print(class_table)
    
    # Retour du dataset
    return(x)
  }
