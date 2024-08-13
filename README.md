# Auxiliaire_enseignement
Ce dépôt contient l'ensemble des documents et fichiers de code R développés durant ma période en tant qu'auxiliaire pour le cours Introduction à l'actuariat II et de Théorie du risque. Il inclut les ressources pédagogiques, les scripts de code et tout autre matériel créé dans le cadre de ce rôle.

# Description des dossiers

## ACT-2001

### Note explicative
Ce dossier contient des fichiers RMarkdown qui expliquent diverses notions du cours à travers plusieurs exemples d'application.
- Inver_Fn_Rep : comment inverser des fonctions de répartition par la méthode par interpolation dans R
- Optimize : comment utiliser la fonction `Optimize` de R
- Simulation : explication de la méthode de simulation Monte-Carlo et ses applications dans R
- TVaR : comment calculer la mesure de risque VaR et TVaR dans R pour différents cas de figure

### Code

Ce dossier contient des exemples de codes R pour diverses notions du cours.
- `Contriubtion.R` : calculer la contribution à la mesure VaR et TVaR selon le principe d'Euler par simulation
- `Simulation_SA.R` : simualtion de v.a. composée (4 techniques)


### Exercices 
Ce dossier contient les solutions de quelques exercices pour les chapitres 1 à 4.

### Examen pratique
Ce dossier contient les solutions pour une partie des examens pratiques.

## ACT-3000

### Note explicative
Ce dossier contient des fichiers RMarkdown qui expliquent diverses notions du cours à travers plusieurs exemples d'application.
- `Methode_Discretisation.Rmd` : exemple des différentes méthode de discrétisation vu dans le cours
- `Methode aggégation.Rmd` : exemple des méthodes d'aggrégations de base

### Code

Ce dossier contient des exemples de codes R pour diverses notions du cours.
- `Code_Base.R` : fonction de base pour le cours
- `Como_Anti` : calcule des fonctions de masse multivariées dans le cas comonotonne et antimonotonne
- `Contribution.R` : calcule exacte de la contribution au mesure de risque VaR et TVaR
- `Discret_bivariee_base.R` : calcule des caractéristiques pour une loi bivariée de base
- `FFT.R` : divers exemples de comment utiliser l'algorithme FFT pour faire de l'aggrégation
- `PoTeicher_comp_simul.R` : simulation de la loi Poisson Teicher composée
- `PoissonTeicher.R` : calcule de la fonction de masse et des caractéristique pour la loi Poisson Teicher
- `PoissonTeicher_Comp.R` : exemple de loi Poisson Teicher composée avec loi exponentielle/Erlang, bêta différent
- `loi_comp_multi.R` : exemple de loi multivarié arbitraire composée avec loi exponentielle/Erlang, bêta identique

### Ressources

Ce dossier contient des documents PDF présentant les démonstrations et des exemples numériques de diverses notions du cours.

- `Convo_Exp_Gam.pdf` : convolution de v.a. de loi exponentielle/Gamma lorsque les paramètres bêta son différent et gamma CRMM
- `Convo_VA_Comp.pdf` : convolution de v.a. composée avec loi exponentielle/Gamma lorsque les paramètres bêta son différent et v.a. multivariée discète
- `PoissonTeicher_ExpoEFGM.pdf` : preuves relative aux lois Poisson Teicher; pmf, espérance conditionelle et expo EFGM; espérance conditionelle


### Exercices 
Ce dossier contient les solutions de quelques exercices pour les chapitres 1 à 4.

### Examen pratique
Ce dossier contient les solutions pour une partie des examens pratiques.
