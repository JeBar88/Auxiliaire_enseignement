# Auxiliaire_enseignement
Ce dépôt contient l'ensemble des documents et fichiers de code R développés durant ma période en tant qu'auxiliaire pour le cours Introduction à l'actuariat II et de Théorie du risque. Il inclut les ressources pédagogiques, les scripts de code et tout autre matériel créé dans le cadre de ce rôle.

# Description des dossiers

## ACT-2001

### Note explicative
Ce dossier contient des fichiers RMarkdown qui expliquent diverses notions du cours à travers plusieurs exemples d'application.
- Inver_Fn_Rep : comment inverser des fonctions de répartition par la méthode par interpolation dans R.
- Optimize : comment utiliser la fonction `Optimize` de R.
- Simulation : explication de la méthode de simulation Monte-Carlo et ses applications dans R.
- TVaR : comment calculer la mesure de risque VaR et TVaR dans R pour différents cas de figure.

### Code
Ce dossier contient des exemples de codes R pour diverses notions du cours.

- `Code_typique.R` : code de base du cours. Ex : espérance, variance, stop-loss, VaR et TVaR.
- `Contribution.R` : calculer les contributions au mesure VaR et TVaR par simulation.
- `Convolution.R` : convolution de v.a. aléatoire discrète.
- `Simulation_IC.R` : calculer des intevalles de confiance pour les approximations faites par simulation.
- `Simulation_SA.R` : simulation de v.a. composée (4 techniques).

### Ressources
Ce dossier contient des documents PDF présentant les démonstrations et des exemples numériques de diverses notions du cours.
- `Convo_PoComp.pdf` : preuve de la loi de la somme de v.a. de loi Poisson composée.
- `Convo_Uniforme_Beta.pdf` : exemple de convolution pour deux loi bêta.
- `Mesure_risque.pdf` : Résumé des propriétés suivies par certaines mesures de risque, avec preuves pour certaines d'entre elles.
- `Preuve_EspT_TVaR.pdf` : Preuves des expressions de l'espérance tronquée, de la VaR et de la TVaR pour plusieurs lois.

### Exercices 
Ce dossier contient les solutions de quelques exercices pour les chapitres 1 à 4.

### Examen pratique
Ce dossier contient les solutions pour une partie des examens pratiques.

## ACT-3000

### Code
Ce dossier contient des exemples de codes R pour diverses notions du cours.
- `Code_Base_FFT.R` : Exemples de bases univariées et exemple d'utilisation de l'algorithme FFT pour l'agrégation.
- `Contribution.R` : Exemples de code pour le calcul exact de la contribution aux mesures de risque VaR et TVaR.
- `Discretisation.R` : Exemples de code pour l'utilisation de la discrétisation.
- `Fonctions_Utiles.R` : Exemples de code pour différentes fonctions utiles dans le cadre du cours.
- `Graph_discret_MV.R` : Exemple de construction de graphiques 3D pour les fonctions de masse de probabilité bivariée.
- `Lois_Discrete_Comp_Multi.R` : Exemples de code pour les v.a. discrètes composées multivariées.
- `Lois_Discrete_Multi.R` : Exemples de code pour les v.a. discrètes multivariées.
- `Lois_Multi_simul.R` : Exemples de code pour la simulation de v.a. multivariées composées.

### Ressources
Ce dossier contient des documents PDF présentant des démonstrations et des exemples numériques de diverses notions du cours.
- `Convo_Exp_Gam.pdf` : Convolution de v.a. de loi exponentielle/Gamma lorsque les paramètres bêta sont différents et Gamma CRMM.
- `Convo_VA_Comp.pdf` : Convolution de v.a. composées avec loi exponentielle/Gamma lorsque les paramètres bêta sont différents ou identiques, et v.a. multivariée discrète.
- `Panjer_DePril.pdf` : Preuves pour les algorithmes de Panjer et DePril.
- `PoissonTeicher_ExpoEFGM.pdf` : Preuves relatives aux lois Poisson-Teicher (pmf, espérance conditionnelle) et exponentielle EFGM (espérance conditionnelle).
- `Tetraedre.pdf` : Représentation géométrique d'une classe de Fréchet lorsque les marginales sont Bernoulli.

### Exercices 
Ce dossier contient les solutions de quelques exercices.

### Examen pratique
Ce dossier contient les solutions de quelques examens pratiques.
