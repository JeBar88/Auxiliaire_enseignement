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

- `Code_typique.R` : code typique à savoir faire impérativement. Ex : espérance, variance, stop-loss, VaR et TVaR.
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

### Note explicative
Ce dossier contient des fichiers RMarkdown qui expliquent diverses notions du cours à travers plusieurs exemples d'application.
- `Methode_Discretisation.Rmd` : exemple des différentes méthodes de discrétisation vu dans le cours.
- `Methode agrégation.Rmd` : exemple des méthodes d'agrégation de base.

### Code
Ce dossier contient des exemples de codes R pour diverses notions du cours.
- `Code_Base.R` : fonction de base pour le cours.
- `Como_Anti` : calcule des fonctions de masse multivariées dans le cas comonotonne et antimonotonne.
- `Contribution.R` : calcule exacte de la contribution aux mesures de risque VaR et TVaR.
- `Discret_bivariee_base.R` : calcule des caractéristiques pour une loi bivariée de base.
- `FFT.R` : divers exemples de comment utiliser l'algorithme FFT pour faire de l'agrégation.
- `Graph_discret_MV.R` : construction de graphique 3D pour les fonctions de masse de probabilité bivariée.
- `PoTeicher_comp_simul.R` : simulation de la loi Poisson Teicher composée.
- `PoissonTeicher.R` : calcule de la fonction de masse et des caractéristiques pour la loi Poisson Teicher.
- `PoissonTeicher_Comp.R` : exemple de loi Poisson Teicher composée avec loi exponentielle/Erlang, bêta différent.
- `loi_comp_multi.R` : exemple de loi multivariée arbitraire composée avec loi exponentielle/Erlang, bêta identique.

### Ressources
Ce dossier contient des documents PDF présentant les démonstrations et des exemples numériques de diverses notions du cours.

- `Convo_Exp_Gam.pdf` : convolution de v.a. de loi exponentielle/Gamma lorsque les paramètres bêta son différent et gamma CRMM.
- `Convo_VA_Comp.pdf` : convolution de v.a. composée avec loi exponentielle/Gamma lorsque les paramètres bêta son différent ou identique et v.a. multivariée discrète.
- `Panjer_DePril.pdf` : preuve pour les algorithmes de Panjer et DePril.
- `PoissonTeicher_ExpoEFGM.pdf` : preuves relatives aux lois Poisson Teicher; pmf, espérance conditionnelle et expo EFGM; espérance conditionnelle.
- `Tetraedre.pdf` : représentation géométrique d'une classe de Fréchet lorsque les marginals sont Bernoulli.


### Exercices 
Ce dossier contient les solutions de quelques exercices.

### Examen pratique
Ce dossier contient les solutions pour une partie des examens pratiques.
