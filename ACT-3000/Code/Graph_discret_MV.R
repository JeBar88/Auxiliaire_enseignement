### Act-3000
## Graphique fonction de masse bivarié
## Jérémie Barde

#### Package ####
library('plot3D') # Il faut intaller un logiciel externe, dumoins sur MacOS (XQuartz)
library('rgl')

#### Foncions utils ####
source(here::here("ACT-3000/Code/Fonctions_Utiles.R"))

### Paramètres
k0 <- 22
m1 <- 0:k0
m2 <- m1
lam <- c(5, 3)
a <- 0:3

### Poisson Teicher
fm12 <- array(numeric(length(m1)*2*length(a)),
              dim = c(length(m1), length(m1), length(a)))
for (i in a) {
  fm12[, , i + 1] <- dPoTeibiv(m1, m2, lam, i)
}

### Cas como et anti
Fm1 <- ppois(m1, lam[1]) 
Fm2 <- ppois(m2, lam[2]) 

fcomo <- dcomo(m1, m2, Fm1, Fm2)
fanti <- dcomo(m1, m2, Fm1, Fm2)

### Graphique pour le cas Poisson Teicher
par(mar = rep(2, 4))
## Changer 1 par 2,3,4 pour avoir les autres cas
hist3D(m1, m2, fm12[, , 1], axes = TRUE, border = 'black',
       alpha = 0.90,
       theta = 160,
       phi = 10,
       contour = TRUE,
       ticktype='detailed')



### Graphique pour les cas como et anti
## Changer fcomo par fanti
hist3D(m1, m2, fcomo, axes = TRUE, border = 'black',
       alpha = 0.9,
       theta = 40,
       phi = 15,
       contour = TRUE,
       ticktype='detailed')




















