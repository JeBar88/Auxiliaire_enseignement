### Act-3000
## Graphique fonction de masse bivarié
## Jérémie Barde

##### Package #####
library('plot3D') # Il faut intaller un logiciel externe, dumoins sur MacOS (XQuartz)
library('rgl')

##### Foncions utils #####
dPoTei <- function(m1, m2, lam1, lam2, a){
  f <- matrix(numeric(length(m1) * length(m2)), ncol = length(m2))
  
  for (i in m1){
    for (j in m2){
      
      k <- 0:min(i, j)
      
      f[i + 1, j + 1] <-  sum(dpois(k, a) * dpois(i - k, lam1 - a) * dpois(j - k, lam2 - a))
    }
  }
  f
}
dcomo <- function(x1, x2, Fx1, Fx2){
  fx <- matrix(numeric(length(x1) * length(x2)), ncol = length(x2))
  pcomo <- function(x1, x2) {
    ifelse(x1 < 0 | x2 < 0, 0, pmin(Fx1[x1 + 1], Fx2[x2 + 1]))
  }
  
  for (i in x1){
    for (j in x2){
      fx[i + 1, j + 1] <-  pcomo(i, j) - pcomo(i-1, j) - pcomo(i, j-1) + pcomo(i-1, j-1)
    }
  }
  fx
}
danti <- function(x1, x2, Fx1, Fx2){
  fx <- matrix(numeric(length(x1) * length(x2)), ncol = length(x2))
  panti <- function(x1, x2) {
    ifelse(x1 < 0 | x2 < 0, 0, pmax(Fx1[x1 + 1] + Fx2[x2 + 1] - 1, 0))
  }
  
  for (i in x1){
    for (j in x2){
      
      fx[i + 1, j + 1] <-  panti(i, j) - panti(i-1, j) - panti(i, j-1) + panti(i-1, j-1)
    }
  }
  fx
}

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
  fm12[, , i + 1] <- dPoTei(m1, m2, lam[1], lam[2], i)
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




















