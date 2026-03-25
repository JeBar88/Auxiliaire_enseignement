### ACT-2001
## Code de base pour v.a. discrètes
## # V.a. discrète, domaine finie
## # V.a. discrète, domaine infinie
## # V.a. discrète, domaine non aritmétique
## Example pour v.a. discrète
## # Multiplication par une constante
## # Franchise et limite
## Jérémie Barde

#### V.a. discrète domaine finie -- M~Bin(n = 10, q = 0.25) ####
### Définir les vairiables
k <- 0:10 # Le domaine de la v.a M est k=1,2,...,10. Si on mets plus les props seront 0.
n <- 10
q <- 0.25
u <- 0.95

### Fonction de masse de probabilité (fmp) -- dloi
fm <- dbinom(k, n, q)
## Pr(M = 5)
fm[5 + 1] # Il faut faire +1 car les vecteur dans R commence à 1, fm[0] retourne rien
cbind(k, fm) # Visuel pour mieux comprendre

### Fonction de répartition (cdf) -- ploi
Fm <- pbinom(k, n, q)
## Pr(M <= 5)
Fm[5 + 1]

### Espérance -- E[M] = \sum_{k = 0}^10 k*Pr(M = k)
Em <- sum(k*fm)
# Verif
cbind("Em_test"=Em, "Em"=n*q) 

### Variance -- \sum_{k = 0}^10 k^2*Pr(M = k) - E[M]^2
Vm <- sum(k^2*fm) - Em^2
cbind("Vm_test"=Vm, "Em"=n*(1-q)*q) 

### Skewness et Kurtosis -- E[((M - E[M])/sigma)^3] et E[((M - E[M])/sigma)^4]
sum(((k - Em)/sqrt(Vm))^3 * fm)
sum(((k - Em)/sqrt(Vm))^4 * fm)

### Stop-Loss -- E[max(M - d, 0)]
d <- 0.5
SLm <- sum(pmax(k - d, 0)*fm)

### Espérance tronquée -- E[M * I(2 < M <= 5)]
EmT <- sum(k*I(k > 2 & k <= 5)*fm)

### Espérance conditionelle -- E[M | M <= d]
d <- 6
Emc <- sum(k*I(k <= 6)*fm)/Fm[6 + 1]

### Fonction inverse ou Value-at-Risk (VaR) -- qloi
VaRm <- qbinom(u, n, q)
Fm[VaRm + 1] # Remarque on ne tombe pas toujours directement sur u avec les v.a. discrètes

### Tail-Value-at-Risk (TVaR)
## Méthode 1 : Formule avec l'espérance tronquée
TvaRm1 <- 1/(1 - u)*(sum(k*I(k > VaRm)*fm) + VaRm*(Fm[VaRm + 1] - u))

## Méthode 2 : Formule avec la Stop-Loss
TvaRm2 <- VaRm + sum(pmax(k - VaRm, 0)*fm)/(1 - u) # Recommandation de Jé!

### Mesure entropique -- 1/rho*ln(fgmM(rho))
rho <- 0.5
Entm <- 1/rho*log(sum(exp(k*rho)*fm))
# Vérif: La mesure entropique tend vers E[M] quand rho -> 0

#### V.a. discrète domaine infinie -- M1~Po(lam1 = 5) et M2~Po(lam2 = 500) ####
### Le domaine est infinie, donc il faut mettre un domaine assez grand pour que la fmp somme à 1
### Par exemple kmax=100
kmax <- 100
k <- 0:kmax
lam <- c(5, 80)

### fmp
fm1 <- dpois(k, lam[1])
fm2 <- dpois(k, lam[2])
sum(fm1) # Good!
sum(fm2) # Not Good!, normale la moyenne est 80 et le domaine à un valeur max de 100

#### Par exemple kmax=1000
kmax <- 1000
k <- 0:kmax
lam <- c(5, 80)

### fmp
fm1 <- dpois(k, lam[1])
fm2 <- dpois(k, lam[2])
sum(fm1) # Good!
sum(fm2) # Good!

### Fonction de répartition
Fm1 <- dpois(k, lam[1])
Fm2 <- dpois(k, lam[2])
#### Les caractéristiques ce calcule exactement de la même façon ensuite

#### V.a. discrète, domaine non aritmétique -- k=(5, 12, 64, 89, 230) #### 
### Il n'est pas possible d'utiliser les fonctions dloi et ploi
p <- c(0.45, 0.1, 0.25, 0.15, 0.05)
vk <- c(5, 12, 64, 89, 230)
u <- 0.9

### Construire un vecteur complet
kmax <- max(vk)
k <- 0:kmax
fm <- replace(numeric(kmax + 1), vk + 1, p) # Les probabilité p on été mis au bonne place

### Foncton de répartition : on ne peux pas utiliser ploi
Fm <- cumsum(fm)
# Pr[M <= 60]
Fm[60 + 1]

### VaR : on ne peux pas utiliser qloi, on utilise la définiton VaR=inf{k;Pr(M >= k) = u}
VaR <- k[min(which(Fm >= u))]
cbind(k, Fm) # La première valeurs de k qui donne un Pr(M >= k) > u est 89

#### Multiplication par une constante -- Y = 10X, X~Po(5) ####
k <- 0:1000
lam <- 5
b <- 10

### Espérance théorique
ExTh <- lam
EyTh <- b*ExTh

### Fmp de la v.a. Y -- Pr(Y = k) = Pr(bX = k) = Pr(X = k/5)
fy <- dpois(k/b, lam) # Warings normale, 

### Espérance
Ey <- sum(k*fy)
cbind(EyTh, Ey)

#### Franchise et limite -- Y = min(max(X - 2, 0), 8) ####
k <- 0:1000
r <- 2.5
q <- 0.7
u <- 0.9

### Fmp de X
fx <- dnbinom(k, r, q)
sum(fx*k)

### Fmp de Y
ky <- pmin(pmax(k - 2, 0), 8) # On transforme le domaine
k <- unique(ky) # On garde seulement les valeurs uniques
fy <- tapply(fx, ky, sum) # On somme les fx selon les valeurs de ky

