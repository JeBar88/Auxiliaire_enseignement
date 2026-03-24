### ACT-2001
## Code de base à maîtriser pour v.a. discrètes
## Code pour v.a. discrète
## Jérémie Barde

#### Exemple typique ####
### X in {0, 5, 10, 15, 20}
x <- seq(0, 20, 5)
fx <- c(0.05, 0.1, 0.5, 0.3, 0.05) # serait donné

## Fonction de répartition
Fx <- cumsum(fx)

## Espérance et variance
Esp <- sum(x * fx)
Var <- sum(x^2 * fx) - Esp^2
cbind(Esp, Var)

## Skewness et Kurtosis
sum(((x - Esp)/sqrt(Var))^3 * fx)
sum(((x - Esp)/sqrt(Var))^4 * fx)

## Fonction stop-loss
d <- 2
SL <- sum(pmax(x - d, 0)*fx)
SL

## E[min(X, 5)]
d <- 5
sum(pmin(x, d)*fx)

## VaR et TVaR
u <- 0.9
VaR <- x[min(which(Fx >= u))]
TVaR <- VaR + sum(pmax(x - VaR, 0)*fx)/(1 - u) # approche par la Stop-loss
cbind(VaR, TVaR)

## Mesure entropique
p <- 0.1 # jouer avec rho pour comprendre la dynamique
ME <- 1/p*log(sum(exp(p*x)*fx))
ME











#### Code pour v.a. discrète -- Calcules des caractéristiques ####
#### Exemple 1 : X~Bin(n = 10, q = 0.25)
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
Fm <- cumsum(fm)
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



