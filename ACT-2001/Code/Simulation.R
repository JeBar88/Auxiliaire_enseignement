### ACT-2001
## Simulation de v.a.
## # Concept de base de la simulation
## # Approximation par simulation v.a. discrète
## # Approximation par simulation v.a. continue
## # Approximation par simulation somme fini de v.a.
## # Approximation par simulation somme fini de v.a. iid
## # Approximation par simulation v.a. composée
## # Approximation par simulation v.a. mélange
## # Approximation par simulation de la contribution
## # Intervalle de confiance pour approximation par simulation
## Jérémie Barde

#### Package utiles ####
library(actuar)

#### Concept de base ####
### On utilise le théorème de la fonction quantile

### Simuler une réalisation de M~Po(5)
## Remarque la réalisation change si le code est relancé
U <- runif(1)
M <- qpois(U, 5)
cbind(U, M)

## On Ajoute un acrange pour avoir des résultats reproductible
set.seed(1943)
U <- runif(1)
M <- qpois(U, 5)
cbind(U, M)

## Simuler 10 milles réalisations de M
m <- 1e4
lam <- 5
U <- runif(m)
M <- qpois(U, lam)
head(cbind(U, M))

## Utiliser les focntion rloi (N'existe pas pour toutes les lois)
m <- 1e4
lam <- 5
M <- rpois(m, lam)
head(M)

### Simuler 10 milles réalisations de M où fx=(0.1, 0.15, 0.3, 0.2, 0.25), k=1,...5
## Il n'est pas possible d'utiliser dloi, on utilise la méthode l'inverse
fx <- c(0.1, 0.15, 0.3, 0.2, 0.25)
k <- 0:5

## Fonction de répartition
Fx <- cumsum(fx)

## Fonction inverse ou VaR
VaRm <- function(u) k[min(which(Fx >= u))]

## Simultion de 10 milles réalisation
m <- 1e4
U <- runif(m)
M <- sapply(U, VaRm)
head(M)

#### V.a. discrète ####
#### M~Bin(10, 0.25)
m <- 1e6
n <- 10
q <- 0.25
u <- 0.95

### Simulation de 1 millions de réalisations de M
M <- rbinom(m, n, q)
# ou
# M <- qbinom(runif(m), n, q)

### Approximation de la fmp -- Pr(M = 5)
mean(M == 5) # Vrai valeur : 0.0583992

### Approximation de la cdf -- Pr(M <= 5)
mean(M <= 5) # Vrai valeur : 0.9802723
# ou
Fn <- ecdf(M)
Fn(5)

### Approximation de l'espérance
Em <- mean(M) 
# verif
cbind(Em, "EmTh"=n*q)

### Approximation de la variance et l'ecart-type
Vm <- var(M)
sdm <- sd(M)
# verif
cbind(Vm, "VmTh"=n*q*(1 - q))
cbind(sdm, "sdmTh"=sqrt(n*q*(1 - q)))

### Approximation du Skewness et du Kurtosis 
mean(((M - Em)/sdm)^3) # Vrai valeur : 0.3651484
mean(((M - Em)/sdm)^4) # Vrai valeur : 2.933333

### Approximation Stop-Loss
d <- 0.5
mean(pmax(M - d, 0)) # VRai valeur : 2.028157

### Approximation espérance tronquée -- E[M * I(2 < M <= 5)]
mean(M*I(M > 2 & M <= 5)) # Vrai valeur : 1.626835

### Approximation espérance conditionelle -- E[M | M <= d]
d <- 6
mean(M * I(M <= 6))/mean(M <= 6) # Vrai valeur : 2.483721

### Approximation de la mesure de risque VaR
VaRm <- quantile(M, u) # Vrai valeur : 5
# ou
# sort(M)[u*m]

### Approximation de la mesure de risque TVaR
VaRm + mean(pmax(M - VaRm, 0))/(1 - u) # Vrai valeur : 5.473595

### Approximation mesure entropique
rho <- 0.5
1/rho*log(mean(exp(M*rho))) # Vrai valeur : 3.005957

#### V.a. continue ####
#### X~Ga(2, 0.1)
### Même chose que pour les lois discrètes
m <- 1e6
a <- 2
b <- 0.1
u <- 0.95

### Simulation de 1 millions de réalisations de X
X <- rgamma(m, a, b)
# ou
# X <- qgamma(runif(m), n, q)

### Approximation de la cdf -- Pr(X <= 5)
mean(X <= 50) # Vrai valeur : 0.9595723

### Approximation de l'espérance
Ex <- mean(X) 
# verif
cbind(Ex, "ExTh"=a/b)

### Approximation de la variance et l'ecart-type
Vx <- var(X)
sdx <- sd(X)
# verif
cbind(Vx, "VxTh"=a/b^2)
cbind(sdm, "sdxTh"=sqrt(a/b^2))

### Approximation du Skewness et du Kurtosis 
mean(((X - Ex)/sdx)^3) # Vrai valeur : 
mean(((X - Ex)/sdx)^4) # Vrai valeur :

### Approximation Stop-Loss
d <- 5
mean(pmax(X - d, 0)) # VRai valeur : 0.9595723

### Approximation espérance tronquée -- E[X * I(X < 20)]
d <- 20
mean(X*I(X < d)) # Vrai valeur : 6.466472

### Approximation espérance conditionelle -- E[X | X < d]
d <- 20
mean(X * I(X < 20))/mean(X < 20) # Vrai valeur : 10.88642

### Approximation de la mesure de risque VaR
VaRx <- quantile(X, u) # Vrai valeur : 47.43865
# ou
# sort(X)[u*m]

### Approximation de la mesure de risque TVaR
mean(X[X > VaRx]) # Si loi continue seulement
VaRx + mean(pmax(X - VaRx, 0))/(1 - u) # Vrai valeur : 59.26517

### Approximation mesure entropique
rho <- 0.05
1/rho*log(mean(exp(X*rho))) # Vrai valeur : 27.72589



#### Somme fini de v.a. ####
#### X1~Ga(a, b) et X2~ln(mu, sig)
m <- 1e6
a <- 5
b <- 0.2
sig <- 0.6
mu <- log(20) - sig^2/2

### Espérance théorique
Ex1Th <- a/b 
Ex2Th <- mlnorm(1, mu, sig)
EsTh <- Ex1Th + Ex2Th

### Simuler 1 million de réalisation de Xi
X1 <- rgamma(m, a, b, )
X2 <- rlnorm(m, mu, sig)

# Vérif
cbind("Ex1"=mean(X1), Ex1Th)
cbind("Ex2"=mean(X2), Ex2Th)

### Simuler 1 million de réalisation de S
S <- X1 + X2

# Vérif
cbind("Es"=mean(S), EsTh)

#### Somme fini de v.a. iid ####
#### Exemple 1 : Mi~M~BinNeg(50, 0.6), i=1,...,5
m <- 1e6
n <- 5
r <- 50
q <- 0.45

### Espérance théorique
Ex <- r*(1 - q)/q
Es <- n*Ex

### Simuler 1 million de réalisation de Mi
X1 <- rnbinom(m, r, q)
X2 <- rnbinom(m, r, q)
X3 <- rnbinom(m, r, q)
X4 <- rnbinom(m, r, q)
X5 <- rnbinom(m, r, q)
mean(X1)

### Simuler 1 million de réalisation de S
S <- X1 + X2 + X2 + X4 + X5
mean(S)

## ATTENTION !!! ## 
# Ceci n'est pas vrai S = nX
# Cette simplification est le résultat d'une relation de dépendance précise (ACT-3000) 
# Exemple
S_test <- n*X1
cbind(Es, "Test"=mean(S_test))

#### Exemple 2 : Mi~M~BinNeg(50, 0.6), i=1,...,25
m <- 1e6
n <- 25
r <- 50
q <- 0.45

### Espérance théorique
Ex <- r*(1 - q)/q
Es <- n*Ex

### Simuler 1 million de réalisation de Mi
X <- replicate(n, rnbinom(m, r, q))
dim(X)
apply(X, 2, mean)

### Simuler 1 million de réalisation de S
S <- rowSums(X)
mean(S)

## ATTENTION !!! ## 
# Ceci n'est pas vrai S = nX
# Cette simplification est le résultat d'une relation de dépendance précise (ACT-3000) 
# Exemple
S_test <- n*X[, 1]
cbind(Es, "Test"=mean(S_test))

#### V.a. composée #####
#### X~PoComp(2, Fb) avec B~Pareto(3, 10)
m <- 1e5
a <- 3
b <- 10
lam <- 2
u <- 0.99

### Méthode 1a
M <- numeric(m)
W <- numeric(m)
for (i in 1:m) {
  M[i] <- qpois(runif(1), lam)
  
  if(M[i] > 0){
    W[i] <- sum(qpareto(runif(M[i]), a, b))
  }
  
}

### Méthode 1b
M <- numeric(m)
X <- numeric(m)
for (i in 1:m) {
  M[i] <- rpois(1, lam)
  
  if(M[i] > 0){
    X[i] <- sum(rpareto(M[i], a, b))
  }
  
}

#### Méthode 2a
M <- rpois(m, lam)
Y <- numeric(m)
for (i in 1:m) {
  Y[i] <- sum(rpareto(M[i], a, b))
}

#### Méthode 2b
M <- rpois(m, lam)
Z <- sapply(1:m, function(i) sum(rpareto(M[i], a, b)))

#### Result
data <- cbind(W, X, Y, Z)
apply(data, 2, mean)
apply(data, 2, quantile, u)












#### V.a. Mélange ####
#### Exemple 1 : X|Y~Po(Y) où fy=(2/5, 3/5) avec k={4, 10}
m <- 1e6
k <- c(4, 10)
fy <- c(2/5, 3/5)

### Espérance théorique
Ex <- sum(k*fy)

### (1) Simuler 1 millions de réalisations de Y
Fy <- cumsum(fy)
Y <- replicate(m, k[min(which(Fy >= runif(1)))])

### (2) Simuler 1 millions de réalisations de X
X <- rpois(m, Y)

# Vérif
cbind("Ex_test"=mean(X), Ex)

#### Exemple 2: X|Y~Exp(Y) où Y~Ga(2, 0.5)
m <- 1e6
a <- 2
b <- 0.5

### Espérance théorique -- E[E[X|Y]]=E[1/Y]
Ex <- b/(a - 1) # Espérance d'une loi Pareto

### (1) Simuler 1 millions de réalisations de Y
Y <- rgamma(m, a, b)

### (2) Simuler 1 millions de réalisations de X
X <- rexp(m, Y)

# Vérif
cbind("Ex_test"=mean(X), Ex)

#### Contribution ####
m <- 1e6
a <- c(2.5, 3)
lam <- c(150, 200)
u <- 0.95

### Espérance Théorique
ExTh <- lam/(a - 1) 
Es <- sum(ExTh)

### Simuler 1 millions de réalisation de Xi
X <- mapply(function(a, lam) rpareto(m, a, lam), a, lam)
apply(X, 2, mean)

### Variable S
S <- rowSums(X)
mean(S)

## VaR
VaRs <- sort(S)[u*m] 
# ou
quantile(S, u, type=1) # type=1 car on veux la valeur exacte dans S

## TVaR
TVaRs <- mean(S[S > VaRs])

## Contribution à la VaR : Pas une très bonne approximation
CVaRx1 <- sum(X[, 1] * I(S == VaRs))
CVaRx2 <- sum(X[, 2] * I(S == VaRs))

# Verif
cbind(CVaRx1 + CVaRx2, VaRs)

## Contribution à la VaR : Pas une très bonne approximation
TVaR <- function(k) mean(S[S > sort(S)[k * m]])
CTVaRx1 <- function(k) mean(X[, 1] * I(S > sort(S)[k*m]))/(1 - k)
CTVaRx2 <- function(k) mean(X[, 2] * I(S > sort(S)[k*m]))/(1 - k)

kappa <- 0.9
cbind(CTVaRx2(kappa) + CTVaRx1(kappa), TVaR(kappa))


















#### Intervalle de confiance ####
m <- 1e4 # modfier pour voir l'impacte de m sur les intervalles
a <- 3
b <- 0.1
lam <- 60

### Espérance théorique
ExTh <- a/b
EyTh <- lam/(a - 1)
EsTh <- ExTh + EyTh

### Simuler des réalisation de X et Y
X <- rgamma(m, 3, 0.1)
Y <- rpareto(m, 3, 60)
S <- X + Y
mean(S)

### Fonction de répartition empirique
Fn <- function(x) mean(S <= x)

## Intervalle de confiance
ES_Fn <- sd(S <= 50)/sqrt(m)
Fn(50) + c(-1, 0, 1) * qnorm(0.975)*ES_Fn

### Moyenne
Es <- mean(S)

## Intervalle de confiance
ES_Es <- sd(S)/sqrt(m)
Es + c(-1, 0, 1) * qnorm(0.975)*ES_Es

### Stop-loss
d <- 20
SLs <- mean(pmax(S - d, 0))

## Intervalle de confiance
ES_SLs <- sd(pmax(S - d, 0))/sqrt(m)
SLs + c(-1, 0, 1) * qnorm(0.975)*ES_SLs

### Mesure VaR
u <- 0.9
VaRs <- sort(S)[u*m]

## Intervalle de confiance
sdR <- sqrt(m*u*(1 - u))
k0 <- floor(qnorm(0.975)*sdR + 0.5)
VaRs_IC <- sort(S)[(m*u) + c(-1, 1) * k0]
c("VaRs"=VaRs, VaRs_IC)

### Mesure TVaR
TVaRs <- mean(S[S > VaRs]) # ou VaR + mean(pmax(S - VaR, 0))/(1 - u)

## Intervalle de confiance
ES_TVaRs <- sqrt(var(pmax(S - VaRs, 0))/((1 - u)^2*m))
TVaRs_IC <- TVaRs + c(-1, 1)*qnorm(0.975)*ES_TVaRs
c("TVaRs"=TVaRs, TVaRs_IC)

### Mesure entropique
p <- 0.001
entros <- 1/p*log(mean(exp(p*S)))

## Intervalle de confiance
# On commence par faire l'intervalle de la fgm
fgms <- mean(exp(p*S))
ES_fgms <- sd(exp(p*S))/sqrt(m)
fgms_IC <- fgms + c(-1, 1)*qnorm(0.975)*ES_fgms

entros_IC <- 1/p*log(fgms_IC)
c("entros"=entros, entros_IC)









