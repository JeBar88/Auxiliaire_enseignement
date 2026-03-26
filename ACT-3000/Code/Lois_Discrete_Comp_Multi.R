### ACT-3000
## Exemple de base v.a multivariée composée
## Loi discrète quelconque composée Expo/Gamma, bêta indentique
## Jérémie Barde

#### Package utile ####
library(actuar)

#### Fonction utile ####
source(here::here("ACT-3000/Code/Fonctions_Utiles.R"))

#### Loi discrète bivariée composée -- Bi~Ga(ai, b) ####
fm12 <- matrix(c(0.30, 0.10, 0.02, 0.03,
                 0.20, 0.10, 0.01, 0.01,
                 0.04, 0.02, 0.05, 0.03,
                 0.01, 0.02, 0.01, 0.05), ncol = 4, byrow =TRUE)
k <- 0:3
a <- c(2.1, 3.4)
b <- 0.1
u <- 0.95

### Marignales
fm1 <- rowSums(fm12)
fm2 <- colSums(fm12)

### Espérance de Xi et S
Emi <- c(sum(k*fm1), sum(k*fm2))
Ex <- a/b*Emi 
Es <- sum(Ex)

### Masse de probabilité pij
pm <- as.vector(fm12)

### Calcules paramètres loi gamma
comb <- arrayInd(seq_along(fm12), dim(fm12)) - 1
am1 <- a[1]*comb[, 1]
am2 <- a[2]*comb[, 2]
am <- am1 + am2

### Fonction de répartition conjointe
Fx12 <- function(x1, x2) pm[1] + sum(pm[-1] * pgamma(x1, am1[-1], b)*pgamma(x2, am2[-1], b))
Fx12(2, 5)

### Fonction de répartition de S
Fs <- function(x) pm[1] + sum(pm[-1] * pgamma(x, am[-1], b))
Fs(500)

### Vérif : Espérance
Es_test <- sum(pm*am/b)
cbind(Es_test, Es)

### Autres caractéristiques
SL <- function(d) sum(pm*(am/b * pgamma(d, am + 1, b, low=FALSE) - d*pgamma(d, am, b, low=FALSE)))
VaRS <- function(k) optimize(function(x) abs(Fs(x) - k), c(0, 500))$min
TVaRS <- function(k){
  VaR <- VaRS(k)
  1/(1 - k)*sum(pm*am/b * pgamma(VaR, am + 1, b, low = FALSE))
}

SL(100)
VaRS(u)
TVaRS(u)

#### Loi Poisson Teicher bivariée composée -- Bi~B~Exp(b) ####
k <- 0:50
lam <- c(3, 2)
a0 <- 1
b <- 0.2

### Fonction de densité de la Teicher
fm12 <- dPoTeibiv(k, k, lam, a0)

### Marignales
fm1 <- rowSums(fm12)
fm2 <- colSums(fm12)

### Espérance de Xi et S
Emi <- c(sum(k*fm1), sum(k*fm2))
Ex <- 1/b*Emi 
Es <- sum(Ex)

### Masse de probabilité pij
pm <- as.vector(fm12)

### Calcules paramètres loi gamma
ami <- expand.grid(k, k)
am <- rowSums(ami)

### Repartition de FS
Fs <- function(x) pm[1] + sum(pm[-1]*pgamma(x, am[-1], b))
Fs(100)

### Vérif : Espérance
Es_test <- sum(pm*am/b)
cbind(Es_test, Es)

#### Loi Poisson Teicher bivariée composée -- Bi~Geo(qi)) k \in N1 ####
nfft <- 2^13
k <- 0:(nfft - 1)
m1 <- 0:50
lam <- c(5, 7)
a0 <- 3
q <- c(0.2, 0.4)

### Fonction de masse de probabilité conjointe du vecteur M
fm12 <- dPoTeibiv(m1, m1, lam, a0)

### Fonction de masse des v.a.s Bi
fbi <- sapply(q, function(i) dgeom(k - 1, i))

### Marignales
fm1 <- rowSums(fm12)
fm2 <- colSums(fm12)

### Espérance de Xi et S
Emi <- lam
Ex <- 1/q*Emi 
Es <- sum(Ex)

### V.a. S
fbit <- mvfft(fbi)
fst <- exp((lam[1] - a0)*(fbit[, 1] - 1))*exp((lam[2] - a0)*(fbit[, 2] - 1))*exp(a0*(fbit[, 1]*fbit[, 2] - 1))
fs <- Re(fft(fst, TRUE))/nfft

### Repartition de FS
Fs <- cumsum(fs)

### Vérif : Espérance
Es_test <- sum(k*fs)
cbind(Es_test, Es)

#### Loi Bernoulli CGMR -- Bi~Ga(ai, b), i=1,...,20 ####
n <- 20
k <- 0:1
i <- 1:n
q <- 0.3+0.015*i
q0 <- 0.2
qt <- 1 - (1 - q)/(1 - q0)
a <- c(2.2, 2.4, 2.6, 2.3, 2.7,
       3.2, 2.1, 2.1, 2.1, 3.3, 
       2.6, 2.7, 3.0, 3.4, 3.2,
       2.9, 2.2, 3.1, 2.2, 3.1)
b <- 0.25

### Fonction de massse de probabilités conjointe
fm <- dBernCGMR(qt, q0)

### Espérance de Xi et S
Emi <- q
Ex <- a/b*Emi 
Es <- sum(Ex)

### Masse de probabilité pm
pm <- as.vector(fm)

### Calcules paramètres loi gamma
comb <- arrayInd(seq_along(fm), dim(fm)) - 1 
ami <- sweep(comb, 2, a, "*") # Clean
# ami <- comb*rep(a, each=nrow(comb)) Beaucoup plus rapide
am <- as.vector(comb %*% a) # ou rowSums(ami) 

### Repartition conjointe (Moins pratique en haute dimension)
Fx <- function(x) {
  gamma_prod <- sapply(seq_along(x), function(i) pgamma(x[i], ami[-1, i], b))
  fm[1] + sum(fm[-1] * apply(gamma_prod, 1, prod))
}

### Repartition de FS
Fs <- function(x) fm[1] + sum(fm[-1]*pgamma(x, am[-1], b))
Fs(100)

### Vérif : Espérance
Es_test <- sum(fm*am/b)
cbind(Es_test, Es)





#### Loi discrète bivariée composée, marginales Bernoulli -- Bi~Exp(bi) ####
fm12 <- matrix(c(0.6, 0.1,
          0.25, 0.05), ncol=2, byrow=TRUE)
k <- 0:1
b <- c(0.20, 0.05)

### Marignales
fm1 <- rowSums(fm12)
fm2 <- colSums(fm12)

### Espérance de Xi et S
Emi <- c(sum(k*fm1), sum(k*fm2))
Ex <- 1/b*Emi 
Es <- sum(Ex)

### Repartition de conjointe
Fx12 <- function(x1, x2) fm12[1, 1] + fm12[2, 1]*pexp(x1, b[1]) +
  fm12[1, 2]*pexp(x2, b[2]) + fm12[2, 2]*pexp(x1, b[1])*pexp(x2, b[2])

### Repartition de FS
Fs <- function(x) fm12[1, 1] + fm12[2, 1]*pexp(x, b[1]) +
  fm12[1, 2]*pexp(x, b[2]) + fm12[2, 2]*perlgen(x, b)
Fs(20)

### Vérif : Espérance
Es_test <- sum(fm*am/b)
cbind(Es_test, Es)

#### Loi Poisson Teicher bivariée composée -- Bi~Ga(ai, bi) ####
# Voir ACT-3000/Ressources/Convo_VA_Comp.pdf (GitHub)
nfft <- 2^10
k <- 0:(nfft - 1)
a <- c(2, 3)
b <- c(0.2, 0.3)
lam <- c(3, 2)
a0 <- 1
q <- b[1]/b[2]

### Esperance
Emi <- lam
Exi <- Emi*a/b
Es <- sum(Exi)

### Variable Ji
fj1t <- dnbinom(k - a[1], a[1], q) |> fft()
fj2t <- replace(numeric(nfft), a[2] + 1, 1) |> fft()
fj0t <- fj1t*fj2t

### Variable Wi
fw1t <- exp((lam[1] - a0)*(fj1t - 1))
fw2t <- exp((lam[2] - a0)*(fj2t - 1))
fw0t <- exp(a0*(fj0t - 1))

### Variable L
flt <- fw1t*fw2t*fw0t
fl <- Re(fft(flt, TRUE))/nfft

### Fonction de répartition de S
Fs <- function(x) fL[1] + sum(fL[-1]*pgamma(x, k[-1], b[2]))

### Verif : Esperance
cbind("Es_test"=sum(k/max(b)*fl), Es)








