### ACT-3000
## Calcule de contriubiton au mesure VaR et TVaR
## # Contribution v.a. discrète
## ## Exemple 1 : Deux risques
## ## Exemple 1.5 -- version d dimension
## ## Exemple 2 -- d risques
## # Contribution v.a. continue
## ## Loi Gamma, bêta identiques
## ## Loi mélange d'Erlang
## Jérémie Barde

#### Fonction utile ####
source(here::here("ACT-3000/Code/Fonctions_Utiles.R"))

#### Contribution v.a. discrète ####
#### Exemple 1 -- Deux risques ####
nfft <- 2^16
k <- 0:(nfft - 1)
n <- c(200, 300)
q <- c(0.15, 1/15)
u <- 0.99

fm1 <- dbinom(k, n[1], q[1])
fm2 <- dbinom(k, n[2], q[2])

## Esperance Mi et S
Emi <- n*q
ES <- sum(Emi)

## Calcule de la v.a. S
fst <- fft(fm1) * fft(fm2)
fs <- Re(fft(fst, TRUE))/nfft
Fs <- cumsum(fs)
# Verif
c(ES, sum(k*fs))

VaRS <- k[min(which(Fs >= u))]
TVaRS <- VaRS + 1/(1 - u)*sum(pmax(k - VaRS, 0) * fs)

## E[Mi|S = k]
gm1 <- k * fm1/Emi[1]
gm2 <- k * fm2/Emi[2]

fam1t <- fft(gm1) * fft(fm2) 
fam1 <- Re(fft(fam1t, TRUE))/nfft

fam2t <- fft(gm2) * fft(fm1)
fam2 <- Re(fft(fam2t, TRUE))/nfft

Eam1 <- Emi[1] * fam1/fs
Eam2 <- Emi[2] * fam2/fs

## Verifier les instabilitees numériques
# Eam1 + Eam2

## Contribution au mesure de risque VaR et TVaR
CVaRm1 <- Eam1[VaRS + 1]
CVaRm2 <- Eam2[VaRS + 1]
cbind(CVaRm1, CVaRm2, "VaRS_test"=CVaRm1 + CVaRm2, VaRS)

# On met une limite supérieure pour ne pas tenir compte des erreurs numériques
# Il peut y avoir des instabilitees numeriques
CTVaRm1 <- (sum(Eam1[k > VaRS & k < 500] * fs[k > VaRS & k < 500]) + CVaRm1 * (Fs[VaRS + 1] - u))/(1 - u) 
CTVaRm2 <- (sum(Eam2[k > VaRS & k <= 500] * fs[k > VaRS & k <= 500]) + CVaRm2 * (Fs[VaRS + 1] - u))/(1 - u)
cbind(CTVaRm1, CTVaRm2, "VaRS_test"=CTVaRm1 + CTVaRm2, TVaRS)

# Façon alternative: pas besoin de limite
Fam1 <- cumsum(fam1)
Fam2 <- cumsum(fam2)
CTVaRm1 <- (Emi[1]*(1 - Fam1[VaRS + 1]) + CVaRm1 * (Fs[VaRS + 1] - u))/(1 - u)
CTVaRm2 <- (Emi[2]*(1 - Fam2[VaRS + 1]) + CVaRm2 * (Fs[VaRS + 1] - u))/(1 - u) 
cbind(CTVaRm1, CTVaRm2, "TVaRS_test"=CTVaRm1 + CTVaRm2, TVaRS)

#### Exemple 1.5 -- version d dimension ####
nfft <- 2^16
k <- 0:(nfft - 1)
d <- 2
n <- c(200, 300)
q <- c(0.15, 1/15)
lam <- 45
u <- 0.9

fm1 <- dbinom(k, n[1], q[1])
fm2 <- dbinom(k, n[2], q[2])
fm3 <- dpois(k, lam)
fmi <- cbind(fm1, fm2)

## Esperance Mi et S
Emi <- c(n*q)
ES <- sum(Emi)

## Calcule de la v.a. S
fmit <- mvfft(fmi)
fst <- apply(fmit, 1, prod)
fs <- Re(fft(fst, TRUE))/nfft
Fs <- cumsum(fs)
# Verif
c(ES, sum(k*fs))

VaRS <- k[min(which(Fs >= u))]
TVaRS <- VaRS + 1/(1 - u)*sum(pmax(k - VaRS, 0) * fs)

## E[Mi|S = k]
gmi <- k*sweep(fmi, 2, Emi, "/") # ou t(t(fmi)/Emi)

# Note: l'argument drop est necessaire en dimension 2
famit <- sapply(1:d, function(i) fft(gmi[, i])*apply(fmit[, -i, drop=FALSE], 1, prod))
fami <- Re(mvfft(famit, TRUE))/nfft

Eami <- sweep(fami, 2, Emi, "*")/fs

## Verifier les instabilitees numériques
# apply(Eami, 1, sum)

## Contribution au mesure de risque VaR et TVaR
CVaRmi <- sapply(1:d, function(i) Eami[VaRS + 1, i])
c(CVaRmi, "VaRS_test"=sum(CVaRmi), "VaRS"=VaRS)

Fami <- apply(fami, 2, cumsum)
CTVaRmi <- sapply(1:d, function(i) 
  (Emi[i]*(1 - Fami[VaRS + 1, i]) + CVaRmi[i] * (Fs[VaRS + 1] - u))/(1 - u))
c(CTVaRmi, "TVaRS_test"=sum(CTVaRmi), "TVaRS"=TVaRS)

#### Exemple 2 -- d risques ####
nfft <- 2^10
k <- 0:(nfft - 1)
d <- 5
i <- 1:d
r <- i/2
q <- (i %% 5 + 1) / 6
u <- 0.99

fmi <- mapply(function(r, q) dnbinom(k, r, q), r, q)
apply(fmi, 2, sum)

## Esperance Mi et S
Emi <- r*(1-q)/q
ES <- sum(Emi)

## Calcule de la v.a. S
fmit <- mvfft(fmi)
fst <- apply(fmit, 1, prod)
fs <- Re(fft(fst, TRUE))/nfft
Fs <- cumsum(fs)
# Verif
c(ES, sum(k*fs))

VaRS <- k[min(which(Fs >= u))]
TVaRS <- VaRS + 1/(1 - u)*sum(pmax(k - VaRS, 0) * fs)

## E[Mi|S = k]
gmi <- k*sweep(fmi, 2, Emi, "/") # ou t(t(fmi)/Emi)

# Note: l'argument drop est necessaire en dimension 2
famit <- sapply(1:d, function(i) fft(gmi[, i])*apply(fmit[, -i, drop=FALSE], 1, prod))
fami <- Re(mvfft(famit, TRUE))/nfft

Eami <- sweep(fami, 2, Emi, "*")/fs

## Verifier les instabilitees numériques
# apply(Eami, 1, sum)

## Contribution au mesure de risque VaR et TVaR
CVaRmi <- sapply(1:d, function(i) Eami[VaRS + 1, i])
c(CVaRmi, "VaRS_test"=sum(CVaRmi), "VaRS"=VaRS)

Fami <- apply(fami, 2, cumsum)
CTVaRmi <- sapply(1:d, function(i) 
  (Emi[i]*(1 - Fami[VaRS + 1, i]) + CVaRmi[i] * (Fs[VaRS + 1] - u))/(1 - u))
c(CTVaRmi, "TVaRS_test"=sum(CTVaRmi), "TVaRS"=TVaRS)









#### ---------------------------------------------------------------------- ####
#### Contribution v.a. continue #### 
#### Loi Gamma, bêta identiques ####
a <- c(2.5, 4, 3)
b <- 0.1
u <- 0.99

### Espérance
Exi <- a/b
Es <- sum(Exi)

### VaR et TVaR de S
at <- sum(a)
VaRs <- qgamma(u, at, b)
TVaRs <- 1/(1 - u)*at/b*pgamma(VaRs, at + 1, b, low=FALSE)

### Contriubtion a la mesure VaR
CVaRxi <- sapply(1:3, function(i) a[i]/b*dgamma(VaRs, at + 1, b))/dgamma(VaRs, at, b)
cbind(sum(CVaRxi), VaRs)

### Contriubtion a la mesure TVaR
CTVaRxi <- sapply(1:3, function(i) 1/(1 - u)*a[i]/b*pgamma(VaRs, at + 1, b, low=FALSE))
cbind(sum(CTVaRxi), TVaRs)

#### Loi mélange d'Erlang, bêta identique ####
#### Ji~Bin(ni, qi), i=1, 2, 3
nfft <- 2^5
k <- 0:(nfft - 1)
n <- c(3, 5, 10)
q <- c(0.2, 0.2, 0.4)
b <- 0.2
u <- 0.99

### Espérance
Eji <- n*q
Ebi <- 1/b
Exi <- Eji*Ebi
Es <- sum(Exi)

### Calculer les poids p de la v.a. S
fji <- mapply(function(n, q) dbinom(k, n, q), n, q)
fjit <- mvfft(fji)
pjt <- apply(fjit, 1, prod)
pj <- Re(fft(pjt, TRUE))/nfft

### Fonction de répartition de S
Fs <- function(x) pj[1] + sum(pj[-1]*pgamma(x, k[-1], b))
Fs(50)

### VaR et TVaR de S
VaRs <- optimize(function(x) abs(Fs(x) - u), c(0, 200))$min
TVaRs <- 1/(1 - u)*sum(pj*k/b*pgamma(VaRs, k + 1, b, low=FALSE))

### Contriubtion a la mesure VaR
CVaR <- function(i){
  poidst <- fft(k*fji[, i])*apply(fjit[, -i], 1, prod)
  poids <- Re(fft(poidst, TRUE))/nfft
  sum(poids/b*dgamma(VaRs, k + 1, b))/sum(pj[-1]*dgamma(VaRs, k[-1], b))
}
CVaRxi <- sapply(1:3, CVaR)
cbind(sum(CVaRxi), VaRs)

### Contriubtion a la mesure TVaR
CTVaR <- function(i){
  poidst <- fft(k*fji[, i])*apply(fjit[, -i], 1, prod)
  poids <- Re(fft(poidst, TRUE))/nfft
  1/(1 - u)*sum(poids/b*pgamma(VaRs, k + 1, b, low=FALSE))
}
CTVaRxi <- sapply(1:3, CTVaR)
cbind(sum(CTVaRxi), TVaRs)

#### Loi Erlang généralisée, 2 risques ####
#### Ji~Bin(ni, qi), i=1, 2, 3
b <- c(0.1, 0.2)
u <- 0.99

### Espérance
Exi <- 1/b
Es <- sum(Exi)

### Fonction de répartition de S
Fs <- function(x) perlgen(x, b)
Fs(50)

### VaR et TVaR de S
VaRs <- optimize(function(x) abs(Fs(x) - u), c(0, 200))$min
TVaRs <- terlgen(VaRs, b)

### Contriubtion a la mesure VaR
CVaR <- function(i, j) {
  A <- VaRs*(exp(-b[i]*VaRs)/(b[j] - b[i]))
  B <- (exp(-b[i]*VaRs) - exp(-b[j]*VaRs))/(b[j] - b[i])^2
  prod(b)*(A - B)/derlgen(VaRs, b)
}
CVaRxi <- mapply(CVaR, c(1, 2), c(2, 1))
cbind(sum(CVaRxi), VaRs)

### Contriubtion a la mesure TVaR
CTVaR <- function(i, j){
  A <- (b[j]*exp(-b[i]*VaRs)*(VaRs + 1/b[1]))/(b[j] - b[i])
  B <- (exp(-b[i]*VaRs) - exp(-b[j]*VaRs))/(b[j] - b[i])^2
  1/(1 - u)*(A - B)
}
CTVaRxi <- mapply(CTVaR, c(1, 2), c(2, 1))
cbind(sum(CTVaRxi), TVaRs)

#### Loi Erlang généralisée, n risques ####
b <- c(0.1, 0.2, 2, 0.5)
u <- 0.99

### Espérance
Exi <- 1/b
Es <- sum(Exi)

### Fonction de répartition de S
Fs <- function(x) perlgen(x, b)
Fs(50)

### VaR et TVaR de S
VaRs <- optimize(function(x) abs(Fs(x) - u), c(0, 200))$min
TVaRs <- terlgen(VaRs, b)

### Contriubtion a la mesure VaR
CVaR <- function(i){
  bk <- b[-i]
  P <- sapply(seq_along(bk), function(i) prod(bk[-i]/(bk[-i] - bk[i])))
  A <- VaRs*(exp(-b[i]*VaRs))/(bk - b[i])
  B <- (exp(-b[i]*VaRs) - exp(-bk*VaRs))/(bk - b[i])^2
  sum(P*b[i]*bk*(A - B))/derlgen(VaRs, b)
}
CVaRxi <- sapply(seq_along(b), CVaR)
cbind(sum(CVaRxi), VaRs)

### Contriubtion a la mesure TVaR
CTVaR <- function(i){
  bk <- b[-i]
  P <- sapply(seq_along(bk), function(i) prod(bk[-i]/(bk[-i] - bk[i])))
  A <- (bk*exp(-b[i]*VaRs)*(VaRs + 1/b[i]))/(bk - b[i])
  B <- (bk*exp(-b[i]*VaRs) - b[i]*exp(-bk*VaRs))/(bk - b[i])^2
  1/(1 - u)*sum(P*(A - B))
}
CTVaRxi <- sapply(seq_along(b), CTVaR)
cbind(sum(CTVaRxi), TVaRs)

