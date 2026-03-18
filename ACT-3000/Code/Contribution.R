### ACT-3000
## Calcule de contriubiton au mesure VaR et TVaR
## Jérémie Barde

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
Eam1 + Eam2

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








