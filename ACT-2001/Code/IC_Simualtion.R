### ACT-2001
## Intervalle de confiance pour les estimations par simulation
## Jérémie Barde
set.seed(1943)
m <- 1e6 # jouer avec m pour voir l'impacte sur les intervalles
X <- rgamma(m, 3, 0.1)
Y <- rpareto(m, 3, 60)
S <- X + Y

### Fonction de répartition empirique
Fn <- function(x) mean(S <= x)
Fn(100)

## Intervalle de confiance
ES_Fn <- sd(S <= 30)/sqrt(m)
Fn(100) + c(-1, 1) * qnorm(0.975)*ES_Fn

### Moyenne
Esp <- mean(S)
Esp

## Intervalle de confiance
ES_Esp <- sd(S)/sqrt(m)
Esp + c(-1, 1) * qnorm(0.975)*ES_Esp

### Stop-loss
d <- 20
SL <- mean(pmax(S - d, 0))
SL

## Intervalle de confiance
ES_SL <- sd(pmax(S - d, 0))/sqrt(m)
SL + c(-1, 1) * qnorm(0.975)*ES_SL

### Mesure VaR
u <- 0.9
VaR <- sort(S)[u*m]
VaR

## Intervalle de confiance
sdR <- sqrt(m*u*(1 - u))
k0 <- floor(qnorm(0.975)*sdR + 0.5)
VaR_IC <- sort(S)[(m*u) + c(-1, 1) * k0]
c(VaR, VaR_IC)

### Mesure TVaR
TVaR <- mean(S[S > VaR]) # ou VaR + mean(pmax(S - VaR, 0))/(1 - u)
TVaR

## Intervalle de confiance
TVaR_ES <- sqrt(var(pmax(S - VaR, 0))/((1 - u)^2*m))
TVaR_IC <- TVaR + c(-1, 1)*qnorm(0.975)*TVaR_ES
c(TVaR, TVaR_IC)

### Mesure entropique
p <- 0.001
entro <- 1/p*log(mean(exp(p*S)))
entro

## Intervalle de confiance
# On commence par faire l'intervalle de la fgm
fgm <- mean(exp(p*S))
ES_fgm <- sd(exp(p*S))/sqrt(m)
fgm_IC <- fgm + c(-1, 1)*qnorm(0.975)*ES_fgm

entro_IC <- 1/p*log(fgm_IC)
c(entro, entro_IC)





