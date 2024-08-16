### ACT-2001
## Code de base à maitriser pour v.a. discrète
## Jérémie Barde

##### Exemple typique #####
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










