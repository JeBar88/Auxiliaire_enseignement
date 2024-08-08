### ACT-2001
## Chapitre 2 informatique : Question 4 
## Jérémie Barde

k <- 0.99
# Hypothèse 1
fM1 <- dpois(0:100, 5)
FM1 <- ppois(0:100, 5)
VaR <- min(which(FM1 >= k)) - 1

# Définition de base
TVaR <- (sum(fM1[-(1:(VaR+1))] * (VaR+1):100) + VaR * (FM1[VaR + 1] - k))/(1-k)
# Avec la stop-loss
TVaRSL <- (sum(pmax(0:100 - VaR, 0) * fM1))/(1-k) + VaR
# Résultats
cbind(VaR, TVaR, TVaRSL)

