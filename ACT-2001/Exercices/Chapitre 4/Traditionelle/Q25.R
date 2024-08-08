### ACT-2001
## Chapitre 4 traditionelle : Question 25
## Jérémie Barde

vk <- 1:1000 # Prendre uns valeur assez élevé comme en réalité ce devrait être l'infini
lam <- 0.2
b <- 1/1000
k <- 0.99

## b)
Fx <- function(x) dpois(0, lam) + sum(dpois(vk, lam) * pgamma(x, vk, b))
1 - Fx(3000)
# Ou directement avec la fonction de survie
Sx <- function(x) sum(dpois(vk, lam) * pgamma(x, vk, b, low = F))
Sx(3000)

## c)
VaR <- optimize(function(x) abs(Fx(x) - k), c(0, 10000))$min
VaR

TVaR <- sum(dpois(vk, lam) * 1/(1-k) * vk/b * pgamma(VaR, vk + 1, b, low=F))
TVaR