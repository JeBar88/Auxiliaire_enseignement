### ACT-2001
## Chapitre 2 informatique : Question 11
## Jérémie Barde

library("actuar")
## Pour X
f <- function(par){
  abs(mpareto(1, par[1], par[2]) - 1000) + abs(qpareto(0.992, par[1], par[2]) - 12000)
}
# Avec constrOptim
constrOptim(c(2, 400), f, grad = NULL, ui = diag(2), ci = c(1, 0))$par

## Pour Y
# On trouve les deux en même temps
g <- function(par){
  abs(mlnorm(1, par[1], par[2]) - 1000) + abs(qlnorm(0.992, par[1], par[2]) - 12000)
}
par_ln <- optim(c(1, 1), g)$par

# On isole un paramètre en fonction d'un autre
mu_fnc <- function(s) log(12000) - s*qnorm(0.992)
f <- function(s){
  exp(mu_fnc(s) + s^2/2)
}

s2 <- optimise(function(x) abs(f(x) - 1000), c(0, 5))$min
mu2 <- mu_fnc(s)
cbind(mu, s)

# Première méthode meilleur
mlnorm(1, mu, s)
mlnorm(1, par_ln[1], par_ln[2])
