###| exercice 2.5 TVaR

###| Trouver les VaR à c(0.01,0.5,0.99) (Jé : BIEN pour la VaR)
qbinom(c(0.01,0.5,0.99), 1000, 0.2)
qpois(c(0.01,0.5,0.99), lambda = 200)

###| Chercher les TVaR (moyenne au dessus de la VaR)
x <- c(1:100)
i <- dpois(x,lambda = 200)
sum(i*x)
min(i)
#test (Jé : Cette méthode peut être utiliser seulement avec la simulatiom)
mean(i>0) # mean(i[i > 0]) = TvaR_0(X) dans le cas de la simualation

###| Ou bien quelque chose comme ca (Jé : nope)
###| 
mean(qbinom(c(0.01,0.5,0.99),1000, 0.2))


###| Chercher les TVaR (moyenne au dessus de la VaR) (Version Jé)
x <- 1:1e4
fx <- dpois(x, lambda = 200)
sum(fx)
sum(fx * x)

# TVaR avec la fonction Stop-Loss selon plus simple
TVaRN <- function(k) qpois(k, 200) + sum(pmax(x - qpois(k, 200), 0) * fx)/(1 - k)
sapply(c(0.01, 0.5, 0.99), TVaRN)

TVaRM <- function(k) qbinom(k, 1e3, 0.2) + sum(pmax(x - qbinom(k, 1e3, 0.2), 0) * fx)/(1 - k)
sapply(c(0.01, 0.5, 0.99), TVaRM)

