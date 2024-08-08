### ACT-2001
## Chapitre 4 informatique : Question 1
## Jérémie Barde

r <- 5
q <- 0.2
a <- 2
b <- 0.0002
vk <- 1:1E4

Fx <- function(x) dnbinom(0, r, q) + sum(dnbinom(vk, r, q) * pgamma(x, a*vk, b))
Fx(6E5)
Fx(5E5)

VaR <- optimize(function(x) abs(Fx(x) - 0.99), c(5E5, 6E5))$min

TVaR <- sum(dnbinom(vk, r, q) * a*vk/b * 1/(1-0.99) * pgamma(VaR, vk*a + 1, b, low = F))
cbind(VaR, TVaR)

# e)
EspM <- r * (1-q)/q
VarM <- r * (1-q)/q^2
EspB <- a/b
VarB <- a/b^2

EspX <- EspM  * EspB
VarX <- EspM * VarB + VarM * EspB^2

EVXM <- EspM * VarB
VEXM <- EspB^2 * VarM

H1 <- EspX
H2 <- EVXM/VarX * (TVaR - EspX)
H3 <- VEXM/VarX * (TVaR - EspX)

cbind(H1, H2, H3)
sum(cbind(H1, H2, H3)) == TVaR














