### ACT-2001
## Chapitre 4 traditionelle : Question 28
## Jérémie Barde

a <- 2
b <- 1/1000
lam <- 0.25
t <- 0.5
k <- 0.99
vk <- 0:100

fm <- (lam*(1-t) * (lam*(1-t) + t*vk)^(vk-1))/factorial(vk) * exp(-(lam*(1-t) + t*vk))
sum(fm)
sum(fm * vk)

## i)
Esp <- lam * a/b
Var <- lam * a/(b^2) + lam/(1-t)^2 * (a/b)^2 
cbind(Esp, Var)

## ii)
v0 <- 1:100
Fx <- function(x) fm[1] + sum(fm[-1] * pgamma(x, v0*a, b))
Fx(2E5)

VaR <- optimize(function(x) abs(Fx(x) - k), c(0, 2E5))$min
TVaR <- 1/(1-k) * sum(fm[-1] * v0 * 2/b * pgamma(VaR, (v0*2) + 1, b, low=F))
cbind(VaR, TVaR)
