### ACT-2001 ###
### Examen informatique 2018 h) 
## iii
b <- c(0.5, 1/6, 1/12)
FX <- function(x){
  prod(b[-1]/(b[-1] - b[1])) * pexp(x, b[1]) + 
    prod(b[-2]/(b[-2] - b[2])) * pexp(x, b[2]) +
    prod(b[-3]/(b[-3] - b[3])) * pexp(x, b[3])
}

FY <- function(y){
  G <- function(i) prod(b[-i]/(b[-i] - b[i])) * pexp(y, b[i])
  sum(sapply(1:length(b), G))
}

sapply(c(50, 80, 100), FX)
sapply(c(50, 80, 100), FY)

## iV
k <- 0.995
V <- optimize(function(x) abs(FY(x) - k), c(0, 1000))$min

# TVaR
H <- function(i) prod(b[-i]/(b[-i] - b[i])) * (V * exp(-b[i] * V) + exp(-b[i] * V)/b[i])
1/(1-k) * sum(sapply(1:length(b), H))

b <- c(0.1, 0.04)

FY(10)

# Stop-Loss
SL <- function(y){
  G <- function(i) prod(b[-i]/(b[-i] - b[i])) * 1/b[i] *  pexp(y, b[i], low = F)
  sum(sapply(1:length(b), G))
}

SL(50)
