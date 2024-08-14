### ACT-2001
## Convolution par la méthode brute
## Jérémie Barde

##### Exemple 1 #####
### X - Po(3) et Y - Po(5), S = X + Y
lam <- c(1, 2)
x <- 0:20

fx <- dpois(x, lam[1])
fy <- dpois(x, lam[2])
apply(cbind(fx, fy), 2, sum) # vérif

fs <- sapply(x, function(k) sum(fx[0:k + 1] * fy[k:0 + 1]))
# Vérif
sum(fs)
cbind(sum(lam), sum(fs * x))
cbind(dpois(5, sum(lam)), fs[6])

##### Exemple 2 #####
### domaine arbitraire
x <- c(5, 10, 15, 20, 50) # serait donné
y <- c(5, 15, 23, 60) # serait donné

px <- c(0.5, 0.20, 0.15, 0.1, 0.05) # serait donné
py <- c(0.55, 0.15, 0.2, 0.1) # serait donné

## Expérance X et Y
EX <- sum(px * x)
EY <- sum(py * y)
cbind(EX, EY)

## S = X + Y
s <- 0:110
fx <- numeric(length(s)) # vecteur de 0
fy <- fx

fx[x + 1] <- px # on remplie au bonne place
fy[y + 1] <- py

fs <- sapply(s, function(k) sum(fx[0:k + 1]*fy[k:0 + 1]))

# Vérif
sum(fs)
cbind(EX + EY, sum(fs*s))

##### Exemple 2 #####
### Xi - Po(i/25) i in {1, 100}, S = X1 + ... + X100
i <- 1:100
lam <- i/25
x <- 0:300

ES <- sum(lam)

fx <- sapply(lam, function(i) dpois(x, i))
apply(fx, 2, sum)

fs <- sapply(x, function(k) sum(fx[0:k + 1, 1]*fx[k:0 + 1, 2]))
for (i in 3:100) {
  fs <- sapply(x, function(k) sum(fs[0:k + 1]*fx[k:0 + 1, i]))
}

# Vérif
sum(fs)
cbind(ES, sum(fs * x))
cbind(dpois(210, ES), fs[211])












