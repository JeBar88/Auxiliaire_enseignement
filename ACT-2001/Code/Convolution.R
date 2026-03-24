### ACT-2001
## Exemple de code pour le produit de convolution
## Jérémie Barde

#### Exemple 1 -- M1~Po(3) et M2~Po(5) ####
lam <- c(1, 2)
k <- 0:20

### Fmp de Mi
fm1 <- dpois(k, lam[1])
fm2 <- dpois(k, lam[2])
# Vérif
apply(cbind(fm1, fm2), 2, sum)

### V.a. S
## Basé sur la formule théorique
fs <- sapply(k, function(k) sum(fm1[0:k + 1] * fm2[k:0 + 1]))

## Code alternatif (Dominique Chevalier)
sums <- outer(k, k, '+')
fm12 <- outer(fm2, fm1)
fs <- sapply(k, function(k) sum(fm12[sums == k]))

# Vérif
sum(fs)
cbind(sum(lam), sum(fs * k))
cbind(dpois(5, sum(lam)), fs[6])

#### Exemple 2 -- Domaine non aritmétique ####
### Domaine et fmp
m1 <- c(5, 10, 15, 20, 50) 
m2 <- c(5, 15, 23, 60) 

pm1 <- c(0.5, 0.20, 0.15, 0.1, 0.05) 
pm2 <- c(0.55, 0.15, 0.2, 0.1)

### Espérance
Em1 <- sum(m1*pm1)
Em2 <- sum(m2*pm2)
Es <- Em1 + Em2

### Construction de vecteur complet
k <- 0:150
fm1 <- replace(numeric(length(k)), m1 + 1, pm1)
fm2 <- replace(numeric(length(k)), m2 + 1, pm2)

### V.a. S 
fs <- sapply(k, function(k) sum(fm1[0:k + 1]*fm2[k:0 + 1]))

# Vérif
sum(fs)
cbind(Es, "E_test"=sum(fs*k))

#### Exemple 2 -- Xi - Po(i/25) i in {1, 100}, S = M1 +...+ M100 ####
n <- 100
i <- 1:n
lam <- i/25
k <- 0:300

### Espérance
Es <- sum(lam)

### Fmp de Mi
fm <- sapply(lam, function(i) dpois(k, i))
apply(fm, 2, sum)

### V.a. S : Méthode brute
fs <- fm[, 1]
for (i in 2:n) {
  fs <- sapply(k, function(k) sum(fs[0:k + 1] * fm[k:0 + 1, i]))
}

# Vérif
sum(fs)
cbind(Es, sum(fs * k))
cbind(dpois(210, Es), fs[211])











