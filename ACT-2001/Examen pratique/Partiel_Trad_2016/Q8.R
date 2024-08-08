### ACT-2001 ###
## Examen partiel trad 2016 : Question 8
dbinom(0, 1, 0.2) * dbinom(0, 2, 0.3) + 
  sum(dbinom(c(1,0,1,0,1), 1, 0.2) * dbinom(c(0,1,1,2,2), 2, 0.3) * pgamma(400, c(3,2,5,4,7), 1/100)) 
