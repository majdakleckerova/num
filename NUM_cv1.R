# Gaussova eliminační metoda

n <- 4                              # počet dimenzí výsledného vektoru
A <- matrix(runif(n*n),n,n)         # tvorba matice o vel. n x n
b <- runif(n)                       # tvorba náhodného vektoru délky n
x <- solve(A,b)


Gauss <- function(A, b){
  # přímá metoda
  Ab <- cbind(A, b)               # tvorba rozšířené matice
  for(k in 1:(n - 1)){            # vynuluje prvky pod hlavní diagonálou
    for (i in (k + 1):n){
      c <- -Ab[i, k] / Ab[k, k]
      j <- (k + 1) : (n + 1)
      Ab[i,j] <- Ab[i, j] + c * Ab[k, j]
    }
  }
  
  # zpětná metoda
  x <- b
  x[n] <- Ab[n, n+1]/Ab[n, n]
  for(i in (n-1):1){
    j <- (i + 1):n
    x[i] <- (Ab[i, n + 1] - sum(Ab[i , j]*x[j]))/Ab[i, i]
  }
  return(x)
}


# GAUSSOVA METODA S PIVOTACÍ
GaussPivot <- function(A, b){
  Ab <- cbind(A, b)
  for(k in 1:(n - 1)){
    pivot <- which.max(abs(Ab[k:n, k])) + k - 1
    if(pivot != k){
      j <- k:(n + 1)
      pom <- Ab[k, j]
    }
    for (i in (k + 1):n){
      c <- -Ab[i, k] / Ab[k, k]
      j <- (k + 1) : (n + 1)
      Ab[i,j] <- Ab[i, j] + c * Ab[k, j]
    }
  }
  
  x <- b
  x[n] <- Ab[n, n+1]/Ab[n, n]
  for(i in (n-1):1){
    j <- (i + 1):n
    x[i] <- (Ab[i, n + 1] - sum(Ab[i , j]*x[j]))/Ab[i, i]
  }
  return(x)
}

print(x)
print(GaussPivot(A, b))
print(Gauss(A,b))



# LU Dekompozice
n <- 4                             
A <- matrix(runif(n*n),n,n) 
# funkce diag(A) ... vypíše diagonální prvky matice

n <- 4                              # počet dimenzí výsledného vektoru
A <- matrix(runif(n*n),n,n) 
LU_Dekompozice <- function(A){
  n <- nrow(A)             
  for(k in 1:(n - 1)){            
    for (i in (k + 1):n){
      c <- -A[i, k] / A[k, k]
      j <- (k+1)/n
      A[i, k] <- c                        # plnění dolního trojúhelníku
      A[i,j] <- A[i, j] + c * A[k, j]     # plnění horního trojúhelníku
    }
  }
  return(A)
}

LU <- LU_Dekompozice(A)
L <- LU
U <- LU
L[upper.tri(LU)] <- 0
diag(L) <- 0
U[lower.tri(LU)] <- 0

A
L %*% U



