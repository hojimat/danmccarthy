rm(list = ls()) # clear cache
setwd("/home/ravshan/Dropbox/project/valuation/")
# input data
df <- read.csv("data/data.csv") #marketing costs
P <- read.csv("data/population.csv")
cpc <- read.csv("data/cpc.csv")$cpc[1] #customer acquisition cost

# public data for model
ADD <- df$marketing/cpc #quarterly customer acquisition number
LOSS <- df$customers - ADD #quarterly customer loss number
END <- df$customers #end of quarter customer number
REV <- END*df$revpc #quarterly revenue data
POP <- P$population #monthly pop data

# automated calculations
Qs <- nrow(df) #number of quarters
Ms <- Qs * 3
q_a <- 0 #data quarter
q_b <- Qs - nrow(na.omit(df)) + 1 #data starting quarter

#quarter dummies
num <- ceiling(Ms/9)
q1 <- rep(c(rep(1,3), rep(0,6), rep(0,3)), num)[1:Ms]
q2 <- rep(c(rep(0,3), rep(1,3), rep(0,6)), num)[1:Ms]
q3 <- rep(c(rep(0,6), rep(1,3), rep(0,3)), num)[1:Ms]



B <- function(m, mp, c, b1, b2, b3){
  ans <- 0
  if(mp > m){
    for(i in (m+1):mp){
      left <- (i-m)^c - (i-m-1)^c
      right <- b1 * q1[i] + b2 * q2[i] + b3 * q3[i]
      ans <- ans + left * exp(right)
    }
  }
  return(ans)
}

Sr <- function(m, mp, c, alpha, r, b1, b2, b3){
  ans <- (alpha / (alpha + B(m, mp-m, c, b1, b2, b3))) ^ r
  return(ans)
}

Fa <- function(m, mp, c, alpha, r, pina, b1, b2, b3){
  ans <- (1 - pina) * (1 - (alpha / (alpha + B(m, mp, c, b1, b2, b3)))^r)
  return(ans)
}



mega <- function(cA, alphaA, rA, pina, cR, alphaR, rR, b1A, b2A, b3A, b1R, b2R, b3R){
  L <- rep(0, Ms)
  M <- rep(0, Ms)
  A <- rep(0, Ms)
  C <- matrix(NA,Ms,Ms)
  
  L[1] <- 0
  M[1] <- POP[1]
  C[1,1] <- 1
  A[1] <- 1
  
  for(m in 2:Ms){
    #solve for M(L-1):
    M[m] <- POP[m] - POP[m-1] + L[m-1]
    #solve for A(M)
    for(i in 1:(m-1)){
      A[m] <- A[m] + M[i] * (Fa(i, m, cA, alphaA, rA, pina, b1A, b2A, b3A) - Fa(i, m-1, cA, alphaA, rA, pina, b1A, b2A, b3A))
    }
    #solve for C(A,m)
    for(i in 1:m){
      if(m == i){
        C[i,m] <- A[m]
      }else{
        C[i,m] <- A[m] * Sr(i, m, cR, alphaR, rR, b1R, b2R, b3R)
      }
    }

    #solve for L(C)
    L[m] <- sum(C[c(1:(m-1)),(m-1)]) - sum(C[c(1:m),m]) + C[m,m]
  }
  
  results <- list("M" = M,"A" = A,"C" = C, "L" = L)
  return(results)
}


SSE <- function(params){
  ans <- 0
  
  cA <- params[1]
  alphaA <- params[2]
  rA <- params[3]
  pina <- params[4]
  cR <- params[5]
  alphaR <- params[6]
  rR <- params[7]
  b1A <- params[8]
  b2A <- params[9]
  b3A <- params[10]
  b1R <- params[11]
  b2R <- params[12]
  b3R <- params[13]
  
  megane <- mega(cA, alphaA, rA, pina, cR, alphaR, rR, b1A, b2A, b3A, b1R, b2R, b3R)
  for(q in 1:Qs){
    ADDhat <- megane$A[3*q - 2] + megane$A[3*q - 1] + megane$A[3*q]
    LOSShat <- megane$L[3*q - 2] + megane$L[3*q - 1] + megane$L[3*q]
    ENDhat <- sum(megane$C[1:(3*q),3*q])
    
    ans <- ans + (ADD[q] - ADDhat[q])^2 + (LOSS[q] - LOSShat[q])^2
  }
  ans <- ans + (END[q] - ENDhat[q])^2
  return(ans)
}

x <- nlm(f = SSE, p=c(1,1,1,1,1,1,1,1,1,1,1,1,1))