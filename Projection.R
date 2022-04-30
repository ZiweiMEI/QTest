
find_mu <- function(Sigma.hat, loading) {
  
  loading <- as.vector(loading)
  # loading.norm <- sqrt(sum(loading^2))
  px <- ncol(Sigma.hat)
  Id <- diag(1, px)
  zeros <- matrix(0,ncol=px,nrow=px)
  
  
  # constraint marix A
  A1 <- cbind(Sigma.hat, -Id, rep(0,px))
  A2 <- cbind(-Sigma.hat, -Id, rep(0,px))
  A3 <- cbind(zeros, Id, rep(-1, px))
  A <- rbind(A1, A2, A3) %>%
    Matrix(sparse=T)
  
  # lower bounds on variables
  blx <- c(rep(-Inf, px), rep(0, px), 0)
  bux <- c(rep(Inf, 2*px), Inf)
  
  objective <- c(rep(0,2*px), 1)
  
  # for ( j in 1:px ){
  blc <- rep(-Inf, 3*px)
  # buc1 <- rep(0, px); buc1[j] <- 1
  # buc2 <- rep(0, px); buc2[j] <- -1
  
  buc1 <- loading 
  buc2 <- -buc1
  buc3 <- rep(0, px)
  buc <- c(buc1, buc2, buc3)
  
  mosek.prob <- list(
    sense="min",
    c=objective,
    A=A,
    bc=rbind(blc, buc),
    bx=rbind(blx, bux)
  )
  mosek.res <- try(mosek(mosek.prob, list(verbose=0, soldetail=1)), silent=TRUE)
  mu <- mosek.res$sol$bas$pobjval
  mu
  
}




.u.hat_CLIME <- function(Sigma.hat, loading, mu) {
  
  # if (Mosek){
  #  loading.norm <- sqrt(sum(loading^2))
  px <- ncol(Sigma.hat)
  Id <- diag(1, px)
  zeros <- matrix(ncol=px,nrow=px); zeros[,] <- 0
  Theta.hat <- matrix(ncol=px,nrow=px)
  
  # constraint matrix
  A1 <- cbind(Id, -1*Id)
  A2 <- cbind(-1*Id, -1*Id)
  A3 <- cbind(Sigma.hat, zeros)
  A4 <- cbind(-1*Sigma.hat, zeros)
  A <- rbind(A1, A2, A3, A4) %>%
    Matrix(sparse=T)
  
  # lower and upper bounds on variables
  blx <- c(rep(-Inf, px), rep(0, px))
  bux <- rep(Inf, 2*px)
  
  objective <- c(rep(0,px), rep(1,px))
  
  # lower and upper bounds on constraints
  blc <- rep(-Inf, 4*px)
  buc3 <- mu + loading
  buc4 <- mu - loading
  buc <- c(rep(0, 2*px), buc3, buc4)
  
  mosek.prob <- list(
    sense="min",
    c=objective,
    A=A,
    bc=rbind(blc, buc),
    bx=rbind(blx, bux)
  )
  mosek.res <- try(mosek(mosek.prob, list(verbose=0)), silent=TRUE)
  u.hat <- mosek.res$sol$bas$xx[1:px]
}