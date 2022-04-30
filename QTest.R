# Packages required: 
#   glmnet, MASS, igraph 
#   Rmosek (Mannualy added on. An instruction for its installation: 
#             https://docs.mosek.com/latest/rmosek/install-interface.html)
# 
# Inputs: 
#    Y: n by 1 vector of the outcome variable.
#    D: n by 1 vector of the endogenous variable.
#    Z: n by pz matrix of instrumental variables. 
#    X: n by px matrix of covariates. 
#    
#    eta: n by 1 vector for calibration. Follow the setting in Fan et al. (2022) by default.
#         Valid choices can be any vectors with $\lfloor n/2 \rfloor$ -1's and $ n \lfloor n/2 \rfloor $ 1's. 
#         n i.i.d. standard normal variables are also feasible. See Fan et al. (2022) for more details. 
#    tau0: The user defined constant for calibration level.
#          Valid only if tau = NULL. 
#    tau: The scalar of calibration. Follow the setting in Fan et al. (2022) by default. 
#    kappa: A user defined tuning parameter for projection vector constructions. Follow the setting in Fan et al. (2022) by default. 
#    sig.level: Significance level. 
#    intercept: Logical. FALSE for a model without intercept 
#    seed: random seed. 
#    
# Outputs: 
#   invalid: Logical. TRUE for rejecting the validity of instruments.  
#   sig.level: Significance level. 
#   pval: P-value of the test. 

QTest <- function(Y, D, Z, X = NULL, eta = NULL, 
                  tau0 = 1, tau = NULL, kappa = 1.2, 
                  sig.level = 0.05,
                  intercept = F, 
                  seed = 2022){
  
  if (intercept){
    if (is.null(X)){
      X = matrix(1,n,1)
    }else{
    X = cbind(1,X)
    }
    intercept = F 
  }
  
  if (is.null(X)){
    px = 0
  }else{
    px = dim(X)[2]
  }
  
  
  W = cbind(X,Z)
  
  n = dim(Z)[1]
  pz = dim(Z)[2]
  p = px + pz 
  
  
  Sigma.hat = t(W)%*%W/n
  
  Sigma.hat <- 1/2 * (Sigma.hat + t(Sigma.hat))
  
  
  set.seed(seed)
  
  ## generate eta following Fan et al. (2022)
  if (is.null(eta)){
    
    eta <- rep(-1,n)
    eta.temp <- eta
    obj <- Inf
    BB = 5e3
    
    for (bb in 1:BB){
      ind.temp <- sample(1:n, size = ceiling(n/2))
      eta.temp[ind.temp] <- 1
      obj.new <- max(abs( colSums(W * eta.temp ) ))
      if (obj.new < obj){
        ind.keep <- ind.temp
        obj <- obj.new
      }
      eta.temp <- eta
    }
    eta[ind.keep] <- 1
  }
  
  ############# Auxilliary parameter estimation 
  
  ## estimate \gamma 
  out.B <- Initialization.step(W, D, lambda = "CV", 
                               intercept = intercept)
  B.hat <- matrix(out.B$lasso.est, ncol = 1)
  e2_h = D - W %*% B.hat 
  gamma.hat <- B.hat[(px+1):p,]
  
  ## \hat{u}_2 
  if (sum(gamma.hat^2) == 0){u2 <- matrix(0,p,1)}else{
    if (p > 0.5*n & p < 1.5*n){
      sub <- sample(1:n,floor(n/2),replace=F)
      W.sub <- W[sub,]
      D.sub <- D[sub,]
      Sigma.hat.sub <- t(W.sub) %*%  (W.sub) / floor(n/2)
      out.B.sub <- Initialization.step(W.sub, D.sub, lambda =  "CV", 
                                       intercept = intercept)
      B.hat.sub <- matrix(out.B.sub$lasso.est, ncol = 1)
      
      gamma.hat.sub <- B.hat.sub[(px+1):p,]
      mu <- find_mu(Sigma.hat.sub, rbind(matrix(0,px,1),matrix(gamma.hat.sub,pz,1))/ sqrt(sum(gamma.hat.sub^2))) * kappa
      mu <- max(mu,1e-12) / sqrt(2)
    }
    else{
      mu <- find_mu(Sigma.hat, rbind(matrix(0,px,1),matrix(gamma.hat,pz,1))/ sqrt(sum(gamma.hat^2))) * kappa
      mu <- max(mu, 1e-12)
    }
    u2 <- .u.hat_CLIME(Sigma.hat, rbind(matrix(0,px,1),matrix(gamma.hat,pz,1))/ sqrt(sum(gamma.hat^2)), mu)
    u2 <- matrix(u2,ncol=1) * sqrt(sum(gamma.hat^2))
  }
  
  ## estimate \|\gamma\|_2^2
  gamma.norm.sq.hat <- sum(gamma.hat^2) + 2/n * sum((W %*% u2) * (D-W%*%B.hat))
  
  
  ## estimate \Gamma 
  out.A <- Initialization.step(W, Y, lambda = "CV", 
                               intercept = intercept)
  
  A.hat <- matrix(out.A$lasso.est, ncol = 1)
  e1_h = Y - W %*% A.hat 
  Gamma.hat <- A.hat[(px+1):p,]
  
  # \hat{u}_1 
  if (sum(Gamma.hat^2) == 0){u1 <- matrix(0,p,1)}else{
    if (p > 0.5*n & p < 1.5*n){
      sub <- sample(1:n,floor(n/2),replace=F)
      W.sub <- W[sub,]
      Y.sub <- Y[sub,]
      Sigma.hat.sub <- t(W.sub) %*%  (W.sub) / floor(n/2)
      
      out.A.sub <- Initialization.step(W.sub, Y.sub, lambda = "CV", intercept = intercept)
      A.hat.sub <- matrix(out.A.sub$lasso.est, ncol = 1)
      
      Gamma.hat.sub <- A.hat.sub[(px+1):p,]
      mu <- find_mu(Sigma.hat.sub, rbind(matrix(0,px,1),matrix(Gamma.hat.sub,pz,1))/ sqrt(sum(Gamma.hat.sub^2))) * kappa
      mu <- max(mu,1e-12) / sqrt(2)
    }
    else{
      mu <- find_mu(Sigma.hat, rbind(matrix(0,px,1),matrix(Gamma.hat,pz,1))/ sqrt(sum(Gamma.hat^2))) * kappa
      mu <- max(mu,1e-12)
    }
    u1 <- .u.hat_CLIME(Sigma.hat, rbind(matrix(0,px,1),matrix(Gamma.hat,pz,1))/ sqrt(sum(Gamma.hat^2)), mu)
    u1 <- matrix(u1,ncol=1) * sqrt(sum(Gamma.hat^2))
    
  }
  
  ## estimate <\Gamma,\gamma>
  GammaTgamma.hat <- sum(Gamma.hat*gamma.hat) + 1/n * sum((W %*% u1) * (D-W%*%B.hat)) +
    1/n * sum((W %*% u2) * (Y-W%*%A.hat))
  
  ## estimate \beta 
  if (gamma.norm.sq.hat <= 0){
    beta.hat = 0
  }else{beta.hat <- GammaTgamma.hat / gamma.norm.sq.hat}
  
  
  ######################## Estimate Q = \|\widetilde{\pi}\|_2^2 
  
  zeta.hat <- Y - beta.hat * D
  
  out.C <- Initialization.step(W, zeta.hat, lambda = "CV",
                               intercept = intercept)
  C.hat <- matrix(out.C$lasso.est, ncol = 1)
  pi.hat <- C.hat[(px+1):p,]
  
  ## \hat{u}_3
  if (sum(pi.hat^2) == 0){u3 <- matrix(0,p,1)}else{
    if (p > 0.5*n & p < 1.5*n){
      sub <- sample(1:n,floor(n/2),replace=F)
      W.sub <- W[sub,]
      zeta.hat.sub <-  zeta.hat[sub,]
      Sigma.hat.sub <- t(W.sub) %*%  (W.sub) / floor(n/2)
      out.C.sub <- Initialization.step(W.sub, zeta.hat.sub, lambda =  "CV", 
                                       intercept = intercept)
      C.hat.sub <- matrix(out.C.sub$lasso.est, ncol = 1)
      
      pi.hat.sub <- C.hat.sub[(px+1):p,]
      mu <- find_mu(Sigma.hat.sub, rbind(matrix(0,px,1),matrix(pi.hat.sub,pz,1))/ sqrt(sum(pi.hat.sub^2))) * kappa
      mu <- max(mu,1e-12) / sqrt(2)
    }
    else{
      mu <- find_mu(Sigma.hat, rbind(matrix(0,px,1),matrix(pi.hat,pz,1)) / sqrt(sum(pi.hat^2))) * kappa
      mu <- max(mu,1e-12)
    }
    u3 <- .u.hat_CLIME(Sigma.hat, rbind(matrix(0,px,1),matrix(pi.hat,pz,1))/ sqrt(sum(pi.hat^2)), mu)
    u3 <- matrix(u3,ncol=1) * sqrt(sum(pi.hat^2))
  }
  
  ## \hat{Q}_0
  pi.norm.sq.hat <- sum(pi.hat^2) + 2/n * sum((W %*% u3) * (zeta.hat-W%*%C.hat))
  
  e_h <- as.numeric(zeta.hat - W %*% C.hat)
  
  # define \tau 
  if (is.null(tau)){
    tt = 1 / (1 + log(log(n*p)) *   sqrt(n) * max(pi.norm.sq.hat, 0)  ) 
    tau = as.numeric(tau0 * tt)
  }
  
  # calibrated estimator of Q 
  pi.norm.sq.hat.cal <- pi.norm.sq.hat +  2 * sqrt(tau)/n * sum((eta) * (zeta.hat-W%*%C.hat)) 
  
  # estimated variance 
  V_h <- 4/n * sum((  sqrt(tau) * eta + as.numeric(W %*% u3) )^2*(e_h)^2)
  sd <- sqrt( V_h ) /sqrt(n)
  
  # the Q test 
  invalid <- (pi.norm.sq.hat.cal > qnorm(1-sig.level) * sd)
  pval = 1 - pnorm(pi.norm.sq.hat.cal / sd) 
  
  return(list(invalid = invalid, sig.level = sig.level, pval = pval))
  
}