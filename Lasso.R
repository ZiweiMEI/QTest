
diagXtX <- function(x, MARGIN = 1, ...) {
  if(MARGIN == 1) {
    # 1 indicates rows
    rowSums(x^2, ...)
  } else {
    # 2 indicates columns
    rowSums(t(x)^2, ...)
  }
}

cv.net <- function(x, y){
  Al <- seq(0,1,0.05)
  mse <- rep(NA,length(Al))
  foldid <- rep(1:10,ceiling(dim(z)[2]/10))[1:dim(z)[1]]
  num <- 0
  for (a in Al){
    num = num + 1
    fit <- cv.glmnet(x, y, foldid = foldid, alpha = a, intercept = F)
    mse[num] <- min(fit$cvm)
  }
  a <- Al[which.min(mse)]
  fit <- cv.glmnet(x, y, foldid = foldid, alpha = a, intercept = F)
  theta <- as.vector(coef(fit, s = fit$lambda.min))
  return(list("Alpha" = a, "htheta" = theta))
}


Lasso <- function(X, y, lambda = NULL, intercept = TRUE) {
  p <- ncol(X)
  n <- nrow(X)
  
  htheta <- if (is.null(lambda)) {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  } else if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "scalreg") {
    Xc <- if (intercept) {
      cbind(rep(1, n), X)
    } else {
      X
    }
    outLas <- scalreg(Xc, y)
    # return object
    if (intercept) {
      outLas$coefficients
    } else {
      # add a coefficient for the (not estimated) intercept b/c of implementation
      c(0, outLas$coefficients)
    }
  } else if (lambda == "net"){
    outLas <- cv.net(X, y)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(outLas$htheta)
  } else if (lambda == "sqrt"){
    y_numeric <- as.numeric(y)
    outLas <- sqrt_lasso(X, y_numeric)
    outLas
  }
  else {
    outLas <- glmnet(X, y, family = "gaussian", alpha = 1,
                     intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = lambda))
  }
  
  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}


Initialization.step <- function(X, y, lambda = NULL, intercept = FALSE) {
  n <- nrow(X)
  col.norm <- 1 / sqrt((1 / n) * diagXtX(X, MARGIN = 2))
  Xnor <- X %*% diag(col.norm)
  
  ### Call Lasso
  htheta <- Lasso(Xnor, y, lambda = lambda, intercept = intercept)
  
  ### Calculate return quantities
  if (intercept == TRUE) {
    Xb <- cbind(rep(1, n), Xnor)
    col.norm <- c(1, col.norm)
  } else {
    Xb <- Xnor
  }
  htheta <- htheta * col.norm
  returnList <- list("lasso.est" = htheta)
  return(returnList)
}


QJ <- function(y,x,z,var){
  
  # lambda = "scalreg"
  n = dim(z)[1]
  K = dim(z)[2]
  col.norm <- 1 / sqrt((1 / n) * diagXtX(z, MARGIN = 2))
  znor <- z %*% diag(col.norm)
  tau = 0.5 
  lam <- sqrt(2.01*log(K)/n) / 3
  Sigma_h = t(z)%*%z/n 
  Gamma_h <- Initialization.step(z, y, lambda = "CV", intercept = F)$lasso.est
  gamma_h <- Initialization.step(z, x, lambda = "CV", intercept = F)$lasso.est
  
  u1 <- Direction(z,loading = gamma_h)
  
  u2 <- Direction(z,loading = Gamma_h)
  
  u1.norm <- u1*sqrt(sum(gamma_h^2))
  u2.norm <- u2*sqrt(sum(Gamma_h^2))
  
  
  # u1.norm <- gamma_h 
  # u2.norm <- Gamma_h 
  
  
  inner_h <- sum(Gamma_h*gamma_h) + 1/n * t(u1.norm) %*% t(z) %*% (y-z%*%Gamma_h) +
    1/n * t(u2.norm) %*% t(z) %*% (x-z%*%gamma_h)
  
  gamma.normsq_h <- sum(gamma_h^2) + 2/n* t(u1.norm) %*% t(z) %*% (x-z%*%gamma_h) 
  Gamma.normsq_h <- sum(Gamma_h^2) + 2/n* t(u2.norm) %*% t(z) %*% (y-z%*%Gamma_h) 
  beta_h <- inner_h/gamma.normsq_h
  # psi_h <- Gamma_h - gamma_h %*% beta_h 
  Q_h <- Gamma.normsq_h + gamma.normsq_h%*%beta_h^2 - 2*inner_h%*%beta_h
  
  
  # u <- Direction(z,loading = psi_h)
  beta_h <- as.numeric(beta_h)
  if (var == "homo"){
    Theta11_h = 1/n*sum((y-znor%*%Gamma_h)^2)
    Theta22_h = 1/n*sum((x-znor%*%gamma_h)^2)
    Theta12_h = 1/n*sum((y-znor%*%Gamma_h)*(x-znor%*%gamma_h))
    V_h = 4/n*(Theta11_h+beta_h^2*Theta22_h-2*beta_h*Theta12_h)*
      t(u2.norm-u1.norm * beta_h)%*%Sigma_h%*%(u2.norm-u1.norm * beta_h) + tau/n
  }else if(var == "hetero"){
    e1_h <- as.numeric(y-z%*%Gamma_h)
    e2_h <- as.numeric(x-z%*%gamma_h)
    var_h <- e1_h^2 + beta_h^2 * e2_h^2 - 2 * beta_h * e1_h*e2_h
    V_h = 4/n^2*t(u2.norm-u1.norm * beta_h)%*%t(z*var_h)%*%
      z%*%(u2.norm-u1.norm * beta_h) + tau/n
  }
  p_val <- 1 - pnorm(Q_h/sqrt(V_h))
  return(c(p_val,beta_h))
}

Regularized <- function(y,x,z,Alpha = seq(0.01,0.5,0.01)){
  n <- dim(z)[1]
  K <- dim(z)[2]
  G <- dim(x)[2]
  mse0 = 1e6
  for (alpha in Alpha){
    P_ah = z%*%solve(t(z)%*%z+alpha*diag(K))%*%t(z) 
    temp <- (P_ah)/(1 - diag(P_ah))
    C_ah <- temp - diag(diag(temp)) 
    H_h <- t(x)%*%t(C_ah)%*%x 
    delta_h <- solve(H_h) %*% (t(x)%*%t(C_ah)%*%y)
    e_h <- y - x%*%delta_h 
    mrss <- sum(e_h^2)/n 
    H_t <- t(x)%*%C_ah
    u_h <- x - P_ah%*%x 
    s_ue <- (sum(u_h*e_h)/n)^2
    mse <- mrss*sum((x%*%rep(1,G)-C_ah%*%x%*%rep(1,G))^2)/n + s_ue*sum(diag(C_ah^2))/n
    if (mse < mse0){
      a = alpha
      mse0 = mse 
    }
  }
  alpha = a
  P_ah = z%*%solve(t(z)%*%z+alpha*diag(K))%*%t(z) 
  temp <- (P_ah)/(1 - diag(P_ah))
  C_ah <- temp - diag(diag(temp)) 
  H_h <- t(x)%*%t(C_ah)%*%x 
  delta_h <- solve(H_h) %*% (t(x)%*%t(C_ah)%*%y)
  e_h <- y - x%*%delta_h 
  V_h <- t(e_h)^2 %*% C_ah^2 %*% e_h^2 /sum(diag(P_ah))
  J_r <- t(e_h) %*% C_ah %*% e_h / sqrt(V_h) + sum(diag(P_ah))
  p_val <- 1 - pchisq(J_r,df=sum(diag(P_ah))-G)
  return(c(p_val,delta_h))
}

