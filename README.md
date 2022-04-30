# The Q Test 

This is the github repo for the Q test proposed by Fan et al. (2022) that tests overidentifying restrictions (or instrumental variable validity) with high-dimensional data, which is robust to heteroskedastic errors. The method is called the Q test since it is based on estimation and inference for quadratic functionals of high-dimensional vectors.

## Required packages 

You can install the required package with 

``` r
install.pacakges(c("glmnet,MASS,igraph"))
```
In addition, users need to install the "Rmosek" package manually. The instructions for installation are available at  [Installation of MOSEK Rmosek package.](https://docs.mosek.com/latest/rmosek/install-interface.html)


"Lasso.R" and "Projection.R" contain necessary R functions for Lasso estimators and projection vectors for bias correction. 

"QTest.R" is the function to implement the Q test.

"example.R" provides an application example that is shown below.  


## Example

This is a basic example that shows you how to use the Q test. 

```{r example}
rm(list = ls())

library(glmnet)
library(Rmosek)
library(MASS)
library(igraph)

source("Lasso.R")
source("Projection.R")
source("QTest.R")

n = 300
px = 150
pz = 100 

p = px + pz
Sigma = diag(p)

phi <- c(seq(0.1,0.5,0.1),rep(0,145))
psi <- c(seq(0.3,0.7,0.1),rep(0,145))

s1 = 10 # num of relevant IV
gamma <- matrix(c(rep(0.5,10),rep(0,pz-10)),ncol = 1)

beta = 1 

## covariates and instruments 
W <- MASS::mvrnorm(n=n, rep(1, p), Sigma)
if (px == 0){
  X <- NULL
  Z <- W
}else{
  X <- W[,1:px]
  Z <- W[,(px+1):p]
}

## error terms 
err <- MASS::mvrnorm(n=n, rep(0, 2), matrix(1.5 * c(1, .5, .5,1),2))
e <- err[,1]
eps2 <- err[,2]
```


Test results with valid instruments 
```{r}
D <-  X%*%psi  + Z %*% gamma + eps2
Y <-   D*beta + X%*%phi + e

QTest(Y,D,Z,X) 
# $invalid
# [1] FALSE
# 
# $sig.level
# [1] 0.05
# 
# $pval
# [1] 0.9317884
```


Test results with invalid instruments 
```{r}
pi <-  c(seq(0.1,1,0.1),rep(0,90))
Y <- D*beta + X%*%phi + Z %*% pi + e

QTest(Y,D,Z,X) 
# $invalid
# [1] TRUE
# 
# $sig.level
# [1] 0.05
# 
# $pval
# [1] 3.160893e-05
```