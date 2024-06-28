#' Data generation
#'
#' This function generates data for the simulation study.
#'
#' @author Soumik Purkayastha, \email{soumikp@@umich.edu}
#'
#' @param x Input for generative model
#' @param type Which generative function should we use?
#'
#' @return A vector of output values from the generative exposure model
#' @examples
#' ## simulation from y = x*x
#' set.seed(1234)
#' x <- runif(1000)
#' y <- data_gen(x, 2)
#' @export

data_gen <- function(x, type = 1){
  if(type == 1){
    y <- sqrt(x)
  }else if(type == 2){
    y <- x^2
  }else if(type == 3){
    y <- (2*x - 1)/(x^2 - x + 5)
  }else if(type == 4){
    y <- log(1/(1 + exp(-x)))
  }else if(type == 5){
    y <- (sin(pi*(x - 0.5)) + exp(x))/(1 + exp(x))
  }else if(type == 6){
    y <- (x^2 + sqrt(x))/(1 + sqrt(1 - x^2))
  }else if(type == 7){
    y <- exp(x)
  }else if(type == 8){
    y <- sin(pi*(x - 0.5))
  }
  return(y)
}

#' Affine transformation to \[0,1]
#'
#' This function affine transforms data to lie in compact \[0,1].
#'
#' @author Soumik Purkayastha, \email{soumikp@@umich.edu}
#'
#' @param x Input for transformation
#'
#' @return A vector of output values from the affine transformation
#' @examples
#' ## simulation from y = x*x
#' set.seed(1234)
#' x <- runif(1000)
#' y <- data_gen(x, 2)
#' z <- affine_unif(y)
#' @export

affine_unif <- function(x){
  (x - min(x))/(max(x) - min(x))
}

#' Differential entropy calculation
#'
#' This function calculates differential entropy given density estimates
#'
#' @author Soumik Purkayastha, \email{soumikp@@umich.edu}
#'
#' @param x Input density estimates for calculation #'
#' @return Differential entropy estimate
#' @examples
#' set.seed(1234)
#' x <- runif(1000)
#' p <- density(x)
#' entropy_calc(p$y)
#' @export

entropy_calc <- function(x){
  x <- x[x>0]
  return(mean(-x*log(x)))
}

#' Conditional differential entropy calculation
#'
#' This function evaluates the difference of conditional differential entropy of \eqn{H(X \mid Z = z) - H(Y \mid Z = z)} in the generative exposure model \eqn{Y = g(X, Z)} over a grid of values of \eqn{Z}, where the grid is controlled by number of grid points = \code{f_eval*n}, \code{n} being the sample size of the data.
#'
#' @author Soumik Purkayastha, \email{soumikp@@umich.edu}
#'
#' @param x x serves as one input of the generative model
#' @param y y serves as the output of the generative model
#' @param z z serves as another input of the generative model, will serve as conditioning variable
#' @param f_eval Evaluation grid will contain \code{frac} fraction of the data points
#' @return Conditional entropy estimates \eqn{H(X \mid Y = y) - H(Z \mid Y = y)} for values of \code{y} in the evaluation grid
#' @examples
#' set.seed(1234)
#' x <- runif(1000)
#' y <- runif(1000)
#' z <- x^2 + y^2
#' cond_dens <- estim_conditional(x, y, z, f_eval = 0.25)
#' @export

estim_conditional <- function(x, y, z, f_eval = 0.25){

  library(hdrcde)

  obj_xz <- hdrcde::cde(x, z, nxmargin = round(length(x)*f_eval))
  obj_yz <- hdrcde::cde(y, z, nxmargin = round(length(x)*f_eval))

  tem <- cbind(obj_xz$x,
               apply(obj_xz$z, 1, entropy_calc) - apply(obj_yz$z, 1, entropy_calc))

  return(list(x = tem[,1],
              y = tem[,2]))
}

#' Nonparametric kernel smoothing (NKS)
#'
#' This function fits a local quadratic regression model \eqn{y = \theta(x)} to the data using a Gaussian kernel and evaluates equilibrium points of the fitted model \eqn{\hat{\theta}^\prime(x)} at the data points.
#'
#' @author Soumik Purkayastha, \email{soumikp@@umich.edu}
#'
#' @param x x serves as the input: x coordinate
#' @param y y serves as the response: y coordinate
#' @param h h serves as the bandwidth used for local fitting
#' @return values of fitted \eqn{\hat{\theta}(x)} and \eqn{\hat{\theta}^\prime(x)}
#' @examples
#' set.seed(1234)
#' x <- runif(1000)
#' y <- runif(1000)
#' z <- x^2 + y^2
#' cond_dens <- estim_conditional(x, y, z, f_eval = 0.25)
#' @export

fit_lq <- function(x, y, h){
  l <- length(x)
  beta0 <- rep(0, l)
  beta1 <- rep(0, l)
  for(i in 1:l) {
    x.reg <- x - x[i]
    w <- dnorm(x - x[i], 0, h)
    fit <- lm(y ~ x.reg + I(x.reg^ 2), weights = w)
    beta0[i] <- fit$coe[1]
    beta1[i] <- fit$coe[2]
  }
  beta <- cbind(beta0, beta1)
  return( beta)
}

#' Adding confidence intervals to Nonparametric kernel smoothing (NKS)
#'
#' This function fits a local quadratic regression model \eqn{y = \theta(x)} to the data using a Gaussian kernel and evaluates equilibrium points of the fitted model \eqn{\hat{\theta}^\prime(x)} at the data points. Finally, estimates and associated CI's are returned
#'
#' @author Soumik Purkayastha, \email{soumikp@@umich.edu}
#'
#' @param x x serves as the input: x coordinate
#' @param y y serves as the response: y coordinate
#' @param h h serves as the bandwidth used for local fitting
#' @param beta beta serves as the values of fitted \eqn{\hat{\theta}(x)} and \eqn{\hat{\theta}^\prime(x)}
#' @return 95% CI's associated with each beta.
#' @examples
#' set.seed(1234)
#' x <- runif(1000)
#' y <- runif(1000)
#' z <- x^2 + y^2
#' cond_dens <- estim_conditional(x, y, z, f_eval = 0.25)
#' @export
#'

add_ci <- function(x, y, h, beta){
  l <- length(x)
  beta0 <- beta[, 1]
  beta1 <- beta[, 2]
  ##Test for equilibruim points
  w <- rep( 0, l)
  diag <- rep( 0, l)

  upperCI_0 <- rep( 0, l)
  lowerCI_0 <- rep( 0, l)

  upperCI_1 <- rep( 0, l)
  lowerCI_1 <- rep( 0, l)

  se_0 <- rep( 0, l)
  se_1 <- rep( 0, l)

  z_0 <- rep( 0, l)
  z_1 <- rep( 0, l)

  p_0 <- rep( 0 ,l)
  p_1 <- rep( 0 ,l)
  options( object.size = 1000000000)
  ##Estimate sigma^2
  for( i in 1:l){
    Xi <- cbind( rep( 1, l), x - x[i], ( x - x[i]) ^ 2)
    ker <- dnorm( Xi[,2], 0, h)
    A <- matrix( 0, ncol = 3, nrow = 3)
    A[1,1] <- sum( ker)
    A[1,2] <- ker %*% Xi[,2]
    A[2,1] <- A[1,2]
    A[1,3] <- ker %*% Xi[,3]
    A[2,2] <- A[1,3]
    A[3,1] <- A[1,3]
    A[2,3] <- ker %*% Xi[,2] ^ 3
    A[3,2] <- A[2,3]
    A[3,3] <- ker %*% Xi[,3] ^ 2
    B <- solve( A)[1,]
    C <- rbind( ker, ker * Xi[,2], ker * Xi[,3])
    wi <- B %*% C
    diag[i] <- wi[i]
    w <- rbind( w, wi)
  }
  w <- w[ 2:( l + 1), ]
  second <- sum( w ^ 2)
  first <- 2 * sum( diag)
  v <- first - second
  vari <- 1 / ( l - v ) * sum( ( y - beta0) ^ 2)
  ##Calculate the 95% confidence band
  for( i in 1:l) {
    X <- cbind( rep( 1, l), x - x[i], ( x - x[i]) ^ 2)
    kernel <- dnorm( X[, 2], 0, h)
    An <- matrix( 0, ncol = 3, nrow = 3)
    Bn <- matrix( 0, ncol = 3, nrow = 3)
    An[1,1] <- sum( kernel) / l
    An[1,2] <- kernel %*% X[,2] / l
    An[2,1] <- An[1,2]
    An[1,3] <- kernel %*% X[,3] / l
    An[2,2] <- An[1,3]
    An[3,1] <- An[1,3]
    An[2,3] <- kernel %*% X[,2] ^ 3 / l
    An[3,2] <- An[2,3]
    An[3,3] <- kernel %*% X[,3] ^ 2 / l
    kernel2 <- kernel ^ 2
    Bn[1,1] <- sum( kernel2) / l / l
    Bn[1,2] <- kernel2 %*% X[,2] / l / l
    Bn[2,1] <- Bn[1,2]
    Bn[1,3] <- kernel2 %*% X[,3] / l / l
    Bn[2,2] <- Bn[1,3]
    Bn[3,1] <- Bn[1,3]
    Bn[2,3] <- kernel2 %*% X[,2] ^ 3 / l / l
    Bn[3,2] <- Bn[2,3]
    Bn[3,3] <- kernel2 %*% X[,3] ^ 2 / l / l
    sol <- solve( An)
    temp <- sol %*% Bn %*% sol

    temp1 <- temp[1,1]
    se_0[i] <- sqrt( vari * temp1)
    z_0[i] <- beta0[i] / se_0[i]
    p_0[i] <- (1 - pnorm(z_0[i]))
    upperCI_0[i] <- beta0[i] + 1.96 * se_0[i]
    lowerCI_0[i] <- beta0[i] - 1.96 * se_0[i]

    temp2 <- temp[2,2]
    se_1[i] <- sqrt( vari * temp2)
    z_1[i] <- abs( beta1[i] / se_1[i])
    p_1[i] <- ( 1 - pnorm( z_1[i])) * 2
    upperCI_1[i] <- beta1[i] + 1.96 * se_1[i]
    lowerCI_1[i] <- beta1[i] - 1.96 * se_1[i]


  }
  upperCI_0 <- round(upperCI_0, 5)
  lowerCI_0 <- round(lowerCI_0, 5)

  upperCI_1 <- round(upperCI_1, 5)
  lowerCI_1 <- round(lowerCI_1, 5)

  p_0 <- round(p_0, 5)
  p_1 <- round(p_1, 5)

  CIp_0 <- cbind( upperCI_0, lowerCI_0, p_0)
  CIp_1 <- cbind( upperCI_1, lowerCI_1, p_1)
  return(cbind(CIp_0, CIp_1))
}

#' Finding bandwidth for Nonparametric kernel smoothing (NKS)
#'
#' This function obtains optimal bandwidth to fit a local quadratic regression model \eqn{y = \theta(x)} to the data using a Gaussian kernel and evaluates equilibrium points of the fitted model \eqn{\hat{\theta}^\prime(x)} at the data points.
#'
#' @author Soumik Purkayastha, \email{soumikp@@umich.edu}
#'
#' @param x x serves as the input: x coordinate
#' @param y y serves as the response: y coordinate
#' @return bandwidth of kernel-based local quadratic fit
#' @examples
#' set.seed(1234)
#' x <- runif(1000)
#' y <- runif(1000)
#' z <- x^2 + y^2
#' cond_dens <- estim_conditional(x, y, z, f_eval = 0.25)
#' @export
#'

find_bw <- function(x, y){
  return(nprobust::lpbwselect(x = x, y = y,
                       eval = x[y == min(y)])$bws[2])
}




