# Part of the varbvs package, https://github.com/pcarbo/varbvs
#
# Copyright (C) 2012-2017, Peter Carbonetto
#
# This program is free software: you can redistribute it under the
# terms of the GNU General Public License; either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANY; without even the implied warranty of
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# Implements the fully-factorized variational approximation for
# Bayesian variable selection in linear regression. It finds the
# "best" fully-factorized variational approximation to the posterior
# distribution of the coefficients for a linear regression model of a
# continuous outcome (quantitiative trait), with spike and slab priors
# on the coefficients. By "best", we mean the approximating
# distribution that locally minimizes the K-L divergence between the
# approximating distribution and the exact posterior.
#
# Input X is an n x p matrix of observations of the variables (or
# features), where n is the number of samples, and p is the number of
# variables. Input y contains samples of the outcome; it is a vector
# of length n.
#
# Inputs sigma, sa and logodds are additional model parameters; sigma
# and sa are scalars. Input sigma specifies the variance of the
# residual, and sa specifies the prior variance of the coefficients
# (scaled by sigma). Input logodds is the prior log-odds of inclusion
# for each variable. Note that the prior log-odds here is defined with
# respect to the *natural* logarithm, whereas in function varbvs the
# prior log-odds is defined with respect to the base-10 logarithm, so
# a scaling factor of log(10) is needed to convert from the latter to
# the former.
#
# Output logw is the variational estimate of the marginal
# log-likelihood given the hyperparameters at each iteration of the
# co-ordinate ascent optimization procedure. Output err is the maximum
# difference between the approximate posterior probabilities (alpha)
# at successive iterations. Outputs alpha, mu and s are the
# parameters of the variational approximation and, equivalently,
# variational estimates of posterior quantites: under the variational
# approximation, the ith regression coefficient is normal with
# probability alpha(i); mu[i] and s(i) are the mean and variance of
# the coefficient given that it is included in the model.
#
# When update.sa = TRUE, there is the additional option of computing
# the maximum a posteriori (MAP) estimate of the prior variance
# parameter (sa), in which sa is drawn from a scaled inverse
# chi-square distribution with scale sa0 and degrees of freedom n0.
# library(BayesLogit)
tgen_probitnorm <- function (X, y, sigma, sa, logodds, alpha, mu, update.order,
                        tol = 1e-4, maxiter = 1e4, verbose = TRUE,
                        outer.iter = NULL, update.sigma = TRUE,
                        update.sa = TRUE, n0 = 10, sa0 = 1, annot,# @@ annot added by W.L. 12/26/2017
                        a0,b0){ # @@ a0 and b0 added by W.L. 1/5/2018

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)
  a_annot <- ncol(annot)
  #cat("!!!!!!!!!!!!!!!!!n:",n,"!!!!!",p,"\n")
  # (1) INITIAL STEPS
  # -----------------
  # Compute a few useful quantities. 
  xy <- c(y %*% X)
  d  <- diagsq(X)
  Xr <- c(X %*% (alpha*mu))
  
  # Calculate the variance of the coefficients.
  s <- sa*sigma/(sa*d + 1)

  # Initialize storage for outputs logw and err.
  logw <- rep(0,maxiter)
  err  <- rep(0,maxiter)
  # logodds_arra <- rep(0,maxiter)
  # eta_arra <- list()
  
  # Initialize parameters for logistic regression in annotation layer @@ added by W.L. 1/5/2018
  aN <- a0
  bN <- b0
  #cauchy <- rep(1,p)
  # initialize the expectation of annotation coefficients for probit layer by W.L. 09/08/2019
  e_w <- rep(0,a_annot)
  t_a_a <- t(annot) %*% annot
  logodds_original <- logodds
  
  # (2) MAIN LOOP
  # -------------
  # Repeat until convergence criterion is met, or until the maximum
  # number of iterations is reached.
  for (iter in 1:maxiter) {

    # Save the current variational and model parameters.
    alpha0 <- alpha
    mu0    <- mu
    s0     <- s
    sigma0 <- sigma
    sa.old <- sa
    
    # (2a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood.
   # if(iter == 1)
    logw0 <- int.linear(Xr,d,y,sigma,alpha,mu,s) +
             int.gamma(logodds,alpha) +
             int.klbeta(alpha,mu,s,sigma*sa) 
     cat(iter,"\t",logw0,"\n")
   #if(iter > 1){
   #   logw0 <- logw0 + 0.5 * t(WN)%*%VN.inverse%*%WN +
   #             0.5 *log(det(VN)) + part.LB(cauchy) -
   #             log(gamma(a0)) + a0*log(b0) -
   #             b0*aN/bN - aN*log(bN) +
   #             log(gamma(aN)) + aN
   #    logw0 <- logw[iter-1]
  
   #}
             

    # (2b) UPDATE VARIATIONAL APPROXIMATION
    # -------------------------------------
    # Run a forward or backward pass of the coordinate ascent updates.
    if (iter %% 2)
      i <- update.order
    else
      i <- rev(update.order)
    out   <- tgen_probitnormupdate(X,sigma,sa,logodds,xy,d,alpha,mu,Xr,i)
    alpha <- out$alpha
    mu    <- out$mu
    Xr    <- out$Xr
    rm(out)

    

    # (2d) UPDATE RESIDUAL VARIANCE
    # -----------------------------
    # Compute the maximum likelihood estimate of sigma, if requested.
    # Note that we must also recalculate the variance of the regression
    # coefficients when this parameter is updated.
    if (update.sigma) {
      # sigma <- (norm2(y - Xr)^2 + dot(d,betavar(alpha,mu,s)) +
                # dot(alpha,(s + mu^2)/sa))/(n + sum(alpha))
      # modified by W.L. 1/7/2018
      sa <- (sa0*n0 + dot(alpha,s + mu^2)/sigma^2)/(n0+2+sum(alpha))
      s     <- sa*sigma/(sa*d + 1)
    }
    
    # (2e) UPDATE PRIOR VARIANCE OF REGRESSION COEFFICIENTS
    # -----------------------------------------------------
    # Compute the maximum a posteriori estimate of sa, if requested.
    # Note that we must also recalculate the variance of the
    # regression coefficients when this parameter is updated.
    if (update.sa) {
      sa <- (sa0*n0 + dot(alpha,s + mu^2))/(n0 + sigma*sum(alpha))
      s  <- sa*sigma/(sa*d + 1)
    }
    
#     # added by W.L. 12/29/2017
#     # (2f) UPDATE LOGODDS
#     # ------------------------
#     # Update LOGODDS using the Bayeslogit package, first deal with updated
#     # alphas, get gammas from alphas
#     cat("when updating LOGODDS, dim of annot,",nrow(annot),"\t",head(annot),"\n")
     gammas <- as.numeric(alpha >= .5)
#     predictor <- as.array(annot)
#     predictor.intcpt <- cbind(predictor,rep(1,nrow(annot)))
#     cat("Before Logit, dimension of predictor.intcpt, head of predictor.intcpt",dim(predictor.intcpt),"\t",head(predictor.intcpt),"\n")
#     model.logit <- BayesLogit::logit(y=gammas,X=predictor.intcpt,samp=1,burn=100)
#     tmp.log <- predictor.intcpt%*%t(model.logit$b)
#     logodds <- exp(tmp.log)/(1+exp(tmp.log))
    
    # added by W.L. 1/5/2018
    # (2f) UPDATE \omega for annotation layer
    # ---------------------------------------  
    #gammas <- as.numeric(alpha >= .5)
    ## expectation of inversed variance of regression coefficients
    #E.eta = as.numeric(aN/bN)
    ## eta_arra[[iter]] <- E.eta
    #lambd.cauchys = c()
    ## update variance matrix for regression coefficients
    #for(ca.i in 1:p){
    #  lambd.cauchys = c(lambd.cauchys,lamd(cauchy[ca.i]))
    #}
    
    cat("aN ",aN," bN ",bN,"!!!\n")
    ## cat("for E.eta, ",E.eta,"\n")
    
    #VN.inverse  = inv.V(E.eta,lambd.cauchys,annot)
    #VN = solve(VN.inverse)
    ## update the expectation vector of regression coefficients
    #WN = Expect.W(VN,gammas,annot)
    
    ## update parameters for variance of omega
    #aN = a0 + ncol(annot)/2
    #bN = b0 + 0.5 * (crossprod(WN) + tr(VN.inverse))
   
    
    ## update cauchies 
    #for(ca.i in 1:p){
    #  tmp = t(annot[ca.i,]) %*% VN %*% annot[ca.i,] + crossprod(t(t(annot[ca.i,]) %*% VN))
    #  cauchy[ca.i] <- sqrt(tmp)
    #}

    #(2f) UPDATE coefficients w for annotation layer
    mu_annot <- annot %*% e_w
    # get the expectation of the hidden factor z_annot
    mu_z <- rep(0,p)
    mu_z[gammas==1] <- mu_annot[gammas==1] + dnorm(-mu_annot[gammas==1])/(1-pnorm(-mu_annot[gammas==1]))
    mu_z[gammas==0] <- mu_annot[gammas==0] - dnorm(-mu_annot[gammas==0])/pnorm(-mu_annot[gammas==0])
    # update the expectation of coefficients
    cat("a_annot ",a_annot,"\n")
    var_w <- chol2inv(chol((c(aN/bN) * diag(1,a_annot) + t_a_a)))
    mu_w <- var_w %*% t(annot) %*% mu_z
    # update the aN and bN
    aN <- a0 + a_annot/2
    bN <- b0 + 1/2*(t(mu_w) %*% mu_w + sum(diag(var_w)))
    WN <- mu_w
    
    # (2g) uPDATE LOGODDS
    logodds = annot %*% WN + logodds_original
    # logodds_arra[iter] <- logodds
    m = dim(annot)[2]
    
    # (2c) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood.
    logw[iter] <- int.linear(Xr,d,y,sigma,alpha,mu,s) +
      int.gamma(logodds,alpha) +
      int.klbeta(alpha,mu,s,sigma*sa) #+
    #  0.5 * t(WN)%*%VN.inverse%*%WN +
    #  0.5 *log(det(VN)) + part.LB(cauchy) -
    #  log(gamma(a0)) + a0*log(b0) -
    #  b0*aN/bN - aN*log(bN) +
    #  log(gamma(aN)) + aN
    # cat("lower bound ",int.linear(Xr,d,y,sigma,alpha,mu,s) +
          # int.gamma(logodds,alpha) +
          # int.klbeta(alpha,mu,s,sigma*sa),"\n")
    # calculate some variables which would be used later
    #E.squares = (crossprod(WN) + tr(VN.inverse))
    # part 1 (log likelihood of part 1)
    # tmp.pt1 <- 0
#     for(ca.i in 1:p){
#       tmp.pt1 <- tmp.pt1 + log(sigm(cauchy[ca.i])) + t(annot[ca.i,])%*%WN%*%(gammas[ca.i] -.5) -
#               lambd.cauchys[ca.i]*t(annot[ca.i,])%*%annot[ca.i,]*E.squares + lambd.cauchys[ca.i] * cauchy[ca.i]^2
#     }
#     
#     tmp.pt2 <- -0.5*m*log(2*pi) - 0.5*aN/bN*E.squares + 
#               0.5*m*(log(bN) - log(aN-1)-1/(2*(aN-2)))
#     tmp.pt3 <- -log(gamma(a0)) -a0*log(b0) -(a0+1)*(log(bN) -log(aN-1) -1/(2*(aN-2))) - a0 -1
#     tmp.pt4 <- -0.5*log(det(VN)) - 0.5*m*(1+log(2*pi))
#     tmp.pt5 <- -log(gamma(aN)) + (aN+1)*log(aN-1) - (aN+1)/(2*(aN-2)) - aN+1 - log(bN)
#     
#     cat("tmp parts ", tmp.pt1,"\t",tmp.pt2,"\t",tmp.pt3,"\t",tmp.pt4,"\t",tmp.pt5,"\n")
#     logw[iter] <- int.linear(Xr,d,y,sigma,alpha,mu,s) +
#       int.gamma(logodds,alpha) +
#       int.klbeta(alpha,mu,s,sigma*sa) + 
#       tmp.pt1 + 
#       tmp.pt2 +
#       tmp.pt3 -
#       tmp.pt4 -
#       tmp.pt5
      
    # cat("after computing updated variational lower bound: ",int.linear(Xr,d,y,sigma,alpha,mu,s) +
          # int.gamma(logodds,alpha) +
          # int.klbeta(alpha,mu,s,sigma*sa),"\n")
    # cat("the second part: ", t(WN),"\n")
    # cat(VN.inverse,"\n")
    # (2h) CHECK CONVERGENCE
    # ----------------------
    # Print the status of the algorithm and check the convergence
    # criterion. Convergence is reached when the maximum difference
    # between the posterior probabilities at two successive iterations
    # is less than the specified tolerance, or when the variational
    # lower bound has decreased.
    err[iter] <- max(abs(alpha - alpha0))
    if (verbose) {
      if (is.null(outer.iter))
        status <- NULL
      else
        status <- sprintf("%05d ",outer.iter)
      progress.str <-
          paste(status,sprintf("%05d %+13.6e %0.1e %06.1f %0.1e %0.1e",
                               iter,logw[iter],err[iter],sum(alpha),
                               sigma,sa),sep="")
      cat(progress.str)
      cat(rep("\r",nchar(progress.str)))
    }
     cat(iter,"!!!!!!hhhh\t",logw0,"!!!!hhh\n")
     cat("~~~",logw[iter],"hhh\n")
    if (logw[iter] < logw0) {
      logw[iter] <- logw0
      err[iter]  <- 0
      sigma      <- sigma0
      sa         <- sa.old
      alpha      <- alpha0
      mu         <- mu0
      s          <- s0
      break
    } else if (err[iter] < tol)
      break
  }
  
  return(list(logw = logw[1:iter],err = err[1:iter],sigma = sigma,sa = sa,
              alpha = alpha,mu = mu,s = s))
}

# ----------------------------------------------------------------------
# Computes an integral that appears in the variational lower bound of
# the marginal log-likelihood. This integral is the expectation of the
# linear regression log-likelihood taken with respect to the1
# variational approximation.
int.linear <- function (Xr, d, y, sigma, alpha, mu, s) {
  n <- length(y)
  return(-n/2*log(2*pi*sigma) - norm2(y - Xr)^2/(2*sigma) 
         - dot(d,betavar(alpha,mu,s))/(2*sigma))
}

# ---------------------------------------------------------------------
# compute the sigmoid function @@ added by W.L. on 1/5/2018
sigm <- function(x){
  temp <- 1+exp(-x)
  return(1/temp)
}

# ---------------------------------------------------------------------
# compute the lambda function for another lower bound in the logistic regression
# added by W.L. on 1/5/2018
lamd <- function(cauchy){
  tmp <- sigm(cauchy)
  return((tmp-.5)/(2*cauchy))
}

# ---------------------------------------------------------------------
# compute updated invered covariance matrix in the logistic regression
# added by W.L. on 1/5/2018
inv.V <- function(E.eta, lambd.cauchys,annotat.array){
  p <- nrow(annotat.array) 
  m <- ncol(annotat.array) # number of annotation types
#   cat("in the inv.v function, dimension of annotat.array ", m,"\n")
#   cat("in the inv.v function, E.eta ", E.eta,"\n")
#   cat("in the inv.v function, class of E.eta ", class(E.eta),"\n")
  tmp1 = diag(E.eta,nrow=m,ncol=m)
  tmp2 = diag(0,nrow=m,ncol=m)
  
  # cat("dim of tmp2 ",dim(tmp2),"\n")
  # cat("lambda.cauchys[i] ", lambd.cauchys[i],"\n")
  # cat("annotat.array dimension ", dim(crossprod(t(annotat.array[i,]))),"\n")
  for(i in 1:p){
    tmp2 = tmp2 + 2 * lambd.cauchys[i] * crossprod(t(annotat.array[i,]))
  }
  
  return(tmp1 + tmp2)
}


# ---------------------------------------------------------------------
# compute updated expectation of regression coefficient vector in the logistic regression
# added by W.L. on 1/5/2018
Expect.W <- function(V,gammas,annot){
  p <- length(gammas)
  m <- ncol(annot)
  tmp <- matrix(0,m,1)
  for(i in 1:p){
    tmp = tmp + (gammas[i] - 0.5)* (annot[i,])
  }
  
  return(V %*% tmp)
}

# ---------------------------------------------------------------------
# compute part of the lower bound for logistic regression in annotation layer
# added by W.L. on 1/5/2018
part.LB <- function(cauchy){
  m = length(cauchy)
  
  result <- 0
  for(i in 1:m){
    tmp = -log(sigm(cauchy[i])) - cauchy[i]/2 + lamd(cauchy[i])*cauchy[i]^2
    result = result + tmp
  }
  
  return(result)
}
