
func_detail <- function(theta, func, grad, hess, eps, ...) {
  func_result <- func(theta, ...)
  attr(func_result, "grad") <- grad(theta, ...)
  
  if (is.null(hess)) {
    hess <- matrix(0, length(theta), length(theta))
    
    for (i in 1:length(theta)){
      theta_1 <- theta
      theta_1[i] <- theta_1[i] + eps
      grad_theta_1 <- grad(theta_1, ...)
      hess[i,] <- (grad_theta_1 - attr(func_result, "grad")) / eps
    }
  }
  
  else hess <- hess(theta, ...)
  attr(func_result, "hess") <- hess
  return(func_result)
}


perturb_hess <- function(func_result, flag_hess = FALSE) {
  multiple <- 1e-6
  hess <- attr(func_result, "hess")
  while (!flag_hess) {
    
    hess <- hess + multiple * norm(hess) * diag(dim(hess)[1])
    
    flag_hess <- TRUE
    
    tryCatch(
      {
        ########################################
        # hess_upper_tri not inverse
        #########################################
        attr(func_result, "hess_inverse") <- chol(hess)
      }, 
      error =  function(e){
        flag_hess <<- FALSE
      })
    
    multiple <- multiple * 10
    
  }
  
  return(func_result)
}

theta_calculate <- function(theta, func, func_result, max.half, ...) {

  delta <- - chol2inv(attr(func_result, "hess_inverse")) %*% attr(func_result, "grad")
  
  flag_delta <- func(theta + delta, ...) < func(theta, ...)
  
  half <- 0
  
  while (flag_delta == FALSE & half <= max.half) {
    half <- half + 1
    delta <- delta / 2
    #############
    
    #############
    #need to check func(theta + delta, ...) is not non-finite
    
    #############
    flag_delta <- func(theta + delta, ...) < func(theta, ...)
  }
  
  if (flag_delta == TRUE) return(theta + delta)
  
  if (half > max.half) {
    stop("The step has failed to improve the objective.")
  }
}


newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  iter <- 0
  for (i in 1:maxit) {
    # func_result <- func_detail(theta, func, grad, hess, eps, ...)
    func_result <- func_detail(theta = theta, func = func, grad = grad, hess = hess, eps = eps, ...)
    
    flag_grad = FALSE ### change
    flag_hess = TRUE

    if (sqrt(sum(attr(func_result, "grad") ^ 2)) < (abs(func_result) * tol + fscale)) flag_grad = TRUE

    
    tryCatch(
      {
        attr(func_result, "hess_inverse") <- chol(attr(func_result, "hess")) ### change
      }, 
      error =  function(e){
        flag_hess <<- FALSE
      })
    
    if (flag_grad == TRUE & flag_hess == TRUE) break
    
    if (flag_hess == FALSE) func_result <- perturb_hess(func_result) 
    
    theta <- theta_calculate(theta = theta, func = func, func_result = func_result, max.half = max.half, ...)
    
    iter <- iter + 1
  }
  answer <- list(f = func_result[1], theta = theta, iter = iter
                 , g = attr(func_result, "grad"), Hi = attr(func_result, "hess_inverse"))
  return(answer)
}



#################### test

rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}




newt(theta = c(1.2,1), func = rb, grad = gb, k=2, tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)


nll <- function(theta,t,y=y) {
  ## -ve log likelihood for AIDS model y_i ~ Poi(alpha*exp(beta*t_i))
  ## theta = (alpha,beta)
  mu <- theta[1] * exp(theta[2] * t) ## mu = E(y)
  -sum(dpois(y,mu,log=TRUE)) ## the negative log likelihood
} ## nll

gll <- function(theta,t,y=y) {
  ## grad of -ve log lik of Poisson AIDS early epidemic model
  alpha <- theta[1];beta <- theta[2] ## enhances readability
  ebt <- exp(beta*t) ## avoid computing twice
  -c(sum(y)/alpha - sum(ebt), ## -dl/dalpha
     sum(y*t) - alpha*sum(t*ebt)) ## -dl/dbeta
} ## gll

hll <- function(theta,t,y) {
  ## Hessian of -ve log lik of Poisson AIDS early epidemic model
  alpha <- theta[1];beta <- theta[2] ## enhances readability
  ebt <- exp(beta*t) ## avoid computing twice
  H <- matrix(0,2,2) ## matrix for Hessian of -ve ll
  H[1,1] <- sum(y)/alpha^2
  H[2,2] <- alpha*sum(t^2*ebt)
  H[1,2] <- H[2,1] <- sum(t*ebt)
  H
} ## hll
th0 <- c(10,.1)
t80 <- 1:13 ## years since 1980
y <- c(12,14,33,50,67,74,123,141,165,204,253,246,240) ## AIDS cases
newt(theta = c(10,.1), func = nll,grad = gll, t = t80,y = y)

gll(th0,t80,y)

