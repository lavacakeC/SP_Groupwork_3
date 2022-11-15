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



func_detail <- function(theta, func, grad, hess=NULL, eps, ...) {
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
    
    func_result <- func_detail(theta, func, grad, hess, eps, ...)
    
    flag_grad = FALSE
    flag_hess = TRUE
    
    grad_zeros <- sum(attr(func_result, "grad") < tol * func_result + fscale)
    if (grad_zeros == length(theta)) flag_grad = TRUE

    
    tryCatch(
      {
        attr(func_result, "hess_inverse") <- chol(attr(func_result, "hess"))
      }, 
      error =  function(e){
        flag_hess <<- FALSE
      })
    
    if (flag_grad == TRUE & flag_hess == TRUE) break
    
    if (flag_hess == FALSE) func_result <- perturb_hess(func_result) 
    
    theta <- theta_calculate(theta, func, func_result, max.half, ...)
    
    iter <- iter + 1
  }
  answer <- list(f = func_result[1], theta = theta, iter = iter
                 , g = attr(func_result, "grad"), Hi = attr(func_result, "hess_inverse"))
  return(answer)
}


newt(theta = c(1.2,1), func = rb, grad = gb, k=2, tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)



