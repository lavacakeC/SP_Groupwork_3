## Group member: 
## Sixiang Cheng, s2437109; Ruishuo Cheng, s2305931; Qinxuan Li, s2299101
## github repo's address:https://github.com/lavacakeC/SP_Groupwork_3.git
## Ruishuo and Qinxuan each undertake about 30% of the work,
## Sixiang undertake about 40% of the work

## In calculus, Newton's method is an iterative method for finding the roots of 
## a differentiable function, as such, Newton's method can be applied to the 
## derivative of a twice-differentiable function to find the roots of the derivative, 
## also known as the critical points. These solutions may be minima, maxima, 
## or saddle points.
## Codes implement Newton's method for minimization of functions.
## Newton's method uses the first few terms of the Taylor series 
## for the function f(x) to iterate and seek the minimum of approximation 
## at each step.
## Step 1: Calculate the objective function, gradient function and 
##         Hessian matrix function
## Step 2: Test whether convergence by seeing whether all elements of
##         the gradient vector have absolute value less than tol times the 
##         absolute value of the objective function plus f-scale. And then test 
##         whether Hessian is positive definite.
## Step 3: If Hessian is not positive definite, perturb it by setting the 
##         multiplier to a small multiple of a matrix norm of the Hessian. If it 
##         still doesn't do it, repeatedly multiply the multiplier until we get 
##         the positive definite Hessian
## Step 4: Solve for the search direction delta
## Step 5: If the value of function on theta plus delta is not less than the 
##         value of function on theta, then we repeatedly halve delta until it is
## Step 6: Iterate theta by theta plus delta, increment by k=1 and return to step 1


func_detail <- function(theta, func, grad, hess, eps, iter, ...) {
  
  ## This function is a function that 
  ## get the result from objective function, 
  ## gradient function and hessian matrix function.
  ## The inputs is a vector of initial values theta 
  ## Function will return the result from objective 
  ## function, gradient vector of the objective 
  ## and hessian matrix of the objective.
  
  func_result <- func(theta, ...) ## get result from objective function 
  
  ## check if the objective or derivatives are not finite at the initial theta
  if (iter == 0 & !is.finite(func_result)) stop("objective is not finite at the initial theta")
  ## use attr() to specify attributes "grad" 
  ## and associate gradient function's results
  attr(func_result, "grad") <- grad(theta, ...)
  
  ## check if the objective or derivatives are not finite at the initial theta
  if (iter == 0 & !all(is.finite(attr(func_result, "grad")))) stop("derivatives are not finite at the initial theta")
  
  if (is.null(hess)) {
    ## if hessian matrix function not supplied, 
    ## set a hessian approximate matrix and calculate each value in the matrix
    hess <- matrix(0, length(theta), length(theta))
    
    for (i in 1:length(theta)){
      theta_1 <- theta
      theta_1[i] <- theta_1[i] + eps ## increase th0[i] by eps
      grad_theta_1 <- grad(theta_1, ...) 
      # approximate -dl/dth[i]
      hess[i,] <- (grad_theta_1 - attr(func_result, "grad")) / eps 
    }
  }
  ## If Hessian matrix function is provide use attr() to specify attributes and 
  ## associate hessian matrix function's results.
  else hess <- hess(theta, ...)
  attr(func_result, "hess") <- hess
  return(func_result)
} ## func_detail


perturb_hess <- function(func_result) {
  
  ## This function calculate perturbed hessian matrix
  ## when the original hessian matrix is not positive definite
  ## The input is results from func_detail function
  ## It will return the same value except that the 
  ## hessian matrix becomes a perturbed hessian matrix
  
  ## flag shows hessian matrix is positive definite or not
  flag = FALSE
  multiple <- 1e-6   ## set a multiple, starting with 1e-6
  hess <- attr(func_result, "hess")  ## get original hessian matrix
  while (!flag) {
    ## false flag means the hessian matrix is not positive definite
    ## multiply the identity matrix by the original Hessian norm, 
    ## and then multiply by multipler
    hess <- hess + multiple * norm(hess) * diag(dim(hess)[1])
    
    tryCatch(
      {
        ## Try to inverse hessian matrix
        attr(func_result, "hess_inverse") <- chol2inv(chol(hess))
        flag <- TRUE  ## change the flag
      }, 
      error =  function(e){
        flag <<- FALSE   ## flag changed if there has error
      })
    
    multiple <- multiple * 10 ## set new multiple, 10 times larger than last one
    
  }
  
  return(func_result)
} ## perturb_hess

theta_calculate <- function(theta, func, func_result, max.half, ...) {
  
  ## This function is a function calculate new theta 
  ## The inputs are vector of initial values
  ## It will return new theta
  
  half <- 0 ## number of times of half delta
  ## negative value of multiply inverse of hessian matrix 
  ## and gradient vector to get delta
  delta <- - attr(func_result, "hess_inverse") %*% attr(func_result, "grad")
  ## Flag is true when new theta smaller than original result
  flag_delta <- func(theta + delta, ...) < func(theta, ...)
  ## set a new delta, it only used when function gets non-finite
  delta_new <- delta 
  
  if (!is.finite(func(theta + delta, ...))) {
    ## if putting the first delta into function gets non-finite, find a new delta
    delta_new <- delta / 2 ## set delta to be hald of original one
    half <- 1 ## set number of times of half delta plus 1
    
    ## Flag is true when new theta smaller than original result
    flag_delta <- func(theta + delta_new, ...) < func(theta, ...)
    
  }
  
  while (flag_delta == FALSE & half <= max.half) {
    ## when delta is greater than original theta 
    ## and half is smaller than or equal to max.half
    half <- half + 1 ## set number of times of half delta plus 1
    delta <- delta_new
    delta_new <- delta_new / 2 ## set delta halve each time
    
    if (!is.finite(func(theta + delta, ...))) { 
      ## if function gets non-finite, we need to find a new delta
      ## which lies between
      ## the old delta and the new one
      delta_new <- (delta + delta_new) / 2
    }
    ## change the flag
    flag_delta <- func(theta + delta_new, ...) < func(theta, ...)
  }
  
  if (flag_delta == TRUE) return(theta + delta_new) ## return new theta
  
  if (half > max.half) {
    stop("The step has failed to improve the objective.")
  }
} ## theta_calculate



newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8
                 ,fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  ## This function is main optimization function 
  
  for (iter in 0:maxit) {
    ## get the function result of the function objective, grad and hess matrix 
    func_result <- func_detail(theta = theta, func = func, grad = grad
                               , hess = hess, eps = eps, iter, ...)
    
    flag_grad <- FALSE ## change
    
    ## check Convergence, if it converged, set flag_grad TRUE
    if (all(abs(attr(func_result, "grad")) < tol * (abs(func_result) + fscale))) flag_grad = TRUE

    tryCatch(
      {
        ## try to inverse the hess matrix use the method,
        ## that only positive definite matrices can succeed
        attr(func_result, "hess_inverse") <- chol2inv(chol(attr(func_result, "hess")))
        ## no error means hess matrix is positive definite, set flag_hess to TRUE
        flag_hess <- TRUE
      }, 
      error =  function(e){ ## error means hess matrix is not positive definite
        flag_hess <<- FALSE ## set flag_hess to FALSE
      })
    
    ## if Convergence and positive definite hess matrix are both matched
    ## so we find the optimization parameters, then break
    ## if we reach maximum number of Newton iterations to try before giving up
    ## we also break
    if (all(flag_grad, flag_hess) | iter == maxit ) break
    
    ## If the Hessian is not positive definite at convergence,
    ## we give a warning about that
    if (flag_grad == TRUE & flag_hess == FALSE) warning("The Hessian matrix is not positive definite at convergence")
    
    ## if hess matrix is not positive definite, we need to perturb it
    if (flag_hess == FALSE) func_result <- perturb_hess(func_result) 
    
    ## calculate the new theta
    theta <- theta_calculate(theta = theta, func = func
                             , func_result = func_result
                             , max.half = max.half, ...)
  }
  ## if we do not find the optimization parameters 
  ## and reach maximum number of Newton iterations
  ## we give a warning about that
  if (iter == maxit & !all(flag_grad, flag_hess)) warning("The maximum number of Newton iterations is reached without convergence")
  
  ## create a list containing all the outputs
  answer <- list(f = func_result[1], theta = theta, iter = iter
                 , g = attr(func_result, "grad")
                 , Hi = attr(func_result, "hess_inverse"))
  return(answer)
} ## newt
