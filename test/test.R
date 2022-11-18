############ test



z <- 2
attr(z, 'grad') <- gb(c(1,2))
z[1]
z + 2
attr(z,'grad')

dim(test)[1]

tryCatch(
  {
    attr(func_result, "hess_inverse") <- chol2inv(chol(test))
  }, 
  error =  function(e){
  print("error")
})


test1 <- function(){
  for (i in 10:11)
    print(i)
}

for (i in 1:5){
  if(i == 5 ){
    warning("stop")
  }
  print(i)
}

if(i == 5 ){
  stop("stop")
}


flag_hess <- FALSE
while (!flag_hess) {
  
  # hess <- hess + multiple * norm(hess) * diag(dim(hess)[1])
  
  flag_hess <- TRUE
  
  tryCatch(
    {
      attr(func_result, "hess_inverse") <- chol(test)
    }, 
    error =  function(e){
      flag_hess <<- FALSE
    })

  # multiple <- multiple * 10
  
}

test <- c(1,2,3)
sum(test^2)

he <- matrix(c(3.035084,669.9264,669.926369,164428.5608),2,2)
he
solve(he)

all(is.finite(as))

is.nan()
test <- TRUE
test1 <- TRUE
all(test,test1)
1 & !0
as <- c(1,2,NA)

log(-1)

###############################################
test <- c(0,0,0)
test2 <- c(1,-1,2)
test < test2
test <- NaN

x <- Hilbert(9)

chol(test)

test_1 <- test + (1e0) * norm(test) * diag(3)

chol(test_1)


norm(test,'I')

class(test)
?norm
test<1

sum(test < 10) == NA

test <- matrix(c(1:9), 3, 3)
sum(test == 0) == length(test)
chol(test)
chol2inv(chol(test + 1e0 * norm(test) * diag(dim(test)[1])))

test1 <- function(a,b,c,...,d){
  a <-  ...
  return(a)
}
test1(a=1,b=2,c=3,4,5,d=6)

?paste
newt()
for( i in 1:maxit){
  func(theta,...)
  grad(theta,...)
  hess(theta,...)
  
  # Flag hess NULL or not
  
  # if grad() < tol * func() + fscale   flag1 = T
  
  # if hess positive definie flag2 = T
  
  # else
  
  # if flag = T 
  
  # if flag 1 = F, flag 2 = T
  
  # solve delta
  
  # if flag 1,2 = F
  
  # while flag3 = F  
  
  # hess = hess + 1en * norm(hess) * diag(length(theta)
  
  # try chol
  
  # delta = - chol2inv(chol(hess_new)) %*% grad
  
  # func(theta + delta) is non-finite
  # while func(...) >= func(theta)
  # func(theta + delta) is non-finite
  # delta = 0.5 delta
  
  # theta_next = theta + delta
  
  
}
func(theta,...)
grad(theta,...)
hess(theta,...)

for (i in 1:9){
  print(i)
}


####### FINAL TEST


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


newt(c(0,0), rb, gb, hb)

newt(theta = c(0,0), func = rb, grad = gb, k=2, tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6)


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
