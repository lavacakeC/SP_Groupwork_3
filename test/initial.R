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
