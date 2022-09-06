
source("objectives.R")
source("aux.R")


### backtracking line search (with default values)
backtracking <- function(f, gradf, x, d, tmax = 1, c = 0.01, rho = 0.5) {
  t <- tmax
  a <- f(x)                      # intercept of linear function
  b <- c * gradf(x) %*% d        # slope of linear function
  while (f(x + t*d) > a + b*t) {
    t <- rho*t
  }
  return(t)
}

### Steepest descent method with backtracking line search
steepest_descent = function(f, gradf, x, t = NA, tol = 1e-10){
  n_iter <- 0                     # iteration counter
  iterates <- list(x)             # save iterates (for visualization)
  linesearch <- is.na(t)          # check if constant step size or linesearch
  while (sqrt(sum(gradf(x)^2)) > tol ) {     # check stopping criterion
    d <- -gradf(x)                           # steepest descent direction
    if (linesearch) {
      t <- backtracking(f, gradf, x, d)      # backtracking line search
    }
    x <- x + t*d                             # steepest descent step
    n_iter <- n_iter+1
    iterates[[n_iter+1]] <- x
  }
  result <- list()
  result$results <- list("solution" = x,        # solution vector
                         "opt_value" = f(x),    # optimal value
                         "n_iter" = n_iter)     # number of iterations
  result$iterates <- iterates                   # all iterates (visualization)
  return(result)   
}

### Newton method
newton = function(f, gradf, hessf, x, t = NA, tol = 1e-10){
  n_iter <- 0                     # iteration counter
  iterates <- list(x)             # save iterates (for visualization)
  linesearch <- is.na(t)          # check if constant step size or linesearch
  while (sqrt(sum(gradf(x)^2)) > tol) {   # check stopping criterion
    d <- solve(hessf(x), -gradf(x))       # Newton-direction as solution of the Newton system
    if (linesearch) {
      t <- backtracking(f, gradf, x, d)     # backtracking line search
    }
    x <- x + t*d              # Newton step
    n_iter <- n_iter+1
    iterates[[n_iter+1]] <- x
  }
  result <- list()
  result$results <- list("solution" = x,        # solution vector
                         "opt_value" = f(x),    # optimal value
                         "n_iter" = n_iter)     # number of iterations
  result$iterates <- iterates                   # all iterates (visualization)
  return(result)
}

#### BFGS-quasi-Newton method
#
inv_bfgs_update = function(gradf, x, x_new, H){             # compute inverse BFGS-update
  s <- x_new - x
  y <- gradf(x_new) - gradf(x)
  rho <- as.numeric( 1/crossprod(s,y) )
  v <- H%*%y                                         # compute matrix-vector product only once
  U <- -rho*( tcrossprod(s,v) + tcrossprod(v,s) ) +
    as.numeric(rho^2*(1/rho + crossprod(y,v)))*tcrossprod(s,s)
  return(U)
}
#
qnewton_bfgs = function(f, gradf, x, t=NA, tol = 1e-10){
  n_iter <- 0                     # iteration counter
  iterates <- list(x)             # save iterates (for visualization)
  linesearch <- is.na(t)          # check if constant step size or linesearch
  H <- diag(length(x))            # initial inverse Hessian approximation
  while (sqrt(sum(gradf(x)^2)) > tol) {     # check stopping criterion
    d <- as.vector(-H %*% gradf(x))                    # quasi-Newton direction
    if (linesearch) {
      t <- backtracking(f, gradf, x, d)    # backtracking line search
    }
    x_new <- x + t*d                            # compute new iterate
    U <- inv_bfgs_update(gradf, x, x_new,  H)   # inverse BFGS-update
    H <- H + U                                  # update quasi-Newton matrix
    x <- x_new                                  # update iterate
    n_iter <- n_iter+1
    iterates[[n_iter+1]] <- x
  }
  result <- list()
  result$results <- list("solution" = x,        # solution vector
                         "opt_value" = f(x),    # optimal value
                         "n_iter" = n_iter)     # number of iterations
  result$iterates <- iterates                   # all iterates (visualization)
  return(result)
}