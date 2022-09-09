###
### Optimization Algorithms
###

### backtracking line search (with default values) and descent direction check
backtracking <- function(f, gradf, x, d, tmax = 1, c = 0.01, rho = 0.5) {
  t <- tmax
  a <- f(x)               # intercept of linear function
  v <- gradf(x) %*% d
  if (v >= 0) stop("d is not a descent direction")
  b <- c * v              # slope of linear function
  while (f(x + t*d) > a + b*t) {
    t <- rho*t
  }
  return(t)
}

### Steepest descent method with backtracking line search
steepest_descent = function(f, gradf, x, t = NA, tol = 1e-10){
  n_iter <- 0
  iterates <- list(x)             # save iterates (for visualization)
  linesearch <- is.na(t)          # check if constant step size or linesearch
  while (sqrt(sum(gradf(x)^2)) > tol ) {
    d <- -gradf(x)                        # steepest descent direction
    if (linesearch) {
      t <- backtracking(f, gradf, x, d)   # backtracking line search
    }
    x <- x + t*d                          # update
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
newton = function(f, gradf, hessf, x, tol = 1e-10){
  n_iter <- 0 
  iterates <- list(x)             # save iterates (for visualization)
  while (sqrt(sum(gradf(x)^2)) > tol) {
    d <- solve(hessf(x), -gradf(x))   # Newton-direction via Newton system
    x <- x + d                        # Newton step with step length 1
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

### Damped Newton
newton_damped = function(f, gradf, hessf, x, tol = 1e-10){
  n_iter <- 0
  iterates <- list(x)             # save iterates (for visualization)
  while (sqrt(sum(gradf(x)^2)) > tol) {
    d <- solve(hessf(x), -gradf(x))       # Newton-direction via Newton system
    t <- backtracking(f, gradf, x, d)     # backtracking line search
    x <- x + t*d                          # update
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

### Inverse BFGS-update for quasi-Newton method
inverse_bfgs_update = function(gradf, x, x_new, H){
  s <- x_new - x
  y <- gradf(x_new) - gradf(x)
  rho <- as.numeric(1/crossprod(s,y))
  v <- H %*% y
  U <- -rho*(tcrossprod(s,v) + tcrossprod(v,s)) +
    as.numeric(rho^2*(1/rho + crossprod(y,v)))*tcrossprod(s,s)
  return(H + U)
}

#### BFGS-quasi-Newton method
quasi_newton_bfgs = function(f, gradf, x, tol = 1e-10){
  n_iter <- 0 
  iterates <- list(x)             # save iterates (for visualization)
  H <- diag(length(x))            # initial inverse Hessian approximation
  while (sqrt(sum(gradf(x)^2)) > tol) {
    d <- as.vector(-H %*% gradf(x))       # quasi-Newton direction
    t <- backtracking(f, gradf, x, d)     # backtracking line search
    x_new <- x + t*d 
    H <- inverse_bfgs_update(gradf, x, x_new, H)   # inverse BFGS-update
    x <- x_new                                     # update iterate
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