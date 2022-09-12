###
### 2d test functions
###

### Rosenbrock function
f_rosen <- function(x) 100*(x[2] - x[1]^2)^2 + (1 - x[1])^2
grad_rosen <- function(x) {
  c(-400*x[1]*(x[2] - x[1]^2) - 2*(1 - x[1]),
    200*(x[2] - x[1]^2))
}
hess_rosen <- function(x) {
  matrix(c(1200*x[1]^2 - 400*x[2]^2 + 2, -400*x[1], -400*x[1], 200),ncol = 2)
}

### Himmelblau function
f_himmel <- function(x) (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
grad_himmel <- function(x) {
  c(4*x[1]*(x[1]^2 + x[2] - 11) + 2*(x[1] + x[2]^2 - 7),
    2*(x[1]^2 + x[2] - 11) + 4*x[2]*(x[1] + x[2]^2- 7))
}
hess_himmel <- function(x) {
  matrix(c(12*x[1]^2 + 4*x[2] - 42, 4*(x[1] + x[2]),
           4*(x[1] + x[2]), 4*x[1] + 12*x[2]^2 - 26), ncol = 2)
}