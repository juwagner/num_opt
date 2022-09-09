###
### Auxiliary functions for visualization
###

### required libraries
library("grid")     # for visualization
library("lattice")  # for visualization

### plot of objective
plot_objective_2d <- function(f, xlim, ylim, res = c(0.01, 0.01)) {
  x <- seq(xlim[1], xlim[2], res[1])
  y <- seq(ylim[1], ylim[2], res[2])
  X <- expand.grid(x, y)
  fx <- sapply(1:nrow(X), function(i) f(as.vector(unlist(X[i,]))))
  levelplot(fx~X[,1]*X[,2],
            xlab="x", ylab="y",
            col.regions = terrain.colors(100))
}

### plot of objective and iterates
plot_iterates_2d <- function(f, iterates, xlim, ylim, res = c(0.01, 0.01)) {
  x <- seq(xlim[1], xlim[2], res[1])
  y <- seq(ylim[1], ylim[2], res[2])
  X <- expand.grid(x, y)
  fx <- sapply(1:nrow(X), function(i) f(as.vector(unlist(X[i,]))))
  x_k <- do.call(cbind, iterates) 
  levelplot(fx~X[,1]*X[,2],
            xlab="x", ylab="y",
            col.regions = terrain.colors(100),
            panel = function(...){
              panel.levelplot(...)
              grid.points(x_k[1,], x_k[2,], pch=3, size=unit(0.2,"char"))
            }
  )
}