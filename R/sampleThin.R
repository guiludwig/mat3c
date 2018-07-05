sampleThin <- function(fmat, GAMMA = gamma[1], KAPPA = kappa[1], SIGMA = sigma[1],
                       WIN = win, R = R_clusters, UP = up){
  # This function creates a ghost of the thinned fingers in the dataset generation
  #!# Why limx, limy and a1-2, b1-2 are not the same?

  fmat.cond <- list()

  for (ij in seq_along(fmat)){

    nx <- max(length(fmat[[ij]]$x), nrow(fmat[[ij]]$fingers))

    #!# BIG CHANGE: I'm only regenerating times if they didn't exist
    if(is.null(fmat[[ij]]$t)){
      # t.cond <- sort(rexp(nx, 2.3)) #!# Why 2.3?
      t.cond <- sort(rbeta(nx, 1, 8)) # f <- function(x) dbeta(x, 1, 8); curve(f, xlim = c(0,1))
      # d.s.t <- as.numeric(t.cond > 1)
      #
      # while (sum(d.s.t) > 0) {
      #   t.cond[d.s.t == 1] <- rexp(sum(d.s.t), 2.3) #!# Why 2.3?
      #   d.s.t <- as.numeric(t.cond > 1) #!# I think this step might be bad
      # }

      ordem <- order(t.cond)

      if(is.null(fmat[[ij]]$fingers)) {
        fmat.cond[[ij]] <- list(centers = c(fmat[[ij]]$centers[1],
                                            fmat[[ij]]$centers[2]),
                                x = fmat[[ij]]$x[ordem],
                                y = fmat[[ij]]$y[ordem],
                                t = sort(t.cond),
                                direcao = UP[ij])
      } else {
        fmat.cond[[ij]] <- list(centers = c(fmat[[ij]]$centers[1],
                                            fmat[[ij]]$centers[2]),
                                x = fmat[[ij]]$fingers[ordem, 1],
                                y = fmat[[ij]]$fingers[ordem, 2],
                                t = sort(t.cond),
                                direcao = UP[ij])
      }
    } else {
      fmat.cond[[ij]] <- fmat[[ij]]
    }

    marked.locations <- geradorVonmises(UP[[ij]], GAMMA,
                                        fmat[[ij]]$centers[1:2],
                                        KAPPA, SIGMA, WIN)

    marked.n <- nrow(marked.locations[, , drop = FALSE])

    if (marked.n > 0) {
      marked.times <- runif(marked.n)
      ordem <- order(marked.times)
      fingers0 <- matrix(c(marked.locations[ordem, 1],
                           marked.locations[ordem, 2],
                           marked.times[ordem]),
                         nrow = marked.n, ncol = 3)

      ind.fingers <- logical(marked.n)

      if (marked.n > 1){
        for (il in 1:marked.n){
          for (jl in seq_along(fmat.cond[[ij]]$x)) {
            dist.l <- (fingers0[il,1] - fmat.cond[[ij]]$x[jl])^2 + (fingers0[il,2] - fmat.cond[[ij]]$y[jl])^2
            ind.fingers[il] <- ind.fingers[il] | (dist.l < R^2 & fingers0[il,3] > fmat.cond[[ij]]$t[jl])
          }
        }
      }

      fmat.cond[[ij]]$thinned <- fingers0[ind.fingers, 1:3, drop = FALSE]

    }
  }

  return(fmat.cond)

}
