resampleTimes <- function(fmat.cond, R = R_clusters){
  # This function creates a new sample of times for the non-thinned fingers, given
  # the results from sampleThin().
  # Rao et al. 2015, p.10: t_min = time of oldest event thinned only by g
  # u_i ~ U[0,t_min]

  for (i in seq_along(fmat.cond)) {

    nx <- length(fmat.cond[[i]]$x) #!# fmat is used only here?
    nxx <- length(fmat.cond[[i]]$thinned[,1])

    exclusively <- logical(nxx)
    # t_min <- numeric(nx)
    t_min <- 1

    if(nxx == 0){
      t_min <- 1 # t_min + 1 # whole vector
    } else {
      for(j in seq_len(nx)){
        exclusively <- ((fmat.cond[[i]]$thinned[,1] - fmat.cond[[i]]$x[j])^2 + (fmat.cond[[i]]$thinned[,2] - fmat.cond[[i]]$y[j])^2 < R^2) & (fmat.cond[[i]]$t[j] < fmat.cond[[i]]$thinned[,3])
        for(j2 in seq_len(nx)[-j]){
          exclusively <- exclusively & (!(((fmat.cond[[i]]$thinned[,1] - fmat.cond[[i]]$x[j2])^2 + (fmat.cond[[i]]$thinned[,2] - fmat.cond[[i]]$y[j2])^2 < R^2) & (fmat.cond[[i]]$t[j2] < fmat.cond[[i]]$thinned[,3])))
        }
        if(any(exclusively)){
          # t_min[j] <- which.min(fmat.cond[[i]]$thinned[exclusively, 3])
          t_min <- min(fmat.cond[[i]]$thinned[exclusively, 3])
        } else {
          # t_min[j] <- 1
          t_min <- 1
        }

        fmat.cond[[i]]$t[j] <- runif(1, 0, t_min)

      }
    }

    #!# Atualizar o tempo aqui
    sample.times <- fmat.cond[[i]]$t # runif(nx, rep(0, nx), t_min) # runif recycles arguments
    fmat.cond[[i]]$t <- sort(sample.times)
    fmat.cond[[i]]$x <- fmat.cond[[i]]$x[order(sample.times)]
    fmat.cond[[i]]$y <- fmat.cond[[i]]$y[order(sample.times)]

    # accept <- min(exp(likelihood_matern(fmat.prop, beta[l-1],phi[l-1],gamma[l-1],sigma[l-1],kappa[l-1])-likelihood_matern(fmat.cond, beta[l-1],phi[l-1],gamma[l-1],sigma[l-1],kappa[l-1])), 1)
    # accept <- 1

    # print(accept)
    # u1 <- runif(1,0,1)
    # if (!is.nan(accept) && u1 < accept){
    #   fmat.cond[[i]] <- fmat.prop[[i]]
    # } else {
    #   fmat.prop[[i]] <- fmat.cond[[i]]
    # }

  }

  return(fmat.cond)

}
