sampleZ1s <- function(BE = beta.p, PH = phi.p, GA = gamma.p,
                     SI = sigma.p, KA = kappa.p, WIN = win, JJ = J,
                     RRc = R_centers, RR = R_clusters,
                     debug = TRUE){

  suff.mcmc <- matrix(0, nrow = JJ, ncol = 2)

  for (j in 1:JJ){

    # dominant.j <- geradorVonmisesMCMC(BE, PH, GA, SI, KA, WIN,
                                      # Rc = RRc, R = RR) #!#

    dominant.j <- rpoispp(BE, win = WIN)

    # final.j <- list()
    # k <- 1
    # for (i in 1:length(dominant.j)) {
    #   if (dominant.j[[i]]$color == "keep") {
    #     temp <- matrix(dominant.j[[i]]$fingers, ncol = 3)[, 1:2, drop = FALSE]
    #     final.j[[k]] <- list(centers = c(dominant.j[[i]]$centers[1],
    #                                      dominant.j[[i]]$centers[2]),
    #                          fingers = temp)
    #     k <- k + 1
    #   }
    # }
    #
    final.j <- list(seq_len(dominant.j$n))
    for (i in seq_len(dominant.j$n)) {
      temp <- geradorVonmises(runif(1, 0, 2*pi), GA,
                              c(dominant.j$x[i], dominant.j$y[i]),
                              KA, SI, WIN)

      marked.n <- nrow(temp)
      marked.times <- runif(marked.n)
      ord <- order(marked.times)
      fingers0 <- cbind(temp[ord, , drop=FALSE], marked.times[ord])

      ind.fingers <- nextFingers(marked.n, fingers0, RR)

      final.j[[i]] <- list(centers = c(dominant.j$x[i],
                                       dominant.j$y[i]),
                           fingers = fingers0[ind.fingers, 1:3, drop = FALSE])
    }

    suff.mcmc[j,] <- sufficientStat(final.j, Rc = RRc)

    if(debug) {
      cat(paste0("\tj = ",j,", number of centers = ", suff.mcmc[j,1], ", number of fingers near centers by ", RRc, ": ", -1*suff.mcmc[j,2], ".\n"))
    }

  }

  return(suff.mcmc)

}
