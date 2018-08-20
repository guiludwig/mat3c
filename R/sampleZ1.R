sampleZ1 <- function(BE = beta.p, PH = phi.p, GA = gamma.p,
                     SI = sigma.p, KA = kappa.p, WIN = win, JJ = J,
                     RRc = R_centers, RR = R_clusters,
                     debug = TRUE){

  suff.mcmc <- matrix(0, nrow = JJ, ncol = 2)

  for (j in 1:JJ){

    dominant.j <- geradorVonmisesMCMC(BE, PH, GA, SI, KA, WIN,
                                      Rc = RRc, R = RR) #!#

    final.j <- list()
    k <- 1
    for (i in 1:length(dominant.j)) {
      if (dominant.j[[i]]$color == "keep") {
        temp <- matrix(dominant.j[[i]]$fingers, ncol = 3)[, 1:2, drop = FALSE]
        final.j[[k]] <- list(centers = c(dominant.j[[i]]$centers[1],
                                         dominant.j[[i]]$centers[2]),
                             fingers = temp)
        k <- k + 1
      }
    }

    suff.mcmc[j,] <- sufficientStat(final.j, Rc = RRc)

    if(debug) {
      cat(paste0("\tj = ",j,", number of centers = ", suff.mcmc[j,1], ", number of fingers near centers by ", RRc, ": ", -1*suff.mcmc[j,2], ".\n"))
    }

  }

  return(suff.mcmc)

}
