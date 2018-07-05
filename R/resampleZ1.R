resampleZ1 <- function(BE = beta.p, PH = phi.p, dataset = final, WIN = win, JJ = J,
                       RRc = R_centers, RR = R_clusters,
                       debug = TRUE, TIME = 30, DEATH = 1){

  # centerlessFingers <- lapply(dataset, function(x) list(x = x$x - x$centers[1], y = x$y - x$centers[2]))
  centerlessFingers <- lapply(dataset, function(x) sweep(x$fingers, 2, x$centers))
  m <- length(centerlessFingers)
  suff.mcmc <- matrix(0, nrow = JJ, ncol = 2)

  for (j in 1:JJ){

    centerP <- rpoispp(BE*TIME, win = WIN)
    n <- centerP$n
    centerP$births <- runif(n, 0, TIME)
    centerP$deaths <- centerP$births + rexp(n, DEATH)
    while(any(duplicated(c(centerP$births, centerP$deaths)))){
      centerP$births <- runif(n)
      centerP$deaths <- centerP$births + rexp(n, DEATH)
    }
    dominant.j <- vector("list", n)
    for(i in seq_len(n)){
      dominant.j[[i]]$centers <- c(centerP$x[i], centerP$y[i])
      dominant.j[[i]]$fingers <- sweep(centerlessFingers[[sample(seq_len(m), 1)]], 2, dominant.j[[i]]$centers, FUN = "+")
      dominant.j[[i]]$color <- "keep"
    }

    # final.j = delete fingers according to birth and death process
    # Maybe remove that portion from geradorVonmisesMCMC code and use it here?

    indexesFingers <- rep(seq_along(dominant.j),
                          sapply(dominant.j, function(x) nrow(x$fingers)))
    ramos <- lapply(dominant.j, function(x) x$fingers)
    ramosTotal <- cbind(unlist(sapply(ramos, function(x) x[,1])),
                        unlist(sapply(ramos, function(x) x[,2])))
    activeSet <- rep(c(TRUE, FALSE),
                     c(nrow(dominant.j[[1]]$fingers),
                       nrow(ramosTotal) - nrow(dominant.j[[1]]$fingers)))

    marks <- c(centerP$births, centerP$deaths)
    ocorr <- c(rep("birth", length(centerP$births)), rep("death", length(centerP$deaths)))
    ocorr <- ocorr[order(marks)]
    marks <- sort(marks)

    for (j2 in 2:(2*length(dominant.j))) {
      if (marks[j2] <= TIME) {
        if (ocorr[j2] == "death") {
          t1 <- which(centerP$deaths == marks[j2])
          if (dominant.j[[t1]]$color == "keep"){
            activeSet <- activeSet & !(indexesFingers == t1)
            dominant.j[[t1]]$color <- "delete"
          }
        }
        if (ocorr[j2] == "birth") {
          t3 <- which(centerP$births == marks[j2])
          if (sum(activeSet) > 0) { #!# SLOW
            num <- sum(rowSums(sweep(ramosTotal[activeSet, 1:2, drop = FALSE], 2, dominant.j[[t3]]$centers)^2) <= RRc^2)
            if (runif(1) > PH^(-num)) { # dominant.j[[t3]]$flag
              dominant.j[[t3]]$color <- "not born"
            } else {
              activeSet <- activeSet | (indexesFingers == t3)
            }
          } else {
            activeSet <- activeSet | (indexesFingers == t3)
          }
        }
      } else { # No need to check further
        break
      }
    }

    final.j <- dominant.j[sapply(dominant.j, function(x) x$color == "keep")]

    #1# Improve

    suff.mcmc[j,] <- sufficientStat(final.j, Rc = RRc)

    if(debug) {
      cat(paste0("\tj = ",j,", number of centers = ", suff.mcmc[j,1], ", number of fingers near centers by ", RRc, ": ", -1*suff.mcmc[j,2], ".\n"))
    }

  }

  return(suff.mcmc)

}
