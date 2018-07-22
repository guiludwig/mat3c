geradorVonmisesMCMC <- function(beta, phi, gamma, sigma, kappa,
                                win, Rc = R_centers, R = R_clusters,
                                TIME = 30, DEATH = 1, completeObs = FALSE){

  pp_germ <- rpoispp(beta*TIME, win = win) # possible location of centers
  while(pp_germ$n == 0){
    pp_germ <- rpoispp(beta*TIME, win = win)
  }

  pp_germ$marks <- runif(pp_germ$n, 0, TIME) # possible birth times

  ordem_c <- order(pp_germ$marks)
  vida <- rexp(pp_germ$n, DEATH) # vida is the duration of the event
  nasc <- pp_germ$marks[ordem_c] # sort(pp_germ$marks)
  flag <- runif(pp_germ$n) #!# Do we need to store flag?
  morte <- nasc + vida
  while(any(duplicated(c(nasc, morte)))){ #!# seed 1234, iteration 2, has -> # 1638 1639
    pp_germ$marks <- runif(pp_germ$n, 0, TIME) # possible birth times
    vida <- rexp(pp_germ$n, DEATH) # vida is the duration of the event
    nasc <- pp_germ$marks[ordem_c] # sort(pp_germ$marks)
    flag <- runif(pp_germ$n)
    morte <- nasc + vida
  }
  u <- runif(pp_germ$n, 0, 2*pi) # sample(angulos, pp_germ$n, repl = TRUE)

  dominant <- vector("list", pp_germ$n) # list()

  for (i in 1:pp_germ$n){
    dominant[[i]] <- list(centers = c(pp_germ$x[ordem_c[i]], pp_germ$y[ordem_c[i]]),
                          nasc = nasc[i],
                          morte = morte[i],
                          flag = flag[i],
                          direction = u[i],
                          color = "keep")
  }

  marcas <- c(nasc, morte)
  ocurr <- c(rep("birth", pp_germ$n), rep("death", pp_germ$n))
  ocorr <- ocurr[order(marcas)] #!# note ocorr != ocurr
  marcas <- sort(marcas)

  xf <- rep(0, length(dominant))
  yf <- rep(0, length(dominant))

  # In this step, generates the fingers for all centers
  for (j2 in seq_along(dominant)) {

    mux <- dominant[[j2]]$centers[1]
    muy <- dominant[[j2]]$centers[2]
    xf[j2] <- mux
    yf[j2] <- muy
    dir <- dominant[[j2]]$direction

    marked.locations <- geradorVonmises(dir, gamma, c(mux, muy), kappa, sigma, win)
    marked.n <- nrow(marked.locations)

    marked.times <- runif(marked.n)
    ordem <- order(marked.times)
    fingers0 <- matrix(c(marked.locations[ordem, 1],
                         marked.locations[ordem, 2],
                         marked.times[ordem]),
                       nrow = marked.n, ncol = 3)

    ind.fingers <- nextFingers(marked.n, fingers0, R)
    dominant[[j2]]$fingers <- fingers0[ind.fingers, 1:3, drop = FALSE]

    if(completeObs) dominant[[j2]]$thinned <- fingers0[!ind.fingers, 1:3, drop = FALSE]

  }

  # Now, we start deleting centers as a function of fingers from older, active centers
  # For each element i in dominant, creates j replicates of i, where j is the number of fingers
  indexesFingers <- rep(seq_along(dominant),
                        sapply(dominant, function(x) nrow(x$fingers)))
  # Gets the fingers' matrix for each element of dominant
  ramos <- lapply(dominant, function(x) x$fingers)
  # This is a nfingers x 4 matrix, each row is a finger, ordered by center in dominant then
  # by time order, columns represent: position x, position y, time t, corresponding center i.
  ramosTotal <- cbind(unlist(sapply(ramos, function(x) x[,1])),
                      unlist(sapply(ramos, function(x) x[,2])))
  activeSet <- rep(c(TRUE, FALSE),
                   c(nrow(dominant[[1]]$fingers),
                     nrow(ramosTotal) - nrow(dominant[[1]]$fingers)))

  for (j2 in 2:(2*length(dominant))) {
    if (marcas[j2] <= TIME) {
      if (ocorr[j2] == "death") {
        t1 <- which(morte == marcas[j2])
        if (dominant[[t1]]$color == "keep"){
          activeSet <- activeSet & !(indexesFingers == t1)
          dominant[[t1]]$color <- "delete"
        }
      }
      if (ocorr[j2] == "birth") {
        t3 <- which(nasc == marcas[j2])
        if (sum(activeSet) > 0) { #!# SLOW
          num <- sum(rowSums(sweep(ramosTotal[activeSet, 1:2, drop = FALSE], 2, dominant[[t3]]$centers)^2) <= Rc^2)
          if (dominant[[t3]]$flag > phi^(-num)) {
            dominant[[t3]]$color <- "not born"
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

  return(dominant)

}
