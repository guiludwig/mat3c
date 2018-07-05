sufficientStatAll <- function(proc.fmat, sigma, kappa,
                              A1 = a1, A2 = a2, B1 = b1, B2 = b2,
                              NN = N, Rc = R_centers, R = R_clusters){ #' R_c
  suff <- rep(0, 5)
  suff[1] <- length(proc.fmat) #' L
  if(suff[1] <= 1){
    return(suff)
  } else {
    for (i in 1:suff[1]){
      #!# Not used anywhere?
      # suff[2] #' sum_i sum_{j \neq i} sum_{z \in \eta_j}
      #!#  for (j in (1:suff[1])[-i]){
      #!#  if (length(proc.final[[j]]$fingers) > 0) {
      #!#  xx <- sweep(proc.final[[j]]$fingers, 2,
      #!#           proc.final[[i]]$centers)
      #!#  suff[2] <- suff[2] - sum(sqrt(apply(xx^2, 1, sum)) <= Rc)
      #!#  } else if (length(proc.final[[j]]$x) > 0) {
      #!#  suff[2] <- suff[2] - sum(((proc.final[[j]]$x - proc.final[[i]]$centers[1])^2 +
      #!#  (proc.final[[j]]$y - proc.final[[i]]$centers[2])^2) <= Rc^2)
      #!#  } else {
      #!#  stop("Cannot evaluate sufficient statistics, does object have $fingers?")
      #!#  }
      # suff[3]
      suff[3] <- suff[3] + length(proc.fmat[[i]]$x) #' sum_i #(eta_{x_i})
      # suff[4]
      mu <- proc.fmat[[i]]$centers
      X.st <- rnorm(NN, 0, 1)
      # note circular:::RvonmisesRad works on [0, 2pi]!
      temp1 <- circular:::RvonmisesRad(NN, proc.fmat[[i]]$direcao, kappa)
      X <- cbind(mu[1] + sigma*abs(X.st)*cos(temp1),
                 mu[2] + sigma*abs(X.st)*sin(temp1))

      d.s.1 <- !(X[, 1] > A1 & X[, 1] < A2 & X[, 2] > B1 & X[, 2] < B2)
      n.s.1 <- sum(d.s.1)
      prob.tot <- 1 - n.s.1/NN

      while(any(d.s.1)){ # While there are fingers outside the box...
        temp1 <- circular:::RvonmisesRad(n.s.1, proc.fmat[[i]]$direcao, kappa)
        temp2 <- rnorm(n.s.1, 0, sigma)
        X[d.s.1, 1] <- mu[1] + abs(temp2)*cos(temp1)
        X[d.s.1, 2] <- mu[2] + abs(temp2)*sin(temp1)
        d.s.1[d.s.1] <- !(X[d.s.1, 1] > A1 & X[d.s.1, 1] < A2 & X[d.s.1, 2] > B1 & X[d.s.1, 2] < B2)
        n.s.1 <- sum(d.s.1) #!# previous line indexing by 0-1s, while checking for logicals?
      }

      #!# suff[4] <- suff[4] + area.final(proc.fmat[[i]], X, sigma, kappa)
      list.fmat <- proc.fmat[[i]]
      suff[4] <- suff[4] + shadowcpp(NN, X, length(list.fmat$x), list.fmat$x, list.fmat$y,
                                     list.fmat$t, sigma, R)

      # suff[5]
      z1 <- proc.fmat[[i]]$x - proc.fmat[[i]]$centers[1]
      z2 <- proc.fmat[[i]]$y - proc.fmat[[i]]$centers[2]
      # radius <- sqrt(z1^2 + z2^2)
      radius <- z1^2 + z2^2
      angle <- atan2(z2, z1)
      angle[angle < 0] <- 2*pi + angle[angle < 0] # atan2 -> [-pi,pi] but we sample u on [0, 2pi]
      #        suff[5] <- suff[5] + sum(dvonmises(angle, proc.fmat[[i]]$direcao, kappa, log = TRUE) +
      #                                  log(2) + dnorm(radius, 0, sigma, log = TRUE) - log(prob.tot))
      #!# Why dnorm(radius)??? #!# p. 6 manuscript
      # suff[5] <- suff[5] + sum(circular:::DvonmisesRad(angle, proc.fmat[[i]]$direcao, kappa, log = TRUE)
      #                          + log(radius) + dnorm(radius, 0, sigma, log = TRUE) - log(prob.tot))
      #!# Should be the same thing tho, right? Constant = proportion fingers inside window
      # Z2 <- sum(outer(circular:::DvonmisesRad(seq(0, 2*pi - pi/20, pi/20), proc.fmat[[i]]$direcao, kappa),
      #                 dnorm(seq(-3,3,.1)*sigma, 0, sigma)))*(pi/20)*.1*sigma
      Z2 <- prob.tot
      suff[5] <- suff[5] + sum(circular:::DvonmisesRad(angle, proc.fmat[[i]]$direcao, kappa, log = TRUE)
                               + 0.5*(log(2/pi) + log(radius)) - log(sigma) - 0.5*radius/(sigma^2) - log(Z2)) #!# formula (7)
      # IDEA: Carry angle, direcao, Z2 in list, do evaluate likelihood manually
    } #!# log(prob.tot) inside sum?
  }
  return(suff)
}
