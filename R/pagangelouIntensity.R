#' Papangelou Intensity
#'
#' Generates Papangelou conditional intensity function for points in a quadrature.
pagangelouIntensity <- function(xi, beta, phi, fmat, R_c = R_centers, Boot = 1000){

  kern1 <- 1 # Sum dist(\xi, existing fingers)
  kern2 <- 1 # average{ Sum_i dist(x_i, bootstrapped fingers on \xi) }
  which.m <- sample(1:length(fmat), Boot, replace = TRUE)
  m <- lapply(which.m, function(k) list(x = fmat[[k]]$x,
                                        y = fmat[[k]]$y))

  return(beta*exp(log(phi)*(-kern1-kern2)))

}
