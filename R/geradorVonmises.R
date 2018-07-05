geradorVonmises <- function(u, gamma, centers, kappa, sigma, win){

  limx <- win$xrange
  limy <- win$yrange
  A1 <- limx[1]
  A2 <- limx[2]
  B1 <- limy[1]
  B2 <- limy[2]

  mux <- centers[1]
  muy <- centers[2]

  n <- rpois(1, gamma)
  while (n == 0) {
    n <- rpois(1, gamma)
  }

  data1 <- circular:::RvonmisesRad(n, u, kappa)
  data2 <- rnorm(n, 0, sigma)

  x <- mux + abs(data2)*cos(data1)
  y <- muy + abs(data2)*sin(data1)

  d.x <- !(x > A1 & x < A2 & y > B1 & y < B2)
  n.s <- sum(d.x)

  while(n.s > 0) {

    # rvonmises(n.s, u, kappa, control.circular = list(modulo = "2pi"))
    data1 <- circular:::RvonmisesRad(n.s, u, kappa)
    data2 <- rnorm(n.s, 0, sigma)

    x[d.x] <- mux + abs(data2)*cos(data1)
    y[d.x] <- muy + abs(data2)*sin(data1)

    d.x[d.x] <- !(x[d.x] > A1 & x[d.x] < A2 & y[d.x] > B1 & y[d.x] < B2)
    n.s <- sum(d.x)

  }

  return(cbind(x, y))

}
