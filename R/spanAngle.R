spanAngle <- function(x){
  z <- x
  if(!is.null(x$fingers)){
    z$x <- x$fingers[,1]
    z$y <- x$fingers[,2]
  }
  z$x <- z$x - x$center[1]
  z$y <- z$y - x$center[2]
  if(is.null(z$x)) return(0)
  if(length(z$x) == 1) return(0)
  theta <- atan2(z$y, z$x) + pi
  if(length(z$x) == 2) return(as.numeric(abs(theta[2]-theta[1])))
  theta <- sort(theta)
  N <- length(theta)
  d <- numeric(N)
  d[1:(N-1)] <- diff(theta)
  d[N] <- 2*pi - theta[N] + theta[1]
  return(2*pi - max(d))
}
