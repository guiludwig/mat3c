minDist <- function(x){
  if(!is.null(x$fingers)){
    d <- dist(x$fingers)
  } else {
    d <- dist(cbind(x$x, x$y))
  }
  if(length(d) == 0){
    return(0)
  } else {
    return(min(d))
  }
}
