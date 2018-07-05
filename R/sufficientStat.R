sufficientStat <- function(proc.final, Rc = R_centers){
  suff <- numeric(2)
  suff[1] <- length(proc.final)
  if(suff[1] <= 1){
    return(suff)
  } else {
    for (i in 1:suff[1]){
      for (j in (1:suff[1])[-i]){
        if (length(proc.final[[j]]$fingers) > 0) {
          xx <- sweep(proc.final[[j]]$fingers, 2,
                      proc.final[[i]]$centers)
          suff[2] <- suff[2] - sum(sqrt(apply(xx^2, 1, sum)) <= Rc)
        } else if (length(proc.final[[j]]$x) > 0) {
          suff[2] <- suff[2] - sum(((proc.final[[j]]$x - proc.final[[i]]$centers[1])^2 +
                                      (proc.final[[j]]$y - proc.final[[i]]$centers[2])^2) <= Rc^2)
        } else {
          stop("Cannot evaluate sufficient statistics, does object have $fingers?")
        }
      }
    }
  }
  return(suff)
}
