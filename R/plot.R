#' @export
plot.mat3 <- function(dominant, scol = "purple", slwd = 1, ...){
  if(any(class(dominant) == "complete")){
    indexes <- sapply(dominant, function(x) x$color == "keep")
  } else {
    indexes <- rep(TRUE, length(dominant))
    tagged <- FALSE
  }
  xf <- sapply(dominant[indexes], function(x) x$centers[1])
  yf <- sapply(dominant[indexes], function(x) x$centers[2])
  plot(xf, yf, ...)

  for (i in 1:length(dominant)) {
    if(tagged){
      if (dominant[[i]]$color == "keep") {
        temp <- rbind(dominant[[i]]$fingers[, 1:2])
        if (length(temp) > 0) {
          # s <- seq(1:(length(temp)/2))
          segments(dominant[[i]]$centers[1], dominant[[i]]$centers[2],
                   temp[, 1], temp[, 2], col = scol, lwd = slwd)
        }
      }
    } else {
      temp <- rbind(dominant[[i]]$fingers[, 1:2])
      if (length(temp) > 0) {
        segments(dominant[[i]]$centers[1], dominant[[i]]$centers[2],
                 temp[, 1], temp[, 2], col = scol, lwd = slwd)
      }
    }
  }
}
# Incorporate this later!
# plot(fmat[[ij]]$centers[1], fmat[[ij]]$centers[2],
#      xlim = range(c(fmat[[ij]]$centers[1], fmat[[ij]]$x)),
#      ylim = range(c(fmat[[ij]]$centers[2], fmat[[ij]]$y)))
# segments(fmat[[ij]]$centers[1], c(fmat[[ij]]$centers[2]),
#          fmat[[ij]]$x, fmat[[ij]]$y)
