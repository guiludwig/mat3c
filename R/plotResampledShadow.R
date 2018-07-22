#' Plot resampled thinned events from shadow
#'
#' This is a debugging function since, by definition, thinned events are not observed in data.
#'
#' @export
#' @examples
#' set.seed(1234)
#' x <- rmat3(70, 2, 5, 0.05, 3, returnTimes = TRUE, returnThin = TRUE)
#' plotResampledShadow(x, 1, 1)
plotResampledShadow <- function(fmat.cond, which.center = 1, which.finger = 1,
                                which.proj = c("x", "y"),
                                R_clusters = 0.005,
                                R_centers = 0.02){

  if(is.null(fmat.cond[[1]]$t)) stop("This function is for debugging purposes, and will only work if an mat3 object has returnTimes = TRUE.")
  which.proj <- match.arg(which.proj)
  probe <- fmat.cond[[which.center]]
  temp <- resampleTimes(fmat.cond, R = R_clusters)[[which.center]] #
  allData <- with(probe, rbind(cbind(x, y, t), thinned))
  n <- length(probe$y)
  if(!n) stop("No fingers in the selected center")
  # layout(matrix(1:2, ncol = 2))
  # for(i in 1:n){ # projetada
  #   index <- which((allData[,1] - allData[i,1])^2 + (allData[,2] - allData[i,2])^2 < Rc^2)
  #   if(i == 1){
  #     switch(which.proj,
  #            x = plot(NA, ylab = "t", xlab = "x", ylim = c(0,1), xlim = range(allData[,1]),
  #                     main = "Projection in y"),
  #            y = plot(NA, ylab = "t", xlab = "y", ylim = c(0,1), xlim = range(allData[,2]),
  #                     main = "Projection in x"))
  #   }
  #   colIndex <- ifelse(which.proj == "x", 1, 2)
  #   points(allData[index, colIndex], allData[index, 3])
  #   points(allData[i, colIndex], allData[i,3], pch = 20)
  #   polygon(c(allData[i,colIndex] - Rc, allData[i,colIndex] - Rc,
  #             allData[i,colIndex] + Rc, allData[i,colIndex] + Rc),
  #           c(1, allData[i,3], allData[i,3], 1))
  #   lines(c(allData[i,colIndex], temp$x[which(temp$x == allData[i,colIndex])]),
  #         c(allData[i,3], temp$t[which(temp$x == allData[i,1])]), lty = 2, col = "Red")
  #   if(which.proj == "x") {
  #     lines(c(temp$x[i] - Rc, temp$x[i] + Rc), c(temp$t[i], temp$t[i]), lty = 2, col = "Red")
  #   } else {
  #     lines(c(temp$y[i] - Rc, temp$y[i] + Rc), c(temp$t[i], temp$t[i]), lty = 2, col = "Red")
  #   }
  # }
  th <- seq(0, 2*pi, length.out=100)
  Rc <- R_centers
  R <- R_clusters
  fingersDeleted <- nrow(allData) > length(probe$x)
  # for(i in 1:n){ # projected
  i <- which.finger
  # Which column
  colIndex <- ifelse(which.proj == "x", 1, 2)
  # which points are near finger i
  index <- which((allData[,1] - allData[i,1])^2 + (allData[,2] - allData[i,2])^2 < R^2)
  # # which fingers have radii intersections with finger i's radius
  # indexCylinder <- ((1:n)[-i])[(allData[i,1] - allData[(1:n)[-i],1])^2 +
  #                                (allData[i,2] - allData[(1:n)[-i],2])^2 < (2*R)^2]
  # which fingers lay near the same line
  indexCylinder <- switch(which.proj,
                          x = which((allData[i,2] - allData[,2])^2 < R^2),
                          y = which((allData[i,1] - allData[,1])^2 < R^2))
  # highlight <- c("black", "darkgrey")[2 - 1:n %in% c(i, indexCylinder)]
  if(fingersDeleted){
    highlight <- c(rep("purple", nrow(allData)), rep("darkgrey", nrow(allData) - length(probe$x)))
  } else {
    highlight <- rep("purple", nrow(allData))
  }

  plot(probe$x, probe$y,
       xlim = range(c(allData[,1] - R,
                      allData[,1] + R,
                      probe$x - Rc,
                      probe$x + Rc)),
       ylim = range(c(allData[,2] - R,
                      allData[,2] + R,
                      probe$y - Rc,
                      probe$y + Rc)),
       pch = 20, xlab = "x", ylab = "y",
       col = highlight)
  for(k in 1:n){
    lines(probe$x[k] + R*cos(th), probe$y[k] + R*sin(th), col = highlight[k])
  }
  switch(which.proj,
         x = abline(h = allData[i,2]),
         y = abline(v = allData[i,1]))
  points(probe$centers[1], probe$centers[2], pch = 20, cex = 1.5, col = "Black")
  lines(probe$centers[1] + Rc*cos(th), probe$centers[2] + Rc*sin(th), col = "Black")
  segments(probe$centers[1], probe$centers[2],
           probe$x, probe$y, "Purple")
  if(fingersDeleted){
    for(k in seq_len(nrow(allData) - length(probe$x))){
      points(allData[k+length(probe$x), 1],
             allData[k+length(probe$x), 2], col = "Darkgrey")
      lines(allData[k+length(probe$x), 1] + R*cos(th),
            allData[k+length(probe$x), 2] + R*sin(th), col = "Darkgrey")
      segments(probe$centers[1], probe$centers[2],
               allData[k+length(probe$x), 1],
               allData[k+length(probe$x), 2], col = "Darkgrey")
    }
  }

  switch(which.proj,
         x = plot(NA, ylab = "time", xlab = "x", ylim = c(0,1),
                  xlim = range(c(allData[,1] - R,
                                 allData[,1] + R,
                                 probe$x - Rc,
                                 probe$x + Rc)),
                  main = "Slice along x"),
         y = plot(NA, ylab = "time", xlab = "y", ylim = c(0,1),
                  xlim = range(c(allData[,2] - R,
                                 allData[,2] + R,
                                 probe$y - Rc,
                                 probe$y + Rc)),
                  main = "Slice along y"))
  points(allData[index, colIndex], allData[index, 3])
  points(allData[i, colIndex], allData[i,3], pch = 20)
  polygon(c(allData[i,colIndex] - R, allData[i,colIndex] - R,
            allData[i,colIndex] + R, allData[i,colIndex] + R),
          c(1, allData[i,3], allData[i,3], 1))
  if(which.proj == "x") {
    newIndex <- which(temp$x == allData[i,colIndex])
    lines(c(allData[i,colIndex], temp$x[newIndex]),
          c(allData[i,3], temp$t[newIndex]), lty = 2, col = "Red")
    lines(c(temp$x[newIndex] - R, temp$x[newIndex] + R),
          c(temp$t[newIndex], temp$t[newIndex]), lty = 2, col = "Red")
  } else {
    newIndex <- which(temp$y == allData[i,colIndex])
    lines(c(allData[i,colIndex], temp$y[newIndex]),
          c(allData[i,3], temp$t[newIndex]), lty = 2, col = "Red")
    lines(c(temp$y[newIndex] - R, temp$y[newIndex] + R),
          c(temp$t[newIndex], temp$t[newIndex]), lty = 2, col = "Red")
  }
  for(j in indexCylinder){
#
#     exclusively <- which((allData[,1] - allData[j,1])^2 + (allData[,2] - allData[j,2])^2 < R^2)
#     exclusively <- exclusively[!exclusively %in% index]
#     points(allData[exclusively, colIndex], allData[exclusively, 3], pch = 4)

    points(allData[j, colIndex], allData[j,3], pch = ifelse(j > length(probe$x), 1, 20))
    # varRad <- sqrt(R^2 - (sqrt(sum((allData[i,1:2] - allData[j,1:2])^2)) - R)^2)
    varRad <- switch(which.proj,
                     x = sqrt(R^2 - (sqrt(sum((allData[i,2] - allData[j,2])^2)) - R)^2),
                     y = sqrt(R^2 - (sqrt(sum((allData[i,1] - allData[j,1])^2)) - R)^2))
    polygon(c(allData[j,colIndex] - varRad, allData[j,colIndex] - varRad,
              allData[j,colIndex] + varRad, allData[j,colIndex] + varRad),
            c(1, allData[j,3], allData[j,3], 1), lty = 2)

    if(which.proj == "x") {
      newIndex <- which(temp$x == allData[j,colIndex])
      lines(c(allData[j,colIndex], temp$x[newIndex]),
            c(allData[j,3], temp$t[newIndex]), lty = 2, col = "Red")
      lines(c(temp$x[newIndex] - varRad, temp$x[newIndex] + varRad),
            c(temp$t[newIndex], temp$t[newIndex]), lty = 2, col = "Red")
    } else {
      newIndex <- which(temp$y == allData[j,colIndex])
      lines(c(allData[j,colIndex], temp$y[newIndex]),
            c(allData[j,3], temp$t[newIndex]), lty = 2, col = "Red")
      lines(c(temp$y[newIndex] - varRad, temp$y[newIndex] + varRad),
            c(temp$t[newIndex], temp$t[newIndex]), lty = 2, col = "Red")
    }

  }
  # }
  invisible()
}
