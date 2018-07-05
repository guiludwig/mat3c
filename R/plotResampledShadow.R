plotResampledShadow <- function(fmat.cond, which.center = 1, which.finger = 1,
                                which.proj = c("x", "y"),
                                Rc = R_clusters,
                                R = R_centers){
  which.proj <- match.arg(which.proj)
  probe <- fmat.cond[[which.center]]
  temp <- resampleTimesDebug(fmat.cond)[[which.center]] #
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
  layout(matrix(1:2, ncol = 2))
  # for(i in 1:n){ # projected
  i <- which.finger
    # Which column
    colIndex <- ifelse(which.proj == "x", 1, 2)
    # which points are near finger i
    index <- which((allData[,1] - allData[i,1])^2 + (allData[,2] - allData[i,2])^2 < Rc^2)
    # which fingers have radii intersections with finger i's radius
    indexCylinder <- ((1:n)[-i])[(allData[i,1] - allData[(1:n)[-i],1])^2 +
                                   (allData[i,2] - allData[(1:n)[-i],2])^2 < (2*Rc)^2]
    highlight <- c("black", "darkgrey")[2 - 1:n %in% c(i, indexCylinder)]
    plot(probe$x, probe$y,
         xlim = range(c(allData[,1] - Rc,
                        allData[,1] + Rc)),
         ylim = range(c(allData[,2] - Rc,
                        allData[,2] + Rc)),
         pch = 20, xlab = "x", ylab = "y",
         col = highlight)
    lines(probe$centers[1] + R*cos(th), probe$centers[2] + R*sin(th), col = "Red")
    for(k in 1:n){
      lines(probe$x[k] + Rc*cos(th), probe$y[k] + Rc*sin(th), col = highlight[k])
    }
    switch(which.proj,
           x = abline(h = allData[i,2]),
           y = abline(v = allData[i,1]))
    points(probe$centers[1], probe$centers[2], pch = 20, cex = 1.5, col = "Red")
    # indexCylinder <- switch(which.proj,
    #                         x = ((1:n)[-i])[abs(allData[i,1] - allData[(1:n)[-i],1]) < 2*Rc],
    #                         y = ((1:n)[-i])[abs(allData[i,2] - allData[(1:n)[-i],2]) < 2*Rc])
    switch(which.proj,
           x = plot(NA, ylab = "t", xlab = "x", ylim = c(0,1),
                    xlim = range(c(allData[,1] - Rc,
                                   allData[,1] + Rc)),
                    main = "Projection in y"),
           y = plot(NA, ylab = "t", xlab = "y", ylim = c(0,1),
                    xlim = range(c(allData[,2] - Rc,
                                   allData[,2] + Rc)),
                    main = "Projection in x"))
    points(allData[index, colIndex], allData[index, 3])
    points(allData[i, colIndex], allData[i,3], pch = 20)
    polygon(c(allData[i,colIndex] - Rc, allData[i,colIndex] - Rc,
              allData[i,colIndex] + Rc, allData[i,colIndex] + Rc),
            c(1, allData[i,3], allData[i,3], 1))
    lines(c(allData[i,colIndex], temp$x[which(temp$x == allData[i,colIndex])]),
          c(allData[i,3], temp$t[which(temp$x == allData[i,1])]), lty = 2, col = "Red")
    if(which.proj == "x") {
      lines(c(temp$x[i] - Rc, temp$x[i] + Rc), c(temp$t[i], temp$t[i]), lty = 2, col = "Red")
    } else {
      lines(c(temp$y[i] - Rc, temp$y[i] + Rc), c(temp$t[i], temp$t[i]), lty = 2, col = "Red")
    }
    for(j in indexCylinder){

      exclusively <- which((allData[,1] - allData[j,1])^2 + (allData[,2] - allData[j,2])^2 < Rc^2)
      exclusively <- exclusively[!exclusively %in% index]
      points(allData[exclusively, colIndex], allData[exclusively, 3], pch = 4)

      # points(allData[j, colIndex], allData[j,3], pch = 20)
      varRad <- sqrt(Rc^2 - (sqrt(sum((allData[i,1:2] - allData[j,1:2])^2)) - Rc)^2)
      polygon(c(allData[j,colIndex] - varRad, allData[j,colIndex] - varRad,
                allData[j,colIndex] + varRad, allData[j,colIndex] + varRad),
              c(1, allData[j,3], allData[j,3], 1), lty = 2)
      lines(c(allData[j,colIndex], temp$x[which(temp$x == allData[j,colIndex])]),
            c(allData[j,3], temp$t[which(temp$x == allData[j,1])]), lty = 2, col = "Red")
      if(which.proj == "x") {
        lines(c(temp$x[j] - varRad, temp$x[j] + varRad), c(temp$t[j], temp$t[j]), lty = 2, col = "Red")
      } else {
        lines(c(temp$y[j] - varRad, temp$y[j] + varRad), c(temp$t[j], temp$t[j]), lty = 2, col = "Red")
      }
    }
    Sys.sleep(1)
  # }
  layout(1)
}
