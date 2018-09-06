#' Sample from Matern-III marked point process
#'
#' Generates a sample from Matern-III marked spatial point process
#'
#' @param beta Intensity of the dominant Poisson process.
#' @param phi Centers' inhibition parameter; rate in which centers are deleted based
#'            on how many fingers from other centers surround it.
#' @param gamma Average number of fingers generated.
#' @param sigma Length (in standard deviations, according to half-Normal distribution)
#'              of fingers.
#' @param kappa Concentration parameter from von Mises distribution.
#' @param win Square window in which the process is to be sampled. See \code{\link{owin}}.
#' @param R_clusters Hardcore radius of inhibition for finger process.
#' @param R_centers Radius of inhibition for center process.
#' @param time Length of time for which a birth-and-death process will be sampled.
#'             Defaults to 30.
#' @param death Rate in which the birth-and-death process centers will die. Defaults
#'              to 1.
#' @param returnTimes Whether the data should be returned with birth times, or just
#'                    the sampled spatial locations. Defaults to FALSE. Birth times
#'                    are assumed unknown in sampled data, so the birth times should
#'                    be used for debugging only.
#' @param formatData Saves the data as text files in the current directory, with the
#'                   corresponding file names. Set to NULL to override.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @import circular spatstat
#' @return \code{rmat3} returns an object of \code{\link{class}} "mat3" containing
#' at least the following components
#'   \item{centers}{The spatial point pattern.}
#'   \item{fingers}{The marks.}
#'
#' @examples
#' set.seed(1234)
#' x <- rmat3(70, 2, 5, 0.05, 3)
#' plot(x)
#' z <- sapply(x, function(x) x$centers)
#' plot(Kest(ppp(z[1,],z[2,])))
#' plot(envelope(ppp(z[1,],z[2,])))
#' # Saves data to  files "centers.txt" and "fingers.txt"
#' \dontrun{set.seed(1234)
#' x <- rmat3(70, 2, 5, 0.05, 3, formatData = c("centers.txt", "fingers.txt"))
#' # When phi < 1, uses MCMC
#' x <- rmat3(70, 1/2, 5, 0.05, 3)
#' }
#'
#' @author Guilherme Ludwig and Nancy Garcia
#'
#' @references
#'
#'   Garcia, N., Guttorp, P. and Ludwig, G. (2018) TBD
#'
#' @seealso \code{\link{fitmat3}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
rmat3 <- function(beta, phi, gamma, sigma, kappa,
                  win = owin(c(0,1), c(0,1)),
                  R_clusters = 0.005, R_centers = 0.02,
                  time = 30, death = 1,
                  returnTimes = FALSE, returnThinned = FALSE,
                  formatData = NULL, ...){

  if(returnThinned) returnTimes = TRUE

  if(phi < 1){

    # Initial configuration:
    attempt  <- rmat3(beta, 1, gamma, sigma, kappa,
                      win, R_clusters, R_centers,
                      time, death,
                      returnTimes, returnThinned)
    # attempt <- runifpoint(1, win = win)

    for(i in 1:max(10*length(attempt), 1000)){ # Illian et al. suggestion: 10*n

      print(i)
      add <- sample(c(TRUE, FALSE), 1)

      if(add){
        x_new <- runifpoint(1, win = win)
        m <- length(attempt)

        # Test inclusion

        # If the point pattern is empty, add automatically
        if(m >= 1){
          nearbyFinger <- 0
          for(i in seq_along(attempt)){
            nearbyFinger <- nearbyFinger + sum(rowSums(sweep(attempt[[i]]$fingers,
                                                             2,
                                                             c(x_new$x, x_new$y))^2) < R_centers)
          }

          rho_birth <- area(win)*beta*exp(-log(phi)*nearbyFinger)/(m+1)
          alpha_birth <- min(rho_birth, 1)
          if(runif(1) > alpha_birth) next() # Skip iteration if failure
        }

        candidate <- list(centers = c(x_new$x, x_new$y))

        # If successful, adds fingers:

        mux <- x_new$x
        muy <- x_new$y
        dir <- runif(1, 0, 2*pi)

        marked.locations <- geradorVonmises(dir, gamma, c(mux, muy), kappa, sigma, win)
        marked.n <- nrow(marked.locations)

        marked.times <- runif(marked.n)
        ordem <- order(marked.times)
        fingers0 <- matrix(c(marked.locations[ordem, 1],
                             marked.locations[ordem, 2],
                             marked.times[ordem]),
                           nrow = marked.n, ncol = 3)

        ind.fingers <- nextFingers(marked.n, fingers0, R_clusters)

        if(returnTimes) {
          candidate$fingers <- fingers0[ind.fingers, 1:3, drop = FALSE]
          if(returnThinned) candidate$thinned <- fingers0[!ind.fingers, 1:3, drop = FALSE]
        } else {
          candidate$fingers <- fingers0[ind.fingers, 1:2, drop = FALSE]
        }


        attempt <- append(attempt, list(candidate))

      } else {
        m <- length(attempt)
        k <- sample(seq_along(attempt), 1)

        # Test deletion

        # If the point pattern is empty, skip iteration
        if(length(attempt) == 0){
          next()
        }
        nearbyFinger <- 0
        for(i in seq_along(attempt)[-k]){
          nearbyFinger <- nearbyFinger + sum(rowSums(sweep(attempt[[i]]$fingers,
                                                           2,
                                                           attempt[[k]]$centers)^2) < R_centers)
        }

        rho_death <- (m-1)*exp(log(phi)*nearbyFinger)/(area(win)*beta)
        alpha_death <- min(rho_death, 1)

        # If successful
        if(runif(1) < alpha_death) {
          attempt[[k]] <- NULL
        }

      }

    }

    class(attempt)  <- "mat3"
    return(attempt)

  }

  limx <- win$xrange
  limy <- win$yrange

  dominant <- geradorVonmisesMCMC(beta, phi, gamma, sigma, kappa, win,
                                  R_centers, R_clusters, time, death,
                                  completeObs = returnThinned)

  if(returnTimes){
    fmat <- vector("list", sum(sapply(dominant, function(x) x$color == "keep")))
  } else {
    final <- vector("list", sum(sapply(dominant, function(x) x$color == "keep")))
  }

  k <- 1
  for(i in 1:length(dominant)) {
    temp <- NULL
    if(dominant[[i]]$color == "keep") {
      if(length(dominant[[i]]$fingers) > 0){
        temp <- matrix(matrix(dominant[[i]]$fingers, ncol = 3)[, 1:2], ncol = 2)
        ordem <- order(matrix(dominant[[i]]$fingers, ncol = 3)[, 3])
        if(returnThinned){
          fmat[[k]] <- list(centers = c(dominant[[i]]$centers[1], dominant[[i]]$centers[2]), x = temp[ordem, 1], y = temp[ordem, 2], t = sort(matrix(dominant[[i]]$fingers, ncol = 3)[, 3]), direcao = dominant[[i]]$direction, thinned = dominant[[i]]$thinned) # complete process (centers, fingers) + directions + times of birth
        } else if(returnTimes){
          fmat[[k]] <- list(centers = c(dominant[[i]]$centers[1], dominant[[i]]$centers[2]), x = temp[ordem, 1], y = temp[ordem, 2], t = sort(matrix(dominant[[i]]$fingers, ncol = 3)[, 3]), direcao = dominant[[i]]$direction) # complete process (centers, fingers) + directions + times of birth
        } else {
          final[[k]] <- list(centers = c(dominant[[i]]$centers[1], dominant[[i]]$centers[2]), fingers = temp) # final process of centers and fingers
        }
        k <- k+1
      }
    }
  }

  if(!is.null(formatData)){

    if(returnTimes){
      stop("set returnTimes to FALSE to save data.")
    }

    data <- NULL
    for (i in 1:length(final)){
      data <- rbind(data, c(i, final[[i]]$centers))
    }
    write.table(data, file = formatData[1], col.names = c("Tree", "X1", "Y1"), quote = FALSE)

    data <- NULL
    for (i in 1:length(final)){
      data <- rbind(data, cbind(i, final[[i]]$fingers))
    }
    write.table(data, file = formatData[2], col.names = c("Tree", "X1", "Y1"), quote = FALSE)
  }

  if(returnTimes){
    class(fmat) <- c("mat3", "complete")
    return(fmat)
  } else {
    class(final) <- c("mat3")
    return(final)
  }

}
