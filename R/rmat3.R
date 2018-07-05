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
#' x <- rmat3(70, 2, 5, 0.05, 3, formatData = c("centers.txt", "fingers.txt"))}
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
                  time = 30, death = 1, returnTimes = FALSE,
                  formatData = NULL, ...){

  limx <- win$xrange
  limy <- win$yrange

  dominant <- geradorVonmisesMCMC(beta, phi, gamma, sigma, kappa, win,
                                  R_centers, R_clusters, time, death)

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
        if(returnTimes){
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
