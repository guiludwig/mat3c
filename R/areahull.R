# https://web.archive.org/web/20100405070507/http://valis.cs.uiuc.edu/~sariel/research/CG/compgeom/msg00831.html
# Let 'vertices' be an array of N pairs (x,y), indexed from 0
# Let 'area' = 0.0
# for i = 0 to N-1, do
# Let j = (i+1) mod N
# Let area = area + vertices[i].x * vertices[j].y
# Let area = area - vertices[i].y * vertices[j].x
# end for
# Return 'area'
areahull <- function(X){
  ix <- chull(X)
  vertices <- X[ix,]
  N <- length(ix)
  if(N <= 2) return(0)
  area <- 0
  for(i in 1:N){
    j <- (i %% N) + 1
    area <- area + vertices[i,1]*vertices[j,2] - vertices[i,2]*vertices[j,1]
  }
  return(abs(area)/2)
}