#' @title Spatio-temporal hot-spots cluster point process model
#' @description Generate a random spatio-temporal point pattern, a simulated realisation of the cluster process.
#' @param lambp Intensity of the Poisson process of cluster centres. A single positive number, a function, or a pixel image.
#' @param r Radius parameter of the clusters
#' @param mu Mean number of points per cluster (a single positive number) or reference intensity for the cluster points (a function or a pixel image).
#' @param s.region A two-column matrix specifying a polygonal region containing all data locations. If \code{s.region} is missing, the Ripley-Rasson estimate convex spatial domain is considered.
#' @param t.region A vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.
#' @details This function generates a realisation of spatio-temporal cluster process, which can be considered as generalisation of the classical Matern cluster process, inside the spatio-temporal window.
#'
#' Consider a Poisson point process in the plane with intensity \eqn{\lambda_{p}} as cluster centres for all times “parent”, as well as a infinite cylinder of radius \eqn{R} around of each Poisson point, orthogonal to the plane. The scatter uniformly in all cylinders of all points which are of the form \eqn{(x,y,z)}, the number of points in each cluster being random with a Poisson (\eqn{\mu}) distribution. The resulting point pattern is a spatio-temporal cluster point process with \eqn{t=z}. This point process has intensity \eqn{\lambda_{p} * \mu}.
#' @return The simulated spatio-temporal point pattern.
#' @author Francisco J. Rodriguez-Cortes <cortesf@@uji.es> \url{https://fjrodriguezcortes.wordpress.com}
#' @references Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.
#' @references Gabriel, E., Rowlingson, B., Diggle P J. (2013) \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software 53, 1-29.
#' @references Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.
#' @examples
#' ## Not run:
#' #################
#'
#' require(scatterplot3d)
#'
#' # homogeneous
#' X <- sthpcpp(lambp=20, r=0.05, mu=100)
#'
#' plot(X$xyt)
#' par(mfrow=c(1,1))
#' scatterplot3d(X$xyt[,1],X$xyt[,2],X$xyt[,3],xlab="x",ylab="y",zlab="t")
#'
#' ## End(Not run)
sthpcpp <- function(lambp, r, mu, s.region, t.region){

  if (missing(s.region)) {s.region <- matrix(c(1,1,0,0,0,1,1,0),ncol=2)}
  if (missing(t.region)) {t.region <- c(0,1)}

  bdry <- owin(poly=list(x=s.region[,1],y=s.region[,2]))
  spp <- rMatClust(lambp, r, mu, win=bdry)
  tpp <- runif(spp$n,t.region[1],t.region[2])
  stc1 <- cbind(spp$x,spp$y,runif(spp$n,t.region[1],t.region[2]))
  stc1 <- as.3dpoints(stc1)

  invisible(return(list(xyt=stc1,s.region=s.region,t.region=t.region)))
}
