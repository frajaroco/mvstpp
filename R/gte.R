#' @title Temporal mark variogram function
#' @description Computes an estimator of the temporal mark variogram function.
#' @param xyt Spatial coordinates and times \eqn{(x,y,t)} of the point pattern.
#' @param dt A vector of times \code{v} at which \eqn{\gamma_{te}(v)} is computed.
#' @param kt A kernel function for the temporal distances. The default is the \code{"box"} kernel. It can also be \code{"epanech"} for the Epanechnikov kernel, or \code{"biweight"}.
#' @param ht A bandwidth of the kernel function \code{kt}.
#' @details \eqn{N = {[x_{1},t_{1}],[x_{2},t_{2}],...,[x_{n},t_{n}]}} with \eqn{x_{i}} in \eqn{W} subset of \eqn{R^2} and \eqn{t_{i}} in \eqn{T} subset of \eqn{R}, the temporal mark variogram \eqn{\gamma_{te}(t)} is based on the stationary one-dimensional point process of the times \eqn{t_{i}}. This approach makes sense if the \eqn{x_{i}} come from a bounded spatial window \eqn{W}. The function \eqn{\gamma_{te}(t)} is given by the mean of the half-squared distances of the locations of points of time lag \eqn{t}, and the estimator is given by \deqn{\widehat{\gamma}_{te}(t) = sum_{i = 1,...,n} sum_{j = 1,...,n; j != i} (1/2 (||x_{1}-x_{2}||)^2 \kappa_{\delta}(|t_{i}-t_{j}|-t) / sum_{i = 1,...,n} sum_{j = 1,...,n; j != i} (\kappa_{\delta} (|t_{i}-t_{j}|-t) ),} where \eqn{\kappa} is a one-dimensional kernel function with bandwidth \eqn{\delta}.
#' @return A list containing:
#' \itemize{
#'   \item \code{gteke}: A vector containing the values of \eqn{\widehat{\gamma}_{te}(t)} estimated.
#'   \item \code{dt}: If \code{dt} is missing, a vector of distances \code{v} at which \eqn{\gamma_{te}(t)} is computed under the restriction \eqn{0<\delta<t}.
#'   \item \code{kernel}: A vector of names and bandwidth of the temporal kernel.
#'   }
#' @author Francisco J. Rodriguez-Cortes <cortesf@@uji.es> \url{https://fjrodriguezcortes.wordpress.com}
#' @references Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.
#' @references Chiu, S. N., Stoyan, D., Kendall, W. S., and Mecke, J. (2013). Stochastic Geometry and its Applications. John Wiley & Sons.
#' @references Gabriel, E., Rowlingson, B., Diggle P J. (2013) \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software 53, 1-29.
#' @references Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.
#' @examples
#' ## Not run:
#' #################
#'
#' # Realisations of the homogeneous Poisson processes
#' hpp <- rpp(lambda = 100, replace = FALSE)$xyt
#'
#' # This function provides an edge-corrected kernel estimator of the temporal mark variogram
#' out <- gte(hpp,ht=0.1/sqrt(100))
#'
#' # R plot - Temporal mark variogram
#' par(mfrow=c(1,1))
#' gte_range <- range(out$gteke,1/6)
#' plot(out$dt,out$gteke,type="l",ylim=gte_range,xlab="t = times",
#' ylab=expression(hat(gamma)[te](t)),main="Temporal mark variogram")
#' lines(out$dt,rep(1/6,length(out$dt)),type="l",col="red")
#'
#' ## End(Not run)
gte <- function(xyt, dt, kt="epanech", ht){

  if (missing(ht)){
    d <- dist(xyt[,3])
    ht <- dpik(d, kernel=kt, range.x=c(min(d),max(d)))}

  if (missing(dt)) {
    tregion <- range(xyt[,3],na.rm=TRUE)
    bsupt <- max(tregion)
    binft <- min(tregion)
    maxt <- (bsupt-binft)/4
    dt <- seq(0, maxt, len=50)
    dt <- sort(dt)}

  pts <- xyt[,1:2]
  xytimes <- xyt[,3]
  ptsx <- pts[,1]
  ptsy <- pts[,2]
  ptst <- xytimes
  nt <- length(ptst)
  ndt <- length(dt)
  gte <- rep(0,ndt)

  kernel=c(kt=kt,ht=ht)

  if (kt=="box"){ kt=1}
  else if (kt=="epanech"){ kt=2}
  else if (kt=="biweight"){ kt=3}

  storage.mode(gte) <- "double"

  gteout <- .Fortran("gtecore",as.double(ptsx),as.double(ptsy),as.double(ptst),as.integer(nt),as.double(dt),as.integer(ndt),as.integer(kt),as.double(ht),(gte))

  gte <- gteout[[9]]

  invisible(return(list(gteke=gte,dt=dt,kernel=kernel)))
}
