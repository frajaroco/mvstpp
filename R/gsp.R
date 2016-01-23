#' @title Spatial mark variogram function
#' @description Computes an estimator of the spatial mark variogram function.
#' @param xyt Spatial coordinates and times \eqn{(x,y,t)} of the point pattern.
#' @param ds A vector of distances \code{u} at which \eqn{\gamma_{sp}(u)} is computed.
#' @param ks A kernel function for the spatial distances. The default is the \code{"box"} kernel. It can also be \code{"epanech"} for the Epanechnikov kernel, or \code{"biweight"}.
#' @param hs A bandwidth of the kernel function \code{ks}.
#' @details \eqn{N = {[x_{1},t_{1}],[x_{2},t_{2}],...,[x_{n},t_{n}]}} with \eqn{x_{i}} in \eqn{W} subset of \eqn{R^2} and \eqn{t_{i}} in \eqn{T} subset of \eqn{R}, the spatial mark variogram \eqn{\gamma_{sp}(r)} considers the times \eqn{t_{i}} in the points \eqn{[x_{i},t_{i}]} as marks of the locations \eqn{x_{i}}. This approach makes sense if, as in statistical applications, the \eqn{t_{i}} come from some time interval \eqn{T}. The function \eqn{\gamma_{sp}(r)} is given by the mean of the half-squared differences of the times of points with spatial distance \eqn{r}, and the estimator is given by \deqn{\widehat{\gamma}_{sp}(r) = sum_{i = 1,...,n} sum_{j = 1,...,n; j! = i}(1/2 (t_{1}-t_{2})^2 \kappa_{\epsilon}(||x_{i}-x_{j}||-r)  / sum_{i = 1,...,n} sum_{j = 1,...,n; j != i} (\kappa_{\epsilon} (||x_{i}-x_{j}||-r)),} where \eqn{\kappa} is a one-dimensional kernel function with bandwidth \eqn{\epsilon}.
#' @return A list containing:
#' \itemize{
#'   \item \code{gspke}: A vector containing the values of \eqn{\widehat{\gamma}_{sp}(r)} estimated.
#'   \item \code{ds}: If \code{ds} is missing, a vector of distances \code{u} at which \eqn{\gamma_{sp}(u)} is computed under the restriction \eqn{0<\epsilon<r}.
#'   \item \code{kernel}: A vector of names and bandwidth of the spatial kernel.
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
#' # R plot
#' plot(hpp)
#'
#' # This function provides an kernel estimator of the spatial mark variogram
#' out <- gsp(hpp,hs=0.1/sqrt(100))
#'
#' # R plot - Spatial mark variogram
#' par(mfrow=c(1,1))
#' gsp_range <- range(out$gspke,1/12)
#' plot(out$ds,out$gspke,type="l",ylim=gsp_range,xlab="r = distance",
#' ylab=expression(hat(gamma)[sp](r)),main="Spatial mark variogram")
#' lines(out$ds,rep(1/12,length(out$ds)),type="l",col="red")
#'
#' ## End(Not run)
gsp <- function(xyt, ds, ks="epanech", hs){

  if (missing(hs)){
    d <- dist(xyt[,1:2])
    hs <- dpik(d,kernel=ks,range.x=c(min(d),max(d)))}

  if (missing(ds)){
    x <- xyt[,1]
    y <- xyt[,2]
    W <- ripras(x,y)
    poly <- W$bdry
    X <- poly[[1]]$x
    Y <- poly[[1]]$y
    sregion <- cbind(X,Y)
    bdry <- owin(poly=list(x=sregion[,1],y=sregion[,2]))
    rect <- as.rectangle(bdry)
    maxd <- min(diff(rect$xrange),diff(rect$yrange))/4
    ds <- seq(0, maxd, len=50)
    ds <- sort(ds)}

  pts <- xyt[,1:2]
  xytimes <- xyt[,3]
  ptsx <- pts[,1]
  ptsy <- pts[,2]
  ptst <- xytimes
  npt <- length(ptsx)
  nds <- length(ds)
  gsps <- rep(0,nds)

  kernel <- c(ks=ks,hs=hs)

  if (ks=="box"){ ks=1}
  else if (ks=="epanech") {ks=2}
  else if (ks=="biweight") {ks=3}

  storage.mode(gsps) <- "double"

  gspout <- .Fortran("gspcore",as.double(ptsx),as.double(ptsy),as.double(ptst),as.integer(npt),as.double(ds),as.integer(nds),as.integer(ks),as.double(hs),(gsps))

  gsps <- gspout[[9]]
  
  dsf <- rep(0,nds+1)
  dsf[2:(nds+1)] <- ds[1:nds]
  
  gsp <- rep(0,nds+1)
  gsp[2] <- gsps[2]
  gsp[3:(nds+1)] <- gsps[2:nds]
  
invisible(return(list(gspke=gsp,ds=dsf,kernel=kernel)))
}
