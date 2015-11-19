#' @title Spatio-temporal double cluster point process model
#' @description Generate a random spatio-temporal point pattern, a simulated realisation of the cluster process.
#' @param lambp Intensity of the Poisson process of cluster centres. A single positive number or function.
#' @param a Length of the semi-axes x of ellipsoid.
#' @param b Length of the semi-axes y of ellipsoid.
#' @param c Length of the semi-axes z of ellipsoid.
#' @param mu Mean number of points per cluster (a single positive number).
#' @param s.region A two-column matrix specifying a polygonal region containing all data locations. If \code{s.region} is missing, the Ripley-Rasson estimate convex spatial domain is considered.
#' @param t.region t.region A vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.
#' @details We consider the straightforward extension of the classical Matern cluster process on the \eqn{R^3} case (with ellipsoid or balls) by considering the \eqn{z}-coordiantes as times.
#'
#' Consider a Poisson point process in the plane with intensity \eqn{\lambda_{p}} as cluster centres for all times “parent”, as well as a ellipsoid (or ball) where the semi-axes are of lengths \eqn{a}, \eqn{b} and \eqn{c}, around of each Poisson point under a random general rotation. The scatter uniformly in all ellipsoid (or ball)  of all points which are of the form \eqn{(x,y,z)}, the number of points in each cluster being random with a Poisson (\eqn{\mu}) distribution. The resulting point pattern is a spatio-temporal cluster point process with \eqn{t=z}. This point process has intensity \eqn{\lambda_{p} * \mu}.
#' @return The simulated spatio-temporal point pattern.
#' @author Francisco J. Rodriguez-Cortes <cortesf@@uji.es> \url{https://fjrodriguezcortes.wordpress.com} and Jonatan A. Gonzalez <jmonsalv@@uji.es>
#' @references Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.
#' @references Gabriel, E., Rowlingson, B., Diggle P J. (2013) \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software 53, 1-29.
#' @references Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.
#' @examples
#' ## Not run:
#' #################
#'
#' require(scatterplot3d)
#'
#' # Ellipsoid
#' Xe <- stdcpp(lambp=20, a=0.12, b=0.09, c=0.07, mu=100)
#'
#' plot(Xe$xyt)
#' par(mfrow=c(1,1))
#' scatterplot3d(Xe$xyt[,1],Xe$xyt[,2],Xe$xyt[,3],xlab="x",ylab="y",zlab="t")
#'
#' # Balls
#' Xb <- stdcpp(lambp=20, a=0.05, b=0.05, c=0.05, mu=100)
#'
#' plot(Xb$xyt)
#' par(mfrow=c(1,1))
#' scatterplot3d(Xb$xyt[,1],Xb$xyt[,2],Xb$xyt[,3],xlab="x",ylab="y",zlab="t")
#'
#' ## End(Not run)
stdcpp <- function(lambp, a, b, c, mu, s.region, t.region){

  if (missing(s.region)) {s.region <- matrix(c(1,1,0,0,0,1,1,0),ncol=2)}
  if (missing(t.region)) {t.region <- c(0,1)}

  stmc <- NULL
  stPoip <- rpp(lambp, s.region=s.region,t.region=t.region)$xyt
  for(i in 1:length(stPoip[,1])){
  PS <- PoiSph(a,b,c,mu=mu,centre=stPoip[i,])
  stmc <- rbind(stmc,PS)}

  ins <- inpip(stmc[,1:2],s.region)
  insw <- stmc[ins,]
  stc3 <- intim(insw,t.region)
  stc3 <- as.3dpoints(stc3)

  invisible(return(list(xyt=stc3,s.region=s.region,t.region=t.region)))
}

PoiSph <- function(a,b,c,mu,centre){
  n <- rpois(1,mu)

  if (a==b & b==c){
    r <- a*((runif(n))^(1/3))

    theta <- 2 * pi * runif(n)
    phi <- acos(2 * runif(n)-1)

    xr <- r * cos(theta)* sin(phi) + centre[1]
    yr <- r * sin(theta)* sin(phi) + centre[2]
    zr <- r * cos(phi) + centre[3]
  }
  else{
  a0 <- a*(runif(n)^(1/3))
  b0 <- b*(runif(n)^(1/3))
  c0 <- c*(runif(n)^(1/3))

  theta <- 2 * pi * runif(n)
  phi <- acos(2 * runif(n)-1)

  x0 <- a0 * cos(theta) * sin(phi)
  y0 <- b0 * sin(theta) * sin(phi)
  z0 <- c0 * cos(phi)

  aa <- 2*pi*runif(1)
  ab <- 2*pi*runif(1)
  ac <- 2*pi*runif(1)

  xr <- (z0*(sin(aa)*sin(ac)+cos(aa)*cos(ac)*sin(ab))-y0*(cos(aa)*sin(ac)-cos(ac)*sin(aa)*sin(ab))+x0*(cos(ac)*cos(ab))) + centre[1]
  yr <- (y0*(cos(aa)*cos(ac)+sin(aa)*sin(ac)*sin(ab))-z0*(cos(ac)*sin(aa)-cos(aa)*sin(ac)*sin(ab))+x0*(cos(ab)*sin(ac))) + centre[2]
  zr <- (z0*(cos(aa)*cos(ab))-x0*(sin(ab))+y0*(cos(ab)*sin(aa))) + centre[3]}

  invisible(return(cbind(xr,yr,zr)))
}

intim <- function(xyt,t.region){
  int <- NULL
  for(i in 1:length(xyt[,1])){
    if (xyt[i,3] > t.region[1] & xyt[i,3] < t.region[2]){
      int <- rbind(int,xyt[i,])}}

  invisible(return(int))
}
