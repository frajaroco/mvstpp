\name{FMvariogram}
\alias{FMvariogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit Gaussian geostatistical variogram to spatial and temporal mark variogram functions 
%%  ~~function to do ... ~~
}
\description{Fit Gaussian geostatistical variogram to the obtained spatial and temporal mark variograms
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FMvariogram(u, gmv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{vector of spatial or temporal distance at which \code{gsp} and \code{gte} was computed and where it will be fit.
  %%     ~~Describe \code{u} here~~
}
  \item{gmv}{vector containing the values of \eqn{\gamma} estimated through \code{gsp} and \code{gte}.
%%     ~~Describe \code{gmv} here~~
}
}
\details{This function fit Gaussian geostatistical variogram to the obtained mark variogram through \code{gsp} and \code{gte} functions using the R package \code{geoR} Ribeiro and Diggle (2001).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{u}{If u is missing, a vector of distances u at which \code{gsp} and \code{gte} is computed from 0 to until quarter of the maximum distance between the points in the pattern.}
 \item{gmv}{The fitted Gaussian geostatistical variogram.}
 \item{Parameters}{The nugget effects, sill and range fitted for a Gaussian geostatistical variogram to the obtained spatial and temporal mark variograms.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Baddeley, A., Rubak, E., Turner, R. (2015). Spatial Point Patterns: Methodology and Applications with R. CRC Press, Boca Raton.

Chiu, S. N., Stoyan, D., Kendall, W. S., and Mecke, J. (2013). Stochastic Geometry and its Applications. John Wiley & Sons.

Gabriel, E., Rowlingson, B., Diggle P J. (2013) \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53:1--29.

Gonzalez, J. A., Rodriguez-Cortes, F. J., Cronie, O. and Mateu, J. (2016). Spatio-temporal point process statistics: a review. Spatial Statiscts. 18:505-544.

Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.

Stoyan, D., Rodriguez-Cortes, F. J., Mateu, J. and Wilfried, G. (2016). Mark variograms for spatio-temporal point processes. Spatial Statistics, 20:125-147.

Ribeiro, Jr., Diggle, P.J. (2001). \code{geoR}: A package for geostatistical analysis. R New 1, 2:14--18.
}
\author{Francisco J. Rodriguez Cortes <cortesf@uji.es> \url{https://fjrodriguezcortes.wordpress.com}}

\examples{
## Not run:
#################

# A realisation of spatio-temporal double-cluster (ellipsoid) point processes
Xe <- stdcpp(lambp=20, a=0.12, b=0.09, c=0.07, mu=100)
plot(Xe$xyt)

# Spatial mark variogram function
outs <- gsp(Xe$xyt)

# R plot - Spatial mark variogram function
par(mfrow=c(1,1))
xl <- c(0,0.25)
yl <- c(min(0,outs$egsp,1/12),max(0,outs$egsp,1/12))
plot(outs$ds,outs$egsp,type="l",xlab="r = distance",ylab=expression(gamma[sp](r)),
     xlim=xl,ylim=yl,col=1,cex.lab=1.5,cex.axis=1.5)
lines(outs$ds,rep(1/12,length(outs$ds)),col=11)
points(0,0,col=11,cex=1.5)

# Model parameters estimated "by least squares fit of empirical variograms"
fitv <- FMvariogram(outs$ds,outs$egsp)
fitv$Parameters
lines(fitv$u,fitv$gv)

# Temporal mark variogram function
outt <- gte(Xe$xyt)

# R plot - Temporal mark variogram function
par(mfrow=c(1,1))
xl <- c(0,0.25)
yl <- c(min(0,outt$egte,1/6),max(0,outt$egte,1/6))
plot(outt$dt,outt$egte,type="l",xlab="t = time",ylab=expression(gamma[te](t)),
     xlim=xl,ylim=yl,col=1,cex.lab=1.5,cex.axis=1.5)
lines(outt$dt,rep(1/6,length(outt$dt)),col=11)
points(0,0,col=11,cex=1.5)

# Model parameters estimated "by least squares fit of empirical variograms"
fitv <- FMvariogram(outt$dt,outt$egte)
fitv$Parameters
lines(fitv$u,fitv$gv)
## End(Not run)
}
