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
