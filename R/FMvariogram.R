FMvariogram <- function(u,gmv){

xyt <- rpp(100)$xyt
n <- length(u)
mx <- u[n]

# Construction of the geodata object for geoR
gd <- list(coords=xyt[,1:2], data=xyt[,3])
oldClass(gd) <- "geodata"
  
# Empirical variograms geoR
sv <- variog(gd, uvec=u[-1], max.dist=mx)

# Empirical variograms mvstpp
sv$u <- u[-1]
sv$v <- gmv[-c(1,2)]

nugg.ef <- sv$v[1]
sill <- sv$v[which.max(sv$v)]
rang <- sv$u[which.max(sv$v)]
Parameters <- c(nugg.e=nugg.ef,sill=sill,rang=rang)

options(warn = -1)
# Model parameters estimated "by least squares fit of empirical variograms"
fit <- variofit(sv,ini.cov.pars=c(sill,rang), nugget=nugg.ef, cov.model="gau", weights="equal")

options(warn = 0)

gv <- v.g(dst=u,nugg.ef=summary(fit)[[5]][[1]],sill=summary(fit)[[3]][[1]],rang=summary(fit)[[3]][[2]])

invisible(return(list(u=u,gmv=gmv,gv=gv$gvm,Parameters=Parameters)))
}

# The Gaussian variogram model
v.g <- function(dst, nugg.ef,sill,rang){
  
  gvm <- nugg.ef+sill*(1-exp(-(dst/rang)^2))
  
  if (dst[1]==0){
      gvm[1]=0
    }
  invisible(return(list(dst=dst,gvm=gvm)))
}
