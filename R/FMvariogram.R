FMvariogram <- function(u,gmv){

xyt <- rpp(100)
n <- length(u)

# Construction of the geodata object for geoR
gd <- list(coords=xyt[,1:2], data=xyt[,3])
oldClass(gd) <- "geodata"
  
# Empirical variograms geoR
sv <- variog(gd, uvec=u[-1], max.dist=u[n])
sv
plot(sv)

# Empirical variograms mvstpp
sv$u <- u[-1]
sv$v <- gmv[-c(1,2)]
sv
plot(sv)

nugg.ef <- sv$v[1]
sill <- sv$v[which.max(sv$v)]
rang <- sv$u[which.max(sv$v)]

# Model parameters estimated "by eye"

# Gaussian cov.model 
lines.variomodel(cov.model = "gau", cov.pars = c(sil,rang), nug = nugg.ef, max.dist = rsup, col="blue")

# Model parameters estimated "by least squares fit of empirical variograms"
rsoutG <- variofit(sv, ini=c(sil,rang), nugget=nugg.ef, cov.model="gau", wei="equal")
rsoutG
lines(rsoutG, col="red")

c1 <- v.g(dst=u,nugg.ef=summary(rsoutG)[[5]][[1]],sill=summary(rsoutG)[[3]][[1]],rang=summary(rsoutG)[[3]][[2]])
lines(c1$dst,c1$gvm,type="l",col="green")

# Plot Spatial mark variogram

domsNe1 <- rep(0,lgrid+1)
domsNe1[2:(lgrid+1)] <- u[1:length(u)]

rasNe1 <- rep(0,lgrid+1)
rasNe1[2] <- gsp1[2]
rasNe1[3:(lgrid+1)] <- gsp1[2:lgrid]

# spatial component fitted variogram
domsNf1 <- rep(0,lgrid+1)
domsNf1[2:(lgrid+1)]<- c1$dst[1:length(c1$dst)]

rasNf1 <- rep(0,lgrid+1)
rasNf1[2] <- c1$gvm[2]
rasNf1[3:(lgrid+1)]<- c1$gvm[2:length(c1$gvm)]

return()
}

# The Gaussian variogram model
v.g <- function(dst, nugg.ef,sill,rang){
  
  gvm <- nugg.ef+sill*(1-exp(-(dst/rang)^2))
  
  if (dst[1]==0){
      gvm[1]=0
    }
  invisible(return(list(dst=dst,gvm=gvm)))
}