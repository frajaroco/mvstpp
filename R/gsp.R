gsp <- function(xyt,s.region,s.lambda,ds,ks="epanech",hs,correction="none",approach="simplified"){
  
  correc <- c("none","isotropic","border","modified.border","translate","setcovf")
  id <- match(correction,correc,nomatch=NA)
  if (any(nbg <- is.na(id))){
    messnbg <- paste("unrecognised correction method:",paste(dQuote(correction[nbg]),collapse=","))
    stop(messnbg,call.=FALSE)
  }
  id <- unique(id)	
  correc2 <- rep(0,6)
  correc2[id] <- 1	
  
  appro <- c("simplified","standardised")
  im <- match(approach,appro,nomatch=NA)
  if (any(nbm <- is.na(im))){
    messnbm <- paste("unrecognised type of estimator:",paste(dQuote(approach[nbm]),collapse=","))
    stop(messnbm,call.=FALSE)
  }
  im <- unique(im)
  appro2 <- rep(0,2)
  appro2[im] <- 1
  
  ker <- c("box","epanech","biweight")
  ik <- match(ks,ker,nomatch=NA)
  if (any(nbk <- is.na(ik))){
    messnbk <- paste("unrecognised kernel function:",paste(dQuote(ks[nbk]),collapse=","))
    stop(messnbk,call.=FALSE)
  }
  ik <- unique(ik)
  ker2 <- rep(0,3)
  ker2[ik] <- 1
  
  dup <- duplicated(data.frame(xyt[,1],xyt[,2],xyt[,3]),fromLast = TRUE)[1]
  if (dup == TRUE){
    messnbd <- paste("spatio-temporal data contain duplicated points")
    warning(messnbd,call.=FALSE)
  }
  
  if (missing(hs)){
    d <- dist(xyt[,1:2])
    hs <- dpik(d,kernel=ks,range.x=c(min(d),max(d)))
  }
  
  if (missing(s.region)){
    if (appro2[1]==1){
      s.region <- sbox(xyt[, 1:2], xfrac = 0.01, yfrac = 0.01)
    } else{
      x <- xyt[,1]
      y <- xyt[,2]
      W <- ripras(x,y)
      poly <- W$bdry
      X <- poly[[1]]$x
      Y <- poly[[1]]$y
      s.region <- cbind(X,Y)
    }
  }
  
  bsw <- owin(poly=list(x=s.region[,1],y=s.region[,2]))
  
  if (missing(ds)){
    rect <- as.rectangle(bsw)
    maxd <- min(diff(rect$xrange),diff(rect$yrange))/4
    ds <- seq(hs, maxd,len=100)[-1]
    ds <- sort(ds)
  }  
  if(ds[1]==0){ds <- ds[-1]
  }
  
  kernel <- c(ks=ks,hs=hs)
  
  pts <- xyt[,1:2]
  xytimes <- xyt[,3]
  ptsx <- pts[,1]
  ptsy <- pts[,2]
  ptst <- xytimes
  npt <- length(ptsx)
  nds <- length(ds)
  area <- area(bsw)
  pert <- perimeter(bsw)
  gsps <- rep(0,nds)
  
  storage.mode(gsps) <- "double"
  
  if (appro2[1]==1){
    gspout <- .Fortran("gspcore",as.double(ptsx),as.double(ptsy),as.double(ptst),
                       as.integer(npt),as.double(ds),as.integer(nds),as.integer(ker2),
                       as.double(hs),(gsps),PACKAGE="msfstpp")
    
    gsps <- gspout[[9]]
    
    dsf <- rep(0,nds+2)
    dsf[3:(nds+2)] <- ds
    ds <- dsf 
    
    egsp <- rep(0,nds+2)
    egsp[2] <- gsps[1]
    egsp[3:(nds+2)] <- gsps
    
    invisible(return(list(egsp=egsp,ds=ds,kernel=kernel,s.region=s.region)))
  } else {
    
    if(missing(s.lambda)){
      misl <- 1
      s.lambda <- rep(npt/area,npt)
    } else {
      misl <- 0
      if (length(s.lambda)==1){
        s.lambda <- rep(s.lambda,npt)
      }
    }
    
    wrs <- array(0,dim=c(npt,npt))
    wts <- array(0,dim=c(npt,npt))
    wbi <- array(0,dim=c(npt,nds))
    wbimod <- array(0,dim=c(npt,nds))
    wss <- rep(0,nds)
    
    options(warn = -1) 
    
    pxy <- ppp(x=ptsx,y=ptsy,window=bsw)  
    
    # correction="isotropic"
    
    if(correction=="isotropic"){
      wisot <- edge.Ripley(pxy,pairdist(pts))
      wrs <- 1/wisot
    }
    
    # correction="translate"
    if(correction=="translate"){
      wtras <- edge.Trans(pxy)
      wts <- 1/wtras
    }
    
    #  correction=="border" or "modified border"
    
    if(any(correction=="border")|any(correction=="modified.border")){
      bi <- bdist.points(pxy)
      for(i in 1:nds) { 
        wbi[,i] <- (bi>ds[i])/sum((bi>ds[i])/s.lambda)
        wbimod[,i] <- (bi>ds[i])/eroded.areas(bsw,ds[i])
      }
      wbi[is.na(wbi)] <- 0
    }
    
    # correction="setcovf"
    
    if(correction=="setcovf"){
      for (i in 1:nds){
        wss[i] <- area-((pert*ds[i])/pi)
      }
      wss <- 1/wss
    }
    
    options(warn = 0)
    
    gspout <- .Fortran("gspcoreinh",as.double(ptsx),as.double(ptsy),as.double(ptst),
                        as.integer(npt),as.double(ds),as.integer(nds),as.double(s.lambda),
                        as.integer(ker2),as.double(hs),as.double(wrs),as.double(wts),
                        as.double(wbi),as.double(wbimod),as.double(wss),as.integer(correc2),
                        (gsps),PACKAGE="msfstpp")
    gsps <- gspout[[16]]
    
    dsf <- rep(0,nds+2)
    dsf[3:(nds+2)] <- ds
    ds <- dsf 
    
    egsp <- rep(0,nds+2)
    egsp[2] <- gsps[1]
    egsp[3:(nds+2)] <- gsps
    
    invisible(return(list(egsp=egsp,ds=ds,kernel=kernel,s.region=s.region,s.lambda=s.lambda)))
  }
}