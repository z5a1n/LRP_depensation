library(rootSolve)

##################################################################################################################
# Survivorship
##################################################################################################################
# default is f = 0, M can be a constant or vector of length n_ages
# n_ages is number of age classes from the age of recruitment to the plus group
# sel = vector for vulnerability-at-age and is of length (n_age) and only required if f > 0
survivorship_F <- function(f=0,M,n_ages,sel,message=T){
  if(n_ages<=2){message(paste0("number of age classes must be greater than 2"))}
  if(length(M)==1){
    if(message==T){message(paste0("Constant M = ",M," used for all ages"))}
    M=rep(M,n_ages)
  }
  if(f==0){
    sel=rep(0,n_ages)
    if(message==T){message("Assumed F = 0, unfished survivorship")}
  }
  if(n_ages != length(sel)){
    message("age classes in vulnerability vector != n_ages")
    return(NA)
  }
  if(n_ages == length(sel)){
    l_age <- rep(NA,n_ages) 
    l_age[1] <- 1
    for(a in 2:(n_ages-1)){
      l_age[a] <- l_age[a-1]* exp(-(M[a-1]+f*sel[a-1]))
    }
    l_age[n_ages] <- l_age[n_ages-1]*exp(-(M[n_ages-1]+f*sel[n_ages-1]))/(1-exp(-(M[n_ages]+f*sel[n_ages])))
    return(l_age)
  }
}
##################################################################################################################
# Ref Pt calcs
##################################################################################################################
RPcalc <- function(M,waa,mat,sel,SRR,SRRpars,Req,Smax){
  
  output <- list()
  
  if(length(M)==1){
    message(paste0("Constant M = ",M," used for all ages"))
    M=rep(M,length(waa))
  }  
  if(length(waa) != length(mat) | length(waa) != length(sel)){
    message("at-age vectors are not the same length")
    return(NA)
  }
  
  f <- seq(0,10,0.001)
  ypr_f <- spr_f <- eq_ssb_f <- eq_rec_f <- phi_f <- rep(NA,length(f))

  for(i in 1:length(f)){
    phi_f[i] <- sum(survivorship_F(f=f[i],M=M,n_ages=length(sel),sel=sel,message=F)*waa*mat)
    ypr_f[i] <- sum(survivorship_F(f=f[i],M=M,n_ages=length(sel),sel=sel,message=F)*waa*
                      (1-exp(-(M+f[i]*sel)))*f[i]*sel/(M+f[i]*sel))
    spr_f[i] <- phi_f[i]/phi_f[1]
  }  
  
  eq_SSB_f_spr_40 <- phi_f[which(abs(spr_f-0.4) == min(abs(spr_f-0.4)))] * Req
  
  if(SRR=="BH_ab"){
    eq_rec_f <- (SRRpars[1]-1/phi_f)/SRRpars[2]
    eq_ssb_f <- eq_rec_f/(SRRpars[1]-eq_rec_f*SRRpars[2])
    eq_f <- f
  }   
  if(SRR=="BH_Rinf"){ 
    eq_rec_f <- SRRpars[1]-SRRpars[2]/phi_f
    eq_ssb_f <- SRRpars[2]/(SRRpars[1]/eq_rec_f-1)
    eq_f <- f
  }
  #############################################################################
  # sigmoidal BH 
  #############################################################################
  if(SRR=="sBH"){
    eq_rec_l <- list()
    eq_ssb_l <- list()
    f_l <- list() # list to accommodate multiple intersection points for each f
    for(i in 1:length(phi_f)){
      phi <- phi_f[i]
      s_sBH <- function(x) (x/phi-(SRRpars[1]/(1+(SRRpars[2]/x)^SRRpars[3]))) #SSB at intersection of phi and SRR
      eq_ssb_l[[i]] <- ssbs <- unique(round(c(uniroot.all(s_sBH,c(0,2*Smax),tol=0.00001), #uniroot was missing some roots. Needed to adjust the max SSB over which to look and then combine. 
                                              uniroot.all(s_sBH,c(0,Smax),tol=0.00001),
                                              uniroot.all(s_sBH,c(0,Smax/2),tol=0.00001)),4))
      eq_rec_l[[i]] <- SRRpars[1]/(1+(SRRpars[2]/eq_ssb_l[[i]])^SRRpars[3])
      f_l[[i]] <- rep(f[i],length(ssbs))
    }
    eq_rec_f <- unlist(eq_rec_l)
    eq_ssb_f <- unlist(eq_ssb_l)
    eq_f <- unlist(f_l)
  }
  #####################################################################################
  # SL
  ######################################################################################
  if(SRR=="SL"){
    eq_rec_l <- list()
    eq_ssb_l <- list()
    f_l <- list()
    for(i in 1:length(phi_f)){
      phi <- phi_f[i]
      k <- SRRpars[1]
      sk <- SRRpars[2]
      c <- SRRpars[3]
      r_SL <- function(x) (x - k*((phi*x)/sk)^c * exp(c*(1-(phi*x)/sk)))
      r_roots <- uniroot.all(r_SL,c(0,SRRpars[1]+1),tol=0.001)
      r_roots <- r_roots[r_roots>0]

      if(identical(r_roots, numeric(0))){next}
      
      thessbs <- list()
      therecs <- list()
      thefs <- list()
      for(j in 1:length(r_roots)){
        therec <- r_roots[j]
        ssb_SL <- function(x) (therec - k*(x/sk)^c * exp(c*(1-x/sk)))
        s_roots <- uniroot.all(ssb_SL,c(0,Smax),tol=0.001)
        thessbs[[j]] <- s_roots[which(abs(s_roots/therec-phi)==min(abs(s_roots/therec-phi)))] #get the SSB on the replacement line and SRR
        therecs[[j]] <- rep(therec,length(unlist(thessbs[[j]])))
        thefs[[j]] <- rep(f[i],length(unlist(thessbs[[j]])))
      } 
      eq_ssb_l[[i]] <- unlist(thessbs)
      eq_rec_l[[i]] <- unlist(therecs)
      f_l[[i]] <- unlist(thefs)
    }
      eq_rec_f <- unlist(eq_rec_l)
      eq_ssb_f <- unlist(eq_ssb_l)
      eq_f <- unlist(f_l)
  }
  ######################################################################################
  eq_y <- rep(NA,length(eq_rec_f))
  for(i in 1:length(eq_rec_f)){
    eq_y[i] <- eq_rec_f[i]*ypr_f[which(f==eq_f[i])]
  }  
  output[[1]] <- f_msy <- eq_f[which(eq_y==max(eq_y))]
  output[[2]] <- msy <- eq_y[which(eq_y==max(eq_y))]
  output[[3]] <- eq_ssb_f[which(eq_y==max(eq_y))] 
  output[[4]] <- eq_ssb_f 
  output[[5]] <- eq_y
  output[[6]] <- eq_f
  output[[7]] <- eq_SSB_f_spr_40
  output[[8]] <- eq_rec_f
  output[[9]] <- ssb0 <- max(eq_ssb_f[which(eq_ssb_f>0 & eq_f == 0)])
  
  names(output) <- c("Fmsy","msy","SSBmsy","eq_SSB","yield","f","eq_SSB_f_spr_40","eq_rec_f","SSB0")
  return(output)
  
}
##################################################################################################################
# SRR fits
##################################################################################################################
BHnll <- function(theta){
  Rinf <- exp(theta[1])
  Sq <- exp(theta[2])
  bhsig <- exp(theta[3])
  log_R <- log(R)
  P_log_R <- log((Rinf/(1+Sq/S)))
  negloglike <- sum(log(bhsig)+(log_R-P_log_R)^2/(2*bhsig^2))
  return (negloglike)
}
sBHnll <- function(theta){
  rinf <- exp(theta[1])
  sbhSq <- exp(theta[2])
  sbhc <- exp(theta[3])
  sbhsig <- exp(theta[4])
  log_R <- log(R)
  P_log_R <- log(rinf/(1+(sbhSq/S)^sbhc))
  negloglike <- sum(log(sbhsig)+(log_R-P_log_R)^2/(2*sbhsig^2))
  return (negloglike)
}
ricknll <- function(theta){
  rk <- exp(theta[1])
  rSk <- exp(theta[2])
  rsig <- exp(theta[3])
  log_R <- log(R)
  P_log_R <- log(((rk*S/rSk)*exp(1-S/rSk)))
  negloglike <- sum(log(rsig)+(log_R-P_log_R)^2/(2*rsig^2))
  return (negloglike)
}
SLnll <- function(theta){
  rk <- exp(theta[1])
  rSk <- exp(theta[2])
  rc <- exp(theta[3])
  rsig <- exp(theta[4])
  log_R <- log(R)
  P_log_R <- log(((rk*(S/rSk)^rc)*exp(rc*(1-S/rSk))))
  negloglike <- -sum(-log(rsig)-(log_R-P_log_R)^2/(2*rsig^2))
  return (negloglike)
}

HSnll <- function(theta){ # Barrowman and Myers 2000
  s_star <- exp(theta[1])
  a_slope <- mean(R/s_star)
  s_min <- S
  s_min[S>s_star] <- s_star
  negloglike <- log(sum((log(R)-(log(a_slope)+log(s_min)))^2))
  return (negloglike)
}

###################################################################
# Extract model parameters and posterior predictions
###################################################################

sBH_posterior_pred <- function(m,ssb){
  RES <- D0 <- data.frame(SSB=0:ssb); if(stock=="NSAS"){RES <- D0 <- data.frame(SSB=10*(0:(ssb/10)))}
  aRinf <- extract(m)$Rinf; aSq <- extract(m)$Sq
  ac <- extract(m)$c; asig <- extract(m)$sig
  for(j in 1:nrow(D0)){ # D0 has predicted recruitment for each value of SSB
    D0[j,1:length(ac)+1] <- aRinf/(1+(aSq/D0$SSB[j])^ac)
  }
  RES$P05 <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.05)
  RES$P95 <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.95)
  RES$M <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.5)
  mylist<-list(); mylist[[1]]<-median(aRinf); mylist[[2]]<-median(aSq)
  mylist[[3]]<-median(ac); mylist[[4]]<-median(asig); mylist[[5]]<-RES
  names(mylist)=c("Rinf","Sq","c","sig","RES")
  return(mylist)
}
SL_posterior_pred <- function(m,ssb){
  RES <- D0 <- data.frame(SSB=0:ssb); if(stock=="NSAS"){RES <- D0 <- data.frame(SSB=10*(0:(ssb/10)))}
  ak <- extract(m)$k; aSk <- extract(m)$Sk
  ac <- extract(m)$c; asig <- extract(m)$sig
  for(j in 1:nrow(D0)){  D0[j,1:length(ac)+1] <- ak*(D0$SSB[j]/aSk)^ac * exp(ac*(1-D0$SSB[j]/aSk)) }
  RES$P05 <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.05)
  RES$P95 <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.95)
  RES$M <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.5)
  mylist<-list(); mylist[[1]]<-median(ak); mylist[[2]]<-median(aSk)
  mylist[[3]]<-median(ac); mylist[[4]]<-median(asig); mylist[[5]]<-RES
  names(mylist)=c("k","Sk","c","sig","RES")
  return(mylist)
}

BH_posterior_pred <- function(m,ssb){
  RES <- D0 <- data.frame(SSB=0:ssb); if(stock=="NSAS"){RES <- D0 <- data.frame(SSB=10*(0:(ssb/10)))}
  aRinf <- extract(m)$Rinf; aSq <- extract(m)$Sq
  asig <- extract(m)$sig
  for(j in 1:nrow(D0)){  D0[j,1:length(aSq)+1] <- aRinf/(1+(aSq/D0$SSB[j])) }
  RES$P05 <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.05)
  RES$P95 <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.95)
  RES$M <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.5)
  mylist<-list(); mylist[[1]]<-median(aRinf); mylist[[2]]<-median(aSq)
  mylist[[3]]<-median(asig); mylist[[4]]<-RES
  names(mylist)=c("Rinf","Sq","sig","RES")
  return(mylist)
}

rick_posterior_pred <- function(m,ssb){
  RES <- D0 <- data.frame(SSB=0:ssb); if(stock=="NSAS"){RES <- D0 <- data.frame(SSB=10*(0:(ssb/10)))}
  ak <- extract(m)$k; aSk <- extract(m)$Sk
  asig <- extract(m)$sig
  for(j in 1:nrow(D0)){  D0[j,1:length(aSk)+1] <- ak*(D0$SSB[j]/aSk) * exp((1-D0$SSB[j]/aSk)) }
  RES$P05 <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.05)
  RES$P95 <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.95)
  RES$M <- apply(D0[,2:ncol(D0)],1,quantile,probs=0.5)
  mylist<-list(); mylist[[1]]<-median(ak); mylist[[2]]<-median(aSk)
  mylist[[3]]<-median(asig); mylist[[4]]<-RES
  names(mylist)=c("k","Sk","sig","RES")
  return(mylist)
}