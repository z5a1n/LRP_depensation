#R version 4.4.3
library(ggplot2)   # v3.5.1
library(rstan)     # v2.32.7
library(rstanarm)  # v2.32.1
library(shinystan) # v2.6.0
library(loo)       # v2.8.0

source("Functions.R")

stock <- "BS"
#stock <- "GoR" 
#stock <- "NSAS" 
#stock <- "WBSS" 

DF <- read.csv("Data/herring_RS_data.csv")
DF$SSB <- DF$SSB/1000 # put in kt
DF$Rec <- DF$Rec/1000000 # put in billions

DF <- DF[DF$Stock==stock,]

R<-DF$Rec; S<-DF$SSB

if(stock=="BS") {maxS <- 800; maxR <- 15; start = 500}  # Used for plotting 
if(stock=="GoR") {maxS <- 150; maxR <- 8; start = 50} 
if(stock=="NSAS") {maxS <- 4000; maxR <- 150; start = 1000} 
if(stock=="WBSS") {maxS <- 400; maxR <- 6; start = 200} 

P <- ggplot() + geom_point(data=DF, aes(y=Rec,x=SSB),colour="red4") + theme_classic() +
  labs(x='SSB (kt)', y="Recruitment (billions)") + 
  scale_x_continuous(limits=c(0,maxS), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0,maxR), expand = c(0, 0)) 

##########################################################################################
nw <- 1000
ni <- 10000

sBH_data <- list(Y=length(S),S=S,R=R,RinfLow=0,RinfUp=max(R),SqLow=0,SqUp=max(S),cLow=0,cUp=5,sigmaLow=0,sigmaUp=2)
 f_sBH <- stan(file="Stan/sBH_c5.stan",data=sBH_data,chains=4,warmup=nw,iter=ni,cores=1,refresh=0,control=list(adapt_delta=0.99,max_treedepth=12))
SL_data <- list(Y=length(S),S=S,R=R,kLow=0,kUp=max(R),SkLow=0,SkUp=max(S),cLow=0,cUp=5,sigmaLow=0,sigmaUp=2)
 f_SL <- stan(file="Stan/SL_c5.stan",data=SL_data,chains=4,warmup=nw,iter=ni,cores=1,refresh=0,control=list(adapt_delta=0.99,max_treedepth=12))
BH_data <- list(Y=length(S),S=S,R=R,RinfLow=0,RinfUp=max(R),SqLow=0,SqUp=max(S),sigmaLow=0,sigmaUp=2)
 f_BH <- stan(file="Stan/BH.stan",data=BH_data,chains=4,warmup=nw,iter=ni,cores=1,refresh=0,control=list(adapt_delta=0.99,max_treedepth=12))
rick_data <- list(Y=length(S),S=S,R=R,kLow=0,kUp=max(R),SkLow=0,SkUp=max(S),sigmaLow=0,sigmaUp=2)
 f_rick <- stan(file="Stan/rick.stan",data=rick_data,chains=4,warmup=nw,iter=ni,cores=1,refresh=0,control=list(adapt_delta=0.99,max_treedepth=12))

##########################################################################################
# HS fit
##########################################################################################
mle_hs <- optim(fn=HSnll,par=log(start),lower=log(min(S)), upper=log(max(S)), hessian = T, method="L-BFGS-B")
if(mle_hs$convergence!=0){print("Did not conv")}
Sad <- S
Sad[S>exp(mle_hs$par)] <- exp(mle_hs$par)
 
# Visual inspection of diagnostics
#launch_shinystan(f_sBH)
#launch_shinystan(f_SL)
#launch_shinystan(f_BH)
#launch_shinystan(f_rick)

# Extract parameter estimates
p_sBH <- sBH_posterior_pred(m=f_sBH,ssb=maxS)
p_SL <- SL_posterior_pred(m=f_SL,ssb=maxS)
p_BH <- BH_posterior_pred(m=f_BH,ssb=maxS)
p_rick <- rick_posterior_pred(m=f_rick,ssb=maxS)
p_hs <- c(exp(mle_hs$par),exp(mean(log(R/Sad))))

################################################################################
#Biological Data and Selectivity for Ref Pt Calcs 
################################################################################
BD <- read.csv("Data/herring_bio_data.csv")
BD <- BD[BD$Stock==stock,1:13] #exclude units and reference columns

if(stock %in% c("GoR")){cols <- 5:12}
if(stock %in% c("BS","WBSS","NSAS")){cols <- 5:13}
M <- colSums(BD[BD$Parameter=="M",cols])/nrow(BD[BD$Parameter=="M",cols]) # 1 to 8+ or 0 to 8+
waa <- colSums(BD[BD$Parameter=="waa",cols])/nrow(BD[BD$Parameter=="waa",cols])
mat <- colSums(BD[BD$Parameter=="mat",cols])/nrow(BD[BD$Parameter=="mat",cols])
sel <- colSums(BD[BD$Parameter=="vul",cols])/max(colSums(BD[BD$Parameter=="vul",cols]))

phi0 <- sum(survivorship_F(M=M,n_ages = length(M))*waa*mat) #g per recruit
################################################################################
fs=14


P + geom_function(fun=function(x) (p_sBH$Rinf/(1+(p_sBH$Sq/x)^p_sBH$c)),colour="black",linetype=2,linewidth=1.2) +
  geom_function(fun=function(x) (p_SL$k*(x/p_SL$Sk)^p_SL$c * exp(p_SL$c*(1-x/p_SL$Sk))),colour="blue",linetype=2,linewidth=1.2) + 
  geom_function(fun=function(x) (p_BH$Rinf/(1+(p_BH$Sq/x))),colour="black",linewidth=1.2) +
  geom_function(fun=function(x) (p_rick$k*x/p_rick$Sk * exp(1-x/p_rick$Sk)),colour="blue",linewidth=1.2) + 
  geom_segment(aes(x = 0, y = 0, xend = p_hs[1], yend =  p_hs[1]*p_hs[2]), colour = "green",linewidth=1.2) +
  geom_segment(aes(x = p_hs[1], y = p_hs[2]*p_hs[1], xend = maxS, yend =  p_hs[1]*p_hs[2]), colour = "green",linewidth=1.2) +
  geom_function(fun=function(x) (1/phi0*x),colour="grey",linewidth=1.2) +
  ggtitle(paste0(stock," herring")) +
  theme(axis.text.x = element_text(size = fs), axis.text.y = element_text(size = fs), axis.title.x = element_text(size = fs), axis.title.y = element_text(size = fs),
          title = element_text(size = fs), plot.margin = margin(10, 20, 10, 10),plot.title = element_text(hjust = 0.5)) 

ggsave(filename=paste0("Figs/",stock,".png"),plot = last_plot(), width=7,height=5,units="in")
################################################################################
ggplot() + geom_point(data=DF,aes(y=Rec,x=SSB)) + theme_classic() + ggtitle(paste0(stock," herring - sBH")) +
  geom_ribbon(data=p_sBH$RES,mapping=aes(ymax=P95,ymin=P05,x=SSB),fill='grey',alpha=0.5) +
  geom_path(data=p_sBH$RES,mapping=aes(y=M,x=SSB)) +
  scale_x_continuous(name="SSB (kt)",limits=c(0,maxS), expand = c(0, 0)) +
  scale_y_continuous(name="Recruitment (billions)", limits=c(0,maxR), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs), axis.text.y = element_text(size = fs), axis.title.x = element_text(size = fs), axis.title.y = element_text(size = fs),
        title = element_text(size = fs), plot.margin = margin(10, 20, 10, 10),plot.title = element_text(hjust = 0.5)) 
  
ggsave(filename=paste0("Figs/figsBH",stock,".png"),plot = last_plot(), width=7,height=5,units="in")
################################################################################
ggplot() + geom_point(data=DF,aes(y=Rec,x=SSB)) + theme_classic() + ggtitle(paste0(stock," herring - SL")) +
  geom_ribbon(data=p_SL$RES,mapping=aes(ymax=P95,ymin=P05,x=SSB),fill='lightblue',alpha=0.5) +
  geom_path(data=p_SL$RES,mapping=aes(y=M,x=SSB),colour='blue') +
  scale_x_continuous(name="SSB (kt)",limits=c(0,maxS), expand = c(0, 0)) +
  scale_y_continuous(name="Recruitment (billions)", limits=c(0,maxR), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs), axis.text.y = element_text(size = fs), axis.title.x = element_text(size = fs), axis.title.y = element_text(size = fs),
      title = element_text(size = fs), plot.margin = margin(10, 20, 10, 10),plot.title = element_text(hjust = 0.5)) 

ggsave(filename=paste0("Figs/figSL",stock,".png"),plot = last_plot(), width=7,height=5,units="in")
################################################################################

TAB <- data.frame(mod=c("BH","sBH","rick","SL","HS"),P1=NA,P2=NA,c=NA,sig=NA,p=NA,looic=NA)

TAB[TAB$mod=="sBH",2:7] <- c(p_sBH$Rinf,p_sBH$Sq,p_sBH$c,p_sBH$sig,length(which(extract(f_sBH)$c>1))/length(extract(f_sBH)$c),loo(extract_log_lik(f_sBH))$estimates[3,1])
TAB[TAB$mod=="BH",2:7] <- c(p_BH$Rinf,p_BH$Sq,NA,p_BH$sig,NA,loo(extract_log_lik(f_BH))$estimates[3,1])
TAB[TAB$mod=="SL",2:7] <- c(p_SL$k,p_SL$Sk,p_SL$c,p_SL$sig,length(which(extract(f_SL)$c>1))/length(extract(f_SL)$c),loo(extract_log_lik(f_SL))$estimates[3,1])
TAB[TAB$mod=="rick",2:7] <- c(p_rick$k,p_rick$Sk,NA,p_rick$sig,NA,loo(extract_log_lik(f_rick))$estimates[3,1])
TAB[TAB$mod=="HS",2:3] <- p_hs

write.csv(TAB,paste0("Tab/",stock,".csv"),row.names = F)

