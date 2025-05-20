#R version 4.4.3
library(ggplot2)   # v3.5.1
library(rstan)     # v2.32.7
library(rstanarm)  # v2.32.1
library(shinystan) # v2.6.0
library(loo)       # v2.8.0

source("Functions.R")

DF <- read.csv("Data/cod_RS_data.csv")
R<-DF$Rec
S<-DF$SSB
maxS <- 200
maxR <- 250

stock<-"cod"

P <- ggplot() + geom_point(data=DF, aes(y=Rec,x=SSB),colour="red4") + theme_classic() +
  labs(x='SSB (kt)', y="Recruitment (billions)") + 
  scale_x_continuous(limits=c(0,maxS), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0,maxR), expand = c(0, 0)) 

##########################################################################################
nw <- 1000
ni <- 10000

sBH_data <- list(Y=length(S),S=S,R=R,RinfLow=0,RinfUp=max(R)*2,SqLow=0,SqUp=max(S)*2,cLow=0,cUp=5,sigmaLow=0,sigmaUp=1)
  f_sBH <- stan(file="Stan/sBH_c5.stan",data=sBH_data,chains=4,warmup=nw,iter=ni,cores=1,refresh=0,control=list(adapt_delta=0.99,max_treedepth=12))
SL_data <- list(Y=length(S),S=S,R=R,kLow=0,kUp=max(R)*2,SkLow=0,SkUp=max(S)*2,cLow=0,cUp=5,sigmaLow=0,sigmaUp=1)
  f_SL <- stan(file="Stan/SL_c5.stan",data=SL_data,chains=4,warmup=nw,iter=ni,cores=1,refresh=0,control=list(adapt_delta=0.99,max_treedepth=12))

# Visual inspection of diagnostics
#launch_shinystan(f_sBH)
#launch_shinystan(f_SL)

# Extract parameter estimates
p_sBH <- sBH_posterior_pred(m=f_sBH,ssb=maxS)
p_SL <- SL_posterior_pred(m=f_SL,ssb=maxS)

################################################################################
#Biological Data and Selectivity for Ref Pt Calcs 
################################################################################
BD <- read.csv("Data/cod_bio_data.csv")
cols<-3:ncol(BD)
M <- colSums(BD[BD$Parameter=="M",cols])/nrow(BD[BD$Parameter=="M",cols]) 
waa <- colSums(BD[BD$Parameter=="waa",cols])/nrow(BD[BD$Parameter=="waa",cols])
mat <- colSums(BD[BD$Parameter=="mat",cols])/nrow(BD[BD$Parameter=="mat",cols])
sel <- colSums(BD[BD$Parameter=="vul",cols])/max(colSums(BD[BD$Parameter=="vul",cols]))

phi0 <- sum(survivorship_F(M=M,n_ages = length(M))*waa*mat) #g per recruit
################################################################################
fs=14


P + geom_function(fun=function(x) (p_sBH$Rinf/(1+(p_sBH$Sq/x)^p_sBH$c)),colour="black",linetype=2,linewidth=1.2) +
  geom_function(fun=function(x) (p_SL$k*(x/p_SL$Sk)^p_SL$c * exp(p_SL$c*(1-x/p_SL$Sk))),colour="blue",linetype=2,linewidth=1.2) + 
  geom_function(fun=function(x) (1/phi0*x),colour="grey",linewidth=1.2) +
  ggtitle("3Ps cod") +
  theme(axis.text.x = element_text(size = fs), axis.text.y = element_text(size = fs), axis.title.x = element_text(size = fs), axis.title.y = element_text(size = fs),
        title = element_text(size = fs), plot.margin = margin(10, 20, 10, 10),plot.title = element_text(hjust = 0.5)) 

ggsave(filename=paste0("Figs/cod.png"),plot = last_plot(), width=7,height=5,units="in")
################################################################################
ggplot() + geom_point(data=DF,aes(y=Rec,x=SSB)) + theme_classic() + ggtitle("3Ps cod - sBH") +
  geom_ribbon(data=p_sBH$RES,mapping=aes(ymax=P95,ymin=P05,x=SSB),fill='grey',alpha=0.5) +
  geom_path(data=p_sBH$RES,mapping=aes(y=M,x=SSB)) +
  scale_x_continuous(name="SSB (kt)",limits=c(0,maxS), expand = c(0, 0)) +
  scale_y_continuous(name="Recruitment (billions)", limits=c(0,maxR), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs), axis.text.y = element_text(size = fs), axis.title.x = element_text(size = fs), axis.title.y = element_text(size = fs),
        title = element_text(size = fs), plot.margin = margin(10, 20, 10, 10),plot.title = element_text(hjust = 0.5)) 

ggsave(filename=paste0("Figs/figsBHcod.png"),plot = last_plot(), width=7,height=5,units="in")
################################################################################
ggplot() + geom_point(data=DF,aes(y=Rec,x=SSB)) + theme_classic() + ggtitle("3Ps cod - SL") +
  geom_ribbon(data=p_SL$RES,mapping=aes(ymax=P95,ymin=P05,x=SSB),fill='lightblue',alpha=0.5) +
  geom_path(data=p_SL$RES,mapping=aes(y=M,x=SSB),colour='blue') +
  scale_x_continuous(name="SSB (kt)",limits=c(0,maxS), expand = c(0, 0)) +
  scale_y_continuous(name="Recruitment (billions)", limits=c(0,maxR), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = fs), axis.text.y = element_text(size = fs), axis.title.x = element_text(size = fs), axis.title.y = element_text(size = fs),
        title = element_text(size = fs), plot.margin = margin(10, 20, 10, 10),plot.title = element_text(hjust = 0.5)) 

ggsave(filename=paste0("Figs/figSLcod.png"),plot = last_plot(), width=7,height=5,units="in")
################################################################################

TAB <- data.frame(mod=c("sBH","SL"),P1=NA,P2=NA,c=NA,sig=NA,p=NA,looic=NA)

TAB[TAB$mod=="sBH",2:7] <- c(p_sBH$Rinf,p_sBH$Sq,p_sBH$c,p_sBH$sig,length(which(extract(f_sBH)$c>1))/length(extract(f_sBH)$c),loo(extract_log_lik(f_sBH))$estimates[3,1])
TAB[TAB$mod=="SL",2:7] <- c(p_SL$k,p_SL$Sk,p_SL$c,p_SL$sig,length(which(extract(f_SL)$c>1))/length(extract(f_SL)$c),loo(extract_log_lik(f_SL))$estimates[3,1])

write.csv(TAB,paste0("Tab/cod.csv"),row.names = F)
