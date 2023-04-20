#
#  This code generates simulated data similar to the  COV-2069 trial.
#  It estimates a 1 parameter model a 2 parameter logistic,
#  a 3 parameter logistic and a log-linear model
#
#  The code takes about 10 minutes to run
#
#  The ests vector has the following value
#ests
#[1] -0.4025000 -0.4612780  0.9181465 -0.7703170 -0.6351498 13.0533910  0.6986539
#[8] -0.6723965
#



set.seed(17843)
p<- .10/240

Y<-NULL;D<-NULL;Z<-NULL; AB<-NULL

# Generate placebo group data, note the antibody values are counterafactual
for (i in 1:800){
    yi<-rbinom(240,1,p)
    ilast<-ifelse(sum(yi)==0,240, which.max(yi))
    Y<-c(Y,yi[1:ilast])
    D<-c(D,seq(1,ilast,1))
    Z<-c(Z,rep(0,ilast))
    Ab0<-runif(1)*.2 + 4
    AB<-c(AB, Ab0 + seq(1,ilast,1)*(-.012))
    }

# Generate mab group data

for (i in 1:800){
    yi1<-rbinom(120,1,p*.20)
    yi2<-rbinom( 20,1,p*.50)
    yi3<-rbinom(100,1,p)
    yi<-c(yi1,yi2,yi3)
    ilast<-ifelse(sum(yi)==0,240, which.max(yi))
    Y<-c(Y,yi[1:ilast])
    D<-c(D,seq(1,ilast,1))
    Z<-c(Z,rep(1,ilast))
    Ab0<-runif(1)*.2+4
    AB<-c(AB, Ab0 + seq(1,ilast,1)*(-0.012))
    }

par(mfrow=c(1,2))
IM<-((Z==1)&(Y==1))
IP<-((Z==0)&(Y==1))
plot(D[IM])
plot(D[IP])


dat_tve<-data.frame(Y,Z,D,AB,rep(1,length(Y)))
colnames(dat_tve)<-c("event","trt","tstart","log10ID50","off")

#-------------------------------------------------------------------------------------------------------
# Fit a 1 parameter model
#-------------------------------------------------------------------------------------------------------
log_1P = function(z, m, a) {
   p=plogis(a)
   z*log(p)
}
wrap_P1P =function(par){
    dat_tve$off<-log_1P(dat_tve$trt, dat_tve$log10ID50, par[1])
    fit=glm(event ~ tstart+offset(off),
              family = poisson(),
              data = dat_tve)
    -logLik(fit)
}
opt_P1P = optim(-0.40, wrap_P1P)


#--------------------------------------------------------------------------------------------------------
# Fit a 2 parameter logistic model:  Note  exp(c) is the logistic intercept and exp(c) x b the logistic slope
#---------------------------------------------------------------------------------------------------------

log_2PL = function(z, m, b,c) {
   z*log(plogis(exp(c)*(1+b*m)))
 }
wrap_P2PL =function(par){
    dat_tve$off<-log_2PL(dat_tve$trt, dat_tve$log10ID50, par[1], par[2])
    fit=glm(event ~ tstart+offset(off),
              family = poisson(),
              data = dat_tve)
    -logLik(fit)
}
opt_P2PL = optim(c(-0.46,0.92), wrap_P2PL)


#--------------------------------------------------------------------------------------------------------
# Fit a 3 parameter logistic model
#---------------------------------------------------------------------------------------------------------
log_3PL = function(z, m, a, b,c) {
   p=plogis(a)
   z*log(p+(1-p)* plogis(exp(c)*(1+b*m)))
 }
wrap_P3PL =function(par){
    dat_tve$off<-log_3PL(dat_tve$trt, dat_tve$log10ID50, par[1], par[2],par[3])
    fit=glm(event ~ tstart+offset(off),
              family = poisson(),
              data = dat_tve)
    -logLik(fit)
}
opt_P3PL = optim(c( -0.77, -0.63, 13.37), wrap_P3PL)


#---------------------------------------------------------------------------------------------------
# Let's fit a 2 parameter log-linear model
#---------------------------------------------------------------------------------------------------
log_P2P = function(z, m, a,b) {
   z*(a+b*m)
}
wrap_P2P =function(par){
    dat_tve$off<-log_P2P(dat_tve$trt, dat_tve$log10ID50, par[1],par[2])
    fit=glm(event ~ tstart+offset(off),
              family = poisson(),
              data = dat_tve)
    -logLik(fit)
}
opt_P2P = optim(c(0.70,-0.67), wrap_P2P)

#---------------------------------------------------------------------------------------------------
# Let's graph the PE from the P3PL fit
#---------------------------------------------------------------------------------------------------

par(mfrow=c(1,1))
m<-seq(0,5,.1)
a<-opt_P3PL$par[1]
b<-opt_P3PL$par[2]
c<-opt_P3PL$par[3]
p<-plogis(a)
PE<-1- (p+(1-p)*plogis(exp(c)*(1+b*m)))
plot(m,PE,ylim=c(0,1),xlim=c(0,5),las=1,type="n")
lines(m,PE)

ests<-c(opt_P1P$par,opt_P2PL$par,opt_P3PL$par,opt_P2P$par)
