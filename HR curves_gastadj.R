library(surrosurv)
library(survival)
library(survminer)
library(survMisc)
library(rstpm2)
library(bpcp)
data("gastadj")
gastadj$trt<-ifelse(gastadj$trt==0.5,1,0)
gastadj$timeS<-gastadj$timeS/365.25
gastadj$timeT<-gastadj$timeT/365.25

index<-0
for(i in 1:length(unique(gastadj$trialref))){
  index<-index+1
  ldf<-gastadj[gastadj$trialref==unique(gastadj$trialref)[i],]
  mod<-paste("mod",index,sep="")
  assign(mod,stpm2(Surv(timeT,statusT==1)~trt, data=ldf, smooth.formula=~ns(log(timeT),df=2)+log(timeT):trt))
}


plot(mod1, newdata = data.frame(trt = 0), type = "hr",
     var = "trt", ci = FALSE, rug = FALSE,
     main = "Overall Survival", xlim=c(0,17), ylim=c(0,4),
     ylab = "Hazard ratio", xlab = "Time (years)",line.col="grey", xaxt="n",yaxt="n"
) #p=0.663

axis(1, at = c(5,10,15))
axis(2, at = c(1,2,3,4))

lines(mod2, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #p=0.464
lines(mod3, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #p=0.742
lines(mod5, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #p=0.765
lines(mod6, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #p=0.440
lines(mod7, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #p=0.902
lines(mod8, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #p=0.370
lines(mod9, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #p=0.780
lines(mod11, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #p=0.953
lines(mod12, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #p=0.425
lines(mod4, newdata = data.frame(trt = 0), type = "hr",
      var = "trt") #p=0.325
lines(mod14, newdata = data.frame(trt = 0), type = "hr",
      var = "trt") #p=0.136
lines(mod10, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="red") #p=0.068
lines(mod13, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="red")#p=0.029



#surrogate

index<-0
for(i in 1:length(unique(gastadj$trialref))){
  index<-index+1
  ldf<-gastadj[gastadj$trialref==unique(gastadj$trialref)[i],]
  mod<-paste("mod",index,sep="")
  assign(mod,stpm2(Surv(timeS,statusS==1)~trt, data=ldf, smooth.formula=~ns(log(timeS),df=2)+log(timeS):trt))
}

plot(mod1, newdata = data.frame(trt = 0), type = "hr",
     var = "trt", ci = FALSE, rug = FALSE,
     main = "Disease-free Survival ", xlim=c(0,17), ylim=c(0,4),
     ylab = "Hazard ratio", xlab = "Time (years)",line.col="grey", xaxt="n",yaxt="n"
) #p=0.253
axis(1, at = c(5,10,15))
axis(2, at = c(1,2,3,4))
lines(mod7, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #0.644
lines(mod8, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#0.855
lines(mod9, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #0.612
lines(mod2, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #0.521
lines(mod3, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #0.978
lines(mod12, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #0.998
lines(mod13, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #0.571
lines(mod5, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey") #0.352
lines(mod4, newdata = data.frame(trt = 0), type = "hr",
      var = "trt") #0.261
lines(mod6, newdata = data.frame(trt = 0), type = "hr",
      var = "trt") #0.108
lines(mod10, newdata = data.frame(trt = 0), type = "hr",
      var = "trt")#0.129
lines(mod11, newdata = data.frame(trt = 0), type = "hr",
      var = "trt") #0.118
lines(mod14, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="red") #0.009


  
  