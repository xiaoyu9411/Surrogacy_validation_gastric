library(surrosurv)
library(survival)
library(survminer)
library(survMisc)
library(rstpm2)
library(bpcp)
data("gastadv")
gastadv$trt<-ifelse(gastadv$trt==0.5,1,0)
gastadv$timeS<-gastadv$timeS/30
gastadv$timeT<-gastadv$timeT/30


index<-0
for(i in 1:length(unique(gastadv$trialref))){
  index<-index+1
  ldf<-gastadv[gastadv$trialref==unique(gastadv$trialref)[i],]
  mod<-paste("mod",index,sep="")
  assign(mod,stpm2(Surv(timeT,statusT==1)~trt, data=ldf, smooth.formula=~ns(log(timeT),df=2)+log(timeT):trt))
}


plot(mod2, newdata = data.frame(trt = 0), type = "hr",
     var = "trt", ci = FALSE, rug = FALSE,
     main = "Overall Survival", xlim=c(0,26), ylim=c(0,3),
     ylab = "Hazard ratio", xlab = "Time (months)",line.col="grey", xaxt="n",yaxt="n"
)#p=0.772
axis(1, at = c(6, 12, 18, 24))
axis(2, at = c(1,2,3))

lines(mod3, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.445
lines(mod5, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.939
lines(mod6, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.461
lines(mod7, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.527
lines(mod10, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.721
lines(mod13, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.945
lines(mod15, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.749
lines(mod16, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.932
lines(mod20, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.143
lines(mod8, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.174
lines(mod12, newdata = data.frame(trt = 0), type = "hr",
      var = "trt")#p=0.161
lines(mod9, newdata = data.frame(trt = 0), type = "hr",
      var = "trt")#p=0.629
lines(mod11, newdata = data.frame(trt = 0), type = "hr",
      var = "trt")#p=0.152
lines(mod14, newdata = data.frame(trt = 0), type = "hr",
      var = "trt")#p=0.107
lines(mod17, newdata = data.frame(trt = 0), type = "hr",
      var = "trt")#p=0.113
lines(mod18, newdata = data.frame(trt = 0), type = "hr",
      var = "trt")#p=0.242
lines(mod19, newdata = data.frame(trt = 0), type = "hr",
      var = "trt")#p=0.247
lines(mod1, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="red")#p=0.065
lines(mod4, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="red")#p=0.042


index<-0
for(i in 1:length(unique(gastadv$trialref))){
  index<-index+1
  ldf<-gastadv[gastadv$trialref==unique(gastadv$trialref)[i],]
  mod<-paste("mod",index,sep="")
  assign(mod,stpm2(Surv(timeS,statusS==1)~trt, data=ldf, smooth.formula=~ns(log(timeS),df=2)+log(timeS):trt))
}

plot(mod3, newdata = data.frame(trt = 0), type = "hr",
     var = "trt", ci = FALSE, rug = FALSE,
     main = "Progression-free Survival", xlim=c(0,26), ylim=c(0,3),
     ylab = "Hazard ratio", xlab = "Time (months)",line.col="grey",xaxt="n",yaxt="n"
)#p=0.880
axis(1, at = c(6, 12, 18, 24))
axis(2, at = c(1,2,3))
lines(mod4, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.260
lines(mod5, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.914
lines(mod6, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.580
lines(mod7, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.273
lines(mod10, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.723
lines(mod11, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.827
lines(mod13, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.791
lines(mod14, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.735
lines(mod16, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.596
lines(mod18, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.285
lines(mod2, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.186
lines(mod17, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="grey")#p=0.192
lines(mod9, newdata = data.frame(trt = 0), type = "hr",
      var = "trt")#p=0.631
lines(mod19, newdata = data.frame(trt = 0), type = "hr",
      var = "trt")#p=0.120
lines(mod20, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="red")#p=0.000
lines(mod1, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="red")#p=0.039
lines(mod8, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="red")#p=0.089
lines(mod12, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="red")#p=0.042
lines(mod15, newdata = data.frame(trt = 0), type = "hr",
      var = "trt",col="red")#p=0.021


