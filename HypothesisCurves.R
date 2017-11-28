#########################
#Figure 1
###################

pdf(file = paste0("/Users/daisy/GoogleDrive/PhD/ENSO/Manuscript/AustralEcology/Figures/Hypothesis",as.Date(Sys.time()),".pdf"),
    width = 4, height = 5.5)
 par(mfrow = c(2,1),
     mar = c(3,4,3,1),
     yaxs="i",
     las=1)


#generate some data
LA<-data.frame(obs=rnorm(2000, 9, 3))
LA$Phase<-"La Nina"
x<-LA$obs
H <- hist(x,breaks=50,freq=TRUE,border="white",axes=FALSE,main="",xlab="",ylim=c(0,150),xlim=c(0,20))
dx <- (H$breaks[2]-H$breaks[1])
m  <- mean(x)
s  <- sd(x)
N <- length(x)
curve(N*dx*dnorm(x, mean=m, sd=s), add=TRUE)
Axis(side=1, labels=c("Jun","Aug","Oct","Nov","Jan"),
     at=c(0,5,10,15,20))
Axis(side=2)
abline(h=0,lwd=1.8)

qE<-quantile(x,probs=c(0.05,.5,.95))
abline(v=qE[1],lty=3)
#abline(v=qE[2],lty=2)
abline(v=qE[3],lty=2)


legend("topleft",legend=c("start","conclusion"),
       lty=c(3,2),bty="n")





#generate some data
LA<-data.frame(obs=rnorm(2000, 9, 3))
LA$Phase<-"La Nina"
x<-LA$obs
H <- hist(x,breaks=50,freq=TRUE,border="white",axes=FALSE,main="",xlab="",ylim=c(0,150),xlim=c(0,20))
dx <- (H$breaks[2]-H$breaks[1])
m  <- mean(x)
s  <- sd(x)
N <- length(x)
curve(N*dx*dnorm(x, mean=m, sd=s), add=TRUE)
Axis(side=1, labels=c("Jun","Aug","Oct","Nov","Jan"),
     at=c(0,5,10,15,20))
Axis(side=2)
abline(h=0,lwd=1.8)

# qE<-quantile(x,probs=c(0.05,.5,.95))
#
# abline(v=qE[1])
# abline(v=qE[2])
# abline(v=qE[3])

El<-data.frame(obs=rnorm(1400, 8.2, 2.2))
El$Phase<-"El Nino"
x<-El$obs
m  <- mean(x)
s  <- sd(x)
N <- length(x)
curve(N*dx*dnorm(x, mean=m, sd=s), add=TRUE,lty=2)


legend("topleft",legend=c(expression(La~Ni*tilde(n)*a),expression(El~Ni*tilde(n)*o)),
       lty=c(1,2),bty="n")


dev.off()

