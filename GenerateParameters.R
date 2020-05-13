###
#new parameters
numruns=10
rzero=c(1.5)# what R0 do you want - can be multiple
totalruns=length(rzero)*numruns # how many runs of each R0

#incu = rlnorm(numruns,meanlog=log(5.2),sdlog=0.35) # incubation period of around 5.2 days
#incu = rnorm(numruns,mean=5.2,sd=0.5)

incu = rep(5,numruns)
meanIP = 7.5-5.2  # mean infecitous period
IP = rnorm(numruns,mean=meanIP,sd=0.1) #
beta=c()
for (i in rzero){beta=c(beta,i/IP)}

IP1=rep(1,length(rzero))
IP2=IP-IP1
incu=rep(incu,length(rzero))
rzero=rep(rzero,each=numruns)

b=-(1/IP + 1/incu)
a=(beta/incu - (1/incu)*(1/IP))
c=-1
dt=log(2)*(-b+sqrt((b^2-4*a*c)))/(2*a)

hosp = rep(0.5,numruns)
crit = rep((1/9.5),numruns)

parameters3=data.frame(beta3=round(beta,3),beta4=round(beta,3),
                       Progress2=round(1/incu,3),
                       Progress3=round(1/IP1,3),
                       Progress4=round(1/IP2,3),
#                       HospitalRate = hosp,
#                       CriticalRate = crit,
                       r0=round(rzero,3),
                       dt=round(dt,3))


#write.table(as.matrix(parameters3[,1:5]),"~/GitHub/MetaWards/Runs/1Apr/ncovparams.csv",row.names=F,col.names =F,sep=",")

#write.csv(parameters3,"~/GitHub/MetaWards/Runs/1Apr/ncovparams_r0dt.csv")