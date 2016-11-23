#KFAS example from:
# https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=24&cad=rja&uact=8&ved=0ahUKEwiLj7rf59_OAhUkIMAKHUUABGg4FBAWCCowAw&url=http%3A%2F%2F54.225.166.221%2Fbehzod%2F182646&usg=AFQjCNEDhl4wNY0gcbmSDYdvn-raTmvhCA


library(Quandl)
library(KFAS)
library(dplyr)
library(ggplot2)
library(lubridate)
library("ggfortify")
library(xts)

cpiauc <- Quandl("FRED/CPIAUCNS", type="zoo", start_date="1955-01-01", end_date="2013-12-31")
str(cpiauc)

lcpiauc <- log(cpiauc)
y <- diff(lcpiauc)
plot(y, xlab="", ylab=expression(""*y[t]*""), main=expression("Month-over-month inflation rate "*y[t]*""))


T <- length(y)
y.model <- SSModel(y~SSMtrend(degree=1,Q=NA)+SSMseasonal(period=12,sea.type="dummy",Q=NA), H=NA)

y.model$T
y.model$R
y.model$Z
y.model$Q
y.model$H

y.updatefn <- function(pars, model) {
  model$Q[,,1] <- diag(exp(pars[1:2]))
  model$H[,,1] <- exp(pars[3])
  model
}

y.fit <- fitSSM(y.model, inits=log(c(0.1,0.01,0.001)), updatefn=y.updatefn, method="BFGS")
y.fit$optim.out$par

y.out <- KFS(y.fit$model)
colnames(y.out$alphahat)

y.KF <- y.out$a %*% as.vector(y.out$model$Z)
y.KS <- y.out$alphahat %*% as.vector(y.out$model$Z)

#find position of level and first seasonal dummy
ilvl <- which(colnames(y.out$a)=="level")
isea <- which(colnames(y.out$a)=="sea_dummy1")

#construct 90% confidence interval for smoothed state
y.KS.CI.lvl <- as.vector(y.out$alphahat[,ilvl]) + qnorm(0.90)* sqrt(cbind(y.out$V[ilvl,ilvl,])) %*% cbind(-1,1)
y.KS.CI.sea <- as.vector(y.out$alphahat[,isea]) + qnorm(0.90)* sqrt(cbind(y.out$V[isea,isea,])) %*% cbind(-1,1)


#log transformed earnings per share - actual, smoothed and seasonal component
plot.ts( cbind(y, y.out$alphahat[,ilvl], y.KS.CI.lvl), col=c("black","red","red","red"), lty=c(1,1,3,3), plot.type="single", xlab="", ylab="", main="actual and smoothed values" )


#seasonal component
plot.ts( cbind(y.out$a[-(T+1),isea], y.KS.CI.sea), col="blue", lty=c(1,3,3), plot.type="single", xlab="", ylab="", main="seasonal component" )
abline(h=0, lty=3)

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
#From Roy's Time-series workshop using KFAS



load("R:/Kate/Data/KFASData.RData")

#--------------------------------------------------------------------
#example of stationary time-series that might look non-stationary
model <- SSModel(matrix(NA, 100, 1) ~ SSMarima(ar=.8, Q=1), H = 1)
set.seed(10000)
sim <- simulateSSM(model, "obs", nsim = 1)
ts.plot(sim)
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#example of a stationary time-series that might look like the variance
# is changing
x = arima.sim(list(order=c(2,0,0), ar=c(1,-.9)), n=2^8)
ts.plot(x)
#---------------------------------------------------------------------

#----------------------------------------------------------------------
#example of local trends

#fit a linear model to the entire time series
autoplot(sst)
sstIndex <- (1:length(sst))
junk <- lm(sst ~ sstIndex)
summary(junk)

#now fit a linear model to an arbitrary 7 year stretch
sstIndex <- (1:length(sst[781:(781+84)]))
junk <- lm(sst[781:(781+84)] ~ sstIndex)
summary(junk)

#add 5 more years to the abitrary 7 year stretch
sstIndex <- (1:length(sst[781:(781+84+60)]))
junk <- lm(sst[781:(781+84+60)] ~ sstIndex)
summary(junk)
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
#working with seasonal data
autoplot(seaIce)

#estimate the quarterly mean
seaIceq <- array(NA_real_, dim=4)
qIndex <- array(NA_real_, dim=c(4,135))
for (iq in 1:4){
  qIndex[iq,]<-seq(from=iq, to = length(seaIce), by = 4 )
  seaIceq[iq] <- mean(seaIce[qIndex[iq,]])
}
#Compare the sea ice quarterly mean with the times series for that quarter.  

#Quarter 1 :
autoplot(seaIce[qIndex[1,]]) + geom_hline(yintercept=seaIceq[1])  

#quarterly mean doesn't give a good estimate of the seasonal cycle because of a pronounced trend in the
# data.

#Let's try with a montly series and a more subtle trend:
sstm <- array(NA_real_, dim=12)
mIndex <- array(NA_real_, dim=c(12,1740))
for (mon in 1:12){
  mIndex[mon,]<-seq(from=mon, to = 1740, by = 12 )
  sstm[mon] <- mean(sst[mIndex[mon,]])
}

#Comparing the April mean with the April time series:
autoplot(sst[mIndex[4,]]) + geom_hline(yintercept=sstm[4])

#Comparing the July mean with the July time series:
autoplot(sst[mIndex[7,]]) + geom_hline(yintercept=sstm[7])

#Particularly for July, we can see that we have not done a good job of estimating seasonal behavior.
# What we would like is a way of estimating the trend and seasonal together, where the trend is not forced to
# be linear, has local properties so it doesn't change dramatically as new data are added to the model, 
# and a seasonal term that is allowed to change over time.  This is what state-space decomposition allows:

#-----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#Setting up state-space model

#the basic local level model for trend is set up like this:
locLevel <- SSMtrend(degree = 1, Q=1., P1inf =1)
str(locLevel)

locLinear <- SSMtrend(degree=2,Q=list(matrix(1.),matrix(0.5)))
str(locLinear)


#illustrate what happens at extreme valus of the variance....when Q=0 this says
# there is no variance in the state equation
nino3Model <- SSModel(nino3 ~ SSMtrend(degree = 1, Q=0., P1inf =1), H = 0.2)
nino3Smooth <- KFS(nino3Model, smoothing = "state")
nino3Level <-signal(nino3Smooth, states = 'level')
nino3Level$signal <- xts(nino3Level$signal, order.by = ninoTimes)
autoplot(cbind(nino3, nino3Level$signal),facets = FALSE)

#Now for very large Q
nino3Model <- SSModel(nino3 ~ SSMtrend(degree = 1, Q=1000., P1inf =1), H = 0.2)
nino3Smooth <- KFS(nino3Model, smoothing = "state")
nino3Level <-signal(nino3Smooth, states = 'level')
nino3Level$signal <- xts(nino3Level$signal, order.by = ninoTimes)
autoplot(cbind(nino3, nino3Level$signal),facets = FALSE)

#in the first case we set up a model with a trend but no variance in the state
# equation...so the smoothed estimates at each observation are just a line...
# then we made the variance huge and the smoothed states are sitting right on top
# of the observations so we can't even see the second line


#Here is another example.  The series is known to be a trend plus error
nino3Model <- SSModel(nino3 ~ SSMtrend(degree = 1, Q=list(matrix(NA)), P1inf =1), H = matrix(NA))
nino3Fit <-  fitSSM(nino3Model, c(log(var(nino3)), log(var(nino3))),method = "BFGS")$model
nino3Smooth <- KFS(nino3Fit, filtering = "state", smoothing = "state")
nino3Level <-signal(nino3Smooth, states = 'level')
nino3Level$signal <- xts(nino3Level$signal, order.by = ninoTimes)
autoplot(cbind(nino3, nino3Level$signal),facets = FALSE)

#with KFAS 'NA' in the model set-up tells us these are parameters to be estimated:
nino3Model <- SSModel(nino3 ~ SSMtrend(degree = 1, Q=list(matrix(NA)), P1inf =1), H = matrix(NA))
str(nino3Model)



#example with some seasonals----------------------------------------
mySeason <- SSMseasonal(period=4, sea.type="dummy", Q = 1.)
str(mySeason)


testSeasonModel <- SSModel(seaIce2 ~ SSMseasonal(period=4, sea.type="dummy", Q = NA), H=NA)
testSeasonFit <- fitSSM(testSeasonModel, inits=c(1., 0.5))$model
testSeasonSmooth <- KFS(testSeasonFit, filtering = "state", smoothing = "state")
testSeason <- signal(testSeasonSmooth, states = 'all')
testSeason$signal <- xts(testSeason$signal, order.by=seaiceTimes)
autoplot(cbind(seaIce2,testSeason$signal), facets = FALSE)
#---------------------------------------------------------------------


#--------------------------------------------------------------------------
#Sea Ice Example with trend and seasonal
seaIceModel<-SSModel(seaIce ~ SSMtrend(degree = 1, Q=list(matrix(NA))) + 
                       SSMseasonal(period=4, sea.type="dummy", Q = NA), H = NA)
str(seaIceModel)

seaIceFit<-fitSSM(seaIceModel,inits=c(0.1,0.05, 0.001),method='BFGS')$model
seaIceSmooth <- KFS(seaIceFit,smooth= c('state', 'mean','disturbance'))
seaIceLevel <-signal(seaIceSmooth, states = 'level')
seaIceLevel$signal <- xts(seaIceLevel$signal, order.by = seaiceTimes)
autoplot(cbind(seaIce, seaIceLevel$signal),facets = FALSE)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################




df.monthly <- readRDS('R:/Kate/code/gf-catchshare-seasonality/data/gf_monthly.RDA')
df.monthly <- tbl_df(df.monthly) %>% filter(year<2016)

trips <- df.monthly$ntrips[df.monthly$area=='south']
#trips <- df.monthly$ntrips[df.monthly$area=='north']

#trips <- df.monthly$tripshare[df.monthly$area=='south']
df.monthly$date <- as.Date(paste(df.monthly$year,"-",df.monthly$month,"-","01",sep=""),format="%Y-%m-%d")
dates <- df.monthly$date[df.monthly$area=='north']
trips <- xts(trips,order.by=dates)
#trips <- xts(trips,order.by=dates)

df.monthly <- df.monthly %>% mutate(july=ifelse(month==7,1,0),oct=ifelse(month==10,1,0))

ggplot(subset(df.monthly,area=='south'),aes(x=date,y=ntrips)) + geom_line() + 
  geom_point(aes(color=factor(oct)))

ggplot(subset(df.monthly,area=='north'),aes(x=date,y=tripshare)) + geom_line() + 
  geom_point(aes(color=factor(july)))

tripsmodel<-SSModel(trips ~ SSMtrend(degree = 1, Q=list(matrix(NA))) + 
                       SSMseasonal(period=12, sea.type="dummy", Q = NA), H = NA)
str(tripsmodel)

tripsFit<-fitSSM(tripsmodel,inits=c(0.1,0.05, 0.001),method='BFGS')$model
tripsSmooth <- KFS(tripsFit,smooth= c('state', 'mean','disturbance'))

tripsLevel <-signal(tripsSmooth, states = 'level')

tripsLevel$signal <- xts(tripsLevel$signal, order.by = dates)

autoplot(cbind(trips, tripsLevel$signal),facets = FALSE)


#plot fitted values with raw data
pred <- data.frame(date=dates,kfas=fitted(tripsSmooth),y=as.numeric(trips))

ggplot(pred,aes(x=date,y=y)) + geom_line() + geom_line(aes(x=date,y=kfas,color='red'),data=pred)


#plot the seasonal
trips.sea <- as.numeric(tripsSmooth$alphahat[,2])
plot.df <- data.frame(date=dates,seasonal=trips.sea)


#add a reference point to the seasonal data frame
plot.df <- plot.df %>% mutate(month=month(date),july=ifelse(month==7,1,0),oct=ifelse(month==10,1,0))
  
ggplot(plot.df,aes(x=date,y=seasonal)) + geom_line() + geom_point(aes(color=factor(oct)))

#compare the seasonal signal from the state-space model to the raw data
df.monthly <- df.monthly %>% mutate(oct=ifelse(month==10,1,0))
ggplot(subset(df.monthly,area=='south'),aes(x=date,y=tripshare)) + geom_line() + 
  geom_point(aes(color=factor(oct)))


#plot level
level <- as.numeric(tripsSmooth$alphahat[,1])
plot.df <- data.frame(date=dates,level=level,z='smoothed')
                
ggplot(plot.df,aes(x=date,y=level)) + geom_line() + 
  geom_line(data=subset(df.monthly,area=='south'),aes(x=date,y=ntrips,color="red"))



#==========================================================================
#get mean absolute error from the KFAS model and compare it to the
# dummy variable model with structural break in 2010.

#the dummy variable model...I know this is clunky but I 
# prefer the traditional dummy variable set up.
chow.df <- df.monthly %>% filter(area=='north') %>%
            mutate(jan=ifelse(month==1,1,0),
                   feb=ifelse(month==2,1,0),
                   march=ifelse(month==3,1,0),
                   april=ifelse(month==4,1,0),
                   may=ifelse(month==5,1,0),
                   june=ifelse(month==6,1,0),
                   july=ifelse(month==7,1,0),
                   aug=ifelse(month==8,1,0),
                   sept=ifelse(month==9,1,0),
                   oct=ifelse(month==10,1,0),
                   nov=ifelse(month==11,1,0),
                   dec=ifelse(month==12,1,0))

model.pre <- lm(ntrips~feb+march+april+may+june+july+aug+
                  sept+oct+nov+dec+factor(year),data=subset(chow.df,year<=2010))
model.post <- lm(ntrips~feb+march+april+may+june+july+aug+
                   sept+oct+nov+dec+factor(year),data=subset(chow.df,year>2010))

yhat <- rbind(
data.frame(date=chow.df$date[chow.df$year<=2010],
              that=predict(model.pre,newdata=chow.df[chow.df$year<=2010,])),
data.frame(date=chow.df$date[chow.df$year>2010],
           that=predict(model.post,newdata=chow.df[chow.df$year>2010,]))
)

yhat$y <- as.numeric(trips)
yhat$model <- 'Seasonal Dummy'

ggplot(yhat,aes(x=date,y=y)) + geom_line() + geom_line(aes(x=date,y=that,color='red'))

# the KFAS model
kfas.pred <- data.frame(date=dates,that=fitted(tripsSmooth),y=as.numeric(trips),model='KFAS')

#combine the two
model.comp <- rbind(yhat,kfas.pred)

model.comp$eps <- model.comp$y-model.comp$that

ggplot(model.comp,aes(x=eps)) + geom_histogram() + facet_wrap(~model)
ggplot(model.comp,aes(x=date,y=eps,color=model)) + geom_line()


#get mean absolute in-sample prediction error:
tbl_df(model.comp) %>% mutate(abs.error=abs(eps)) %>% group_by(model) %>%
              summarise(mean(abs.error,na.rm=T))

#==========================================================================




model.pre <- lm(ntrips~feb+march+april+may+june+july+aug+
                  sept+oct+nov+dec+factor(year),data=subset(chow.df,year<=2009))
model.post <- lm(ntrips~feb+march+april+may+june+july+aug+
                   sept+oct+nov+dec+factor(year),data=subset(chow.df,year>2009))

yhat <- rbind(
  data.frame(date=chow.df$date[chow.df$year<=2009],
             that=predict(model.pre,newdata=chow.df[chow.df$year<=2009,])),
  data.frame(date=chow.df$date[chow.df$year>2009],
             that=predict(model.post,newdata=chow.df[chow.df$year>2009,]))
)

yhat$y <- as.numeric(trips)
yhat$model <- 'Seasonal Dummy'

model.comp <- rbind(yhat,kfas.pred)

model.comp$eps <- model.comp$y-model.comp$that

tbl_df(model.comp) %>% mutate(abs.error=abs(eps)) %>% group_by(model) %>%
  summarise(mean(abs.error,na.rm=T))




##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
#another simulation experiment

z <- data.frame(date=dates,level=tripsSmooth$alphahat[,1],state=tripsSmooth$alphahat[,2],pred=fitted(tripsSmooth))
z$level_plus_state<-z$level+z$state

dates <- ss.data$date
trips <- xts(ss.data$trips,order.by=dates)

tripsModel<-SSModel(trips ~ SSMseasonal(period=12, sea.type="dummy", Q = NA), H = NA)

tripsFit<-fitSSM(tripsModel,inits=c(0.1,0.001),method='BFGS')$model
tripsSmooth <- KFS(tripsFit,filtering='state',smoothing='state')
testtrips <- signal(tripsSmooth,state="all")
testtrips$signal <- xts(testtrips$signal,order.by=dates)
autoplot(cbind(trips,testtrips$signal),facets=FALSE)

#plot the seasonal
trips.sea <- as.numeric(tripsSmooth$alphahat[,2])
plot.df <- data.frame(date=dates,seasonal=trips.sea)
#add a reference point to the seasonal data frame
plot.df <- plot.df %>% mutate(month=month(date),july=ifelse(month==7,1,0))

ggplot(plot.df,aes(x=date,y=seasonal)) + geom_line() + geom_point(aes(color=factor(july)))

fitted <- tripsSmooth$alphahat[,1]+tripsSmooth$alphahat[,2]
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################













