##################################################
##################################################
##################################################
##################################################
#Script to analyze structural breaks with seasonal data
# under different conditions

library(strucchange)
library(dplyr)
library(ggplot2)
library(lubridate)
library(forecast)
##################################################
##################################################
##################################################
##################################################
##################################################

#first simulate a monthly time-series with a 
# constant level and structural break in year 
# 5 of a ten  year time series

monthly.df <- tbl_df(data.frame(year=rep(2000:2010,12),month=c(rep(1,11),rep(2,11),rep(3,11),
                                                               rep(4,11),rep(5,11),rep(6,11),
                                                               rep(7,11),rep(8,11),rep(9,11),
                                                               rep(10,11),rep(11,11),rep(12,11)))) %>%
  mutate(jan=ifelse(month==1,1,0),feb=ifelse(month==2,1,0),march=ifelse(month==3,1,0),
         april=ifelse(month==4,1,0),may=ifelse(month==5,1,0),june=ifelse(month==6,1,0),
         july=ifelse(month==7,1,0),aug=ifelse(month==8,1,0),sept=ifelse(month==9,1,0),
         oct=ifelse(month==10,1,0),nov=ifelse(month==11,1,0),dec=ifelse(month==12,1,0)) %>%
  mutate(beta1=ifelse(year<2005,10,20),
         beta2=ifelse(year<2005,20,40),
         beta3=ifelse(year<2005,30,60),
         beta4=ifelse(year<2005,40,80),
         beta5=ifelse(year<2005,50,100),
         beta6=ifelse(year<2005,60,120),
         beta7=ifelse(year<2005,70,140),
         beta8=ifelse(year<2005,80,160),
         beta9=ifelse(year<2005,90,180),
         beta10=ifelse(year<2005,100,200),
         beta11=ifelse(year<2005,25,50),
         beta12=ifelse(year<2005,15,30),e=rnorm(132,10,20)) %>%
  mutate(trips=(jan*beta1)+(feb*beta2)+(march*beta3)+(april*beta4)+(may*beta5)+
           (june*beta6)+(july*beta7)+(aug*beta8)+(sept*beta9)+
           (oct*beta10)+(nov*beta11)+(dec*beta12)+e) %>%
  mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d")) %>%
  arrange(year,month)

ggplot(monthly.df,aes(x=date,y=trips)) + geom_line() + geom_point(aes(color=factor(oct)))

f_statistics <- Fstats(ts(trips,frequency=12, start=c(2000, 1)) ~ factor(month) , data = monthly.df)
plot(f_statistics)
summary(breakpoints(ts(trips,frequency=12, start=c(2000, 1)) ~ factor(month), data = monthly.df))


#same pattern different intensity

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################

#same pattern with trend and no structural break

# this is a really simple pattern. there is no break, the coefficients are the same in every
# year but there is an annual time trend

monthly.df <- tbl_df(data.frame(year=rep(2000:2010,12),month=c(rep(1,11),rep(2,11),rep(3,11),
                                                               rep(4,11),rep(5,11),rep(6,11),
                                                               rep(7,11),rep(8,11),rep(9,11),
                                                               rep(10,11),rep(11,11),rep(12,11)))) %>%
  mutate(jan=1,feb=ifelse(month==2,1,0),march=ifelse(month==3,1,0),
         april=ifelse(month==4,1,0),may=ifelse(month==5,1,0),june=ifelse(month==6,1,0),
         july=ifelse(month==7,1,0),aug=ifelse(month==8,1,0),sept=ifelse(month==9,1,0),
         oct=ifelse(month==10,1,0),nov=ifelse(month==11,1,0),dec=ifelse(month==12,1,0)) %>%
  mutate(int=100,
         t=year-1999,
         beta1=2000-20*t,
         beta2=20,
         beta3=30,
         beta4=40,
         beta5=50,
         beta6=60,
         beta7=70,
         beta8=80,
         beta9=90,
         beta10=100,
         beta11=25,
         beta12=15,e=rnorm(132,20,20)) %>%
  mutate(trips=(jan*beta1)+(feb*beta2)+(march*beta3)+(april*beta4)+(may*beta5)+
           (june*beta6)+(july*beta7)+(aug*beta8)+(sept*beta9)+
           (oct*beta10)+(nov*beta11)+(dec*beta12)+e) %>%
  mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d")) %>%
  arrange(year,month)

ggplot(monthly.df,aes(x=date,y=trips)) + geom_line() + geom_point(aes(color=factor(oct)))

f_statistics <- Fstats(ts(trips,frequency=12, start=c(2000, 1)) ~ factor(month) , data = monthly.df)
plot(f_statistics)
summary(breakpoints(ts(trips,frequency=12, start=c(2000, 1)) ~ factor(month), data = monthly.df))

#breakpoint analysis of a model with trend
f_statistics <- Fstats(ts(trips,frequency=12, start=c(2000, 1)) ~ factor(month) + t, data = monthly.df)
plot(f_statistics)
summary(breakpoints(ts(trips,frequency=12, start=c(2000, 1)) ~ factor(month) + t, data = monthly.df))

#----------------------------------------------------------------------------------------------------
#breakpoint analysis of a model where observations are 1st differenced
monthly.1d <- monthly.df %>% ungroup() %>% arrange(year,month) %>% 
              mutate(first.diff=trips-lag(trips))

#is a seasonal model even still appropriate?
tsdisplay(diff(ts(monthly.1d$trips,frequency=12,start=c(2000,1))))
#note: looks like acf and pacf still have spikes around 12 even after first differencing so
# a seasonal model appears to still be appropriate

ggplot(monthly.1d,aes(x=date,y=first.diff)) + geom_line() + geom_point(aes(color=factor(oct)))


monthly.1d <- monthly.1d %>% filter(year>2000)
f_statistics <- Fstats(ts(first.diff,frequency=12, start=c(2001, 1)) ~ factor(month), data = monthly.1d)
plot(f_statistics)
summary(breakpoints(ts(first.diff,frequency=12, start=c(2001, 1)) ~ factor(month), data = monthly.1d))
#----------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------
# breakpoint analysis with seasonal difference
monthly.dsea <- monthly.df %>% ungroup() %>% arrange(year,month) %>% 
  mutate(sea.diff=trips-lag(trips,12))
#is a seasonal model even still appropriate?
tsdisplay(ts(monthly.dsea$sea.diff,frequency=12,start=c(2000,1)))
#note: looks like acf and pacf still have spikes around 12 even after first differencing so
# a seasonal model appears to still be appropriate

ggplot(monthly.dsea,aes(x=date,y=sea.diff)) + geom_line() + geom_point(aes(color=factor(oct)))


monthly.dsea <- monthly.dsea %>% filter(year>2000)

f_statistics <- Fstats(ts(sea.diff,frequency=12, start=c(2001, 1)) ~ factor(month), data = monthly.dsea)
plot(f_statistics)
summary(breakpoints(ts(sea.diff,frequency=12, start=c(2001, 1)) ~ factor(month), data = monthly.dsea))

#----------------------------------------------------------------------------------------------------


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
#   same pattern with trend and structural break

monthly.df <- tbl_df(data.frame(year=rep(2000:2010,12),month=c(rep(1,11),rep(2,11),rep(3,11),
                                                               rep(4,11),rep(5,11),rep(6,11),
                                                               rep(7,11),rep(8,11),rep(9,11),
                                                               rep(10,11),rep(11,11),rep(12,11)))) %>%
  mutate(jan=1,feb=ifelse(month==2,1,0),march=ifelse(month==3,1,0),
         april=ifelse(month==4,1,0),may=ifelse(month==5,1,0),june=ifelse(month==6,1,0),
         july=ifelse(month==7,1,0),aug=ifelse(month==8,1,0),sept=ifelse(month==9,1,0),
         oct=ifelse(month==10,1,0),nov=ifelse(month==11,1,0),dec=ifelse(month==12,1,0)) %>%
  mutate(int=100,
         t=year-1999,
         beta1=2000-20*t,
         beta2=ifelse(year<2005,20,40),
         beta3=ifelse(year<2005,30,60),
         beta4=ifelse(year<2005,40,80),
         beta5=ifelse(year<2005,50,100),
         beta6=ifelse(year<2005,60,120),
         beta7=ifelse(year<2005,70,140),
         beta8=ifelse(year<2005,80,160),
         beta9=ifelse(year<2005,90,180),
         beta10=ifelse(year<2005,100,200),
         beta11=ifelse(year<2005,25,50),
         beta12=ifelse(year<2005,15,30),e=rnorm(132,20,20)) %>%
  mutate(trips=(jan*beta1)+(feb*beta2)+(march*beta3)+(april*beta4)+(may*beta5)+
           (june*beta6)+(july*beta7)+(aug*beta8)+(sept*beta9)+
           (oct*beta10)+(nov*beta11)+(dec*beta12)+e) %>%
  mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d")) %>%
  arrange(year,month)

ggplot(monthly.df,aes(x=date,y=trips)) + geom_line() + geom_point(aes(color=factor(oct)))

#preliminary breakpoint analysis
f_statistics <- Fstats(ts(trips,frequency=12, start=c(2000, 1)) ~ factor(month), data = monthly.df)
plot(f_statistics)
summary(breakpoints(ts(trips,frequency=12, start=c(2000, 1)) ~ factor(month), data = monthly.df))

#note: a trending seasonal series with a structural break...breakpoint analysis has a hard time
#    with this one.  If we run a sequence of Chow tests we get the highest F-stat around late 2006


#------------------------------------------------------------------------------------------------
# breakpoint analysis controlling for trend by first differencing

monthly.1d <- monthly.df %>% ungroup() %>% arrange(year,month) %>% 
  mutate(first.diff=trips-lag(trips))

#is a seasonal model even still appropriate?
tsdisplay(diff(ts(monthly.1d$trips,frequency=12,start=c(2000,1))))
#note: looks like acf and pacf still have spikes around 12 even after first differencing so
# a seasonal model appears to still be appropriate

ggplot(monthly.1d,aes(x=date,y=first.diff)) + geom_line() + geom_point(aes(color=factor(oct)))


monthly.1d <- monthly.1d %>% filter(year>2000)
f_statistics <- Fstats(ts(first.diff,frequency=12, start=c(2001, 1)) ~ factor(month), data = monthly.1d)
plot(f_statistics)
summary(breakpoints(ts(first.diff,frequency=12, start=c(2001, 1)) ~ factor(month), data = monthly.1d))

#In this case, controlling for the trend finds the true breakpoint...if we use the F-stat...BIC is biased
# toward 0 breakpoints.
#---------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------
# breakpoint analysis using anova to find the right model
monthly.df <- monthly.df %>% ungroup() %>% arrange(year,month) %>% 
                mutate(trips.L1=lag(trips,1),trips.L2=lag(trips,2),
                       trips.L12=lag(trips,12))

monthly.df <- monthly.df %>% filter(year>2000)
#models
lm_s1<-lm(trips~factor(month),data=monthly.df)
lm_s2<-lm(trips~factor(month)+trips.L1,data=monthly.df)
lm_s3<-lm(trips~factor(month)+trips.L1+trips.L2,data=monthly.df)
lm_s4<-lm(trips~factor(month)+trips.L1+trips.L2+trips.L12,data=monthly.df)

anova(lm_s1,lm_s4)
anova(lm_s2,lm_s4)
anova(lm_s3,lm_s4)

f_statistics <- Fstats(ts(trips,frequency=12, start=c(2001, 1)) ~ factor(month)+trips.L1+trips.L2+trips.L12, data = monthly.df)
plot(f_statistics)
summary(breakpoints(ts(trips,frequency=12, start=c(2001, 1)) ~ factor(month)+trips.L1+trips.L2+trips.L12, data = monthly.df))

#It is encouraging here that the differenced model seems to identify a break around the right time.
#-------------------------------------------------------------------------------------------------




####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

#now we want to see how the breakpoint models treat a seasonal unit root:
library(tseries)

#---------------------------------------------------------------------------------------
#first let's get some intuition by looking at how the breakpoint analysis treats a plain 
# vanilla unit root with no seasonality:

#random walk
x <- rnorm(100)
rw  <- tbl_df(data.frame(y=cumsum(x),t=seq(1:length(x)))) %>% mutate(y.L1=lag(y),ydiff.1=y-y.L1)

ggplot(rw,aes(x=t,y=y)) + geom_line()

#Dickey Fuller Test of Stationarity
adf.test(rw$y)

#OLS based test for unit root
summary(lm(y~t+y.L1,data=rw))

#breakpoint analysis
f_statistics <- Fstats(y ~ 1, data = rw)
plot(f_statistics)
summary(breakpoints(y ~ 1, data = rw))

#model has a peak in the F-stats but BIC is smoothly decreasing as the number of possible 
# breakpoints increases...this is good news.  It seems to suggest that if we could add more
# breakpoints the BIC would continue to go down. This is what we expect because the random walk 
# has what look like breaks in a lot of spots.


#breakpoint analysis if we first difference
ggplot(rw,aes(x=t,y=ydiff.1)) + geom_line()
f_statistics <- Fstats(ydiff.1 ~ 1, data = rw)
plot(f_statistics)
summary(breakpoints(ydiff.1 ~ 1, data = rw))
#---------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#now a random walk with drift
drift=1
variance=40
rw <- rw %>% mutate(ydrift=cumsum(rnorm(n=100,mean=drift,sd=sqrt(variance)))) %>% 
            mutate(ydrift.L1=lag(ydrift))

ggplot(rw,aes(x=t,y=ydrift)) + geom_line()

adf.test(rw$ydrift)

summary(lm(ydrift~t+ydrift.L1,data=rw))

f_statistics <- Fstats(ydrift ~ 1, data = rw)
plot(f_statistics)
summary(breakpoints(ydrift ~ 1, data = rw))

#the BIC criteria likes a model with 2 breaks even though the data are actually generated from a 
# random walk with drift 
#----------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
#seasonal unit root
# from http://freakonometrics.hypotheses.org/13380

montreal=read.table("http://freakonometrics.free.fr/temp-montreal-monthly.txt")
 M=as.matrix(montreal[,2:13])
 X=as.numeric(t(M))
 tsm=ts(X,start=1948,freq=12)
 plot(tsm)


 Xp1=Xp2=as.numeric(t(M))
 for(t in 13:length(M)){
  Xp1[t]=Xp1[t-12]
   Xp2[t]=Xp2[t-12]+rnorm(1,0,2)
   }
 Xp1=Xp1+rnorm(length(Xp1),0,.02)
 tsp1=ts(Xp1,start=1948,freq=12)
 tsp2=ts(Xp2,start=1948,freq=12)

 
  par(mfrow=c(2,1))
 plot(tsp1)
 plot(tsp2)
#-------------------------------------------------------------------------------------


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################






#   same pattern with trend and structural break
#same pattern with seasonal unit root and no break
#   same pattern with seasonal unit root and structural break


#different pattern
#different pattern with trend

