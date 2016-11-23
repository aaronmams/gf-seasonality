library(dplyr)
library(forecast)
library(WaveletComp)
library(lubridate)
library(ggplot2)

#Aggregate data by day

daily <- tbl_df(readRDS("R:/Kate/Data/gf_main.RDA")) %>% filter(area!='vancouver' & year< 2016) %>%
  mutate(date=as.Date(date,format="%Y-%m-%d")) %>%
  mutate(trip_id=paste(vessel_id,"_",date,sep="")) %>%
  group_by(area,date) %>%
  summarise(ntrips=n_distinct(trip_id),total_weight=sum(weight),total_val=sum(val),nvessel=n_distinct(vessel_id))

all.days <- tbl_df(rbind(data.frame(date=seq(as.Date("1994-01-01"),as.Date("2015-12-31"),by="days"),
                                    trips=0,wt=0,val=0,boats=0,area="south"),
                         data.frame(date=seq(as.Date("1994-01-01"),as.Date("2015-12-31"),by="days"),
                                    trips=0,wt=0,val=0,boats=0,area="north")))

all.days <- all.days %>% left_join(daily,by=c('date','area')) %>% 
  mutate(total.trips=ifelse(is.na(ntrips),trips,ntrips),
         total.wt=ifelse(is.na(total_weight),wt,total_weight),
         total.val=ifelse(is.na(total_val),val,total_val),
         total.boats = ifelse(is.na(nvessel),boats,nvessel)) %>%
  mutate(year=year(date),month=month(date))
df.daily <- all.days %>% select(date,area,total.wt,total.val,total.trips,total.boats,year,month)

df.monthly <- df.daily %>% group_by(area,year,month) %>% 
  summarise(total.wt=sum(total.wt),
            total.val=sum(total.val),
            total.trips=sum(total.trips),
            total.boats=sum(total.boats)
  )

########################################Calculate period using spectral density from forecast package####################################################

dailyperiods<-df.daily %>%
  group_by(area) %>%
  summarise_each(funs(findfrequency))

monthlyperiods<-df.monthly %>% 
  group_by(area) %>% 
  summarise_each(funs(findfrequency))

dailyperiods_preIFQ<-df.daily %>%
  filter(date<as.Date("2011/1/1")) %>% 
  group_by(area) %>%
  summarise_each(funs(findfrequency))

monthlyperiods_preIFQ<-df.monthly %>% 
  filter(year<2011) %>% 
  group_by(area) %>% 
  summarise_each(funs(findfrequency)) #might not have a long enough time series pre-IFQ to do a good job of identifying the period?

dailyperiods_postIFQ<-df.daily %>%
  filter(date>as.Date("2011/1/1")) %>% 
  group_by(area) %>%
  summarise_each(funs(findfrequency))

monthlyperiods_postIFQ<-df.monthly %>% 
  filter(year>2011) %>% 
  group_by(area) %>% 
  summarise_each(funs(findfrequency))

########################################Find period using wavelets####################################################
df.monthly2 <- df2 %>% filter(area!='vancouver' & year>2002 & year< 2016) %>% #monthly df with date (as first of month) instead of year and month (better for plotting)
  mutate(date=as.Date(paste0("01-",substr(date,4,9)), format='%d-%b-%y')) %>% 
  group_by(date,area) %>%
  summarise(total_weight=sum(weight),total_val=sum(val),
            ntrips=n_distinct(ftid),nvessel=n_distinct(vessel_id)) 



northdailytrips.w<-analyze.wavelet(filter(df.daily,area=="north"), "ntrips",
                                   loess.span = 0.75, #not sure if this is the right value to detrend data, but changing it doesn't seem to change main resilts much
                                   dt = 1, dj = 1/250,
                                   lowerPeriod = 7,
                                   upperPeriod = 500,
                                   make.pval = T, n.sim = 10)

wt.image(northdailytrips.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),show.date = T,label.time.axis=T)

southdailytrips.w<-analyze.wavelet(filter(df.daily,area=="south"), "ntrips",
                                   loess.span = .75, #not sure if this is the right value to detrend data, but changing it doesn't seem to change main results much
                                   dt = 1, dj = 1/250,
                                   lowerPeriod = 7,
                                   upperPeriod = 500,
                                   make.pval = T, n.sim = 10)

wt.image(southdailytrips.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),show.date = T,label.time.axis=T)

#monthly
northmonthlytrips.w<-analyze.wavelet(filter(df.monthly2,area=="north"), "ntrips",
                                     loess.span = 0.75, #not sure if this is the right value to detrend data, but changing it doesn't seem to change main resilts much
                                     dt = 1, dj = 1/250,
                                     lowerPeriod = 1,
                                     upperPeriod = 24,
                                     make.pval = T, n.sim = 10)

wt.image(northmonthlytrips.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),show.date = T,label.time.axis=T)

southmonthlytrips.w<-analyze.wavelet(filter(df.monthly2,area=="south"), "ntrips",
                                     loess.span = 0.75, #not sure if this is the right value to detrend data, but changing it doesn't seem to change main resilts much
                                     dt = 1, dj = 1/250,
                                     lowerPeriod = 1,
                                     upperPeriod = 24,
                                     make.pval = T, n.sim = 10)

wt.image(southmonthlytrips.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),show.date = T,label.time.axis=T)


par(mfrow=c(1,1))
########################################Moving loess window STL plots####################################################

#get monthly trip shares to mess around with
df.monthly <- df.monthly %>% group_by(area,year) %>% mutate(yr.total=sum(total.trips),
                                                            mon.share=total.trips/yr.total)

plot.stl.daily<-function(Area,variable,Window){
  decomp<-stl(ts(unlist(filter(df.daily,area==Area)[,variable]),frequency=365, start=c(2003, 1)), s.window=Window) 
  plot(decomp)
}

plot.stl.monthly<-function(Area,variable,Window){
  decomp<-stl(ts(unlist(filter(df.monthly,area==Area)[,variable]),frequency=12, start=c(2003, 1)), s.window=Window) 
  plot(decomp)
}

plot.stl.monthly("south","ntrips",7)
plot.stl.monthly("north","ntrips",7)
plot.stl.monthly("north","mon.share",7)

#compare with unchanging seasonality
plot.stl.monthly("south","ntrips","periodic")
plot.stl.monthly("north","ntrips","periodic")

plot.stl.daily("south","total.trips",4)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

#Plot monthly seasonality but only for vessels that fishing before and after catch shares:


daily <- tbl_df(readRDS("R:/Kate/Data/gf_main.RDA")) %>% filter(area!='vancouver' & year< 2016) %>%
  mutate(date=as.Date(date,format="%Y-%m-%d"),month=month(date)) %>%
  mutate(trip_id=paste(vessel_id,"_",date,sep="")) %>%
  group_by(vessel_id,area,year,month) %>%
  summarise(ntrips=n_distinct(ftid),wt=sum(weight),val=sum(val))
 
#get the list of vessels in each area that fished after catch shares
catchshare.boats <- unique(daily$vessel_id[daily$year>=2011])

#create a data frame with only these vessels
noexit.trips <- daily %>% filter(vessel_id %in% catchshare.boats) %>%
                group_by(year,month,area) %>% summarise(ntrips=sum(ntrips)) %>%
                mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d"),
                       oct=ifelse(month==10,1,0),group='noexit') 

exit.trips <- daily %>% filter(!vessel_id %in% catchshare.boats) %>%
  group_by(year,month,area) %>% summarise(ntrips=sum(ntrips)) %>%
  mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d"),
         oct=ifelse(month==10,1,0),group='exit')

tmp <- rbind(noexit.trips,exit.trips)

ggplot(subset(tmp,area=='south'),aes(x=date,y=log(ntrips))) + geom_line(aes(color=factor(group))) + geom_point(aes(color=factor(oct)))


ggplot(subset(noexit.trips,area=='south'),aes(x=date,y=log(ntrips))) + geom_line() + geom_point(aes(color=factor(oct)))





############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
#plot cdf for trips in the south

df <- readRDS("R:/Kate/code/gf-catchshare-seasonality/data/gf_monthly.RDA")

df <- df %>% filter(year<2015) %>%
          group_by(area,year) %>% mutate(total.trips = sum(ntrips)) %>%
          arrange(area,year,month) %>%
      group_by(area,year) %>% mutate(px=cumsum(tripshare)) %>% 
      mutate(marker=ifelse(month==5,5,
                           ifelse(month==6,6,
                                  ifelse(month==7,7,0))))
             
          
ggplot(subset(df,area=='north'),aes(x=month,y=px)) + geom_line() + geom_point(aes(color=factor(marker))) + geom_hline(yintercept=0.5) + 
    facet_wrap(~year) + scale_x_continuous(breaks=c(1:12))


  
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################




#look at the monthly rankings from 2003 - 2010

trips.rank <- df.monthly %>% filter(year<2015) %>% group_by(area,year) %>% arrange(area,year,month) %>% 
                mutate(rank=dense_rank(desc(total.trips)))

trips.rank.q <- df.monthly %>% filter(year<2015) %>% mutate(q=ifelse(month %in% c(1,2,3),1,
                                               ifelse(month %in% c(4,5,6),2,
                                                      ifelse(month %in% c(7,8,9),3,4)))) %>%
              group_by(area,year,q) %>% summarise(total.trips=sum(total.trips)) %>%
              group_by(area,year) %>% arrange(area,year,q) %>% 
              mutate(rank=dense_rank(desc(total.trips)))


ggplot(subset(trips.rank,area=='south'),aes(x=year,y=rank)) + geom_line() + geom_point() + 
  geom_vline(xintercept=c(2003,2010)) + 
      facet_wrap(~month) + scale_x_continuous(breaks=c(1994:2014)) + scale_y_continuous(breaks=c(1:12)) +
  theme(axis.text.x=element_text(angle=90))

#plot quarterly
ggplot(subset(trips.rank.q,area=='north'),aes(x=year,y=rank)) + geom_line() + geom_point() + 
  geom_vline(xintercept=c(2003,2011)) + 
  facet_wrap(~q) + scale_x_continuous(breaks=c(1994:2014)) + scale_y_continuous(breaks=c(1:12)) 


#test rank correlation using spearman:
a <- trips.rank$total.trips[trips.rank$year==1994 & trips.rank$area=='north']
b <- trips.rank$total.trips[trips.rank$year==1995 & trips.rank$area=='north']
cor.test(a,b,method='spearman')


#plot monthly shares as deviations from a constant mean
ggplot(subset(df.monthly,area=='south'),aes(x=year,y=tripshare))+ geom_line() + geom_point() + facet_wrap(~month) + 
    geom_vline(xintercept=c(2003,2010)) + geom_hline(yintercept=0.1)



#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################


library(strucchange)
library(ggplot2)
library(RODBC)
library(dplyr)
library(lubridate)
require(data.table)
#source('data_file.R')
login <- read.csv('R:/Kate/code/login.csv')
pw <- as.character(login$pacfinpw)
uid <- as.character(login$pacfinuid)

#source('data_file.R')
df.monthly <- readRDS('data/gf_monthly.RDA')

#----------------------------------------------------------------------------------
#I'm making a small change to Kate's code because dplyr has a lag function built in

#add columns for levels lagged by 1 or 2 months (ie level 1 or 2 months before current month)
#then remove 2003, make a monthly dummy variable column
df.monthly<-df.monthly%>% filter(year<2016) %>%
  group_by(area)%>%
  mutate(total_weight1=lag(total_weight,1),
         total_weight2=lag(total_weight,2),
         total_val1=lag(total_val,1),
         total_val2=lag(total_val,2),
         ntrips1=lag(ntrips,1),
         ntrips2=lag(ntrips,2),
         lnwt=log(total_weight),
         lnval=log(total_val),
         lntrips=log(ntrips),
         dlogwt=lnwt-lag(lnwt),
         dlogwt1=lag(dlogwt,1),
         dlogwt2=lag(dlogwt,2),
         dlogval=lnval-lag(lnval),
         dlogval1=lag(dlogval,1),
         dlogval2=lag(dlogval,2),
         dlogtrip=lntrips-lag(lntrips),
         dlogtrip1=lag(dlogtrip,1),
         dlogtrip2=lag(dlogtrip,2),
         wtshare1=lag(wtshare,1),
         wtshare2=lag(wtshare,2),
         valshare1=lag(valshare,1),
         valshare2=lag(valshare,2),
         tripshare1=lag(tripshare,1),
         tripshare2=lag(tripshare,2))%>%
  filter(year>2003)%>%
  mutate(fmonth=factor(month))

daily <- tbl_df(readRDS("R:/Kate/Data/gf_main.RDA")) %>% filter(area!='vancouver') %>%
            mutate(date=as.Date(date,format="%Y-%m-%d")) %>%
            mutate(trip_id=paste(vessel_id,"_",date,sep="")) %>%
            group_by(area,date) %>%
            summarise(ntrips=n_distinct(trip_id),weight=sum(weight),value=sum(val))
all.days <- tbl_df(rbind(data.frame(date=seq(as.Date("2001-01-01"),as.Date("2016-06-23"),by="days"),
                              trips=0,wt=0,val=0,area="south"),
                         data.frame(date=seq(as.Date("2001-01-01"),as.Date("2016-06-23"),by="days"),
                                    trips=0,wt=0,val=0,area="north")))

all.days <- all.days %>% left_join(daily,by=c('date','area')) %>% 
  mutate(total.trips=ifelse(is.na(ntrips),trips,ntrips),
         total.wt=ifelse(is.na(weight),wt,weight),
         total.val=ifelse(is.na(value),val,value)) %>%
  mutate(year=year(date),month=month(date))
daily <- all.days

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
#basic summary stats
nvess <- tbl_df(readRDS("R:/Kate/Data/gf_main.RDA")) %>% group_by(year,area) %>% 
          summarise(nvess=n_distinct(vessel_id)) %>% filter(area != 'vancouver')

ggplot(subset(nvess,year<2011),aes(x=year,y=nvess)) + geom_point() + 
        geom_smooth(method="lm") + facet_wrap(~area)
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
# first we want to visualize what the seasonality is

#boxplot of daily data broken into monthly bins
boxplot.raw <- daily %>% filter(year>=2004 & year <=2010) %>%
                mutate(month=month(date))
ggplot(boxplot.raw,aes(x=factor(month),y=total.val)) + geom_boxplot() + 
    facet_wrap(~area,scales='free')

line.raw <- daily %>% filter(year>=2004 & year <=2014) %>%
  mutate(month=month(date))
ggplot(line.raw,aes(x=date,y=total.trips)) + geom_line() + facet_wrap(~area,scales="free") 
  

#now detrend the daily data and 



#line plot of monthly values to see if trend exists
df.monthly <- df.monthly %>% mutate(date=as.Date(paste(year,"-",month,"-","01",sep="")))
ggplot(subset(df.monthly,year<2015),aes(x=date,y=total_weight)) + geom_line()  +
      facet_wrap(~area,scales="free")

#both cases seem to have a little trend, let's remove it by plotting only 
# ols residuals
df.pre <- df.monthly %>% filter(year >= 2004 & year<=2010)
df.pre.south <- df.pre %>% filter(area=='south')
m.south <- lm(total_val~year,data=df.pre.south)
df.pre.south$e <- resid(m.south)
ggplot(df.pre.south,aes(x=date,y=e)) + geom_line() + geom_point() + geom_smooth(method="lm")

df.pre <- df.monthly %>% filter(year >= 2004 & year<=2010)
df.pre.north <- df.pre %>% filter(area=='north')
m.north <- lm(total_val~year,data=df.pre.north)
df.pre.north$e <- resid(m.north)
ggplot(df.pre.north,aes(x=date,y=e)) + geom_line() + geom_point() + geom_smooth(method="lm")


ggplot(subset(df.pre),aes(x=date,y=resid)) + geom_line() + geom_smooth(method="lm") + 
    facet_wrap(~area,scales="free")


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

#difference data from monthly means
df.test <- daily %>% filter(year>2006 & year<2015) %>%
            mutate(itq=ifelse(year>2010,1,0)) %>%
            group_by(area,month,itq) %>%
            mutate(month.mean=mean(total.trips),mean.ratio=total.trips/month.mean)
ggplot(df.test,aes(x=date,y=mean.ratio)) + geom_line() + facet_wrap(~area)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

#make a CDF for leavers and stayers

exit.df <- tbl_df(readRDS("R:/Kate/Data/gf_main.RDA")) %>% filter(area!='vancouver') %>%
      mutate(trip_id=paste(vessel_id,"_",date,sep="")) %>%      
      group_by(vessel_id,year,month,area) %>% 
      summarise(trips=n_distinct(trip_id)) 

#find the leavers and stayers
vessels <- tbl_df(readRDS("R:/Kate/Data/gf_main.RDA")) %>% filter(area!='vancouver') 

v.exit <- data.frame(rbindlist(lapply(unique(vessels$vessel_id),function(x){
  count.pre <- 0
  count.post <- 0
  npre <- 0
  npost <- 0
  count.pre <- length(unique(vessels$ftid[vessels$vessel_id==x & vessels$year<2011]))
  count.post <- length(unique(vessels$ftid[vessels$vessel_id==x & vessels$year>2010]))
  if(count.pre>npre){npre=count.pre}
  if(count.post>npost){npost=count.post}
  return(data.frame(vessel_id=x,npre=npre,npost=npost))
})))  
  
exit.df <- exit.df %>% left_join(v.exit,by='vessel_id') %>% 
              mutate(v.type = ifelse(npost==0,'exit','fish')) 

#
cdf <- exit.df %>% filter(year==2005) %>%
            group_by(v.type,area) %>% mutate(yr.total=sum(trips,na.rm=T)) %>%
            group_by(area,v.type,month) %>% arrange(area,v.type,month) %>% 
            summarise(ntrips=sum(trips),yr.total=max(yr.total),share=ntrips/yr.total) %>%
            group_by(area,v.type) %>% arrange(area,v.type,month) %>% mutate(px=cumsum(share))
            
ggplot(cdf,aes(x=month,y=px,color=v.type)) + geom_line() + geom_point() + 
  facet_wrap(~area) + 
    scale_x_continuous(breaks=c(0:12))

#cdf of trips/month for non-exiting boats only
cdf.fish <- exit.df %>%
  group_by(v.type,area,ifq) %>% mutate(yr.total=sum(trips,na.rm=T)) %>%
  group_by(area,ifq,v.type,month) %>% arrange(area,v.type,month) %>% 
  summarise(ntrips=sum(trips),yr.total=max(yr.total),share=ntrips/yr.total) %>%
  group_by(area,ifq,v.type) %>% arrange(area,ifq,v.type,month) %>% mutate(px=cumsum(share))

ggplot(subset(cdf.fish,v.type=="fish"),aes(x=month,y=px,color=ifq)) + geom_line() + geom_point() + 
  facet_wrap(~area) + 
  scale_x_continuous(breaks=c(0:12))

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#Final tables for the paper:

#test the daily data for seasonality
channel<- odbcConnect(dsn="pacfin",uid=paste(uid),pw=paste(pw),believeNRows=FALSE)
daily.df <- sqlQuery(channel,paste("select VESSEL_NUM, LANDING_DATE, PACFIN_PORT_CODE, INPFC_AREA_TYPE_CODE", 
                             "from PACFIN_MARTS.COMPREHENSIVE_FT", 
                             "where LANDING_YEAR > 2003 AND LANDING_YEAR < 2011 AND FLEET_CODE='LE'"))
close(channel)

#==========================================================================
#try the periodogram for a single year
daily <- tbl_df(readRDS("R:/Kate/Data/gf_main.RDA")) %>% filter(area!='vancouver') %>%
          mutate(trip_id=paste(date,"_",vessel_id,sep=""),
              date=as.Date(date,format="%Y-%m-%d")) %>%
          group_by(date,area) %>%
          summarise(ntrips=n_distinct(trip_id))
all.days <- tbl_df(data.frame(date=seq(as.Date("2001-01-01"),as.Date("2016-06-23"),by="days"),trips=0))
all.days <- all.days %>% left_join(daily,by='date') %>% mutate(total.trips=ifelse(is.na(ntrips),trips,ntrips)) %>%
              mutate(year=year(date))

#plot daily trips in the north and south from 2004 - 2010
rects1 <- data.frame(xstart = as.Date('2007-03-01'), 
                    xend = as.Date('2007-03-30'))
#plot daily trips in the north and south from 2004 - 2010
rects2 <- data.frame(xstart = as.Date('2007-07-01'), 
                    xend = as.Date('2007-07-30'))

ggplot() + 
  geom_rect(data = rects1, aes(xmin = xstart, xmax = xend, 
                              ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_rect(data= rects2, aes(xmin = xstart, xmax = xend, 
                              ymin = -Inf, ymax = Inf), alpha = 0.4) + 
  geom_line(data = subset(all.days,area=="north" & date>="2007-01-01" & date<="2007-12-31"),
            aes(x=date,y=total.trips))

ggplot(subset(all.days,area=="south" & date>="2008-01-01" & date<="2014-12-31"),
       aes(x=date,y=total.trips)) + geom_line() 


#-------------------------------------------------------------------
#difference data from the seasonal mean
mean.diff <- all.days %>% mutate(itq=ifelse(year<2011,"pre","post"),month=month(date)) %>%
              group_by(area,month,itq) %>%
              mutate(seasonal.mean=mean(total.trips),
                     trip.md=total.trips-seasonal.mean) 
ggplot(subset(mean.diff,area=="south" & date>="2008-01-01" & date <= "2014-12-31"),
              aes(x=date,y=trip.md)) + geom_line()

#-------------------------------------------------------------------


#--------------------------------------------------------------------
#A KS test type illustration
ks.test <- all.days %>% filter(area=="south" & year %in% c(2004,2005,2006,2011)) %>%
            mutate(month=month(date)) %>%
            group_by(year,month) %>% 
            summarise(nt=sum(total.trips)) %>%
            group_by(year) %>%
            mutate(yr.total=sum(nt),frac=nt/yr.total,cumsum=cumsum(frac))

ggplot(ks.test,aes(x=month,y=cumsum,color=factor(year)))+ geom_line()

#--------------------------------------------------------------------

daily <- daily %>% filter(area!='vancouver' & year>2006 & year< 2008) %>% 
          mutate(trip_id=paste(date,"_",vessel_id,sep=""),
                 date=as.Date(date,format="%Y-%m-%d")) %>%
          group_by(date,area) %>%
          summarise(ntrips=n_distinct(trip_id)) %>% filter(area=="south")

#generate a list of all dates between 2004-01-01 and 2010-12-31
days <- tbl_df(data.frame(date=seq(as.Date("2007-01-01"),as.Date("2007-12-31"),by="days"),trips=0))
days <- days %>% left_join(daily,by='date') %>% mutate(total.trips=ifelse(is.na(ntrips),trips,ntrips))
ggplot(days,aes(x=date,y=total.trips)) + geom_line()
#use Welch method to compute periodogram
w <- pwelch(days$total.trips,fs=length(days$total),plot=T)



trips <- ts(data=days$ntrips,start=c(2004,1,1),frequency=365.25)
plot.ts(trips)
#noting no clear evidence of trend in the series let's plot the periodogram
trips.pre <- spec.pgram(trips,log="yes")
1/trips.pre$freq[trips.pre$spec==max(trips.pre$spec)]


#let's plot the periodogram to see what (if any) type of seaonal frequency
# we are dealing with
trips.pre <- spec.pgram(df.monthly$ntrips[df.monthly$area=='south' & df.monthly$year<2011],log='yes')
ggplot(data.frame(freq=trips.pre$freq,spec=trips.pre$spec),aes(x=freq,y=spec)) + geom_line() + theme_bw()
1/trips.pre$freq[trips.pre$spec==max(trips.pre$spec)]

trips.pre <- spec.pgram(df.monthly$ntrips[df.monthly$area=='north' & df.monthly$year<2011],log='yes')
ggplot(data.frame(freq=trips.pre$freq,spec=trips.pre$spec),aes(x=freq,y=spec)) + geom_line() + theme_bw()
1/trips.pre$freq[trips.pre$spec==max(trips.pre$spec)]
#====================================================================================================

#--------------------------------------------------------------
#test if monthly dummy variable coefficients are significant in
# the raw data before 2011

raw.trips <- tbl_df(readRDS("R:/Kate/Data/gf_main.RDA")) %>% 
              mutate(trip_id=paste(date,"_",vessel_id,sep="")) %>% 
              group_by(area,year,date) %>% summarise(ntrip=n_distinct(trip_id)) %>%
              mutate(month=month(date)) %>% mutate(date=as.Date(date))

#merge these data with data containing 0's for days where nobody fished
all.days <- tbl_df(data.frame(date=seq(as.Date("2001-01-01"),as.Date("2016-06-23"),by="days"),trips=0))
all.days <- all.days %>% left_join(raw.trips,by='date') %>% mutate(total.trips=ifelse(is.na(ntrip),trips,ntrip)) %>%
              filter(area!='vancouver' & year > 2003 & year < 2011)

summary(lm(total.trips~factor(month)+factor(year)+factor(area),data=all.days))  

#--------------------------------------------------------------


#see if aggregating the data quarterly would produce a quarterly frequency
df.quarter <- df.monthly %>% filter(area=='south' & year < 2011) %>%
              mutate(Q=ifelse(month %in% c(1,2,3),1,
                              ifelse(month %in% c(4,5,6),2,
                                     ifelse(month %in% c(7,8,9),3,4)))) %>%
              group_by(year,Q) %>%
              summarise(ntrips=sum(ntrips))

trips.pre <- spec.pgram(df.quarter$ntrips[df.quarter$year<2011],log='yes')
ggplot(data.frame(freq=trips.pre$freq,spec=trips.pre$spec),aes(x=freq,y=spec)) + geom_line() + theme_bw()
1/trips.pre$freq[trips.pre$spec==max(trips.pre$spec)]


#plot monthly trip shares with a horizontal line for the preITQ mean

preITQmeans <- df.monthly %>% group_by(area,month) %>% filter(year < 2011) %>% summarise(preITQmean=mean(tripshare))
df.monthly <- df.monthly %>% left_join(preITQmeans,by=c('area','month')) %>% mutate(itq=ifelse(year<2011,'preITQ','postITQ'))

ggplot(subset(df.monthly,area=='south'),aes(x=year,y=tripshare,shape=itq)) + geom_point() + 
    geom_line(aes(x=year,y=preITQmean)) + geom_vline(xintercept=2011) + facet_wrap(~month)

ggplot(subset(df.monthly,area=='south'),aes(x=year,y=lntrips,color=itq)) + geom_point() + geom_vline(xintercept=2011) +
facet_wrap(~month)

#plot log trips with the mean value superimposed...use the pre-itq mean and post-itq mean
itqmeans <- df.monthly  %>% mutate(itq=ifelse(year<2011,'preITQ','postITQ')) %>%
            group_by(area,month,itq) %>% summarise(meanvals = mean(lntrips,na.rm=T))

plot.df <- df.monthly %>% select(year,month,area,lntrips,itq) %>% left_join(itqmeans,by=c('area','month','itq'))
ggplot(subset(plot.df,area=='north'),aes(x=year,y=lntrips,shape=itq)) + geom_point() + 
  geom_line(aes(x=year,y=meanvals)) + facet_wrap(~month)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#SOME EXPLORATORY PLOTS USING PACFIN FISH TICKETS

#-------------------------------------------------------------------
#make a table of the number of groundfish boats in our study and 
# their participation in other fisheries

login<-read.csv('R:/Kate/code/login.csv')
uid <- as.character(login$pacfinuid)
pw <- as.character(login$pacfinpw)

channel<- odbcConnect(dsn="pacfin",uid=paste(uid),pw=paste(pw),believeNRows=FALSE)
#odbcGetInfo(channel)
t <- Sys.time()
gf.allfishing <- sqlQuery(channel,paste("select LANDING_YEAR, VESSEL_ID, MANAGEMENT_GROUP_CODE, sum(LANDED_WEIGHT_LBS)",
                             "from PACFIN_MARTS.COMPREHENSIVE_FT",
                             "where LANDING_YEAR > 1993",
                             "GROUP BY VESSEL_ID, MANAGEMENT_GROUP_CODE, LANDING_YEAR"))

Sys.time() - t
close(channel)
names(gf.allfishing) <- c('year','vessel_id','fishery','wt')

boats <- tbl_df(readRDS(file="R:/Kate/Data/gf_main.RDA")) %>% filter(area!='vancouver') %>%
          group_by(year,vessel_id,area) %>%
          summarise(count.gf=n_distinct(ftid))
          

gf.allfishing <- tbl_df(gf.allfishing) %>% inner_join(boats,by=c('vessel_id','year'))

fishery.tbl <- gf.allfishing %>% group_by(year,fishery,area) %>%
                summarise(count=n_distinct(vessel_id))

#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#plot crab landings by month for groundfish vessels
channel<- odbcConnect(dsn="pacfin",uid=paste(uid),pw=paste(pw),believeNRows=FALSE)
#odbcGetInfo(channel)
t <- Sys.time()
crab.df <- sqlQuery(channel,paste("select LANDING_YEAR, LANDING_MONTH, VESSEL_ID, sum(LANDED_WEIGHT_LBS)",
                                        "from PACFIN_MARTS.COMPREHENSIVE_FT",
                                        "where LANDING_YEAR > 1993",
                                        "and",
                                        "MANAGEMENT_GROUP_CODE = 'CRAB'",
                                        "GROUP BY VESSEL_ID, LANDING_YEAR, LANDING_MONTH"))

Sys.time() - t
close(channel)

names(crab.df) <- c('year','month','vessel_id','crab_wt')
boats <- tbl_df(readRDS(file="R:/Kate/Data/gf_main.RDA")) %>% filter(area!='vancouver') %>%
            mutate(month=month(date)) %>% 
            group_by(year,month,vessel_id,area) %>% 
            summarise(gf_wt=sum(weight)) %>%
            inner_join(crab.df,by=c('year','month','vessel_id')) %>%
            group_by(year,month,area) %>%
            summarise(gf_wt=sum(gf_wt),crab_wt=sum(crab_wt)) %>%
            mutate(date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d"))


#plot monthly crab landings
rects = data.frame(xstart=c(as.Date('2000-11-01'),
                            as.Date('2001-11-01'),
                            as.Date('2002-11-01'),
                            as.Date('2003-11-01'),
                            as.Date('2004-11-01'),
                              as.Date('2005-11-01'),
                            as.Date('2006-11-01'),
                               as.Date('2007-11-01'),
                            as.Date('2008-11-01'),
                                as.Date('2009-11-01'),
                            as.Date('2010-11-01')),
                   xend=c(as.Date('2001-02-01'),
                          as.Date('2002-02-01'),
                          as.Date('2003-02-01'),
                          as.Date('2004-02-01'),
                           as.Date('2005-02-01'),
                            as.Date('2006-02-01'),
                          as.Date('2007-02-01'),
                             as.Date('2008-02-01'),
                          as.Date('2009-02-01'),
                              as.Date('2010-02-01'),
                               as.Date('2011-02-01')))
ggplot(subset(boats,area=='south' & year <= 2011 & year >= 2000)) +
      geom_line(aes(x=date,y=crab_wt))  + geom_point(aes(x=date,y=crab_wt)) + 
      geom_rect(data=rects,aes(xmin=xstart,xmax=xend,ymin=-Inf,ymax=Inf),fill='blue',alpha=0.4)
  
  (subset(boats,area=='north' & year <= 2011 & year >= 2000),aes(x=date,y=crab_wt)) + geom_line() + geom_point() +
  scale_x_date(date_breaks="6 month") + 
  theme(axis.text.x=element_text(angle=90))











#--------------------------------------------------------------------
#vessels fishing south of 40'10 revenue shares for GF and other stuff
channel<- odbcConnect(dsn="pacfin",uid=paste(uid),pw=paste(pw),believeNRows=FALSE)
odbcGetInfo(channel)

df <- sqlQuery(channel,paste("select VESSEL_ID, LANDING_YEAR, THOMSON_FISHERY_CODE, PACFIN_PORT_CODE, sum(EXVESSEL_REVENUE)", 
              "from PACFIN_MARTS.COMPREHENSIVE_FT", 
              "where LANDING_YEAR > 2002",
              "group by VESSEL_ID, LANDING_YEAR, THOMSON_FISHERY_CODE, PACFIN_PORT_CODE"))

vessels <- sqlQuery(channel,paste("select distinct VESSEL_ID", 
                                  "from PACFIN_MARTS.COMPREHENSIVE_FT", 
                                  "where LANDING_YEAR > 2010 and FLEET_CODE='LE'"))

  
close(channel)
names(df) <- c('vessel','year','fcode','port','rev')


homeport <- tbl_df(df) %>% group_by(vessel,port) %>% summarise(portcount=n()) %>%
            arrange(vessel,-portcount) %>% group_by(vessel) %>% filter(row_number()==1)
names(homeport) <- c('vessel','homeport','count')

df <- tbl_df(df) %>% filter(vessel %in% vessels$VESSEL_ID) %>%
        mutate(fishery=ifelse(fcode %in% c(6,7,8,9),'gf',
                              ifelse(fcode %in% c(1),'crab',
                                     ifelse(fcode %in% c(3,4), 'shrimp','other')))) %>%
        inner_join(homeport,by='vessel') %>%
        filter(homeport %in% c('ERK','BRG','SF','PRN','MNT','MOS','MRO','AVL','SB','VEN')) %>%
        group_by(vessel,year,fishery) %>% 
        summarise(rev=sum(rev)) %>%
        group_by(vessel,year) %>% mutate(total=sum(rev),share=rev/total)

ggplot(subset(df,fishery=='gf' & year<2015),aes(x=factor(year),y=share)) + geom_boxplot()

#--------------------------------------------------------------------


#----------------------------------------------------------------------
#get a list of all vessels landing IFQ groundfish South of 40 10' and
# use that list to plot crab landings by month and year
channel<- odbcConnect(dsn="pacfin",uid=uid,pw=pw,believeNRows=FALSE)

df <- sqlQuery(channel,paste("select VESSEL_ID, LANDING_YEAR, LANDING_MONTH, LANDED_WEIGHT_LBS", 
                             "from PACFIN_MARTS.COMPREHENSIVE_FT", 
                             "where LANDING_YEAR > 2004 and LANDING_YEAR<2014 and PACFIN_SPECIES_CODE='DCRB'"))

vessels <- sqlQuery(channel,paste("select distinct VESSEL_ID", 
                                  "from PACFIN_MARTS.COMPREHENSIVE_FT", 
                                  "where LANDING_YEAR > 2010 and FLEET_CODE='LE' and PACFIN_PORT_CODE in ('ERK','BRG','SF','PRN','MOS','MNT','MRO','SB','OXN','VEN')"))

close(channel)

names(df) <- c('vessel','year','month','lbs')
names(vessels) <- c('vessel')

crab <- tbl_df(df) %>% inner_join(vessels,by='vessel') %>% group_by(year,month) %>%
        summarise(lbs=sum(lbs)) %>% group_by(year) %>% mutate(yrtotal=sum(lbs),share=lbs/yrtotal,season=ifelse(month%in%c(1,2,3),'winter',
                                                        ifelse(month%in%c(4,5,6),'spring',
                                                        ifelse(month%in%c(7,8,9),'summer',
                                                        ifelse(month%in%c(10,11,12),'fall',NA))))) 

        
ggplot(crab,aes(x=year,y=share,color=factor(month))) + geom_line() + facet_wrap(~season)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#species seasonality...look at monthly distributions of sablefish,
# dover, thornyheads, petrale, and chilipepper

channel<- odbcConnect(dsn="pacfin",uid=uid,pw=pw,believeNRows=FALSE)

sabl <- sqlQuery(channel,paste("select LANDING_YEAR, LANDING_MONTH, PACFIN_SPECIES_CODE, sum(LANDED_WEIGHT_LBS)", 
                             "from PACFIN_MARTS.COMPREHENSIVE_FT", 
                             "where LANDING_YEAR > 2007 and LANDING_YEAR<2014 and PACFIN_SPECIES_CODE in ('PTRL','DOVR','SABL','CLPR') and",
                             "PACFIN_PORT_CODE in ('ERK','BRG','SF','PRN','MNT','MOS','MRO','SB','OXN','VEN') and FLEET_CODE='LE'",
                             "group by LANDING_YEAR, LANDING_MONTH, PACFIN_SPECIES_CODE"))
close(channel)

names(sabl) <- c('year','month','sps','lbs')

sabl <- tbl_df(sabl) %>% mutate(era=ifelse(year<2011,'pre','post')) %>%
        group_by(era,month,sps) %>% summarise(lbs=sum(lbs)) %>% 
        group_by(era,sps) %>% mutate(total=sum(lbs),share=lbs/total)
        

ggplot(subset(sabl,sps=='SABL'),aes(x=month,y=share,fill=factor(era))) + geom_bar(stat='identity',position='dodge') + 
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12))

#----------------------------------------------------------------------

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
#Did the use of fixed gear in the IFQ fishery increase? And is this change a possible 
# explanation for changing seasonal patterns

#by definition the use of fixed gear to harvest IFQ groundfish increase because what is now the 
# IFQ fishery was exclusively the trawl fishery before 2011...therefore any use of fixed gear 
# to harvest IFQ groundfish from 2011-2014 is an increase....the question is whether this is a
# significant increase and whether the popular times for fixed gear usage correspond to
# times that may have gained effort share after catch shares.


#---------------------------------------------------------------------------------------------------------
#first summarize the use of fixed gear vs. trawl gear from 2011-2014
channel<- odbcConnect(dsn="pacfin",uid=uid,pw=pw,believeNRows=FALSE)
#odbcGetInfo(channel)

df <- sqlQuery(channel,paste("select FTID, LANDING_YEAR, LANDING_MONTH, IS_IFQ_LANDING, PACFIN_PORT_CODE, PACFIN_GROUP_GEAR_CODE, DAHL_GROUNDFISH_CODE, FLEET_CODE", 
                             "from PACFIN_MARTS.COMPREHENSIVE_FT", 
                             "where LANDING_YEAR > 2007 and LANDING_YEAR<2015 and MANAGEMENT_GROUP_CODE='GRND'"))


close(channel)

names(df) <- c('FTID','year','month','IFQ','port','gear','dahl','fleet')

gear.df <- tbl_df(df) %>% filter(IFQ=='TRUE' & year %in% c(2011,2012,2013)) %>% group_by(year,port,gear) %>% summarise(tix=n_distinct(FTID)) %>% 
  mutate(state=ifelse(port%in%c('BRG','ERK','MNT','MOS','MRO','PRN','SF','OXN','SB'),'CA','OTH')) %>%
  group_by(state,gear) %>% summarise(ntix=sum(tix))

ggplot(subset(gear.df,state=='CA' & year<2015),aes(x=year,y=ntix,fill=gear)) + geom_bar(stat='identity') 


#plot the seasonality in trawl gear only 
twl <- df %>% filter(fleet=='LE' & gear=='TWL' & port %in% c('ERK','BRG','SF','PRN','MOS','MNT','MRO','SB','OXN','VEN')) %>%
        mutate(era=ifelse(year<2011,'pre','post')) %>%
        group_by(era,month) %>% summarise(ntix=n_distinct(FTID))
ggplot(twl,aes(x=month,y=ntix,fill=era)) + geom_bar(stat='identity',position='dodge') + 
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12))


#now plot the monthly intensity of effort with trawl gear from 2008-2010 & 2011-2013 and compare
# to the monthly intensity of pot and HKL gear in the groundfish fishery from 2008-2010
months.df <- tbl_df(df) %>% filter(dahl %in% c('04','07','08','09','10','20') & year < 2015) %>% #filter out the nearshore and weird efp stuff
              mutate(keepvar=ifelse(year<2011|(year>2010&IFQ=='TRUE'),1,0))  %>%
              filter(keepvar==1)%>%
    
            mutate(era=ifelse(year<2011,'pre','post'))  %>%
             filter(port %in% c('ERK','BRG','SF','PRN','MRO','MNT','MOS','SB','OXN','VEN')) %>%   #the fixed gear stuff south of 40'10 really only matters in MRO and PRN
             group_by(gear,month,era) %>% summarise(ntix=n_distinct(FTID))

ggplot(subset(months.df,gear!='TLS'),aes(x=month,y=ntix,fill=era)) + geom_bar(stat='identity',position='dodge') + 
  facet_wrap(~gear) + scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12))

#------------------------------------------------------------------------------------------------------



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################




