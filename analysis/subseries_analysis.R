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
         lntrips1=lag(lntrips,1),
         lntrips2=lag(lntrips,2),
         lntrips12=lag(lntrips,12),
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
  mutate(fmonth=factor(month),date=as.Date(paste(year,"-",month,"-","01",sep=""),format="%Y-%m-%d")) %>%
  mutate(oct=ifelse(month==10,1,0)) %>%
  group_by(area,year) %>% mutate(yr.avg=mean(lntrips),mean.diff=lntrips-yr.avg,mean.diff1=lag(mean.diff,1),mean.diff2=lag(mean.diff,2))


df.bimonthly <- df.monthly %>% mutate(month2=ifelse(month %in% c(1,2),1,
                                                    ifelse(month %in% c(3,4),2,
                                                           ifelse(month %in% c(5,6),3,
                                                                  ifelse(month %in% c(7,8),4,
                                                                         ifelse(month %in% c(9,10),5,
                                                                                ifelse(month %in% c(11,12),6,NA))))))) %>%
  group_by(area,year,month2) %>%
  summarise(total_weight=sum(total_weight),
            total_val=sum(total_val),
            ntrips=sum(ntrips),
            nvessel=sum(nvessel),
            wt.yrs=max(wt.yrs),
            val.yrs=max(val.yrs),
            trips.yrs=max(trips.yrs)) %>% 
  mutate(tripshare=ntrips/trips.yrs,valshare=total_val/val.yrs,wtshare=total_weight/wt.yrs) %>%
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
         lntrips1=lag(lntrips,1),
         lntrips2=lag(lntrips,2),
         
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
  mutate(fmonth2=factor(month2),date=factor(paste(year,"-",month2,sep=""))) %>%
  mutate(oct=ifelse(month2==5,1,0))%>%
  group_by(area,year) %>% mutate(yr.avg=mean(lntrips),mean.diff=lntrips-yr.avg,mean.diff1=lag(mean.diff,1),mean.diff2=lag(mean.diff,2))


#coefficent of variation for trips
cov <- df.monthly %>% group_by(area,year) %>% summarise(annual.mean=mean(ntrips),annual.sd=sd(ntrips),
                                                        annual.logmean=mean(lntrips),annual.logsd=sd(lntrips)) %>%
  mutate(cov=annual.sd/annual.mean,logcov=annual.logsd/annual.logmean)

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

#Look at the individual monthly series for trips and test them for breaks:

ggplot(subset(df.monthly,area=='north'),aes(x=year,y=lntrips)) + geom_line() + geom_point() + 
  facet_wrap(~month)

#test individual subseries for breaks:
z <- df.monthly[df.monthly$area=='north' & df.monthly$month==12,]
f_statistics <- Fstats(ts(lntrips,frequency=1, start=c(1994, 1)) ~ 1 , data = z)
plot(f_statistics)
breaks <- breakpoints(ts(lntrips,frequency=1, start=c(1994, 1)) ~ 1, data = z)
summary(breaks)
break.dates <- breakdates(breaks)

#----------------------------------------------------------
#for the south, the sup-F test was pretty uninamos in its 
# suggestion of a break in 2010.  Let's test the coefficients
# to see if relative magnitudes change

reg.df <- df.monthly %>% filter(area=='south') %>% mutate(break2010=ifelse(year>2010,1,0))
summary(lm(lntrips~factor(month)+factor(month)*break2010,data=reg.df))



#----------------------------------------------------------


