library(strucchange)
library(dplyr)
library(forecast)
library(broom)
library(stargazer)
library(ggplot2)

#useful function for rounding data frame
round_df <- function(dataframe, digits) {
  nums <- vapply(dataframe, is.numeric, FUN.VALUE = logical(1))
  
  dataframe[,nums] <- round(dataframe[,nums], digits = digits)
  
  (dataframe)
}
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
#ORGANIZE DATA 



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
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
#run some basic Chow tests just to see whether significant breakpoints are detected:

#North-----------------------------

df.n <- df.monthly[df.monthly$area=='north',]

#function to compute a chow test for a break at a particular month and year
chow.fn <- function(a,m,yr){
  df <- df.monthly[df.monthly$area==a,]
  idx <- which(df$month==m & df$year==yr)
  df1 <- df[1:idx,]
  df2 <- df[idx+1:nrow(df),]
  r.reg <- lm(tripshare~factor(month),data=df)
  ur.reg1 <- lm(tripshare~factor(month),data=df1)
  ur.reg2 <- lm(tripshare~factor(month),data=df2)
  SSR=NULL
  SSR$r=r.reg$residuals^2
  SSR$ur1=ur.reg1$residuals^2
  SSR$ur2=ur.reg2$residuals^2
  
  K = r.reg$rank
  
  Chow = ((sum(SSR$r) - (sum(SSR$ur1) + sum(SSR$ur2)))/K)/((sum(SSR$ur1)+sum(SSR$ur2))/(nrow(df) - 2*K))
  p <- 1-pf(Chow,K,(nrow(df)-2*K))  
return(data.frame(Fstat=Chow,p=p))
}


#----------------------------------


#South-----------------------------



#----------------------------------

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
#quick and dirty breakpoint analysis using OLS residuals to detrend the series

#south only to start
df <- df.monthly %>% filter(area=='south')

reg <- lm(lntrips~factor(year)-1,data=df)

df <- df %>% ungroup() %>% mutate(ols.resid=residuals(reg))
ggplot(df,aes(x=date,y=ols.resid)) + geom_line() + geom_point()



df <- df %>% filter(year>1994)

f_statistics <- Fstats(ts(ols.resid,frequency=12, start=c(1995, 1)) ~ factor(month) + lntrips12, data = df)
plot(f_statistics)
breaks <- breakpoints(ts(ols.resid,frequency=12, start=c(1995, 1)) ~ factor(month) + lntrips12, data = df)
break.dates <- breakdates(breaks)
summary(breaks)

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



#To avoid repetitive code, make a loop to get results for all log-differenced and share variables
#Levels will have their own loop for now, see below.
tomodel<-expand.grid(c("ntrips","lntrips","dlogwt","dlogval","dlogtrip","wtshare","valshare","tripshare","mean.diff"),c("north","south"))
names(tomodel)<-c("variable","area")
tomodel <- data.frame(lapply(tomodel, as.character), stringsAsFactors=FALSE)

#let's also keep track of the number of lag terms, break date predicted by fstats, whether it's significant, 
#number of breakdates predicted by B&P
lagterms<-numeric(nrow(tomodel))
fbreakdate<-numeric(nrow(tomodel))
fbreakp<-numeric(nrow(tomodel))
n_bp_breaks<-numeric(nrow(tomodel))

for(i in 1:nrow(tomodel))
#for(i in 1:2)
{#start of var/area loop

#because of the lagged variables we need to include only the years 1995-2014
  df.monthly <- df.monthly %>% filter(year>1994)
  
  df.bimonthly <- df.bimonthly %>% filter(year>1994)
  
varofinterest<-tomodel$variable[i]
Area<-tomodel$area[i]

mydf<-df.monthly[,c("year","month","area",varofinterest,paste0(varofinterest,1),paste0(varofinterest,2),"fmonth","date")]
mydf<-filter(mydf,area==Area)
names(mydf)<-c("year","month","area","var","varlag1","varlag2","fmonth","date")

mydf2<-df.bimonthly[,c("year","month2","area",varofinterest,paste0(varofinterest,1),paste0(varofinterest,2),"fmonth2","date")]
mydf2<-filter(mydf2,area==Area)
names(mydf2)<-c("year","month2","area","var","varlag1","varlag2","fmonth2","date")

#models
lm_s1<-lm(var~fmonth,data=mydf)
lm_s2<-lm(var~fmonth+varlag1,data=mydf)
lm_s3<-lm(var~fmonth+varlag1+varlag2,data=mydf)

a<-anova(lm_s1,lm_s2,lm_s3)

#Identify best model according to ANOVA
best<-tail(which(a[[6]]<0.05),n=1)
if(length(best)==0){best<-1}

if (best==1){
  #Get F statistics, B&P breakpoints
  f_statistics <- Fstats(ts(var,frequency=12, start=c(1995, 1)) ~ factor(month) , data = mydf)
  breaks <- breakpoints(ts(var,frequency=12, start=c(1995, 1)) ~ factor(month), data = mydf)
  break.dates <- breakdates(breaks)

  #Repeat for the bimonthly model
  f_statistics2 <- Fstats(ts(var,frequency=6, start=c(1995, 1)) ~ factor(month2) , data = mydf2)
  breaks2 <- breakpoints(ts(var,frequency=6, start=c(1995, 1)) ~ factor(month2), data = mydf2)
  break.dates2 <- breakdates(breaks)
  
  
    
  #Make plots and save them
  pdf(paste("results/",paste0(Area,"_",varofinterest,".pdf"),sep=""),height=5,width=10)
  par(mfrow=c(1,3))
  plot(f_statistics,main=paste("F statistics",Area,varofinterest),alpha = 0.01)
  lines(breakpoints(f_statistics))
  
  plot(breaks,main=paste("BIC and RSS",Area,varofinterest))
  
  plot(ts(mydf$var,frequency=12, start=c(1995, 1)),type="l",main=paste(Area,varofinterest,"with Bai & Perron breakpoints"),ylab=varofinterest)
  lines(breaks)
  dev.off()
  
  supftest<-sctest(f_statistics, type="supF")
  
  #Create linear model with single breakpoint identified by f-statistics, if supF test is significant
  if(supftest$p.value<0.05){
    
    #Make linear model with break
    mydf$decimaldate<-mydf$year+(mydf$month-1)/12 #need decimal date to id break
    breakpt<-breakdates(breakpoints(f_statistics))
    
    mydf<-mydf%>%
      mutate(era=cut(decimaldate,c(-Inf,breakpt,Inf),labels=c("beforebreak","afterbreak"))) #create dummy for before/after
    
    themodel<-lm(var~fmonth*era,data=mydf)

  }else{
    #Make linear model without break
    themodel<-lm(var~fmonth,data=mydf)
  }

}

if (best==2){
  #Get F statistics, B&P breakpoints
  f_statistics <- Fstats(ts(var,frequency=12, start=c(1995, 1)) ~ factor(month)+varlag1, data = mydf)
  breaks <- breakpoints(ts(var,frequency=12, start=c(1995, 1)) ~ factor(month)+varlag1, data = mydf)
  break.dates <- breakdates(breaks)
  
  #Make plots and save them
  pdf(paste("results/",paste0(Area,"_",varofinterest,".pdf"),sep=""),height=5,width=10)
  par(mfrow=c(1,3))
  plot(f_statistics,main=paste("F statistics",Area,varofinterest),alpha = 0.01)
  lines(breakpoints(f_statistics))
  
  plot(breaks,main=paste("BIC and RSS",Area,varofinterest))
  
  plot(ts(mydf$var,frequency=12, start=c(1995, 1)),type="l",main=paste(Area,varofinterest,"with Bai & Perron breakpoints"),ylab=varofinterest)
  lines(breaks)
  dev.off()
  
  supftest<-sctest(f_statistics, type="supF")
  
  #Create linear model with single breakpoint identified by f-statistics, if supF test is significant
  if(supftest$p.value<0.05){
    
    #Make linear model with break
    mydf$decimaldate<-mydf$year+(mydf$month-1)/12 #need decimal date to id break
    breakpt<-breakdates(breakpoints(f_statistics))
    
    mydf<-mydf%>%
      mutate(era=cut(decimaldate,c(-Inf,breakpt,Inf),labels=c("beforebreak","afterbreak"))) #create dummy for before/after
    
    themodel<-lm(var~fmonth*era+varlag1*era,data=mydf)
    
  }else{
    #Make linear model without break
    themodel<-lm(var~fmonth+varlag1,data=mydf)
  }
}

if (best==3){
  #Get F statistics, B&P breakpoints
  f_statistics <- Fstats(ts(var,frequency=12, start=c(1995, 1)) ~ factor(month)+varlag1+varlag2, data = mydf)
  breaks <- breakpoints(ts(var,frequency=12, start=c(1995, 1)) ~ factor(month)++varlag1+varlag2, data = mydf)
  break.dates <- breakdates(breaks)
  
  #Make plots and save them
#  pdf(paste0(Area,"_",varofinterest,".pdf"),height=5,width=10)

pdf(paste("results/",paste0(Area,"_",varofinterest,".pdf"),sep=""),height=5,width=10)
  par(mfrow=c(1,3))
  plot(f_statistics,main=paste("F statistics",Area,varofinterest),alpha = 0.01)
  lines(breakpoints(f_statistics))
  
  plot(breaks,main=paste("BIC and RSS",Area,varofinterest))
  
  plot(ts(mydf$var,frequency=12, start=c(1995, 1)),type="l",main=paste(Area,varofinterest,"with Bai & Perron breakpoints"),ylab=varofinterest)
  lines(breaks)
  dev.off()
  
  supftest<-sctest(f_statistics, type="supF")
  
  #Create linear model with single breakpoint identified by f-statistics, if supF test is significant
  if(supftest$p.value<0.05){
    
    #Make linear model with break
    mydf$decimaldate<-mydf$year+(mydf$month-1)/12 #need decimal date to id break
    breakpt<-breakdates(breakpoints(f_statistics))
    
    mydf<-mydf%>%
      mutate(era=cut(decimaldate,c(-Inf,breakpt,Inf),labels=c("beforebreak","afterbreak"))) #create dummy for before/after
    
    themodel<-lm(var~fmonth*era+varlag1*era+varlag2*era,data=mydf)
    
  }else{
    #Make linear model without break
    themodel<-lm(var~fmonth+varlag1+varlag2,data=mydf)
  }
  
}

#Make text file of results. Note that tex tables made by stargazer are too long (not compatible with longtable?) so stick with text for now
text_results<-stargazer(themodel, title=paste(Area,varofinterest,"linear model results"), type="text")
write(text_results,file=paste("results/",paste0(Area," ",varofinterest," ","linear model results",".txt"),sep=""))
#We can also make CSVs of the results, easier to paste into word tables
tidy_results<-round_df(tidy(themodel),3)
#write.csv(tidy_results,file=paste0(Area," ",varofinterest," ","linear model coefs",".csv")) #uncomment if you want this
glance_results<-round_df(glance(themodel),3)
#write.csv(glance_results,file=paste0(Area," ",varofinterest," ","linear model summary stats",".csv"))

lagterms[i]<-best-1 #how many lagged predictors are we including?
fbreakdate[i]<-breakpt
fbreakp[i]<-round(supftest$p.value,3)
n_bp_breaks[i]<-length(which(!is.na(breaks$breakpoints)))

}#end of var/area loop

modeled<-cbind(tomodel,lagterms,fbreakdate,fbreakp,n_bp_breaks) #summary of number of lag terms used, break date (according to f test), significance of supF test, number of Bai and Perron breaks

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
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

#Models for levels. These linear models include a time term since there might be a trend

tomodel2<-expand.grid(c("total_weight","total_val","ntrips"),c("north","south"))
names(tomodel2)<-c("variable","area")
tomodel2 <- data.frame(lapply(tomodel2, as.character), stringsAsFactors=FALSE)

#let's also keep track of the number of lag terms, break date predicted by fstats, whether it's significant, number of breakdates predicted by B&P
lagterms2<-numeric(nrow(tomodel2))
fbreakdate2<-numeric(nrow(tomodel2))
fbreakp2<-numeric(nrow(tomodel2))
n_bp_breaks2<-numeric(nrow(tomodel2))

for(i in 1:nrow(tomodel2))
  #for(i in 1:2)
{#start of var/area loop
  
  varofinterest<-tomodel2$variable[i]
  Area<-tomodel2$area[i]
  mydf<-df.monthly[,c("year","month","area",varofinterest,paste0(varofinterest,1),paste0(varofinterest,2),"fmonth")]
  mydf<-filter(mydf,area==Area)
  names(mydf)<-c("year","month","area","var","varlag1","varlag2","fmonth")
  mydf$decimaldate<-mydf$year+(mydf$month-1)/12
  
  #models
  lm_s1<-lm(var~fmonth+decimaldate,data=mydf)
  lm_s2<-lm(var~fmonth+decimaldate+varlag1,data=mydf)
  lm_s3<-lm(var~fmonth+decimaldate+varlag1+varlag2,data=mydf)
  
  a<-anova(lm_s1,lm_s2,lm_s3)
  
  #Identify best model according to ANOVA
  best<-tail(which(a[[6]]<0.05),n=1)
  if(length(best)==0){best<-1}
  
  if (best==1){
    #Get F statistics, B&P breakpoints
    f_statistics <- Fstats(ts(var,frequency=12, start=c(2004, 1)) ~ factor(month)+decimaldate, data = mydf)
    breaks <- breakpoints(ts(var,frequency=12, start=c(2004, 1)) ~ factor(month)+decimaldate, data = mydf)
    break.dates <- breakdates(breaks)
    
    #Make plots and save them
    pdf(paste("results/",paste0(Area,"_",varofinterest,".pdf"),sep=""),height=5,width=10)
    par(mfrow=c(1,3))
    plot(f_statistics,main=paste("F statistics",Area,varofinterest),alpha = 0.01)
    lines(breakpoints(f_statistics))
    
    plot(breaks,main=paste("BIC and RSS",Area,varofinterest))
    
    plot(ts(mydf$var,frequency=12, start=c(2004, 1)),type="l",main=paste(Area,varofinterest,"with Bai & Perron breakpoints"),ylab=varofinterest)
    lines(breaks)
    dev.off()
    
    supftest<-sctest(f_statistics, type="supF")
    
    #Create linear model with single breakpoint identified by f-statistics, if supF test is significant
    if(supftest$p.value<0.05){
      
      #Make linear model with break
      #mydf$decimaldate<-mydf$year+(mydf$month-1)/12 #need decimal date to id break
      breakpt<-breakdates(breakpoints(f_statistics))
      
      mydf<-mydf%>%
        mutate(era=cut(decimaldate,c(-Inf,breakpt,Inf),labels=c("beforebreak","afterbreak"))) #create dummy for before/after
      
      themodel<-lm(var~fmonth*era+decimaldate*era,data=mydf)
      
    }else{
      #Make linear model without break
      themodel<-lm(var~fmonth+decimaldate,data=mydf)
    }
    
  }
  
  if (best==2){
    #Get F statistics, B&P breakpoints
    f_statistics <- Fstats(ts(var,frequency=12, start=c(2004, 1)) ~ factor(month)+decimaldate+varlag1, data = mydf)
    breaks <- breakpoints(ts(var,frequency=12, start=c(2004, 1)) ~ factor(month)+decimaldate+varlag1, data = mydf)
    break.dates <- breakdates(breaks)
    
    #Make plots and save them
    pdf(paste("results/",paste0(Area,"_",varofinterest,".pdf"),sep=""),height=5,width=10)
    par(mfrow=c(1,3))
    plot(f_statistics,main=paste("F statistics",Area,varofinterest),alpha = 0.01)
    lines(breakpoints(f_statistics))
    
    plot(breaks,main=paste("BIC and RSS",Area,varofinterest))
    
    plot(ts(mydf$var,frequency=12, start=c(2004, 1)),type="l",main=paste(Area,varofinterest,"with Bai & Perron breakpoints"),ylab=varofinterest)
    lines(breaks)
    dev.off()
    
    supftest<-sctest(f_statistics, type="supF")
    
    #Create linear model with single breakpoint identified by f-statistics, if supF test is significant
    if(supftest$p.value<0.05){
      
      #Make linear model with break
      #mydf$decimaldate<-mydf$year+(mydf$month-1)/12 #need decimal date to id break
      breakpt<-breakdates(breakpoints(f_statistics))
      
      mydf<-mydf%>%
        mutate(era=cut(decimaldate,c(-Inf,breakpt,Inf),labels=c("beforebreak","afterbreak"))) #create dummy for before/after
      
      themodel<-lm(var~fmonth*era+varlag1*era+decimaldate*era,data=mydf)
      
    }else{
      #Make linear model without break
      themodel<-lm(var~fmonth+varlag1+decimaldate,data=mydf)
    }
  }
  
  if (best==3){
    #Get F statistics, B&P breakpoints
    f_statistics <- Fstats(ts(var,frequency=12, start=c(2004, 1)) ~ factor(month)+decimaldate+varlag1+varlag2, data = mydf)
    breaks <- breakpoints(ts(var,frequency=12, start=c(2004, 1)) ~ factor(month)+decimaldate+varlag1+varlag2, data = mydf)
    break.dates <- breakdates(breaks)
    
    #Make plots and save them
    pdf(paste("results/",paste0(Area,"_",varofinterest,".pdf"),sep=""),height=5,width=10)
    par(mfrow=c(1,3))
    plot(f_statistics,main=paste("F statistics",Area,varofinterest),alpha = 0.01)
    lines(breakpoints(f_statistics))
    
    plot(breaks,main=paste("BIC and RSS",Area,varofinterest))
    
    plot(ts(mydf$var,frequency=12, start=c(2004, 1)),type="l",main=paste(Area,varofinterest,"with Bai & Perron breakpoints"),ylab=varofinterest)
    lines(breaks)
    dev.off()
    
    supftest<-sctest(f_statistics, type="supF")
    
    #Create linear model with single breakpoint identified by f-statistics, if supF test is significant
    if(supftest$p.value<0.05){
      
      #Make linear model with break
      #mydf$decimaldate<-mydf$year+(mydf$month-1)/12 #need decimal date to id break
      breakpt<-breakdates(breakpoints(f_statistics))
      
      mydf<-mydf%>%
        mutate(era=cut(decimaldate,c(-Inf,breakpt,Inf),labels=c("beforebreak","afterbreak"))) #create dummy for before/after
      
      themodel<-lm(var~fmonth*era+varlag1*era+varlag2*era+decimaldate*era,data=mydf)
      
    }else{
      #Make linear model without break
      themodel<-lm(var~fmonth+varlag1+varlag2+decimaldate,data=mydf)
    }
    
  }
  
  #Make text file of results. Note that tex tables made by stargazer are too long (not compatible with longtable?) so stick with text for now
  text_results<-stargazer(themodel, title=paste(Area,gsub("_", ".", varofinterest),"linear model results"), type="text")#apparently can't have underscore in title?
  write(text_results,file=paste("results/",paste0(Area," ",varofinterest," ","linear model results",".txt"),sep=""))
  #We can also make CSVs of the results, easier to paste into word tables
  tidy_results<-round_df(tidy(themodel),3)
  #write.csv(tidy_results,file=paste0(Area," ",varofinterest," ","linear model coefs",".csv")) #uncomment if you want this
  glance_results<-round_df(glance(themodel),3)
  #write.csv(glance_results,file=paste0(Area," ",varofinterest," ","linear model summary stats",".csv"))
  
  lagterms2[i]<-best-1 #how many lagged predictors are we including?
  fbreakdate2[i]<-breakpt
  fbreakp2[i]<-round(supftest$p.value,3)
  n_bp_breaks2[i]<-length(which(!is.na(breaks$breakpoints)))
  
}#end of var/area loop

modeled2<-cbind(tomodel2,lagterms2,fbreakdate2,fbreakp2,n_bp_breaks2) #summary of number of lag terms used, break date (according to f test), number of Bai and Perron breaks


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
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


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


# 
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# #Here we implement the same procedure separately on each of the linear models
# # specified above.  The procedure works as follows:
# 
# # 1. run the F-test on all potential breakpoints and plot
# # 2. create a single plot that has the F-test results along with actual data
# #     for the north and south models and include actual data on the plot
# # 3. run auxillary F-tests (expF, aveF, supF, etc.)
# 
# #============================================================================
# #Step 1: get F-stats 
# 
# #F-test for south region dlog trip model
# fs_dlogtrips_s <- Fstats(ts(dlogtrip,frequency=12, start=c(2003, 1)) ~ factor(month), data = bymonth_s)
# plot(fs_dlogtrips_s)
# lines(breakpoints(fs_dlogtrips_s))
# 
# #F-test for north region dlog trip model
# fs_dlogtrips_n<-Fstats(ts(dlogtrip,frequency=12, start=c(2003, 1)) ~ factor(month), data = bymonth_n)
# plot(fs_dlogtrips_n)
# lines(breakpoints(fs_dlogtrips_n))
# 
# #=======================================================================
# #Step 2: Bai and Perron procedure
# breaks_dlogtrips_s <- breakpoints(ts(dlogtrip,frequency=12, start=c(2003, 1)) ~ factor(month), data = bymonth_s)
# breakdates_dlogtrips_s <- breakdates(breaks_dlogtrips_s)
# #plot(ts(bymonth_s$dlogtrip,frequency=12, start=c(2003, 1)),type="l",main="North log-differenced trips")
# #lines(breaks_dlogtrips_s)
# 
# 
# breaks_dlogtrips_n <- breakpoints(ts(dlogtrip,frequency=12, start=c(2003, 1)) ~ factor(month), data = bymonth_n)
# breakdates_dlogtrips_n <- breakdates(breaks_dlogtrips_n)
# #===================================================================================
# #Step 3: create a combined plot
# par(mfrow=c(3,2))
# 
# plot(fs_dlogtrips_n,main="supF test, north log-differenced trips",alpha = 0.01)
# lines(breakpoints(fs_dlogtrips_n))
# plot(fs_dlogtrips_s,main="supF test, south log-differenced trips",alpha = 0.01)
# lines(breakpoints(fs_dlogtrips_s))
# 
# plot(breaks_dlogtrips_n,main="BIC and RSS, north log-difference trips")
# plot(breaks_dlogtrips_s,main="BIC and RSS, south log-difference trips")
# 
# plot(ts(bymonth_n$dlogtrip,frequency=12, start=c(2003, 1)),type="l",main="North log-differenced trips")
# lines(breaks_dlogtrips_n)
# plot(ts(bymonth_s$dlogtrip,frequency=12, start=c(2003, 1)),type="l",main="South log-differenced trips")
# lines(breaks_dlogtrips_s)
# #====================================================================
# 
# #====================================================================
# #Step 4: 
# #package strucchange performs expF, supF, and aveF test...run these
# # on both the north and south models
# sctest(fs_dlogtrips_n, type="expF")
# sctest(fs_dlogtrips_n, type="supF")
# sctest(fs_dlogtrips_n, type="aveF")
# sctest(fs_dlogtrips_s, type="expF")
# sctest(fs_dlogtrips_s, type="supF")
# sctest(fs_dlogtrips_s, type="aveF")
# #all significant, but does that matter? supF is the only one informative about location anyway
# #=======================================================================
# 
# #====================================================================



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





#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#POSSIBLE THROW-AWAY CODE





# #####################log-differenced time series#####################
# ######################################################Number of trips
# ##################south
# 
# ##########################################################weight
# ####south
# lm_wt_s1<-lm(dlogwt~month,data=bymonth_s)
# lm_wt_s2<-lm(dlogwt~month+dlogwt1,data=bymonth_s)
# lm_wt_s3<-lm(dlogwt~month+dlogwt1+dlogwt2,data=bymonth_s)
# anova(lm_wt_s1,lm_wt_s2,lm_wt_s3)
# #use third model
# 
# fs_dlogwt_s<-Fstats(ts(dlogwt,frequency=12, start=c(2003, 1)) ~ month+dlogwt1+dlogwt2, data = bymonth_s)
# plot(fs_dlogwt_s)
# lines(breakpoints(fs_dlogwt_s))
# 
# ###north
# lm_wt_n1<-lm(dlogwt~month,data=bymonth_n)
# lm_wt_n2<-lm(dlogwt~month+dlogwt1,data=bymonth_n)
# lm_wt_n3<-lm(dlogwt~month+dlogwt1+dlogwt2,data=bymonth_n)
# anova(lm_wt_n1,lm_wt_n2,lm_wt_n3)
# #use third model
# 
# fs_dlogwt_n<-Fstats(ts(dlogwt,frequency=12, start=c(2003, 1)) ~ month+dlogwt1+dlogwt2, data = bymonth_n)
# plot(fs_dlogwt_n)
# lines(breakpoints(fs_dlogwt_n))
# 
# #plot all
# par(mfrow=c(2,2))
# plot(fs_dlogwt_n,main="supF test, north log-differenced weight",alpha = 0.01)
# lines(breakpoints(fs_dlogwt_n))
# plot(fs_dlogwt_s,main="supF test, south log-differenced weight",alpha = 0.01)
# lines(breakpoints(fs_dlogwt_s))
# plot(bymonth_n$yrmonth,bymonth_n$dlogwt,type="l",main="North log-differenced weight")
# plot(bymonth_s$yrmonth,bymonth_s$dlogwt,type="l",main="South log-differenced weight")
# 
# sctest(fs_dlogwt_n, type="expF")
# sctest(fs_dlogwt_n, type="supF")
# sctest(fs_dlogwt_n, type="aveF")
# sctest(fs_dlogwt_s, type="expF")
# sctest(fs_dlogwt_s, type="supF")
# sctest(fs_dlogwt_s, type="aveF")
# 
# ##########################################################################value
# ###south
# lm_val_s1<-lm(dlogval~month,data=bymonth_s)
# lm_val_s2<-lm(dlogval~month+dlogval1,data=bymonth_s)
# lm_val_s3<-lm(dlogval~month+dlogval1+dlogval2,data=bymonth_s)
# anova(lm_val_s1,lm_val_s2,lm_val_s3)
# #use second model
# 
# fs_dlogval_s<-Fstats(ts(dlogval,frequency=12, start=c(2003, 1)) ~ month+dlogval1, data = bymonth_s)
# plot(fs_dlogval_s)
# lines(breakpoints(fs_dlogval_s))
# 
# ####north
# lm_val_n1<-lm(dlogval~month,data=bymonth_n)
# lm_val_n2<-lm(dlogval~month+dlogval1,data=bymonth_n)
# lm_val_n3<-lm(dlogval~month+dlogval1+dlogval2,data=bymonth_n)
# anova(lm_val_n1,lm_val_n2,lm_val_n3)
# #use second model
# 
# fs_dlogval_n<-Fstats(ts(dlogval,frequency=12, start=c(2003, 1)) ~ month+dlogval1, data = bymonth_n)
# plot(fs_dlogval_n)
# lines(breakpoints(fs_dlogval_n))
# 
# #plot
# par(mfrow=c(2,2))
# plot(fs_dlogval_n,main="supF test, north log-differenced value",alpha = 0.01)
# lines(breakpoints(fs_dlogval_n))
# plot(fs_dlogval_s,main="supF test, south log-differenced value",alpha = 0.01)
# lines(breakpoints(fs_dlogval_s))
# plot(bymonth_n$yrmonth,bymonth_n$dlogval,type="l",main="North log-differenced value")
# plot(bymonth_s$yrmonth,bymonth_s$dlogval,type="l",main="South log-differenced value")
# 
# sctest(fs_dlogval_n, type="expF")
# sctest(fs_dlogval_n, type="supF")
# sctest(fs_dlogval_n, type="aveF")
# sctest(fs_dlogval_s, type="expF")
# sctest(fs_dlogval_s, type="supF")
# sctest(fs_dlogval_s, type="aveF")
# 
# ##################################log time series################################
# #add time column
# bymonth_n$t<-c(1:dim(bymonth_n)[1])
# bymonth_s$t<-c(1:dim(bymonth_s)[1])
# 
# #######################################################Number of trips
# #south
# lm_logtrips_s1<-lm(logtrip~month+t,data=bymonth_s)
# lm_logtrips_s2<-lm(logtrip~month+logtrip1+t,data=bymonth_s)
# lm_logtrips_s3<-lm(logtrip~month+logtrip1+logtrip2+t,data=bymonth_s)
# anova(lm_logtrips_s1,lm_logtrips_s2,lm_logtrips_s3)
# #use second model
# 
# fs_logtrips_s <- Fstats(ts(logtrip,frequency=12, start=c(2003, 1)) ~ month+logtrip1+t, data = bymonth_s)
# plot(fs_logtrips_s,alpha=0.01)
# lines(breakpoints(fs_logtrips_s))
# 
# #north
# lm_logtrips_n1<-lm(logtrip~month+t,data=bymonth_n)
# lm_logtrips_n2<-lm(logtrip~month+logtrip1+t,data=bymonth_n)
# lm_logtrips_n3<-lm(logtrip~month+logtrip1+logtrip2+t,data=bymonth_n)
# anova(lm_logtrips_n1,lm_logtrips_n2,lm_logtrips_n3)
# #use second model
# 
# fs_logtrips_n <- Fstats(ts(logtrip,frequency=12, start=c(2003, 1)) ~ month+logtrip1+t, data = bymonth_n)
# plot(fs_logtrips_n,alpha=0.01)
# lines(breakpoints(fs_logtrips_n))
# 
# #plot
# par(mfrow=c(2,2))
# plot(fs_logtrips_n,main="supF test, north log trips",alpha = 0.05)
# lines(breakpoints(fs_logtrips_n))
# plot(fs_logtrips_s,main="supF test, south log trips",alpha = 0.05)
# lines(breakpoints(fs_logtrips_s))
# plot(bymonth_n$yrmonth,bymonth_n$logtrip,type="l",main="North log trips")
# plot(bymonth_s$yrmonth,bymonth_s$logtrip,type="l",main="South log trips")
# 
# sctest(fs_logtrips_n, type="expF")
# sctest(fs_logtrips_n, type="supF")
# sctest(fs_logtrips_n, type="aveF")
# sctest(fs_logtrips_s, type="expF")
# sctest(fs_logtrips_s, type="supF")
# sctest(fs_logtrips_s, type="aveF")
# 
# #######################################################Value
# #south
# lm_logval_s1<-lm(logval~month+t,data=bymonth_s)
# lm_logval_s2<-lm(logval~month+logval1+t,data=bymonth_s)
# lm_logval_s3<-lm(logval~month+logval1+logval2+t,data=bymonth_s)
# anova(lm_logval_s1,lm_logval_s2,lm_logval_s3)
# #use second model
# 
# fs_logval_s <- Fstats(ts(logval,frequency=12, start=c(2003, 1)) ~ month+logval1+t, data = bymonth_s)
# plot(fs_logval_s,alpha=0.01)
# lines(breakpoints(fs_logval_s))
# 
# #north
# lm_logval_n1<-lm(logval~month+t,data=bymonth_n)
# lm_logval_n2<-lm(logval~month+logval1+t,data=bymonth_n)
# lm_logval_n3<-lm(logval~month+logval1+logval2+t,data=bymonth_n)
# anova(lm_logval_n1,lm_logval_n2,lm_logval_n3)
# #use third model
# 
# fs_logval_n <- Fstats(ts(logval,frequency=12, start=c(2003, 1)) ~ month+logval1+logval2+t, data = bymonth_n)
# plot(fs_logval_n,alpha=0.01)
# lines(breakpoints(fs_logval_n))
# 
# par(mfrow=c(2,2))
# plot(fs_logval_n,main="supF test, north log value",alpha = 0.01)
# lines(breakpoints(fs_logval_n))
# plot(fs_logval_s,main="supF test, south log value",alpha = 0.01)
# lines(breakpoints(fs_logval_s))
# plot(bymonth_n$yrmonth,bymonth_n$logval,type="l",main="North log value")
# plot(bymonth_s$yrmonth,bymonth_s$logval,type="l",main="South log value")
# 
# sctest(fs_logval_n, type="expF")
# sctest(fs_logval_n, type="supF")
# sctest(fs_logval_n, type="aveF")
# sctest(fs_logval_s, type="expF")
# sctest(fs_logval_s, type="supF")
# sctest(fs_logtrips_s, type="aveF")
# 
# 
# 
# #######################################################Weight
# #south
# lm_logwt_s1<-lm(logwt~month+t,data=bymonth_s)
# lm_logwt_s2<-lm(logwt~month+logwt1+t,data=bymonth_s)
# lm_logwt_s3<-lm(logwt~month+logwt1+logwt2+t,data=bymonth_s)
# anova(lm_logwt_s1,lm_logwt_s2,lm_logwt_s3)
# #use second model
# 
# fs_logwt_s <- Fstats(ts(logwt,frequency=12, start=c(2003, 1)) ~ month+logwt1+t, data = bymonth_s)
# plot(fs_logwt_s,alpha=0.01)
# lines(breakpoints(fs_logwt_s))
# 
# #north
# lm_logwt_n1<-lm(logwt~month+t,data=bymonth_n)
# lm_logwt_n2<-lm(logwt~month+logwt1+t,data=bymonth_n)
# lm_logwt_n3<-lm(logwt~month+logwt1+logwt2+t,data=bymonth_n)
# anova(lm_logwt_n1,lm_logwt_n2,lm_logwt_n3)
# #use third model
# 
# fs_logwt_n <- Fstats(ts(logwt,frequency=12, start=c(2003, 1)) ~ month+logwt1+logwt2+t, data = bymonth_n)
# plot(fs_logwt_n,alpha=0.01)
# lines(breakpoints(fs_logwt_n))
# 
# #plot
# par(mfrow=c(2,2))
# plot(fs_logwt_n,main="supF test, north log weight",alpha = 0.05)
# lines(breakpoints(fs_logwt_n))
# plot(fs_logwt_s,main="supF test, south log weight",alpha = 0.05)
# lines(breakpoints(fs_logwt_s))
# plot(bymonth_n$yrmonth,bymonth_n$logwt,type="l",main="North log weight")
# plot(bymonth_s$yrmonth,bymonth_s$logwt,type="l",main="South log weight")
# 
# sctest(fs_logval_n, type="expF")
# sctest(fs_logval_n, type="supF")
# sctest(fs_logval_n, type="aveF")
# sctest(fs_logval_s, type="expF")
# sctest(fs_logval_s, type="supF")
# sctest(fs_logval_s, type="aveF")
# 
# ##################################untransformed time series################################
# 
# 
# #######################################################Number of trips
# #south
# lm_trips_s1<-lm(ntrips~month+t,data=bymonth_s)
# lm_trips_s2<-lm(ntrips~month+ntrips1+t,data=bymonth_s)
# lm_trips_s3<-lm(ntrips~month+ntrips1+ntrips2+t,data=bymonth_s)
# anova(lm_trips_s1,lm_trips_s2,lm_trips_s3)
# #use third model
# 
# fs_trips_s <- Fstats(ts(ntrips,frequency=12, start=c(2003, 1)) ~ month+ntrips1+ntrips2+t, data = bymonth_s)
# plot(fs_trips_s,alpha=0.01)
# lines(breakpoints(fs_trips_s))
# 
# #north
# lm_trips_n1<-lm(ntrips~month+t,data=bymonth_n)
# lm_trips_n2<-lm(ntrips~month+ntrips1+t,data=bymonth_n)
# lm_trips_n3<-lm(ntrips~month+ntrips1+ntrips2+t,data=bymonth_n)
# anova(lm_trips_n1,lm_trips_n2,lm_trips_n3)
# #use second model
# 
# fs_trips_n <- Fstats(ts(ntrips,frequency=12, start=c(2003, 1)) ~ month+ntrips1+t, data = bymonth_n)
# plot(fs_trips_n,alpha=0.01)
# lines(breakpoints(fs_trips_n))
# 
# #plot
# par(mfrow=c(2,2))
# plot(fs_trips_n,main="supF test, north trips",alpha = 0.05)
# lines(breakpoints(fs_trips_n))
# plot(fs_trips_s,main="supF test, south trips",alpha = 0.05)
# lines(breakpoints(fs_trips_s))
# plot(bymonth_n$yrmonth,bymonth_n$ntrips,type="l",main="North trips")
# plot(bymonth_s$yrmonth,bymonth_s$ntrips,type="l",main="South trips")
# 
# sctest(fs_logtrips_n, type="expF")
# sctest(fs_logtrips_n, type="supF")
# sctest(fs_logtrips_n, type="aveF")
# sctest(fs_logtrips_s, type="expF")
# sctest(fs_logtrips_s, type="supF")
# sctest(fs_logtrips_s, type="aveF")
# 
# #######################################################Value
# #south
# lm_val_s1<-lm(total_val~month+t,data=bymonth_s)
# lm_val_s2<-lm(total_val~month+total_val1+t,data=bymonth_s)
# lm_val_s3<-lm(total_val~month+total_val1+total_val2+t,data=bymonth_s)
# anova(lm_val_s1,lm_val_s2,lm_val_s3)
# #use second model
# 
# fs_val_s <- Fstats(ts(total_val,frequency=12, start=c(2003, 1)) ~ month+total_val1+t, data = bymonth_s)
# plot(fs_val_s,alpha=0.05)
# lines(breakpoints(fs_val_s))
# 
# #north
# lm_val_n1<-lm(total_val~month+t,data=bymonth_n)
# lm_val_n2<-lm(total_val~month+total_val1+t,data=bymonth_n)
# lm_val_n3<-lm(total_val~month+total_val1+logval2+t,data=bymonth_n)
# anova(lm_val_n1,lm_val_n2,lm_val_n3)
# #use second model
# 
# fs_val_n <- Fstats(ts(total_val,frequency=12, start=c(2003, 1)) ~ month+total_val1+t, data = bymonth_n)
# plot(fs_val_n,alpha=0.01)
# lines(breakpoints(fs_val_n))
# 
# par(mfrow=c(2,2))
# plot(fs_val_n,main="supF test, north value",alpha = 0.05)
# lines(breakpoints(fs_val_n))
# plot(fs_val_s,main="supF test, south value",alpha = 0.05)
# lines(breakpoints(fs_val_s))
# plot(bymonth_n$yrmonth,bymonth_n$val,type="l",main="North value")
# plot(bymonth_s$yrmonth,bymonth_s$val,type="l",main="South value")
# 
# sctest(fs_logval_n, type="expF")
# sctest(fs_logval_n, type="supF")
# sctest(fs_logval_n, type="aveF")
# sctest(fs_logval_s, type="expF")
# sctest(fs_logval_s, type="supF")
# sctest(fs_logtrips_s, type="aveF")
# 
# 
# 
# #######################################################Weight
# #south
# lm_wt_s1<-lm(total_weight~month+t,data=bymonth_s)
# lm_wt_s2<-lm(total_weight~month+total_weight1+t,data=bymonth_s)
# lm_wt_s3<-lm(total_weight~month+total_weight1+total_weight2+t,data=bymonth_s)
# anova(lm_wt_s1,lm_wt_s2,lm_wt_s3)
# #use third? model
# 
# fs_wt_s <- Fstats(ts(total_weight,frequency=12, start=c(2003, 1)) ~ month+total_weight1+total_weight2+t, data = bymonth_s)
# plot(fs_wt_s,alpha=0.05)
# lines(breakpoints(fs_wt_s))
# 
# #north
# lm_wt_n1<-lm(total_weight~month+t,data=bymonth_n)
# lm_wt_n2<-lm(total_weight~month+total_weight1+t,data=bymonth_n)
# lm_wt_n3<-lm(total_weight~month+total_weight1+total_weight2+t,data=bymonth_n)
# anova(lm_wt_n1,lm_wt_n2,lm_wt_n3)
# #use second model
# 
# fs_wt_n <- Fstats(ts(total_weight,frequency=12, start=c(2003, 1)) ~ month+total_weight1+t, data = bymonth_n)
# plot(fs_wt_n,alpha=0.01)
# lines(breakpoints(fs_wt_n))
# 
# #plot
# par(mfrow=c(2,2))
# plot(fs_wt_n,main="supF test, north weight",alpha = 0.05)
# lines(breakpoints(fs_wt_n))
# plot(fs_wt_s,main="supF test, south weight",alpha = 0.05)
# lines(breakpoints(fs_wt_s))
# plot(bymonth_n$yrmonth,bymonth_n$total_weight,type="l",main="North weight")
# plot(bymonth_s$yrmonth,bymonth_s$total_weight,type="l",main="South weight")
# 
# sctest(fs_logval_n, type="expF")
# sctest(fs_logval_n, type="supF")
# sctest(fs_logval_n, type="aveF")
# sctest(fs_logval_s, type="expF")
# sctest(fs_logval_s, type="supF")
# sctest(fs_logval_s, type="aveF")
# 
# 
# 
# #######################################also try changepoint library to id change in mean/variance?
# #doesn't really work
# library(changepoint)
# data("discoveries", package = "datasets")
# dis.pelt <- cpt.meanvar(discoveries, test.stat = "Poisson",method = "PELT")
# plot(dis.pelt, cpt.width = 3)
# 
# strips.pelt <- cpt.meanvar(bymonth_s$ntrips, test.stat = "Poisson",method = "PELT")
# plot(strips.pelt, cpt.width = 3)
# #Might work better on seasonally-adjusted data...
# 
# 
# ########################################What about looking at change in ANNUAL mean and variance?
# #find breakpoints, following nile example from strucchange? or ar1 model?
# 
# #function to create lagged colum
# lagpad <- function(x, k) {
#   c(rep(NA, k), x)[1 : length(x)] 
# }
# 
# byyear<-bymonth%>%
#   group_by(year,area)%>%
#   summarise(total_wt=sum(total_weight),var_weight=var(total_weight),total_vl=sum(total_val),var_val=var(total_val),total_trips=sum(ntrips),var_trips=var(ntrips),total_vessels=sum(nvessels),var_vessels=var(nvessels))%>%
#   group_by(area)%>%
#   mutate(total_wt1=lagpad(total_wt,1),total_vl1=lagpad(total_vl,1),total_trips1=lagpad(total_trips,1),total_vessels1=lagpad(total_vessels,1))
# byyear<-filter(byyear,year>2002)
# #first lagged values are wrong...fix later but for now just make NAs
# byyear$total_wt1[1:2]=NA
# byyear$total_vl1[1:2]=NA
# byyear$total_trips1[1:2]=NA
# byyear$total_vessels1[1:2]=NA
# #also might be interesting to look at mean trips, value, weight per vessel (or mean value, weight per trip)
# 
# 
# 
# byyear_s<-filter(byyear,area=="south")
# byyear_n<-filter(byyear,area=="north")
# 
# ####AR model not working right, try again later
# ###########################################number of trips
# fs_antrips_s <- Fstats(ts(total_trips,frequency=1, start=c(2003)) ~ 1, data = byyear_s)
# plot(fs_antrips_s,alpha=0.05)
# lines(breakpoints(fs_antrips_s))
# 
# fs_antrips_n <- Fstats(ts(total_trips,frequency=1, start=c(2003)) ~ 1, data = byyear_n)
# plot(fs_antrips_n,alpha=0.05)
# lines(breakpoints(fs_antrips_n))
# 
# #plot
# par(mfrow=c(2,2))
# plot(fs_antrips_n,main="supF test, north total annual trips",alpha = 0.05)
# lines(breakpoints(fs_antrips_n))
# plot(fs_antrips_s,main="supF test, south total annual trips",alpha = 0.05)
# lines(breakpoints(fs_antrips_s))
# plot(byyear_n$year,byyear_n$total_trips,type="l",main="North total annual trips")
# plot(byyear_s$year,byyear_s$total_trips,type="l",main="South total annual trips")
# 
# ###############################################value
# fs_anval_s <- Fstats(ts(total_vl,frequency=1, start=c(2003)) ~ 1, data = byyear_s)
# plot(fs_anval_s,alpha=0.05)
# lines(breakpoints(fs_anval_s))
# 
# fs_anval_n <- Fstats(ts(total_vl,frequency=1, start=c(2003)) ~ 1, data = byyear_n)
# plot(fs_anval_n,alpha=0.05)
# lines(breakpoints(fs_anval_n))
# 
# par(mfrow=c(2,2))
# plot(fs_anval_n,main="supF test, north total annual value",alpha = 0.05)
# lines(breakpoints(fs_anval_n))
# plot(fs_anval_s,main="supF test, south total annual value",alpha = 0.05)
# lines(breakpoints(fs_anval_s))
# plot(byyear_n$year,byyear_n$total_vl,type="l",main="North total annual value")
# plot(byyear_s$year,byyear_s$total_vl,type="l",main="South total annual value")
# 
# ###############################################weight
# fs_anwt_s <- Fstats(ts(total_wt,frequency=1, start=c(2003)) ~ 1, data = byyear_s)
# plot(fs_anwt_s,alpha=0.05)
# lines(breakpoints(fs_anwt_s))
# 
# fs_anwt_n <- Fstats(ts(total_wt,frequency=1, start=c(2003)) ~ 1, data = byyear_n)
# plot(fs_anwt_n,alpha=0.05)
# lines(breakpoints(fs_anwt_n))
# 
# par(mfrow=c(2,2))
# plot(fs_anwt_n,main="supF test, north total annual weight",alpha = 0.05)
# lines(breakpoints(fs_anwt_n))
# plot(fs_anwt_s,main="supF test, south total annual weight",alpha = 0.05)
# lines(breakpoints(fs_anwt_s))
# plot(byyear_n$year,byyear_n$total_wt,type="l",main="North total annual weight")
# plot(byyear_s$year,byyear_s$total_wt,type="l",main="South total annual weight")
# 
# ################
# #try breakout detection?
# library(BreakoutDetection)
# 
# #########################################Next: run lms with breakpoints IDed by supF test