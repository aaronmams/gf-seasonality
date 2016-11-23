
require(RODBC)
require(dplyr)
require(lubridate)

login<-read.csv('R:/Kate/code/login.csv')
uid <- as.character(login$pacfinuid)
pw <- as.character(login$pacfinpw)

#some test queries
channel<- odbcConnect(dsn="pacfin",uid=paste(uid),pw=paste(pw),believeNRows=FALSE,rows_at_time=1)
test <- sqlQuery(channel,paste("select YEAR,FTID,GRADE,AGID",
                             "from PACFIN.FTL",
                             "where YEAR > 2015"))
close(channel)



channel<- odbcConnect(dsn="pacfin",uid=paste(uid),pw=paste(pw),believeNRows=FALSE)
#odbcGetInfo(channel)
t <- Sys.time()
df <- sqlQuery(channel,paste("select LANDING_YEAR, LANDING_DATE, FTID, AGENCY_CODE, PARTICIPATION_GROUP_CODE",",",
                             "FLEET_CODE, VESSEL_ID, PACFIN_GROUP_GEAR_CODE, PACFIN_PORT_CODE, PACFIN_SPECIES_CODE",",",
                             "LANDED_WEIGHT_LBS, EXVESSEL_REVENUE, IS_IFQ_LANDING",
                             "from PACFIN_MARTS.COMPREHENSIVE_FT",
                             "where LANDING_YEAR > 1993",
                             "AND",
                             "THOMSON_FISHERY_CODE IN ('01','05','06','07','08','09','10','11','30')",
                             "AND",
                             "PACFIN_SPECIES_CODE<>'PWHT'", 
                             "AND",
                             "MANAGEMENT_GROUP_CODE='GRND'"))
Sys.time() - t
close(channel)


#df<-read.csv("groundfish.csv",header=TRUE,stringsAsFactors = FALSE)

#---------------------------------------------------------------------------------------------------
#Note from KR: I used SQLdeveloper to pull directly from pacfin since we don't have a local database. Here's the query I used:
#select * from PACFIN_MARTS.COMPREHENSIVE_FT where LANDING_YEAR > 2000 AND 
#THOMSON_FISHERY_CODE IN ('01','05','06','07','08','09','10','11','30') AND
#PACFIN_SPECIES_CODE<>'PWHT' AND
#MANAGEMENT_GROUP_CODE='GRND'

#The result from this query was then read into R as the object "df"
#----------------------------------------------------------------------------------------------------
#keep observations if they were from the limited entry fishery before 2011 & if they had IFQ landings
# after 2011

df<-tbl_df(df)%>%
  filter((LANDING_YEAR<2011& FLEET_CODE=="LE")|(LANDING_YEAR>2010&IS_IFQ_LANDING=="TRUE"))

#-----------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------
#Bin ports by north, south, vancouver
ports<-read.csv("data/lbk_port_codes.csv",header=TRUE) 
#note some ports don't have lat/lon, but in our data they 
# all appear to be south of the 40.10 line (ie OMD, OSF, OSM)

north <- ports$PCID[which(ports$LAT>40.1 & ports$INPFC_AREA!='VANCOUVER')]
vanports <- ports$PCID[which(ports$INPFC_AREA=="VANCOUVER")]

df <- df %>% mutate(area=ifelse(PACFIN_PORT_CODE %in% north,'north',
                                ifelse(PACFIN_PORT_CODE %in% vanports,'vancouver','south')))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
nomsps<-read.csv("data/nominal species names.csv",header=TRUE,stringsAsFactors=FALSE) #has species codes for nominal species (a few are NA)
nomsps <- nomsps[,c('Nominal.species.name','Species.name')]
#this gets a little funky, we are actually going to match the nominal name from our list with the 
# species name in the data frame.  using a left join we will keep all the observations and only
# match those that are recorded with nominal names
names(nomsps) <- c('PACFIN_SPECIES_CODE','species')

df <- df %>% left_join(nomsps,by="PACFIN_SPECIES_CODE") %>%
  mutate(SPS=ifelse(is.na(species),PACFIN_SPECIES_CODE,species))

#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#final step: change the column names to be consistent with what Kate has referenced in the rest of the 
#    code

#Just take the columns that might be useful at some point, somewhat different than the columns pulled from FTL_NEW
#df2<-df %>%
#  select(LANDING_YEAR,LANDING_DATE,FTID,AGENCY_CODE,PARTICIPATION_GROUP_CODE,FLEET_CODE,VESSEL_ID,
#         PACFIN_GROUP_GEAR_CODE,PACFIN_PORT_CODE,PACFIN_SPECIES_CODE,LANDED_WEIGHT_LBS,PRICE_PER_POUND, 
#         EXVESSEL_REVENUE,IS_IFQ_LANDING, area, SPS)

df2 <- df 
names(df2) <- c("year","date","ftid","agency","participation_grp",
                "fleet_code","vessel_id","group_gear","port_code",
                "species_code","weight", "revenue","ifq", "area", "species","SPS")


#Account for inflation--not sure if this is important, but might as well. This puts everything in 2005 dollars (arbitrary year).
deflate<-read.csv("data/deflators.csv",stringsAsFactors = FALSE)
names(deflate)<-c("year","deflator")
df2<-left_join(df2,deflate)%>%mutate(val=revenue/deflator) 
#-----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#filter out Vancouver obs, use obs from 2003-2015, and roll up to monthly observations

df.monthly <- df2 %>% filter(area!='vancouver' & year<= 2015) %>% mutate(month=month(date)) %>%
  group_by(year,month,area) %>%
  summarise(total_weight=sum(weight),total_val=sum(val),
            ntrips=n_distinct(ftid),nvessel=n_distinct(vessel_id)) %>%
  group_by(area,year) %>%
  mutate(wt.yrs=sum(total_weight),
         val.yrs=sum(total_val),
         trips.yrs=sum(ntrips),
         wtshare=total_weight/wt.yrs,
         valshare=total_val/val.yrs,
         tripshare=ntrips/trips.yrs)

#To calculate shares
df.yearly<-df2 %>% 
  filter(area!='vancouver' & year<= 2015) %>%
  group_by(year,area) %>%
  summarise(total_annual_weight=sum(weight),
            total_annual_val=sum(val),
            n_annual_trips=n_distinct(ftid))



#----------------------------------------------------------------------------------------------------------

saveRDS(df2,file="R:/Kate/Data/gf_main.RDA")
saveRDS(df.monthly,file="data/gf_monthly.RDA")
saveRDS(df.yearly,file="data/gf_yearly.RDA")
