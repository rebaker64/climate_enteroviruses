require("lubridate")
require("doBy")
require("pals")
library(circular)
require("tidyr")
library("dplyr")

# make center of mass and intensity for EVA71 and CVA16
setwd("~/Dropbox/Enterovirus/Scripts/forGithub")

# load and generate data variables
disease <- read.csv("china_entero.csv")
cid <- read.csv("china_id.csv")
did <- merge(disease, cid, by.x=c("PROV_ID"),by.y=c("ID"))
did$Prov <- toupper(did$Prov)
did$date <- as.Date(paste0(1,"/",did$WEEK,"/",did$YEAR), format = c("%w/%W/%Y"))
did$month <- month(did$date)
did <- did[!is.na(did$date),]

#### get intensity and center of mass
listprov <- unique(did$Prov)
datout <- NULL
for(i in 1:length(listprov)){
  sub <- did[did$Prov==listprov[i],]
  sub <- sub[order(sub$date),]
  sub$week <- week(sub$date)
  if(listprov[i]=="TIBET"){ # outlier in tibet data may be creating spurious peak time
    sub$EVA71[sub$EVA71>30] <- NA
  }
  
  # calculate average cases over whole time series
  subwEVA71 <- summaryBy(EVA71 ~ week, data = sub, na.rm=T)
  subwCVA16 <- summaryBy(CVA16 ~ week, data = sub, na.rm=T) 
  
  # get center of mass (mean timing of cases)
 comEVA71<- as.numeric(subwEVA71 %>% summarize(
    cog=weighted.mean.circular(circular(week/52*360, type="angles", units="degrees"), EVA71.mean)
  )  %>%
    mutate(
      cog=ifelse(cog < 0, cog + 360, cog)/360*52))

  comCVA16<- as.numeric(subwCVA16 %>% summarize(
   cog=weighted.mean.circular(circular(week/52*360, type="angles", units="degrees"), CVA16.mean)
 )  %>%
   mutate(
     cog=ifelse(cog < 0, cog + 360, cog)/360*52)) 
  
  # calculate intensity
  pseva71 <- subwEVA71$EVA71.mean[1:52]/sum(subwEVA71$EVA71.mean[1:52])
  pscva16 <- subwCVA16$CVA16.mean[1:52]/sum(subwCVA16$CVA16.mean[1:52])
  pseva71 <- pseva71[pseva71 != 0]
  pscva16 <- pscva16[pscva16 != 0]
  intEVA71 <-1/(-sum(pseva71*log(pseva71)))
  intCVA16 <- 1/(-sum(pscva16*log(pscva16)))  
 # output data
  datsub <- data.frame(province = listprov[i], comEVA71 = comEVA71, comCVA16 = comCVA16, intEVA71 = intEVA71, intCVA16 = intCVA16)
  datout <-rbind(datout, datsub)
 
}

rm(list=setdiff(ls(), "datout"))

# get climate data
load(file="china_adm1_levelwt.RData")

# merge
mg_ch <- merge(datout, tdatach, by = c("province"))

# plot china data - tmean and center of mass (mean timing of cases)
plot(mg_ch$tmean, mg_ch$comEVA71, ylab="",xlab="", pch=16,
     bty = "n", col=rgb(0,0,0,0), main = "EVA71 China 2009-2013", xlim = c(0,25), ylim = c(15,35))
#polygon(c(seq(5,22,1),rev(seq(5,22,1))), c(rep(0,18),rev(rep(40,18))), border=NA, col=rgb(0,0.5,1,0.2))
points(mg_ch$tmean, mg_ch$comEVA71, ylab="",xlab="", pch=16, bty = "n", col=rgb(0,0,0,0.8))
title(xlab="Temperature (\u00B0C)", line = 2)
title(ylab="Mean timing of cases (wk)", line = 2)
mg_ch <- mg_ch[order(mg_ch$tmean),]
mg_ch <- mg_ch[!is.na(mg_ch$comCVA16),]
plx<-predict(loess(mg_ch$comEVA71 ~ mg_ch$tmean), se=T)
lines(mg_ch$tmean,plx$fit)
lines(mg_ch$tmean,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(mg_ch$tmean,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)

plot(mg_ch$tmean, mg_ch$comCVA16, ylab="",xlab="", 
     pch=16, bty = "n", col=rgb(0,0,0,0), main ="CVA16 China 2009-2013" , xlim = c(0,25), ylim = c(15,35))
points(mg_ch$tmean, mg_ch$comCVA16, ylab="",xlab="", pch=16, bty = "n", col=rgb(0,0,0,0.8))
title(xlab="Temperature (\u00B0C)", line = 2)
title(ylab="Mean timing of cases (wk)", line = 2)
mg_ch <- mg_ch[order(mg_ch$tmean),]
plx<-predict(loess(mg_ch$comCVA16 ~ mg_ch$tmean), se=T)
lines(mg_ch$tmean,plx$fit)
lines(mg_ch$tmean,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(mg_ch$tmean,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)


#normalize intensity between 0 and 1
mg_ch <- mg_ch[!is.nan(mg_ch$intEVA71),]
mg_ch$intCVA16_norm <- (mg_ch$intCVA16 - min(mg_ch$intCVA16))/(max(mg_ch$intCVA16) - min(mg_ch$intCVA16))
mg_ch$intEVA71_norm <- (mg_ch$intEVA71 - min(mg_ch$intEVA71))/(max(mg_ch$intEVA71) - min(mg_ch$intEVA71))

# plot china data - trange and intensity
plot(mg_ch$trange, mg_ch$intEVA71_norm, ylab="",
     xlab="", pch=16, bty = "n", col=rgb(0,0,0,0), main = "EVA71 China 2009-2013", xlim = c(10,40))
points(mg_ch$trange, mg_ch$intEVA71_norm, ylab="",xlab="", pch=16, bty = "n", col=rgb(0,0,0,0.8))
title(xlab="Temperature range (\u00B0C)", line = 2)
title(ylab="Epidemic intensity", line = 2)
mg_ch <- mg_ch[order(mg_ch$trange),]
plx<-predict(loess(mg_ch$intEVA71_norm ~ mg_ch$trange), se=T)
lines(mg_ch$trange,plx$fit)
lines(mg_ch$trange,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(mg_ch$trange,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)

plot(mg_ch$trange, mg_ch$intCVA16_norm, ylab="",
     xlab="", pch=16, bty = "n", col=rgb(0,0,0,0), main ="CVA16 China 2009-2013", xlim = c(10,40))
points(mg_ch$trange, mg_ch$intCVA16_norm, ylab="",xlab="", pch=16, bty = "n", col=rgb(0,0,0,0.8))
title(xlab="Temperature range (\u00B0C)", line = 2)
title(ylab="Epidemic intensity", line = 2)
mg_ch <- mg_ch[order(mg_ch$trange),]
plx<-predict(loess(mg_ch$intCVA16_norm ~ mg_ch$trange), se=T)
lines(mg_ch$trange,plx$fit)
lines(mg_ch$trange,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(mg_ch$trange,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)

