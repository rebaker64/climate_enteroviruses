require("lubridate")
require("doBy")
require("pals")
library("sp")
library("rgdal")
require("tsiR")
require("dplyr")
require("tidyr")


# load enterovirus data
disease <- read.csv("china_entero.csv")
cid <- read.csv("china_id.csv")
did <- merge(disease, cid, by.x=c("PROV_ID"),by.y=c("ID"))
did$Prov <- toupper(did$Prov)
did$date <- as.Date(paste0(1,"/",did$WEEK,"/",did$YEAR), format = c("%w/%W/%Y"))
did$month <- month(did$date)
did <- did[!is.na(did$date),]

# load demographic data
demog <- read.csv("china_demog.csv")
demog$Pop <- demog$Pop*10000
demog$births <- (demog$CBR*demog$Pop/1000)/52
mg <- merge(did, demog, by.x=c("PROV_ID","YEAR"), by.y=c("id","Year"), all.x=T)


# make empirical betas, fitting model to all provinces
provlist <- unique(mg$province)
outallemp <- NULL
storebetas <- NULL
for(i in 1:length(provlist)){
  sub <- mg[mg$province==provlist[i],]
  sub <- sub[order(sub$date),]
  sub <- sub[sub$WEEK !=53, ]
  subtsir <- sub[c("date","EVA71","births","Pop")]
  names(subtsir) <- c("time","cases","births","pop")
  plot(sub$date,sub$EVA71,type="l", main = sub$province[1])
  Prov.Params <- estpars(data = subtsir, IP = 1, alpha = 0.975, sbar = NULL, regtype = 'lm')
  tsir <- simulatetsir(data=subtsir, IP = 1, parms=Prov.Params, nsim=100, inits.fit = TRUE, method="negbin")
  
  if(sum(tsir$res$mean==0) >50){
    tsir2 <- simulatetsir(data=subtsir, IP = 1, parms=Prov.Params, method="negbin")  # re-run without inits fit if selecting on run with a lot of zeros
    if(sum(tsir$res$mean==0) > sum(tsir2$res$mean==0)){tsir <- tsir2}
    if(sum(tsir$res$mean==0) >50){
      tsir3 <- simulatetsir(data=subtsir, IP = 1, parms=Prov.Params, nsim=100, epidemic="break", threshold = 3, method="negbin")
      if(sum(tsir$res$mean==0) > sum(tsir3$res$mean==0)){tsir <- tsir3} # some low observations at the start of time series causing issues, use epidemics = break
    }
  }
  lines(sub$date,tsir$res$mean, lwd = 2, col="red")
  
  It <- sub$EVA71*tsir$rho
  It_lead <- lead(It)
  sub$lcases <- lead(sub$EVA71)
  sub$empbeta <- (sub$lcases*sub$Pop)/((sub$EVA71^(0.975))*(tsir$sbar + tsir$Z))
  outallemp <- rbind(outallemp,sub)
  df <- data.frame(week = seq(1,52,1), prov = rep(sub$province[1],times = 52), beta = tsir$beta)
  storebetas <- rbind(storebetas,df)
}

# merge with climate data
shum <- read.csv("era5.q2m.daily.China.subregion.2009-2014.kg_per_kg.csv")
shum <- shum %>% gather(prov_combined, shum, -time)
shum$province<- gsub("([A-Za-z]+).*", "\\1", shum$prov_combined)
shum$week <- week(shum$time)
shum$year <- year(shum$time)
shum <- summaryBy(shum ~ province + year + week, data =shum)
shum <- shum %>%
  group_by(province) %>%
  mutate(lgshum = dplyr::lag(shum.mean, n = 1, default = NA))
shum$province[shum$province=="Inner"] <- "Inner Mongolia"

temp <- read.csv("era5.t2m.daily.China.subregion.2009-2014.K.csv")
temp <- temp %>% gather(prov_combined, temp, -time)
temp$province<- gsub("([A-Za-z]+).*", "\\1", temp$prov_combined)
temp$week <- week(temp$time)
temp$year <- year(temp$time)
temp <- summaryBy(temp ~ province + year + week, data =temp)
temp <- temp %>%
  group_by(province) %>%
  mutate(lgtemp = dplyr::lag(temp.mean, n = 1, default = NA))
temp$province[temp$province=="Inner"] <- "Inner Mongolia"

prec <- read.csv("chirps.precip.daily.China.subregion.2009-2014.mm_per_day.csv")
prec <- prec %>% gather(prov_combined, prec, -time)
prec$province<- gsub("([A-Za-z]+).*", "\\1", prec$prov_combined)
prec$week <- week(prec$time)
prec$year <- year(prec$time)
prec <- summaryBy(prec ~ province + year + week, data =prec, FUN = sum)
prec <- prec %>%
  group_by(province) %>%
  mutate(lgprec = dplyr::lag(prec.sum, n = 1, default = NA))
prec$province[prec$province=="Inner"] <- "Inner Mongolia"

outallemp$week <- week(outallemp$date)
outallemp$year <- year(outallemp$date)
outall <- merge(outallemp,shum,by=c("year","week","province"))
outall <- merge(outall,temp,by=c("year","week","province"))
outall <- merge(outall,prec,by=c("year","week","province"))

outall$empbeta[is.nan(outall$empbeta)] <- NA
outall$empbeta[is.infinite(outall$empbeta)] <- NA
finaldat <- outall[!names(outall) %in% c("PROV_ID","Prov","CBR","WEEK")]
finaldat$empbeta[finaldat$EVA71==0] <- NA

finaldat$cases <- finaldat$EVA71
finaldat$shum <- finaldat$shum.mean
finaldat$temp <- finaldat$temp.mean
finaldat$prec <- finaldat$prec.sum

finaldat$week <- week(finaldat$date)
finaldat$schlterm <- 1
finaldat$schlterm[finaldat$week >= 5 & finaldat$week <= 8] <- 0
finaldat$schlterm[finaldat$week >= 27 & finaldat$week <= 34] <- 0

cva <- finaldat[c("year","week","province","cases","temp","prec","shum","schlterm","empbeta","Pop","births")]

save(cva, file="EVA71.RData")


