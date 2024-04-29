library("doBy")
library("tsiR")
library("fields")
library("pals")
require("lomb")
require("dplyr")
require("lubridate")
require("scales")
require("tidyr")
require("here")
here::i_am("projections_analysis.R")
includevar =FALSE # whether to include variability

# load longer time series of observational data to use for variability analysis
obs <- read.csv("era5.t2m.daily.China.subregion.1990-2020.K.csv")
obs <- obs %>% gather(prov_combined, temp, -time)
obs$province<- gsub("([A-Za-z]+).*", "\\1", obs$prov_combined)
obs$week <- week(obs$time)
obs$year <- year(obs$time)
obs <- summaryBy(temp ~ province + year + week, data =obs)
obs$province[obs$province=="Inner"] <- "Inner Mongolia"
obs$temp <- obs$temp.mean - 273.15

# get climatology data
historic <- read.csv("historic_daily_proj_China.csv")
historic$temp_historic <- historic$temp
ssp585 <- read.csv("ssp585_daily_proj_China.csv")
ssp585$temp_s585 <- ssp585$temp
temp <- merge(historic, ssp585, by = c("province","day","modname"))
temp$tempchange <- temp$temp_s585 - temp$temp_historic
temp <- temp[c("province","day","modname","tempchange")]
temp <- temp[order(temp$province, temp$modname, temp$day),]
temp <- temp[temp$day != 366,] 
weekid <- floor(seq(1,53,length = 366)[1:365])
temp$week <- rep(weekid, length = 1000000)[1:nrow(temp)]
temp <- temp %>% group_by(province, week, modname) %>% summarize(tempchange = mean(tempchange))

# load empirical beta dataset just to get province list and demographic data
load(file="eva71.RData")
prov_list_alphabetical <- unique(eva$province)
load(file = "China_latlons.RData")
latlons <- latlons[rev(order(latlons$Lat.mean)),]
prov_list <- latlons$province

#load regression results
load(file = "jp_ch_cva_tempschooling.RData")
# for EVA load(file = "jp_ch_eva_tempschooling.RData")
int <- lmalllong$coefficients
int <- as.data.frame(int)
int$province <- "years"
intlist <- c(c("int","temp","schl"), prov_list_alphabetical[2:15],"Japan",prov_list_alphabetical[16:31])
int$province[1:length(intlist)] <- intlist
coeff_temperature <- int$int[int$province=="temp"]
coeff_school <- int$int[int$province=="schl"]

outallmod <- NULL
outallcv <- NULL
outallmaxi <- NULL
outallchangetemp <- NULL
outallchangebeta <- NULL

ls_modnames <- unique(temp$modname)
for(i in 1:length(ls_modnames)){
  # select climate model
  sub_temp <- temp[temp$modname==ls_modnames[i],]
  df1 <- merge(eva, sub_temp, by =c("province","week"))
  
  changeInI <- rep(NA, length = length(prov_list))
  cvRel <- rep(NA, length = length(prov_list))
  changeInMaxI <- rep(NA, length = length(prov_list))
  changetemp <- rep(NA, length = length(prov_list))
  changebeta <- rep(NA, length = length(prov_list))
  
for(j in 1:length(prov_list)){

# load province level data
sub_prov <- df1[df1$province==prov_list[j],]
sub_prov <- sub_prov[order(sub_prov$year, sub_prov$week),]
sub_prov$temp<- sub_prov$temp - 273.15
sumclim <- summaryBy(temp + tempchange + prec ~ week, data = sub_prov, na.rm=T)
sumclim$schlterm <- 1
sumclim$schlterm[sumclim$week >= 5 & sumclim$week <= 8] <- 0
sumclim$schlterm[sumclim$week >= 27 & sumclim$week <= 34] <- 0

# get intercept
intuse <- int$int[1] + max(int$int[int$province==prov_list[j]],0)

# get historic and future climatology using delta change method
temphistoric <- rep(sumclim$temp.mean, times = 100)
tempfuture <- rep((sumclim$temp.mean + sumclim$tempchange.mean), times = 100)
range_temp_change <-  diff(range(tempfuture)) - diff(range(temphistoric))

# get variability - historic and future variabilty is the same
varobs <- obs[obs$province==prov_list[j], ]
lm1 <- lm(temp ~ year + as.factor(week), data = varobs)
var_sub <- residuals(lm1)
varhist <- rep(var_sub, length = 5200)
varfuture <- rep(var_sub, length = 5200)

# choose whether to include variability
if(includevar == TRUE){
  temphistoricV <- temphistoric + varhist
  tempfutureV <- tempfuture + varfuture
}
if(includevar == FALSE){
  temphistoricV <- temphistoric 
  tempfutureV <- tempfuture 
  
}

# make transmission terms 
week_beta <- c(0,int$int[72:96])
beta <- exp(temphistoricV*coeff_temperature + sumclim$schlterm*coeff_school + intuse ) - 1
beta2 <- exp(tempfutureV*coeff_temperature + sumclim$schlterm*coeff_school + intuse ) - 1
betaconst <- exp(temphistoric*coeff_temperature + sumclim$schlterm*coeff_school + intuse ) - 1
betaconst2 <- exp(tempfuture*coeff_temperature + sumclim$schlterm*coeff_school + intuse ) - 1
range_beta_change <-  diff(range(betaconst2)) - diff(range(betaconst))

# projections using historic climate
pred1 <- predicttsir(times = seq(1,5200),births = mean(sub_prov$births[sub_prov$year==2013]), 
                     beta = rep(beta, times = 100)/floor(mean(sub_prov$Pop[sub_prov$year==2013])), 
                     alpha = 0.975, nsim=10,
                     S0 = max(floor(mean(sub_prov$Pop[sub_prov$year==2013]))) - 1, I0 = 1, stochastic = F)
Ilong <- pred1$I$mean[3640:5200]

# projections with climate change
pred2 <- predicttsir(times = seq(1,5200),births = mean(sub_prov$births[sub_prov$year==2013]), 
                     beta = rep(beta2, times = 100)/floor(mean(sub_prov$Pop[sub_prov$year==2013])), 
                     alpha = 0.975, nsim=10,
                     S0 = max(floor(mean(sub_prov$Pop[sub_prov$year==2013]))) - 1, I0 = 1, stochastic = F)
Ilong2 <- pred2$I$mean[3640:5200]

# make dataset for output by calculating annual peak size for 30 years
df1f <- data.frame(time = rep(seq(1,52,1), length = length(Ilong)), year = rep(seq(1,31,1), each = 52)[1:length(Ilong)], I = Ilong)
df1maxs <-df1f %>% group_by(year) %>% summarize(maxI = max(I))
df1maxs <- df1maxs[-31,]
df2f <- data.frame(time = rep(seq(1,52,1), length = length(Ilong)), year = rep(seq(1,31,1), each = 52)[1:length(Ilong)], I = Ilong2)
df2maxs <-df2f %>% group_by(year) %>% summarize(maxI = max(I))
df2maxs <- df2maxs[-31,]
df1maxs$cc <- "Baseline"
df2maxs$cc <- "RCP8.5"

# calculate key statistics
changeInI[j] <- (mean(df2maxs$maxI) - mean(df1maxs$maxI))/mean(df1maxs$maxI)
cvRel[j] <- (sd(df2maxs$maxI)/mean(df2maxs$maxI) - sd(df1maxs$maxI)/mean(df1maxs$maxI))/((sd(df1maxs$maxI)/mean(df1maxs$maxI)))
changeInMaxI[j] <- (max(Ilong2) - max(Ilong))/max(Ilong)
changetemp[j] <- range_temp_change
changebeta[j] <- range_beta_change
}
changeInI <- rev(changeInI)
cvRel <- rev(cvRel)
changeInMaxI <- rev(changeInMaxI)
changetemp  <- rev(changetemp)
changebeta  <- rev(changebeta) # for plotting

outallmod <- rbind(outallmod,changeInI)
outallcv <- rbind(outallcv,cvRel)
outallmaxi <- rbind(outallmaxi, changeInMaxI)
outallchangetemp <- rbind(outallchangetemp, changetemp)
outallchangebeta <- rbind(outallchangebeta, changebeta)

}



# PLOT MEAN CHANGES
# remove outliers for visualization
outallmod[outallmod > 0.1] <- 0.1
outallmod[outallmod < -0.1] <- -0.1
pdf("cva71_runProjAnalysis_meanchanges.pdf",height=8,width=6)
par(mar=c(8,7,1,1))
image.plot(x = seq(1,length(ls_modnames)),
           y =  seq(1,length(prov_list),1), 
           z= outallmod, col = coolwarm(15),ylab="",yaxt="n",xlab=" ",xaxt="n",zlim=c(-0.1,0.1))
axis(2,at = seq(1,length(prov_list),1), rev(prov_list), las=2, cex =0.9 )
axis(1,at=seq(1,length(ls_modnames),1), ls_modnames, srt = 45, las =2)
dev.off()

# PLOT CHANGES IN COEFFICIENT OF VARIATION
# remove outliers for visualization
outallcv[ outallcv > 0.3] <- 0.3
outallcv[ outallcv < -0.3] <- -0.3
pdf("cva71_runProjAnalysis_CVchanges.pdf",height=8,width=6)
par(mar=c(8,7,1,1))
image.plot(x = seq(1,length(ls_modnames)),
           y =  seq(1,length(prov_list),1), 
           z= outallcv, col =  coolwarm(15),ylab="",yaxt="n",xlab=" ",xaxt="n", zlim=c(-0.3,0.3))
axis(2,at = seq(1,length(prov_list),1), rev(prov_list), las=2, cex =0.9)
axis(1,at=seq(1,length(ls_modnames),1), ls_modnames, srt = 45, las =2)
dev.off()

# PLOT CHANGES IN MAX I (ONLY WHEN includevar=TRUE)
# remove outliers for visualization
outallmaxi[ outallmaxi > 0.2] <- 0.2
outallmaxi[ outallmaxi < -0.2] <- -0.2
pdf("cva71_runProjAnalysis_MaxIchanges.pdf",height=8,width=6)
par(mar=c(8,7,1,1))
image.plot(x = seq(1,length(ls_modnames)),
           y =  seq(1,length(prov_list),1), 
           z= outallmaxi, col =  coolwarm(15),ylab="",yaxt="n",xlab=" ",xaxt="n", zlim=c(-0.2,0.2))
axis(2,at = seq(1,length(prov_list),1), rev(prov_list), las=2, cex =0.9)
axis(1,at=seq(1,length(ls_modnames),1), ls_modnames, srt = 45, las =2)
dev.off()


pdf("temperature_range_change.pdf",height=8,width=6)
par(mar=c(8,7,1,1))
image.plot(x = seq(1,length(ls_modnames)),
           y =  seq(1,length(prov_list),1), 
           z= outallchangetemp, col = coolwarm(15),ylab="",yaxt="n",xlab=" ",xaxt="n",zlim=c(-4.25, 4.25))
axis(2,at = seq(1,length(prov_list),1), rev(prov_list), las=2, cex =0.9 )
axis(1,at=seq(1,length(ls_modnames),1), ls_modnames, srt = 45, las =2)
dev.off()

outallchangebeta[outallchangebeta > 1] <- 1
outallchangebeta[outallchangebeta < -1] <- -1
pdf("beta_range_change.pdf",height=8,width=6)
par(mar=c(8,7,1,1))
image.plot(x = seq(1,length(ls_modnames)),
           y =  seq(1,length(prov_list),1), 
           z= outallchangebeta, col = coolwarm(15),ylab="",yaxt="n",xlab=" ",xaxt="n",zlim=c(-1,1))
axis(2,at = seq(1,length(prov_list),1), rev(prov_list), las=2, cex =0.9 )
axis(1,at=seq(1,length(ls_modnames),1), ls_modnames, srt = 45, las =2)
dev.off()
