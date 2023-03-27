require("lfe")

### POLIO
load(file = "polioeb.RData")

lm1po <- felm(log(empbeta + 1) ~ temp | state | 0 | state , data =df)
sumpo1 <- summary(lm1po)

lm2po <- felm(log(empbeta + 1) ~ temp | state + year| 0 | state, data =df)
sumpo2 <- summary(lm2po)

lm3po <- felm(log(empbeta + 1) ~   temp | state + year + week| 0| state, data =df)
sumpo3 <- summary(lm3po)


### CVA16
#China
load(file="cva16.RData")

# Japan
load(file="cvajp.RData")
 
# both 
names(cvajp) <- names(cva)
cva <- rbind(cva, cvajp)
cva$temp <- cva$temp - 273.15

# run regressions
lm1cva <- felm(log(empbeta + 1) ~ temp | province | 0 | province, data =cva)
sumcva1 <- summary(lm1cva)

lm2cva <- felm(log(empbeta + 1) ~ temp | province +  year | 0 | province , data =cva)
sumcva2 <- summary(lm2cva)

lm3cva <- felm(log(empbeta + 1) ~   temp + schlterm | province + year| 0| province, data =cva)
sumcva3 <- summary(lm3cva)


### EVA71
# China
load(file="eva71.RData")

## Japan
load(file="evajp.RData")

# both 
names(evajp) <- names(eva)
eva <- rbind(eva, evajp)
eva$temp <- eva$temp - 273.15

# run regressions
lm1eva <- felm(log(empbeta + 1) ~ temp | province | 0 | province, data =eva)
sumeva1 <- summary(lm1eva)

lm2eva <- felm(log(empbeta + 1) ~ temp | province + year | 0 | province, data =eva)
sumeva2 <- summary(lm2eva)

lm3eva <- felm(log(empbeta + 1) ~   temp + schlterm | province + year | 0| province, data =eva)
sumeva3 <- summary(lm3eva)


coeffdat <- data.frame(model = rep(c("Polio","EVA71","CVA16"), each = 3), 
                       meaneffect= c(  sumpo1$coefficients[1,1], 
                                       sumpo2$coefficients[1,1],
                                       sumpo3$coefficients[1,1],
                                       sumeva1$coefficients[1,1], 
                                       sumeva2$coefficients[1,1],
                                       sumeva3$coefficients[1,1],
                                       sumcva1$coefficients[1,1], 
                                       sumcva2$coefficients[1,1],
                                       sumcva3$coefficients[1,1]),
                       se= c(  sumpo1$coefficients[1,2], 
                                       sumpo2$coefficients[1,2],
                                       sumpo3$coefficients[1,2],
                                       sumeva1$coefficients[1,2], 
                                       sumeva2$coefficients[1,2],
                                       sumeva3$coefficients[1,2],
                                       sumcva1$coefficients[1,2], 
                                       sumcva2$coefficients[1,2],
                                       sumcva3$coefficients[1,2]), 
                       controls = c("location","+ year","+ week", "location","+ year","+ schooling","location",
                                    "+ year","+ schooling"),
                       id = c(9,8,7,6,5,4,3,2,1))
 coeffdat$maxci <- coeffdat$meaneffect + 1.96*coeffdat$se                      
 coeffdat$minci <- coeffdat$meaneffect - 1.96*coeffdat$se   
 
 #### plot
 
 require("pals")
 palpar <- rev(rep(parula(5)[1:3], each = 3))
 
 coeffdat <- coeffdat[order(coeffdat$id),]
 
 pdf("coeffs.pdf", width = 5, height = 4)
 par(mar=c(3,3,3,5))
 plot(coeffdat$meaneffect, seq(1,9),xlim=c(0,max(coeffdat$maxci)),ylab="", yaxt="n", pch = 16
      ,col=palpar, xlab="")
 abline(v = 0, col="grey64", lwd = 2, lty = 3)
 title(xlab="Estimated effect", line = 2)
for(i in 1:nrow(coeffdat)){
  segments(y0 = i,y1 = i, x0 = coeffdat$minci[i], x1 = coeffdat$maxci[i], col=palpar[i])
}
axis(4,at = seq(1,9,1), coeffdat$controls, las = 1, tcl = 0, cex.axis = 0.8)
abline(h = 3.5, col="grey64") 
abline(h = 6.5, col="grey64") 
text(0,8.8,"Polio", pos = 4,col=palpar[8])
text(0,6,"EVA71", pos = 4, col=palpar[4])
text(0,3,"CVA16", pos = 4, col=palpar[1])
 
dev.off()
 
 