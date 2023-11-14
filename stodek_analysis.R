########################################################################
# Updated version of "cjs - stodek.R" 
# No major changes, just restructured so all of your updates to "covars"
# are saved in "stodek.rdata"

## BACKGROUND ##
# We are examining spatiotemporal patterns of dekay's browsnake (Storeria dekayi)
# seasonal migration and road mortality. This pencil-sized snake is a "leaf-litter"
# species: considered highly abundant within its range, but its cryptic life-history
# makes the species difficult to (re)capture. Our study system is notable for hosting
# bi-annual mass migrations of brownsnakes across a 1.7-km road in a state park in IL.
# The road runs along an ecotone separating upland ridge/valley habitat (where
# snakes hibernate in winter) and lowland floodplains (where snakes feed during the
# summer). 

# We are using GIS to map out movement/mortality hotspots along the roadway,
# but we are also interested in how body size (SVL), sex, temperature/rainfall,
# and search effort (survey-days/week) influence capture rates.

## NOTES ##
# The "stodek.rdata" file contains separate capture histories parsed into either daily, weekly,
# or seasonal bins. For our analysis, we settled on the weekly resolution. The object "ch2"
# below contains the capture history from "weeklyCaps", and individual size/sex covariates from
# "dat" (the original dataset, post-QA/QC). Environmental covariates are located in "covars3".

########################################################################
######     ############################################################
####   ##   ##########################################################
########   ##########       Cormack-Jolly-Seber model        ########
######   ###########################################################
####        #######################################################
##################################################################
# Clear the environment
rm(list = ls())

#setwd('C:/Users/Andrew/Documents/EIU/Project_StoreriaFRSP/from Iwo')
setwd("~/Dropbox/stodek_road_mortality/")

load("stodek.rdata")
load("stodek_cjs_models_plots.rdata")

# Install and/or load necessary packages
packages = c("dplyr", "ggplot2", "RMark", "tidyr", "parallel", "R2ucare", "MuMIn")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) { 
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

## TODO: research zero-inflated methods

# Below commented out by Iwo Gross for now
#c <- read.csv('captures.csv')
#c2 <- c[rowSums(c[,-1]) < 5,] # max recaps is 4, min row sum for batch mark is 9
## batch marks are 314, 313, 315, 114, 115, 214, 215, 515, 715
#
#c3 <- c2[,-1]

# Convert weeklyCaps to capture history (ch) file
# First convert counts >1 (more than one capture of a snake in a week) to 1's
weeklyCaps2 <- replace(weeklyCaps, weeklyCaps >1, 1)
ch <- weeklyCaps2 %>% unite(weeklyCaps2, sep="") 
colnames(ch) <- 'ch'

## add sex back to ch
ch$sex <- dat[match(row.names(ch),dat$snakeID),"sex"]
ch$svl <- dat[match(row.names(ch),dat$snakeID),"svl"]
ch2 <- ch[is.na(ch$sex)==F,]

stodek.processed = process.data(ch2,model="CJS", groups = "sex")
stodek.ddl = make.design.data(stodek.processed)

stodek.ddl$Phi = merge_design.covariates(stodek.ddl$Phi,covars3)
#stodek.ddl$Phi
stodek.ddl$p = merge_design.covariates(stodek.ddl$p,covars3)
#stodek.ddl$p

Phi.dot=list(formula=~1)
Phi.sex=list(formula=~sex)
Phi.svl=list(formula=~svl)
Phi.sex.svl=list(formula=~sex+svl)
p.dot=list(formula=~1)
p.sex=list(formula=~sex)
p.svl=list(formula=~svl)
p.sex.svl=list(formula=~sex+svl)
p.precip=list(formula=~total.precip)
p.effort=list(formula=~visPerWeek)
p.tmin=list(formula=~mean.tmin)
p.tmax=list(formula=~mean.tmax)
p.season.tmin=list(formula=~season+mean.tmin)
p.season.tmax=list(formula=~season+mean.tmax)
p.season.precip=list(formula=~season+total.precip)
p.season.visPerWeek=list(formula=~season+visPerWeek)


#mod1=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.dot))
#mod2=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.sex))
#mod3=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.svl))
#mod4=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.sex.svl))
#mod5=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.sex,p=p.dot))
#mod6=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.svl,p=p.dot))
#mod7=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.sex.svl,p=p.dot))
#mod8=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.precip))
#mod9=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.effort))
#mod10=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.tmin))
#mod11=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.tmax))
#mod12=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.season.tmin))
#mod13=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.season.tmax))
#mod14=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.season.precip))
#mod15=mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.season.visPerWeek))
#
#results<-collect.models()
#results

#globalmodel <- mark(stodek.processed,stodek.ddl,model.parameters=list(Phi = list(formula = ~ 1), p = list(formula = ~ (season + total.precip + visPerWeek + mean.tmin + sex + svl)^2 + season + total.precip + visPerWeek + mean.tmin + sex + svl)))

globalmodel <- mark(stodek.processed,stodek.ddl,model.parameters=list(Phi = list(formula = ~ (sex)), p = list(formula = ~ (season + total.precip + visPerWeek + mean.tmin))))
model_dredge <- dredge(globalmodel)#; View(occ_dredge)
model_dredge_delta <- get.models(model_dredge, subset = delta <= 2.5)
model_avg <- model.avg(model_dredge_delta, vcv = TRUE, fit = TRUE)
coef(model_avg)
model_avg$estimates

globalmodel2 <- mark(stodek.processed,stodek.ddl,model.parameters=list(Phi = list(formula = ~ (sex)), p = list(formula = ~ (total.precip + visPerWeek + mean.tmin))))
model_dredge2 <- dredge(globalmodel2)#; View(occ_dredge)
model_dredge_delta2 <- get.models(model_dredge2, subset = delta <= 2.5)
model_avg2 <- model.avg(model_dredge_delta2, vcv = TRUE, fit = TRUE)
coef(model_avg2)
model_avg2$estimates

globalmodel3 <- mark(stodek.processed,stodek.ddl,model.parameters=list(Phi = list(formula = ~ (sex)), p = list(formula = ~ (season*mean.tmin + season*total.precip + visPerWeek + mean.tmin + total.precip + season))))
model_dredge3 <- dredge(globalmodel3)#; View(occ_dredge)
model_dredge_delta3 <- get.models(model_dredge3, subset = delta <= 2.5)
model_avg3 <- model.avg(model_dredge_delta3, vcv = TRUE, fit = TRUE)
coef(model_avg3)
model_avg3$estimates

globalmodel4 <- mark(stodek.processed,stodek.ddl,model.parameters=list(Phi = list(formula = ~ (sex)), p = list(formula = ~ (season*mean.tmin + mean.tmin + season))))
model_dredge4 <- dredge(globalmodel4)#; View(occ_dredge)
model_dredge_delta4 <- get.models(model_dredge4, subset = delta <= 2.5)
model_avg4 <- model.avg(model_dredge_delta4, vcv = TRUE, fit = TRUE)
coef(model_avg4)
model_avg4$estimates

phi.dot.p.tmin <- mark(stodek.processed,stodek.ddl,model.parameters=list(Phi=Phi.dot,p=p.tmin))

#Back-transforming for logit link blah blah
#(exp(2.75228))/(1+exp(2.75228))

# Summary of parameter estimates and standard errors
#coefTable(model_avg2)

#head(get.real(phi.dot.p.tmin, "Phi", se = TRUE))
#head(get.real(phi.dot.p.tmin, "p", se = TRUE))

# plot covariate relationship: MEAN.TMIN
df <- data.frame(mean.tmin=seq(min(covars3$mean.tmin,na.rm=T),max(covars3$mean.tmin,na.rm=T),5))

logit.values=phi.dot.p.tmin$results$beta$estimate[2]+(df$mean.tmin)*phi.dot.p.tmin$results$beta$estimate[3]
deriv=matrix(1,ncol=2,nrow=nrow(df))
deriv[,2]=df$mean.tmin
std.errors=sqrt(diag(deriv%*%phi.dot.p.tmin$results$beta.vcv[2:3,2:3]%*%t(deriv)))
lcl.logit=logit.values-1.96*std.errors
ucl.logit=logit.values+1.96*std.errors
plot(df$mean.tmin,plogis(logit.values),type="l",xlab="mean.tmin",ylab="p")
lines(df$mean.tmin,plogis(lcl.logit),lty=2)
lines(df$mean.tmin,plogis(ucl.logit),lty=2)



# plot covariate relationship: TOTAL.PRECIP
#df <- data.frame(total.precip=seq(min(covars3$total.precip,na.rm=T),max(covars3$total.precip,na.rm=T),0.25))

# For fall?
#logit.values=mod15$results$beta$estimate[2] + (df$total.precip)*mod15$results$beta$estimate[5] 

# Spring?
#logit.values=mod15$results$beta$estimate[2] + (df$total.precip)*mod15$results$beta$estimate[3] +
#                                              (df$total.precip)*mod15$results$beta$estimate[5]
# Summer?
#logit.values=mod15$results$beta$estimate[2] + (df$total.precip)*mod15$results$beta$estimate[4] +
#                                              (df$total.precip)*mod15$results$beta$estimate[5]
   
deriv=matrix(1,ncol=2,nrow=nrow(df))
deriv[,2]=df$total.precip
std.errors=sqrt(diag(deriv%*%mod15$results$beta.vcv[2:3,2:3]%*%t(deriv)))
lcl.logit=logit.values-1.96*std.errors
ucl.logit=logit.values+1.96*std.errors
plot(df$total.precip,plogis(logit.values),type="l",xlab="total.precip",ylab="p")
lines(df$total.precip,plogis(lcl.logit),lty=2)
lines(df$total.precip,plogis(ucl.logit),lty=2)

# Calculate proportion of 0's and 1's in "ch"
sum(apply(weeklyCaps, 1, function(x) sum(x == 0)))
sum(apply(weeklyCaps, 1, function(x) sum(x > 0)))

options(future.globals.maxSize = 4000 * 1024^5)

# Model averaging
#stodek.mod.avg <- model.average(results, vcv=TRUE)

# Goodness-of-fit

ch3 <- data.frame(ch2$ch)
ch3$freq <- 1
write.table(ch3, file = "ch3.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)

ch3.import <- import.chdata(paste0(getwd(),"/ch3.txt"),
                              field.names=c("ch", "freq"),
                              header=FALSE)

ch3.import$freq <- as.numeric(ch3.import$freq)

dip.hist <- matrix(as.numeric(unlist(strsplit(as.character(ch3.import$ch),""))),nrow=length(ch3.import$ch),byrow=T)
dip.freq <- ch3.import$freq
#dip.group <- ch3.import$sex

overall_CJS(dip.hist, dip.freq)
#                          chi2 degree_of_freedom p_value
#Gof test for CJS model: 11.068                27   0.997

#Summary table of model-averaged data
#subset analysis with season interaction and additive
#annual or monthly survivorship (check Gitzen notes for relevant code)

save(list = c("globalmodel", "model_dredge", "model_dredge_delta", "model_dredge2", "model_dredge_delta2", "model_dredge3", "model_dredge_delta3", "model_dredge4", "model_dredge_delta4"), file = "stodek_cjs_models_plots.rdata")