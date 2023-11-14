#######################################################################
# Hotspot analysis and population modeling of midland brownsnakes
# (Storeria dekayi) at Fox Ridge State Park (Coles County, IL)
# Date created: 6 August 2020
#######################################################################
#####   ##############################################################
###     #############################################################
#####   ##########           Data preparation           ############
#####   ###########################################################
###       ########################################################
#################################################################
# Clear the environment
rm(list = ls())

setwd("~/Dropbox/stodek_road_mortality/")

# Install and/or load necessary packages
packages = c("anchors", "dplyr", "ggplot2", "RMark", "tidyr")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) { 
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Import data, classify "." == NA
datum <-read.csv("stodek-frsp.csv", na.strings = ".")
weatherOG <- read.csv("Charleston_weather_raw.csv", na.strings = ".")
interpol <- read.csv('Oct2013 interpolated weather data.csv')

# Exclude captures from 2010-2011
# Also exclude extraneous columns from this first run
# Remember to revert to all captures for the hotspot analysis
dat <- datum[rownames(subset(datum, year > 2011)), -c(6, 10, 25:31)]
# Exclude batch snakeIDs for now
batches <- c("314", "115", "315", "214", "313", "114", "215", "715", "515")
dat <- dat[!as.factor(dat$snakeID) %in% batches, ]

# Format your variables...
dat$snakeID <- as.factor(dat$snakeID)
dat$alive <- as.factor(dat$alive)
dat$orientation <- as.factor(dat$orientation)
dat$sex <- as.factor(dat$sex)
dat$gravid <- as.factor(dat$gravid)
dat$gut <- as.factor(dat$gut)
dat$svl <- as.numeric(dat$svl)
dat$tl <- as.numeric(dat$tl)
dat$bci <- as.numeric(dat$bci)
dat$recap <- as.factor(dat$recap)

# Form a single "date" variable
dat$date <- as.Date(with(dat, paste(year, month, day, sep= "-")), "%Y-%m-%d")

# Create separate columns for julian week and season based on specified row date
dat$week <- NA
dat$season <- NA

# Populating NA vectors with data
# Assign julian week ID
dat$week <- strftime(dat[,"date"], format = "%V")
dat$week <- as.numeric(dat$week)
dat$season[which(dat$week > 9 & dat$week < 24)] <- "spring"
dat$season[which(dat$week > 23 & dat$week < 36)] <- "summer"
dat$season[which(dat$week > 35 & dat$week < 47)] <- "autumn"
dat$season <- as.factor(dat$season)

# relevel "season" 
dat$season <- factor(dat$season, levels = c("spring", "summer", "autumn"))

# Examine the data
#str(dat)
#summary(dat)

# Produce chronological week and season columns
dat$weekYr <- with(dat, paste(year, week, sep = "-"))
dat$seasonYr <- with(dat, paste(year, season, sep = "-"))

# Produce a table ticking off captures per day for all individuals (and batches)
# Produce additional tables that bin the daily captures by week and by season
# The "seasonalCaps" table is complete, but "captures" and "weekJustCaps" need to have no-snake periods incorporated
captures <- as.data.frame.matrix(table(dat$snakeID,dat$date))
weeklyJustCaps <- as.data.frame.matrix(table(dat$snakeID,dat$weekYr))
seasonalCaps <- as.data.frame.matrix(table(dat$snakeID,dat$seasonYr))

# Next, build a similarly structured table encompassing all survey days and weeks
# Function prints all dates within a specific interval
itemizeDates <- function(startDate, endDate, 
                         format="%Y-%m-%d") {
  out <- seq(as.Date(startDate, format=format), 
             as.Date(endDate, format=format), by="days")  
  format(out, format)
}

# Run the above function for year-specific date intervals
# Changed the start and end dates so they correspond to julian weeks 10 and 46, respectively
YR12a <- itemizeDates(startDate = "2012-03-05", endDate = "2012-06-24")
YR12b <- itemizeDates(startDate = "2012-09-29", endDate = "2012-11-16")
YR13 <- itemizeDates(startDate = "2013-03-05", endDate = "2013-11-16")
YR14 <- itemizeDates(startDate = "2014-03-05", endDate = "2014-11-16")
YR15 <- itemizeDates(startDate = "2015-03-05", endDate = "2015-11-15")

### Capture tables ###
# Merge all years (for daily and weekly tables)
allSurveyDays <- c(YR12a, YR12b, YR13, YR14, YR15)
allSurveyWeeks <- unique(strftime(allSurveyDays, format = "%Y-%V"))

# Construct the tables containing all survey days and weeks
nocaptures <- data.frame(matrix(nrow = nrow(captures), ncol = length(allSurveyDays)))
colnames(nocaptures) <- allSurveyDays
weeklyNoCaps <- data.frame(matrix(nrow = nrow(captures), ncol = length(allSurveyWeeks)))
colnames(weeklyNoCaps) <- allSurveyWeeks

# Populate with "0" captures
nocaptures[] <- 0
weeklyNoCaps[] <- 0

# Bind snake-day-only and all-day tables
captures <- cbind(captures, nocaptures)
weeklyCaps <- cbind(weeklyJustCaps, weeklyNoCaps)

# Remove duplicates. The "captures" columns (snake days or weeks) are unaffected, but any duplicated dates
## from "nocaptures" or "weeklyNoCaps" are removed, leaving behind all zero-snake days/weeks
captures <- captures[,!duplicated(names(captures))]
weeklyCaps <- weeklyCaps[,!duplicated(names(weeklyCaps))]

# So "captures" should now have the same dimensions as "nocaptures"
dim(captures) == dim(nocaptures)
#[1] TRUE TRUE # Good!
dim(weeklyCaps) == dim(weeklyNoCaps)
#[1] TRUE TRUE # Good!

# Order date columns chronologically
dailyCaps <- captures[,order(colnames(captures))]
weeklyCaps <- weeklyCaps[,order(colnames(weeklyCaps))]

### FINAL tables: daily, weekly, seasonal ###
# dailyCaps, weeklyCaps, seasonalCaps

### Covariates ###
# Organizing day, week, and season columns
covars <- data.frame("date" = allSurveyDays)
covars$date <- as.Date(covars$date, "%Y-%m-%d")
covars$week <- NA

# Populating "week" column
covars$week <- strftime(covars[,"date"], format = "%V")
covars$week <- as.numeric(covars$week)
covars$season <- NA

# Populating "season" column
# If there's any easier way to command this, or ANY of this code, please let me know
covars$season[which(covars$week > 9 & covars$week < 24)] <- paste(strftime(covars$date[which(covars$week > 9 & covars$week < 24)], format = "%Y"), "spring" , sep = "-")
covars$season[which(covars$week > 23 & covars$week < 36)] <- paste(strftime(covars$date[which(covars$week > 23 & covars$week < 36)], format = "%Y"), "summer" , sep = "-")
covars$season[which(covars$week > 35 & covars$week < 47)] <- paste(strftime(covars$date[which(covars$week > 35 & covars$week < 47)], format = "%Y"), "autumn" , sep = "-")

# Adding a year value to "week" so they aren't binned across years
covars$week[] <- strftime(covars[,"date"], format = "%Y-%V")

# Weather
# Removed all but "DATE", "PRCP", "TMAX", "TMIN"
# Only one station (USC00111436; 39.4762N	-88.1652W) collected reliable temperature and precip data from the Charleston area
weather <- weatherOG[weatherOG$STATION == "USC00111436", c(6, 9, 12, 13)]
# Reformat and rename the "DATE" variable
weather$DATE <- as.Date(weather$DATE, "%m/%d/%y")
weather$date <- weather$DATE
# Reorder by "date"
weather <- weather[order(weather$date), ]
# Select only those rows in "weather" that corresponde to date in "covars"
weather <- weather[weather$date %in% covars$date, ]
# Merge "covars" and "weather" by "date"
bonk <-merge(covars, weather, by=c("date"))
# Remove redundant column and rename others
covars <- bonk[, -c(4)]
colnames(covars)[4:6] <- c("precip", "tmax", "tmin")

# Calculates number of visits made during each week and each season

covars$visPerWeek <- NA
covars$visPerSeason <- NA

# Provides counts of each unique week and date value
weekTab <- as.data.frame((table(covars[,"week"])))
seasonTab <- as.data.frame((table(covars[,"season"])))

# Pulls the count designated for each week or date value and applies it to the appropriate column
for (i in 1 : nrow(covars)) {
covars$visPerWeek[i] <- weekTab[which(weekTab$Var1 %in% covars$week[i]),2]
covars$visPerSeason[i] <- seasonTab[which(seasonTab$Var1 %in% covars$season[i]),2]
}

######################################
## Extracted from "convert covars.R" 
######################################

covars$yday <- as.numeric(format(covars$date, "%j"))

covars2 <- covars %>%
  group_by(week) %>%
  summarize(mean.tmin = mean(tmin,na.rm = T),
            sd.tmin = sd(tmin,na.rm=T),
            mean.tmax = mean(tmax,na.rm = T),
            sd.tmax = sd(tmax,na.rm=T),
            total.precip = sum(precip),
            visPerWeek=mean(visPerWeek),
            min.yday = min(yday))

setdiff(names(weeklyCaps),covars2$week)
setdiff(covars2$week,names(weeklyCaps))

october <- rbind(covars2[grep(39,covars2$week),],covars2[grep(40,covars2$week),],
                 covars2[grep(41,covars2$week),],covars2[grep(42,covars2$week),],
                 covars2[grep(43,covars2$week),],covars2[grep(44,covars2$week),],
                 covars2[grep(45,covars2$week),])
ggplot(october,aes(x=min.yday,y=mean.tmin))+geom_point()+
  geom_smooth(method = 'lm')+
  geom_vline(xintercept = c(285,292,299),linetype=3)
dev.off()

covars3 <- rbind(covars2[1:55,],interpol,covars2[56:132,])
covars3$total.precip[is.na(covars3$total.precip)==T] <- mean(covars3$total.precip,na.rm=T)
covars3$time <- seq(1,135,1)
covars3$season <- NA
covars3$season[which(as.numeric(substr(covars3$week,6,7)) > 9 & as.numeric(substr(covars3$week,6,7)) < 24)] <- "spring"
covars3$season[which(as.numeric(substr(covars3$week,6,7)) > 23 & as.numeric(substr(covars3$week,6,7)) < 36)] <- "summer"
covars3$season[which(as.numeric(substr(covars3$week,6,7)) > 35 & as.numeric(substr(covars3$week,6,7)) < 47)] <-  "autumn"

## covars4: for regression julian week counts against temperature fluctuations
bonk <- do.call(rbind, strsplit(covars3$week, '-'))
colnames(bonk) <- c("year", "week")
covars4 <- cbind(bonk,covars3[,2:10])





save(list = c("dailyCaps", "weeklyCaps", "seasonalCaps", "covars", "covars2", "covars3", "interpol", "weatherOG", "datum", "dat"), file = "stodek.rdata")