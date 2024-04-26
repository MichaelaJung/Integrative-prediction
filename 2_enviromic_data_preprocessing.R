#############
### Setup ###
#############

# load libraries
library(readxl)
library(dplyr)

# set quantile defining periods of time
myquant <- 0.90

# set number of days before flowering to consider for period 1
ndays <- 80


#################
### Load data ###
#################

# load hourly weather data
weather <- as.data.frame(read_xlsx("./Input/Weather_raw.xlsx"))

# load soil data
soil <- as.data.frame(read_xlsx("./Input/Soil_raw.xlsx"))

# load phenotypes
pheno <- as.data.frame(read_xlsx("./Output/Pheno_processed.xlsx"))


######################################
### Hourly to daily weather values ###
######################################

# add identifier for day
weather$Day <- as.POSIXct(round.POSIXt(weather$Date,"days"))

# estimate daily temperature means, daily humidity means, and daily radiation sums 
d1 <- aggregate(weather[,c("Temperature","Humidity")], FUN=mean,
                by=list(weather$Location, weather$Day), na.rm=TRUE)
d2 <- aggregate(weather[,"Radiation"], FUN=sum,
                by=list(weather$Location, weather$Day), na.rm=TRUE)
colnames(d1)[3:ncol(d1)] <- paste0(colnames(d1)[3:ncol(d1)], "_Dmean")
colnames(d2)[3] <- "Radiation_Dsum"
d <- Reduce(dplyr::full_join, list(d1,d2))
colnames(d)[1:2] <- c("Location", "Date")


###########################
### Define time periods ###
###########################

# subset phenotypic data for floral emergence and harvest date
pheno <- pheno[which(pheno$Trait %in% c("Flowering_begin", "Harvest_date")),]

# convert to wide format
pheno <- reshape2::dcast(pheno, Genotype + Envir ~ Trait, value.var="Value")

# obtain quantiles (in days of year)
dates <- pheno %>%
  group_by(Envir) %>%
  summarize(FD = quantile(Flowering_begin, probs = myquant, na.rm=T),
            HD = quantile(Harvest_date, probs = myquant, na.rm=T))
dates$FD <- round(dates$FD)
dates$HD <- round(dates$HD)
dates <- as.data.frame(dates)

# replace missing values of quantiles for flowering date
# Spain 2020 (low heritability):
# derived manually as 0.9 quantile before removal of the environment
dates[which(dates$Envir %in% "ESP.2020"),"FD"] <- 91
# Spain 2022 (not measured):
# estimated for other varieties in the same location (on 6.4.2022)
dates[which(dates$Envir %in% "ESP.2022"),"FD"] <- 96

# convert numeric quantiles into dates and years
dates[,"FDdate"] <- as.POSIXct(NA)
dates[,"HDdate"] <- as.POSIXct(NA)
dates$Year <- as.numeric(stringr::str_split(dates$Envir, pattern="\\.", simplify = T)[,2])
for (i in 1:nrow(dates)){
  dates$FDdate[i] <- eseis::time_convert(input=dates$FD[i],
                                       output="POSIXct",
                                       timezone = "UTC",
                                       year=dates$Year[i])
}
for (i in 1:nrow(dates)){
  dates$HDdate[i] <- eseis::time_convert(input=dates$HD[i],
                                       output="POSIXct",
                                       timezone = "UTC",
                                       year=dates$Year[i])
}

# remove time component from date-time
dates$FDdate <- format(as.POSIXct(dates$FDdate), format = "%Y-%m-%d")
dates$HDdate <- format(as.POSIXct(dates$HDdate), format = "%Y-%m-%d")


#################################################
### Sum weather covariables over time periods ###
#################################################

# assign time periods to the entries in the environmental covariables file
d$Period <- NA
d$Date_char <- format(as.POSIXct(d$Date), format = "%Y-%m-%d")
d$Envir <- paste0(d$Location, ".",substr(d$Date_char, 0, 4))
for (i in 1:nrow(dates)) {
  env <- dates$Envir[i]
  date1 <- format(eseis::time_convert(input=dates$FD[i] - ndays,
                               output="POSIXct",
                               timezone = "UTC",
                               year=dates$Year[i]),format = "%Y-%m-%d")
  date2 <- dates$FDdate[i]
  date3 <- dates$HDdate[i]
  d[which(d$Envir %in% env & d$Date_char >= date1 & d$Date_char <= date2),"Period"] <- "P1"
  d[which(d$Envir %in% env & d$Date_char > date2 & d$Date_char <= date3),"Period"] <- "P2"
}

# save intermediate output
save(d, file="./Output/Weather_daily.Rdata")

# delete entries with no time period assigned (time after harvest/before flowering)
d <- d[-which(is.na(d$Period)),]

# obtain sums over time periods
out <- NULL
for (i in 3:5) {
  x <- aggregate(x=d[,i], by=list(d$Period, d$Envir), FUN=sum, na.rm=TRUE)
  x$Variable <- colnames(d)[i]
  out <- rbind(out, x)
}

# format output
out$Variable <- paste0(out$Variable, "_", out$Group.1, "sum")
out$Group.1 <- NULL
colnames(out)[1:2] <- c("Envir", "Value")


############################
### Add soil covariables ###
############################

# replicate soil covariables for each environment
out2 <- NULL
for (i in unique(substr(d$Date_char, start=0, stop=4))) {
  x <- soil
  x$Envir <- paste0(x$Group.1, ".", i)
  out2 <- rbind(out2, x)
}

# format output
out2$Group.1 <- NULL
out2$Group.2 <- NULL
colnames(out2)[1] <- "Value"


#########################################
## Merge soil and weather covariables ###
#########################################

# merge
out <- bind_rows(out, out2)

# transform to wide format
out <- reshape2::dcast(out, Envir ~ Variable, value.var="Value")

# convert to a numeric matrix
rownames(out) <- out$Envir
out$Envir <- NULL
W <- data.matrix(out)

# save
save(W, file="./Output/W.Rdata")
