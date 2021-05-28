# you probably want to change this
setwd("/home/ben/onedrive/Exeter/Research/covid/")

# Regions, e.g. South West, North East, ...
regions <- paste("E1200000", 1:9, sep = "")
# files with deaths for each regions by age bands
fl0 <- "https://api.coronavirus.data.gov.uk/v2/data?areaType=region&areaCode=E12000001&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv"
## put into list data for each region
 dat <- lapply(regions, function(x) read.csv(gsub(regions[1], x, fl0)))
 names(dat) <- regions
# saveRDS(dat, "RegionDeaths0_20210510.rds")
#dat <- readRDS("RegionDeaths0_20210510.rds")

# bind list into single data.frame
dat <- do.call(rbind, dat)

# ages we want (to avoid overlapping age bands)
ages <- c("00_04", "05_09", "10_14", "15_19", "20_24", "25_29", "30_34", 
  "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69", 
  "70_74", "75_79", "80_84", "85_89", "90+")
dat <- subset(dat, age %in% ages)
# dates
dates <- sort(unique(dat$date))

# set as factors, to make life easy for putting into an array
dat$age <- factor(dat$age, levels = ages)
dat$date <- factor(dat$date, levels = dates)
dat$areaCode <- factor(dat$areaCode, levels = regions)

# then put into array in form dates x regions x ages
# (maybe not the most logical order: use aperm to change)
dimnames <- list(date = dates, regin = regions, age = ages)
arr <- array(0L, dim = sapply(dimnames, length), dimnames = dimnames)
ind <- cbind(as.integer(dat$date), as.integer(dat$areaCode), as.integer(dat$age))
arr[ind] <- dat$deaths
# and then save
saveRDS(arr, "RegionDeaths_20210510.rds")

# load in population by region and age data
# from here
# https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland
# and then butchered into a csv
pop0 <- read.csv("ukmidyearestimates20192020ladcodes_v2.csv")
# get just region data
pop0 <- subset(pop0, Geography == "Region")
# put into age bands above (could be done less manually)
band <- rep(seq_along(ages[-1]), each = 5)
band <- c(band, max(band) + 1)
# aggregate over bands
pop <- as.data.frame(t(apply(pop0[, substr(names(pop0), 1, 1) == "X"], 1, tapply, band, sum)))
names(pop) <- ages
# and then save
saveRDS(pop, "RegionPop_20210510.rds")
