library(tidyverse)
library(dbplyr)
library(pracma)
library(tm)

design <- read_csv("inputs.csv")
names(design)
design <- select(design, -lock_1_restrict, - lock_2_release)
design
con <- DBI::dbConnect(RSQLite::SQLite(), "summaries_DW.db")
compact <- tbl(con, "compact")
compact
HDeaths <- dplyr::select(compact, age, week, trustId, HD, output, replicate) %>%
  collect()
HDeaths
TrustLookup <- read_csv("trust19Lookup.csv")
TrustLookup <- select(TrustLookup, trustId)
TrustLookup
WD11ToAcuteTrustIncWalesHB <- read_csv("~/Dropbox/BayesExeter/NewEnsemble/WD11ToAcuteTrustIncWalesHB.csv")
WD11ToAcuteTrustIncWalesHB <- select(WD11ToAcuteTrustIncWalesHB, trustId, trustName)
WD11ToAcuteTrustIncWalesHB
Trusts19wNames <- inner_join(TrustLookup, WD11ToAcuteTrustIncWalesHB, by= "trustId") %>%
  distinct()
Trusts19wNames
TrustNames <- tools::toTitleCase(tolower(Trusts19wNames$trustName))
removewords1 <- c("Nhs", "Trust", "Foundation", "University Hospitals of ", "University Hospital of", "and District", "University Local Board")
removewords2 <- c("University Hospitals", "University Hospital", "Teaching Hospitals", "Teaching Hospitals")
removewords3 <- c("Hospitals", "Hospital", "Teaching", "District", "Health", "Healthcare", "Services", "General", "Countess of ", "The", "Group", "Acute")
TrustNames <- strTrim(stripWhitespace(removeWords(TrustNames,removewords1)))
TrustNames <- strTrim(stripWhitespace(removeWords(TrustNames,removewords2)))
TrustNames <- strTrim(stripWhitespace(removeWords(TrustNames,removewords3)))
TrustNames
Trusts19wNames <- mutate(Trusts19wNames, TrustNames) %>%
  select(-trustName)
Trusts19wNames
HDeaths <- inner_join(HDeaths, Trusts19wNames, by="trustId")
HDeaths
#### Note these steps removed about 20% of the data. We must have data at trusts with codes but no name. Resolve with Rob.

####
#Interesting plot
SWtrusts <- TrustNames[c(51,78,81,87,89,90,98,101,103,106)]
SWplot <- ggplot(filter(HDeaths, age==8 &TrustNames%in%SWtrusts)) +
  geom_line(aes(x=week, y= HD, colour=output,lty=as.factor(replicate))) +
  facet_wrap(~TrustNames)+
  theme(legend.position = "none")+
  ylab("Hospital Deaths")#+
#ylim(c(0,100))
SWplot

####
#Join with design
DeathsToEmulate <- inner_join(HDeaths, design, by="output")
Lewisham <- filter(DeathsToEmulate, trustId == "RJ2", age==8, week==12) %>%
  select(HD, R0:eta) %>%
  gather(key, value, R0:eta)
LewishamPlot <- ggplot(Lewisham)+
  geom_point(aes(x=value, y=HD, colour=key))+
  facet_wrap(~key, scales="free_x")
LewishamPlot
LewishamPlot+ylim(c(0,185))

RXR <- filter(DeathsToEmulate, trustId=="RXR")
RXR[which.max(RXR$HD),]$output

Replicates <- filter(DeathsToEmulate, repeats>1,week==12,age==8)
Replicates
View(Replicates)
Replicates[which.max(Replicates$HD),]$output
ggplot(filter(Replicates,trustId=="RXR"),aes(HD, R0))+
  geom_boxplot(aes(group=HD))