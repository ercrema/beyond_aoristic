library(here)
library(rnaturalearth)
library(dplyr)
library(sf)

# Read Data ----
d  <- read.csv(here('data','site_raw.csv'))

# New Column Names ----
colnames(d)  <- c('ID','Period','Chronology','Phase','SiteName','Address','SiteType','Latitude','Longitude','Elevation')

# Identify Prefectures ----
prefTrans  <- read.csv(here('data','prefectures_translations.csv'))
d$Prefecture  <- NA

for(i in 1:nrow(prefTrans))
{
	d$Prefecture[grep(prefTrans$JpNames[i],d$Address)]  <- prefTrans$Translation[i]
}

# Eliminate cases without Address and Lat/Lon
d  <- d[- which(is.na(d$Address)|d$Address == '_'|is.na(d$Latitude)|is.na(d$Longitude))]

# Estimate prefecture of remaining sites based on Lat/Lon
ii  <- which(is.na(d$Prefecture) & !is.na(d$Latitude) & !is.na(d$Longitude)) 
check.sites  <- d[ii,] |> st_as_sf(coords=c('Longitude','Latitude'),crs=4326)
jp  <- ne_states(country = 'japan',returnclass = 'sf')
d$Prefecture[ii]  <- jp$name[as.numeric(st_intersects(check.sites,jp))]

# Eliminate all remaining cases
d <- subset(d,!is.na(Prefecture))

# Assign Region
d  <- left_join(d, prefTrans, by=c('Prefecture'='Translation'))


# Site Name Extraction ----

# Assign unique SiteID ----
#1. Similar site name, same prefecture, and distance within X km


# Handle Chronology ----
# 1. Remove all cases without chronology
# 2. Separate entries per SiteID
# 3. Record in dunif2 format


		 
