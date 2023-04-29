library(here)
library(rnaturalearth)
library(dplyr)
library(sf)
library(spdep)
library(rcarbon)

# Read Data ----
d  <- read.csv(here('data','site_raw.csv'))

# New Column Names ----
colnames(d)  <- c('ID','Period','Chronology','Phase','SiteName','Address','SiteType','Latitude','Longitude','Elevation')

# Define Chronology  ----
# Remove cases without chronology:
d  <- subset(d,!Chronology%in%c("","?"))

# Remove cases with only start or end
ages=unique(d$Chronology)
partial = ages[which(unlist(lapply(strsplit(ages,"-"),length))==1)]
tages = ages[which(unlist(lapply(strsplit(ages,"-"),length))==2)]
partial= c(partial,tages[which(unlist(lapply(strsplit(tages,"-"),function(x){x[1]==""})))])
d  <- subset(d,!Chronology%in%partial)

# Remove cases with characters
ages.tmp  <- ages
ages.tmp=gsub(pattern='B',replacement='',ages.tmp)
ages.tmp=gsub(pattern='Ｂ',replacement='',ages.tmp)
ages.tmp=gsub(pattern='１',replacement='1',ages.tmp)
ages.tmp=gsub(pattern='\\?',replacement='',ages.tmp)
ages.tmp=gsub(pattern='or',replacement=',',ages.tmp)
ages.tmp=gsub(pattern=' ',replacement='',ages.tmp)
ages.tmp=gsub(pattern='　' ,replacement='',ages.tmp)
ages.tmp=gsub(pattern='-',replacement='',ages.tmp)
ages.tmp=gsub(pattern=',',replacement='',ages.tmp)

texts  <- ages[grep(pattern="[^0-9.-]", ages.tmp)]
d  <- subset(d,!Chronology%in%texts)

# Manual Fixes:
ages=unique(d$Chronology)
problems = c("B2500-B2000,B2500-B20002,B2000-B1500",
	     "14000-500",
	     "14000-B450,B50-1",
	     "9500-B1250",
	     "9500-B2000,B1800-B1250",
	     "2500-B2000",
	     "B2000-B1500,2500-B2000,B2000-B1500",
	     "B100-B180?",
	     "240-700,200-240,B380-B500",
	     "100-50",
	     "B100-B１50",
	     "B9500-B2500,1350-700",
	     "650-1")


fixes = c("B2500-B2000,B2500-B2000,B2000-B1500",
	     "B14000-500",
	     "B14000-B450,B50-1",
	     "B9500-B1250",
	     "B9500-B2000,B1800-B1250",
	     "B2500-B2000",
	     "B2000-B1500,B2500-B2000,B2000-B1500",
	     "100-180?",
	     "240-700,200-240,380-500",
	     "B100-B50",
	     "B100-B150",
	     "B9500-B2500,B1350-B700",
	     "B650-1")


manualfix  <- data.frame(problems,fixes)

for (i in 1:nrow(manualfix))
{
	ages[which(ages==manualfix$problems[i])] = manualfix$fixes[i]
	d$Chronology[which(d$Chronology==manualfix$problems[i])] = manualfix$fixes[i]
}

# Extract Dates:
chronos  <- data.frame(entry=ages)
chronos$start  <- chronos$end  <- chronos$multiple  <- chronos$gap  <- chronos$notes  <- NA
chronos$exclude  <- FALSE


for (i in 1:length(ages))
{
	tmp = ages[i]
	chronos$multiple[i]  <- any(grepl(",",tmp))
	starts <- numeric()
	ends  <- numeric()

	if (chronos$multiple[i])
	{
		tmp  <- strsplit(tmp,",") |> unlist()
	}

	for (j in 1:length(tmp))
	{
		notes  <- ""
		tmp2  <- tmp[j]
		tmp2=unlist(strsplit(tmp2,"\\, |\\,|\\-|\\or"))
		tmp2=gsub(pattern='B',replacement='-',tmp2)
		tmp2=gsub(pattern='Ｂ',replacement='-',tmp2)
		tmp2=gsub(pattern='１',replacement='1',tmp2)
		tmp2=gsub(pattern='\\?',replacement='',tmp2)
		tmp2=gsub(pattern='or',replacement=',',tmp2)
		tmp2=gsub(pattern=' ',replacement='',tmp2)
		tmp2=gsub(pattern='　' ,replacement='',tmp2)
		tmp2=gsub(pattern='\\?',replacement='',tmp2)
		textcheck=!grepl(pattern="[^0-9.-]", tmp2)
		if (length(tmp2)==1){notes  <- paste0(notes,"only one age")}
		if (length(tmp2)==2 & tmp2[1]!="" & tmp2[2]!="" & as.numeric(tmp2)[1]>as.numeric(tmp2)[2]){notes  <- paste0(notes,"age discrepancy")}
		if (notes=="")
		{
			if(length(tmp2>1)&all(tmp2!='')&length(tmp2)%%2 == 0 & all(textcheck)){
				starts = c(starts,as.numeric(tmp2)[1])
				ends = c(ends,as.numeric(tmp2)[2]) 
			}
		} else if (notes!="")
		{
			chronos$exclude[i] = TRUE
			chronos$notes[i] = notes

		}
	}

	if (!chronos$exclude[i])
	{
	starts  <- BCADtoBP(starts)
	ends  <- BCADtoBP(ends)
	chronos$start[i]  <- max(starts)
	chronos$end[i]  <- min(ends)
	if (chronos$multiple[i])
	{
		dfse  <- data.frame(starts=starts,ends=ends)
		mat  <- matrix(NA,nrow=length(starts),ncol=length(starts))
		for (k1 in 1:length(starts))
		{
			for (k2 in 1:length(starts))
			{
				mat[k1,k2] = (dfse$ends[k1] >= dfse$starts[k2] & dfse$starts[k1]<=dfse$ends[k2])|(dfse$starts[k1]<=dfse$ends[k2] & dfse$starts[k1] >= dfse$starts[k2])
			}
			mat[k1,k1]=NA
		}
		if(any(apply(mat,1,sum,na.rm=T)==0)) {chronos$gap[i] = TRUE}
	}
	}
}

chronos  <- subset(chronos,exclude==FALSE) |> unique()
chronos$s  <- abs(chronos$start - chronos$end)/2
chronos$m  <- chronos$start - chronos$s

# Several Entries hav e a chronology dated to 100~100 - this should be assigned to a small margin of error (15 years)
i  <- which(chronos$s==0)
chronos$s[i] = 15
# Merge back to original data
d <- left_join(d,select(chronos,entry,m,s,start,end),by=c("Chronology"="entry"))

# Identify Prefectures and Neighbourhood for ICAR----
win <- ne_states(country = "japan",returnclass = 'sf') |> subset(!name_vi %in% c("Okinawa","Hokkaido"))
Npref <- nrow(win)
win$ID  <-  1:Npref
W_nb <- poly2nb(win, row.names =  win$ID)
nbInfo <- nb2WB(W_nb)


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
d <- subset(d,!is.na(Prefecture)&!Prefecture%in%c("Okinawa","Hokkaido"))

#Eliminate sites with problematic chronologies
d  <- subset(d,!is.na(m))

# Narrow to Yayoi period 
foo  <- function(x1,y1,x2,y2)
{
	d  <- x1 - y1
	p  <- 1/d
	if (x1<=x2 & y1>=y2) {return(1)}
	if (x1>x2 & y1>x2) {return(0)}
	if (x1<y2 & y1<y2) {return(0)}
	if (x1>=x2 & y1>y2) {return(abs(x2-y1)*p)}
	if (x1<x2 & y1<y2) {return(abs(x1-y2)*p)}
	if (x1>=x2 & y1<=y2) {return(abs(x2-y2)*p)}
}

# Compute probabilities within window
d$prob.within  <- unlist(apply(cbind(d$start,d$end,2900,1700),1,function(x){foo(x[1],x[2],x[3],x[4])}))
d <- subset(d,prob.within==1)

# Save Output
sitedb  <- d
save(win,nbInfo,sitedb,file=here('data','sitedata.RData'))

