#title: "Data processing of SEAMAP dataset"
#author: "Nathan Baker"
#date: "November 10, 2016"
#output: html_document

# Make 2 datasets (one with fish and abundance, the other with environmental factors)
# Find out number of nets per trawl
# Subset data to only include fish
# Run simple analyses of 1990 vs 2015 data
    # SAD distribution

#import data
library(readr)
dat = read_csv('./data/bakernj.Coastal Survey.ABUNDANCEBIOMASS.2016-11-03T11.33.55.zip')
##determining number of unique nets (collection number) for each trawl (event name)
n = NULL
uni_event = unique(dat$EVENTNAME)
for (i in 1:length(uni_event)) {
    dat_sub = subset(dat, EVENTNAME == uni_event[i])
    n[i] = length(unique(dat_sub$COLLECTIONNUMBER))
}
n
sum(n)
table(n)

fish_species = read.csv('./data/fish_species.csv', header=FALSE, colClasses='character')
names(fish_species) = 'species'


#library(rfishbase)
#?validate_names
##ignore fishbase for now
##THIS PART REMOVES THE =\" FROM ENTRY NAMES"
dat = as.data.frame(dat)
for (i in 1:ncol(dat)) {
    dat[ , i] = sub('\"', '', sub('=\"', '', dat[ , i]))
}


uni_sp = unique(dat$SPECIESCOMMONNAME)
# spcies in the sp list
gd_common_names = uni_sp[uni_sp %in% fish_species$species]
# species not in the sp list
uni_sp[!(uni_sp %in% fish_species$species)]

gd_sci_names = unique(dat$SPECIESSCIENTIFICNAME[dat$SPECIESCOMMONNAME %in% gd_common_names])

##manual_sci_names = c(" " , " ", " ")

#write.csv(gd_sci_names, file='./data/gd_sci_names.csv', row.names=F)

#manually enter in small subset of species that lack common names

gd_sci_names = read.csv('./data/gd_sci_names.csv')
names(gd_sci_names) = 'species'

dat_sub = subset(dat, dat$SPECIESSCIENTIFICNAME %in% gd_sci_names$species)
dim(dat_sub)
head(dat_sub)

#output names that we have filtered out of the dataset
uni_sci_sp = unique(dat$SPECIESSCIENTIFICNAME)
##write.csv(uni_sci_sp[!(uni_sci_sp %in% gd_sci_names$species)], './data/sp_names_filtered_out.csv', row.names=F)

# subset fish columns
fish_cols = c('EVENTNAME','COLLECTIONNUMBER','SPECIESSCIENTIFICNAME',
              'SPECIESCOMMONNAME','NUMBERTOTAL','EFFORT','LOCATION',
              'REGION','DEPTHZONE','STATIONCODE')
dat_sub[1:5 , fish_cols]

# run check that only one date applies to each event name

out = NULL
for (i in 1:length(dat$EVENTNAME)){
  uni_dates = unique(dat$DATE[dat$EVENTNAME == dat$EVENTNAME[i]])
  out[i] = length(uni_dates)
}

sum(out!=2)
##there was only one date per event name

# altnerative way to 
#sapply(1:nrow(dat_sub), function(i) 
  #length(unique(dat_sub$DATE[dat_sub$EVENTNAME == dat_sub$EVENTNAME[i]]))) 


# create site x sp matrix (sites as rows and species as columns)
dat_sub$NUMBERTOTAL = as.integer(dat_sub$NUMBERTOTAL)

sitexsp = tapply(dat_sub$NUMBERTOTAL,
                 list(dat_sub$EVENTNAME, 
                      dat_sub$SPECIESSCIENTIFICNAME), 
                 sum)
sitexsp = ifelse(is.na(sitexsp), 0, sitexsp)
##write.csv(sitexsp,file = "./data/sitexsp_eventnames.csv", row.names=TRUE)
sitexsp = read.csv('./data/sitexsp_eventnames.csv')
names(sitexsp) = c('EVENTNAME', names(sitexsp)[-1])

env = dat_sub[match(sitexsp$EVENTNAME, dat_sub$EVENTNAME), 
              c('REGION', 'LONGITUDESTART', 'LATITUDESTART')]
names(env) = c('REGION', 'x', 'y')

install.packages('devtools')
library(devtools)
install_github('MoBiodiv/mobr')
library(mobr)

sitexsp2 = subset(sitexsp, env$REGION %in% c('SOUTH CAROLINA', 'GEORGIA'))
env_scga = subset(env, REGION %in% c('SOUTH CAROLINA', 'GEORGIA'))

mob_in = make_mob_in(sitexsp2, env_scga)
deltaS = get_delta_stats(mob_in, 'REGION', ref_group='SOUTH CAROLINA')




