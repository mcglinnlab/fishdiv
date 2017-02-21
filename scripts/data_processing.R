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
              'REGION','DEPTHZONE','STATIONCODE', )
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
sitexsp = read.csv('./data/sitexsp_eventnames.csv', row.names = 1)

env = dat_sub[match(row.names(sitexsp), dat_sub$EVENTNAME), 
              c('EVENTNAME', 'REGION', 'LONGITUDESTART', 'LATITUDESTART')]
names(env) = c('EVENTNAME', 'REGION', 'x', 'y')
#write.csv(env, file="./data/env.csv", row.names=FALSE)
env = read.csv('./data/env.csv')

library(vegan)
fish_ca = cca(sitexsp)
##integer overflow
fish_dca = decorana(sitexsp)
summary(fish_dca)

fish_rda = rda(sitexsp ~ env$REGION + env$x+ env$y)
fish_rda
plot(fish_rda, type = 'n', scaling = 1)
orditorp(fish_rda, display = 'species', col = 'blue')

# test significance of overall model
anova(fish_rda)
##Permutation test for rda under reduced model
##Permutation: free
##Number of permutations: 999

##Model: rda(formula = sitexsp ~ env$REGION + env$x + env$y)
##Df Variance      F Pr(>F)    
##Model    2456 24603869 5.4236  0.001 ***
##Residual 2542  4695327                  
##Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# tests partial effects of each variable in model
anova(fish_rda, by='margin')

summary(env)

fish1990 = subset(sitexsp, env$EVENTNAME %in% 1990001:1990547)
env1990 = subset(env, env$EVENTNAME %in% 1990001:1990547)
fish2015 = subset(sitexsp[ EVENTNAME = 2015001:2015657])



library(devtools)
install_github('MoBiodiv/mobr')
library(mobr)

sitexscga = subset(fish1990, env1990$REGION %in% c('SOUTH CAROLINA', 'GEORGIA'))
env_scga = subset(env1990, REGION %in% c('SOUTH CAROLINA', 'GEORGIA'))

mob_in = make_mob_in(sitexscga, env_scga)

mob_stats = get_mob_stats(mob_in, 'REGION')

deltaS = get_delta_stats(mob_in, 'REGION', ref_group='SOUTH CAROLINA', log_scale = T,
                         nperm=5)
plot(deltaS, 'GEORGIA', 'SOUTH CAROLINA', same_scale=T)

# add Time Period variable to env

env$period = ifelse(env$EVENTNAME %in% 1990001:1995563,
                    'past', ifelse(env$EVENTNAME %in% 2010001:2015657,
                                   'modern', NA))

##historic and modern data for 1990-1995 and 2010-2015 ----
fish_historic = subset(sitexsp, !is.na(env$period))
env_historic = subset(env, !is.na(env$period))

historic_mob_in = make_mob_in(fish_historic, env_historic)
historic_mob_stats = get_mob_stats(historic_mob_in, 'period', nperm=200)
plot(historic_mob_stats, multipanel=T)
historic_deltaS = get_delta_stats(historic_mob_in, 'period', ref_group = 'past',
                                  log_scale = T, nperm = 10)
plot(historic_deltaS, 'modern', 'past', same_scale = T)



##historic analysis for all 5 regions -----
historic_mob_in = make_mob_in(fish_historic, env_historic)
historic_mob_stats = get_mob_stats(historic_mob_in, 'REGION')
historic_deltaS = get_delta_stats(historic_mob_in, 'REGION', ref_group = 'SOUTH CAROLINA', log_scale = T,
                                  nperm = 10)
plot(historic_deltaS, 'GEORGIA', 'SOUTH CAROLINA', same_scale = T)

##historic analysis for FL and SC
hist_sitexscfl = subset(fish_historic, env_historic$REGION %in% c('SOUTH CAROLINA', 'FLORIDA'))
hist_envscfl = subset(env_historic, REGION %in% c('SOUTH CAROLINA', 'FLORIDA'))
hist_mob_inscfl = make_mob_in(hist_sitexscfl, hist_envscfl)
mob_stats = get_mob_stats(hist_mob_inscfl, 'REGION')
hist_deltaS = get_delta_stats(hist_mob_inscfl, 'REGION', ref_group = 'SOUTH CAROLINA', log_scale = T,
                               nperm = 50)
plot(deltaS, 'FLORIDA', 'SOUTH CAROLINA', same_scale = T)

##modern analysis for all 5 regions
modern_mob_in = make_mob_in(fish_modern, env_modern)
modern_mob_stats = get_mob_stats(modern_mob_in, 'REGION')
modern_deltaS = get_delta_stats(modern_mob_in, 'REGION', ref_group = 'SOUTH CAROLINA', log_scale = T,
                                nperm = 10)
plot(modern_deltaS, 'FLORIDA', 'SOUTH CAROLINA', same_scale = T)

##modern analysis for FL and SC
mod_sitexscfl = subset(fish_modern, env_modern$REGION %in% c('SOUTH CAROLINA', 'FLORIDA'))
mod_envscfl = subset(env_modern, REGION %in% c('SOUTH CAROLINA', 'FLORIDA'))
mod_mob_inscfl = make_mob_in(mod_sitexscfl, mod_envscfl)
mob_stats = get_mob_stats(mod_mob_inscfl, 'REGION')
mod_deltaS = get_delta_stats(mod_mob_inscfl, 'REGION', ref_group = 'SOUTH CAROLINA', log_scale = T,
                              nperm = 50)
plot(deltaS, 'FLORIDA', 'SOUTH CAROLINA', same_scale = T)
