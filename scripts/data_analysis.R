sitexsp = read.csv('./data/sitexsp_eventnames.csv', row.names = 1)
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

##this was a test

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
                                  log_scale = T, inds = 10, nperm = 1)
##run 200 perms
plot(historic_deltaS, 'modern', 'past', same_scale = T)
stack_effects(historic_deltaS, 'modern')
stack_effects(historic_deltaS, 'modern', prop = T)
##start looking at compositional effects (constrained ordination)
?"mobr"


##historic analysis for all 5 regions -----
##historic_mob_in = make_mob_in(fish_historic, env_historic)
##historic_mob_stats = get_mob_stats(historic_mob_in, 'REGION')
##historic_deltaS = get_delta_stats(historic_mob_in, 'REGION', ref_group = 'SOUTH CAROLINA', log_scale = T,
##nperm = 10)
##plot(historic_deltaS, 'GEORGIA', 'SOUTH CAROLINA', same_scale = T)

##historic analysis for FL and SC
##hist_sitexscfl = subset(fish_historic, env_historic$REGION %in% c('SOUTH CAROLINA', 'FLORIDA'))
##hist_envscfl = subset(env_historic, REGION %in% c('SOUTH CAROLINA', 'FLORIDA'))
##hist_mob_inscfl = make_mob_in(hist_sitexscfl, hist_envscfl)
##mob_stats = get_mob_stats(hist_mob_inscfl, 'REGION')
##hist_deltaS = get_delta_stats(hist_mob_inscfl, 'REGION', ref_group = 'SOUTH CAROLINA', log_scale = T,
##nperm = 50)
##plot(deltaS, 'FLORIDA', 'SOUTH CAROLINA', same_scale = T)

##modern analysis for all 5 regions
##modern_mob_in = make_mob_in(fish_modern, env_modern)
##modern_mob_stats = get_mob_stats(modern_mob_in, 'REGION')
##modern_deltaS = get_delta_stats(modern_mob_in, 'REGION', ref_group = 'SOUTH CAROLINA', log_scale = T,
##nperm = 10)
##plot(modern_deltaS, 'FLORIDA', 'SOUTH CAROLINA', same_scale = T)

##modern analysis for FL and SC
##mod_sitexscfl = subset(fish_modern, env_modern$REGION %in% c('SOUTH CAROLINA', 'FLORIDA'))
##mod_envscfl = subset(env_modern, REGION %in% c('SOUTH CAROLINA', 'FLORIDA'))
##mod_mob_inscfl = make_mob_in(mod_sitexscfl, mod_envscfl)
##mob_stats = get_mob_stats(mod_mob_inscfl, 'REGION')
##mod_deltaS = get_delta_stats(mod_mob_inscfl, 'REGION', ref_group = 'SOUTH CAROLINA', log_scale = T,
##nperm = 50)
##plot(deltaS, 'FLORIDA', 'SOUTH CAROLINA', same_scale = T)

library(maps)
library(sp)
library(maptools)
library(rgdal)
library(lattice)
library(classInt)

counties =  c("north carolina,dare", "north carolina,pasquotank", "north
              carolina,hyde", "north carolina,pamlico", "north carolina,beaufort", "north
              carolina,onslow", "north carolina,carteret", "north carolina,pender", "north
              carolina,new hanover", "north carolina,brunswick", "north carolina,craven",
              "north carolina,perquimans", "north carolina,camden", "north
              carolina,currituck:knotts", "north carolina,currituck:knotts", "north
              carolina,currituck:spit", "north carolina,tyrell", "north
              carolina,washington", "north carolina,bertie", "north carolina,hertford",
              "north carolina,chowan", "north carolina,jones", "north carolina,columbus",
              "south carolina,horry", "south carolina,marion", "south carolina,georgetown",
              "south carolina,charleston", "south carolina,berkeley", "south
              carolina,dorchester", "south carolina,colleton", "south carolina,beaufort",
              "south carolina,jasper", "south carolina,hampton", "georgia,effingham",
              "georgia,chatham", "georgia,bryan", "georgia,liberty", "georgia,mcintosh",
              "georgia,long", "georgia,glynn", "georgia,wayne", "georgia,camden",
              "georgia,brantley", "georgia,chariton", "florida,nassau", "florida,duval",
              "florida,saint johns", "florida,clay", "florida,putnam", "florida,flagler",
              "florida,volusia", "florida,brevard", "florida,indian river",
              "florida,orange", "florida,osceola", "florida,seminole")

southeast = map(database = "state", 
                regions = c("north carolina", "south carolina", "georgia", "florida"))
se_coast = map(database = "county", regions = counties, fill = T, plot = F)
secoast_sp = map2SpatialPolygons(se_coast, se_coast$names, CRS("+proj=longlat"))
plot(secoast_sp, axes = T)
fish_pastpop = subset(sitexsp %in% env$period == "past")
se_spdf = SpatialPolygonsDataFrame(secoast_sp,fish_historic)
spplot(se_spdf, '')


?vegan::diversity



env_historic$year = as.numeric(substr(row.names(fish_historic), 1, 4))

# ignore aggregation until site id's have been rectified in a meaningful way
tst = aggregate(fish_historic, list(env_historic$period), FUN = 'sum')

env_historic$S = rowSums(fish_historic > 0)

S_breaks = hist(env_historic$S, plot=FALSE)$breaks
S_bins = as.integer(cut(env_historic$S, S_breaks))

pdf('sr_map.pdf')
par(mfrow=c(1,2))
map(database = "county", regions = counties)
points(env_historic$x[env_historic$period == 'past'],
       env_historic$y[env_historic$period == 'past'],
       col=terrain.colors(15)[S_bins[env_historic$period == 'past']],
       pch=19, cex=.5)
map(database = "county", regions = counties)
points(env_historic$x[env_historic$period == 'modern'],
       env_historic$y[env_historic$period == 'modern'],
       col=terrain.colors(15)[S_bins[env_historic$period == 'modern']],
       pch=19, cex=.5)
dev.off()

S_modern = subset(fish_historic, env$period == "modern")
lat_modern = subset(env$x, env$period == "modern")
long_modern = subset(env$y, env$period == "modern")
h_modern = hist(S_modern)
S_binsmod = as.integer(cut(S_modern, h_modern$breaks))


map(database = "county", regions = counties)
points(long, lat, col=terrain.colors(10)[S_bins])