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
                    'historic', ifelse(env$EVENTNAME %in% 2010001:2015657,
                                   'modern', 'other'))

##historic and modern data for 1990-1995 and 2010-2015 ----
fish_historic = subset(sitexsp, !is.na(env$period))
env_historic = subset(env, !is.na(env$period))

historic_mob_in = make_mob_in(fish_historic, env_historic)
historic_mob_stats = get_mob_stats(historic_mob_in, 'period', ref_group = 'historic',
                                   index = c("N","S","S_rare","S_asymp","ENS_PIE"),
                                             n_perm=200)
historic_deltaS = get_delta_stats(historic_mob_in, 'period', ref_group = 'historic',
                                  log_scale = T, nperm = 200)
historic_deltaS

pdf("./figs/mob_stats_boxplots_updates.pdf")
plot(historic_mob_stats, index = c("N", "S", "S_rare", "S_asymp", "ENS_PIE"), 
     ref_group = 'historic', multipanel=T)
dev.off()

##creating rarefaction curves and SAD curve
plot_rarefaction(historic_mob_in, 'period', 'historic', 'indiv', pooled=F, lwd=2,
                 leg_loc='topright')
plot_rarefaction(historic_mob_in, 'period', 'historic', 'indiv', pooled=T, lwd=4,
                 leg_loc='topright')
plot_rarefaction(historic_mob_in, 'period', 'historic', 'spat', 
                 xy_coords = historic_mob_in$spat, lwd = 4, leg_loc = 'topright')
plot_abu(historic_mob_in, 'period', 'historic', type='rad', pooled=F, log='x',
         leg_loc = 'topright')
plot_abu(historic_mob_in, 'period', 'historic', type='rad', pooled=T, log='x')


##run 200 perms
pdf("./figs/deltaS_results.pdf")
plot(historic_deltaS, 'modern', 'historic', same_scale = T)
overlap_effects(historic_deltaS, 'modern', display = "raw")
overlap_effects(historic_deltaS, 'modern', display = "stacked", prop = T)
dev.off()
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

##SR map
env_historic$S = rowSums(fish_historic > 0)

S_breaks = hist(env_historic$S, plot=FALSE)$breaks
S_bins = as.integer(cut(env_historic$S, S_breaks))

pdf('./figs/sr_map.pdf')
par(mfrow=c(1,2))
map(database = "county", regions = counties)
points(env_historic$x[env_historic$period == 'historic'],
       env_historic$y[env_historic$period == 'historic'],
       col=terrain.colors(15)[S_bins[env_historic$period == 'historic']],
       pch=19, cex=.5)
map(database = "county", regions = counties)
points(env_historic$x[env_historic$period == 'modern'],
       env_historic$y[env_historic$period == 'modern'],
       col=terrain.colors(15)[S_bins[env_historic$period == 'modern']],
       pch=19, cex=.5)
dev.off()

##N map
env_historic$N = rowSums(fish_historic)
N_breaks = hist(env_historic$N, plot=FALSE)$breaks
N_bins = as.integer(cut(env_historic$N, N_breaks))


pdf('./figs/N_map.pdf')
par(mfrow=c(1,2))
map(database = "county", regions = counties)
points(env_historic$x[env_historic$period == 'historic'],
       env_historic$y[env_historic$period == 'historic'],
       col=terrain.colors(30)[N_bins[env_historic$period == 'historic']],
       pch=19, cex=.5)
map(database = "county", regions = counties)
points(env_historic$x[env_historic$period == 'modern'],
       env_historic$y[env_historic$period == 'modern'],
       col=terrain.colors(30)[N_bins[env_historic$period == 'modern']],
       pch=19, cex=.5)
dev.off()


##tests?
S_modern = subset(fish_historic, env$period == "modern")
lat_modern = subset(env$x, env$period == "modern")
long_modern = subset(env$y, env$period == "modern")
h_modern = hist(S_modern)
S_binsmod = as.integer(cut(S_modern, h_modern$breaks))


map(database = "county", regions = counties)
points(long, lat, col=terrain.colors(10)[S_bins])



##analysis of log(S)~year
env$S = rowSums(sitexsp > 0)
pdf('./figs/SR~year.pdf')
plot((log(env$S))~jitter(env$YEAR))
lines(lowess(as.numeric(env$YEAR), log(env$S)), col='blue', lwd=5)
mod = lm(log(env$S) ~ as.numeric(env$YEAR))
abline(mod, col='red', lwd = 4)
dev.off()

##analysis of log(S)~trawl/eventname
plot((log(env$S))~env$EVENTNAME)
lines(lowess(as.numeric(env$EVENTNAME), log(env$S)), col='blue', lwd=2)
mod2 = lm(log(env$S) ~ as.numeric(env$EVENTNAME))
abline(mod2, col='red')

##analysis of log(N)~year
env$N = rowSums(sitexsp)
pdf('./figs/N~year.pdf')
plot((log(env$N))~jitter(env$YEAR))
lines(lowess(as.numeric(env$YEAR), log(env$N)), col='blue',lwd = 5)
mod3 = lm(log(env$N) ~ as.numeric(env$YEAR))
abline(mod3, col = 'red', lwd = 4)
dev.off()


##adjusting fishbase compatibility
library("rfishbase")
library("lettercase")
splist = names(sitexsp)
name_split = strsplit(tolower(splist), '.', fixed = T)
genera = sapply(name_split, function(x) str_ucfirst(x[[1]]))
sp = sapply(name_split, function(x) x[2])
names(sitexsp) = paste(genera, sp)
##write.csv(sitexsp, file = "./data/sitexsp_fishbase.csv", row.names=TRUE)
sitexsp_fishbase = read.csv("./data/sitexsp_fishbase.csv", row.names = 1)


##seeing what fish are different between the two communities
sp_historic = colnames(subset(sitexsp, env$period == "historic"))
counts_historic = colSums(subset(sitexsp, env$period == "historic"))
historic_counts = data.frame(sp_historic, counts_historic)
historic_counts$rank = rank(counts_historic, na.last = T)
historic_counts

sp_modern = colnames(subset(sitexsp, env$period == "modern"))
counts_modern = colSums(subset(sitexsp, env$period == "modern"))
modern_counts = data.frame(sp_modern, counts_modern)
modern_counts$rank = rank(counts_modern, na.last = T)
modern_counts


fish_rank = matrix(c(historic_counts$rank, modern_counts$rank), nrow = 2,
                   ncol = 208, byrow = T, dimnames = list(c("historic", "modern"),
                                                          sp_modern))
rank_plot = barplot(fish_rank, beside = T, 
                    col = c("lightblue", "lightcyan"),
                    main = "Fish Abundance Ranks in Southeast Atlantic")
modern_absent = subset(modern_counts, modern_counts$counts_modern == "0")
historic_absent = subset(historic_counts, historic_counts$counts_historic == "0") 
diff = rownames(modern_absent)[!(rownames(modern_absent) %in% rownames(historic_absent))]
diff

qts_historic = ecdf(counts_historic)(counts_historic)
qts_modern = ecdf(counts_modern)(counts_modern)
plot(0, 0, ylim=c(0,1), xlim=c(1,2), frame.plot=F, type='n', axes =F, xlab='', ylab='quantile')
axis(side=2)
for(i in seq_along(qts_historic))
points(c(1,2), c(qts_historic[i], qts_modern[i]), type='l')


##Regional Analysis
regfish_modern = subset(sitexsp, env$period == 'modern')
regenv_modern = subset(env, env$period == 'modern')
regmodern_mob_in = make_mob_in(regfish_modern, regenv_modern)
regmodern_mob_stats = get_mob_stats(regmodern_mob_in, 'REGION',
                                    ref_group = 'SOUTH CAROLINA',
                                   index = c("N","S","S_rare","S_asymp","ENS_PIE"),
                                   n_perm=200)
plot(regmodern_mob_stats, index = c("N", "S", "S_rare", "S_asymp", "ENS_PIE"), 
     ref_group = 'SOUTH CAROLINA', multipanel=T)
regmodern_deltaS = get_delta_stats(regmodern_mob_in, 'REGION', ref_group = 'SOUTH CAROLINA',
                                  log_scale = T, n_perm = 10)
regmodern_deltaS
plot_rarefaction(regmodern_mob_in, 'REGION', 'SOUTH CAROLINA', 'indiv', pooled=F, lwd=2,
                 leg_loc='topright')
plot_rarefaction(regmodern_mob_in, 'REGION', 'SOUTH CAROLINA', 'indiv', pooled=T, lwd=4,
                 leg_loc='topright')
plot_rarefaction(regmodern_mob_in, 'REGION', 'SOUTH CAROLINA', 'spat', 
                 xy_coords = regmodern_mob_in$spat, lwd = 4, leg_loc = 'topright')
plot_abu(regmodern_mob_in, 'REGION', 'SOUTH CAROLINA', type ='rad', pooled=F, log='x',
         leg_loc = 'topright')
plot_abu(regmodern_mob_in, 'REGION', 'SOUTH CAROLINA', type='rad', pooled=T, log='x')


regfish_historic = subset(sitexsp, env$period == 'historic')
regenv_historic = subset(env, env$period == 'historic')
reghistoric_mob_in = make_mob_in(regfish_historic, regenv_historic)
reghistoric_mob_stats = get_mob_stats(regmodern_mob_in, 'REGION',
                                       ref_group = 'SOUTH CAROLINA',
                                       index = c("N","S","S_rare","S_asymp","ENS_PIE"),
                                       n_perm=200)
plot(reghistoric_mob_stats, index = c("N", "S", "S_rare", "S_asymp", "ENS_PIE"), 
     ref_group = 'SOUTH CAROLINA', multipanel=T)
reghistoric_deltaS = get_delta_stats(reghistoric_mob_in, 'REGION', ref_group = 'SOUTH CAROLINA',
                                     log_scale = T, n_perm = 10)
plot_rarefaction(reghistoric_mob_in, 'REGION', 'SOUTH CAROLINA', 'indiv', pooled=F, lwd=2,
                 leg_loc='topright')
plot_rarefaction(reghistoric_mob_in, 'REGION', 'SOUTH CAROLINA', 'indiv', pooled=T, lwd=4,
                 leg_loc='topright')
plot_rarefaction(reghistoric_mob_in, 'REGION', 'SOUTH CAROLINA', 'spat', 
                 xy_coords = regmodern_mob_in$spat, lwd = 4, leg_loc = 'topright')
plot_abu(reghistoric_mob_in, 'REGION', 'SOUTH CAROLINA', type ='rad', pooled=F, log='x',
         leg_loc = 'topright')
plot_abu(reghistoric_mob_in, 'REGION', 'SOUTH CAROLINA', type='rad', pooled=T, log='x')
