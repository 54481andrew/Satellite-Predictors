## First generate the datasets
## Load packages for graphing
library(ggplot2)
library(reshape)
library(reshape2)
library(plotrix)
library(plotmo)
library(parallel)
library(plyr)
library(geosphere)
library(raster)
library(lubridate) ## Date handling
library(zoo) ## Date handline
library(gridExtra)
library(ggcorrplot) # Correlation matrix plot
library(MASS) # stepAIC
library(dismo) # Boosted classification trees
library(lme4) # GLMM (binomial distributed random effects models)
library(glmnet) # Elastic net regression

## My helper functions
source('Tools/Functions.r')

## - Specify paths to EFC and PREEMPT datasets
efc.data.path = '../Trap_Data/Data/Clean_EFC_Data_By_Trap.csv'
preempt.data.path = '../Trap_Data/Data/Clean_PRE_Data_By_Trap.csv'

## --- Load PREEMPT data with satellite predictors
pre.dataset <- read.csv(preempt.data.path)

## --- Load EFC data; keep only those sites with transect-level data
efc.dataset <- read.csv(efc.data.path)
efc.dataset <- subset(efc.dataset, Site %in% c('Bantou', 'Tanganya', 'Gania'))

## --- Join EFC and PREEMPT data
rodent.data <- rbind.fill(efc.dataset, pre.dataset)

## --- Choose any predictor names to omit
omit.names <- c('MODIS')

rodent.data <- rodent.data[,!(names(rodent.data) %in% omit.names)]

## Note that we omit cloud here
sat.vars <- c('Frac_bare', 'Frac_grass', 'Frac_tree', 
                     'Frac_burn', 'Frac_rice', 'Frac_water', 
                     'Frac_mound', 
                     'Density_Buildings', 'Density_Moderns',
                     'Density_Traditionals')
full.sat.vars <- get.features(rvec = c(25, 50,100,200,1000), names = sat.vars)$names
modis.vars <- c('Elev', 'Tmu', 'Pmu', 'Nmu','Pcv','Ncv','Pc','Pm',
                'Nc','Nm','Pmin','Pmax','Nmin','Nmax','Pdur',
                'Ndur', 'Pop')
lc.vars <- c('Water_Bodies','Evergreen_Needleleaf_Forest',
             'Evergreen_Broadleaf_Forest','Deciduous_Needleleaf_Forest',
             'Deciduous_Broadleaf_Forest','Mixed_Forest','Closed_Shrubland','Open_Shrubland',
             'Woody_Savanna','Savannas','Grasslands','Permanent_Wetlands','Croplands','Urban_BuiltUp',
             'Cropland_Natural_Mosaic', 'nBuilds','nModBuild', 'nHuts', 'fHuts')
plag.vars <- paste('P', 1:12, sep = '')

focal.variables <- c(full.sat.vars,
                     modis.vars, lc.vars, plag.vars)

scaled.df <- rodent.data
for(r in c(25, 50,100,200,1000)){
    select.variables <- get.features(rvec = r, names = sat.vars)$names
    norm.out <- norm.preds(df = rodent.data, foc.vars = select.variables)
    scaled.df[,select.variables] = norm.out[[1]]
}
nonsat.vars <- c(modis.vars, lc.vars, plag.vars)
norm.out <- norm.preds(df = rodent.data, foc.vars = nonsat.vars)
scaled.df[,nonsat.vars] <- norm.out[[1]]

## Initialize new dataframe, scaled.df, with norm'd data and Site
scaled.data <- scaled.df
scaled.data$Mna <- factor(1*(scaled.data$Mna), levels = c('0','1'), ordered = TRUE)
scaled.data$ID <- 1:nrow(scaled.data)

##
out = data.table::setDT(scaled.data)[,.(TotCapture = sum(Trap.weight[Mna==1]),
                                        TotTraps = sum(Trap.weight)),
                   by = list(Site)]
out$TS <- with(out, TotCapture / TotTraps)

site.order = out[order(out$TS, decreasing = TRUE),c('Site', 'TotCapture', 'TS')]

site.order[,.(Site, TotCapture, TS = round(TS, 2))]

## By Mn Captures:
##          Site TotCapture   TS
##  1:    Bantou        321 0.05
##  2:  Tanganya        238 0.04
##  3:   Bafodia        230 0.08
##  4:     Largo         15 0.01
##  5:   Benduma         11 0.01
##  6:    Barlie          9 0.01
##  7:    Talama          8 0.01
##  8:   Kapethe          7 0.01
##  9:  Naiawama          7 0.03
## 10: Gbenikoro          6 0.03
## 11:     Guala          5 0.00
## 12:    Makump          5 0.01
## 13:    Badala          3 0.00
## 14:  Njaguima          2 0.00
## 15:   Yekeyor          2 0.00
## 16:     Gania          1 0.00
## 17:    Makuna          1 0.00
## 18:   Mokorie          1 0.00
## 19: Gbainkfay          0 0.00
## 20:    Petema          0 0.00

##          Site TotCapture   TS
##  1:   Bafodia        230 0.08
##  2:    Bantou        321 0.05
##  3:  Tanganya        238 0.04
##  4:  Naiawama          7 0.03
##  5: Gbenikoro          6 0.03
##  6:    Talama          8 0.01
##  7:   Benduma         11 0.01
##  8:    Barlie          9 0.01
##  9:     Largo         15 0.01
## 10:   Kapethe          7 0.01
## 11:    Makump          5 0.01
## 12:  Njaguima          2 0.00
## 13:   Mokorie          1 0.00
## 14:   Yekeyor          2 0.00
## 15:     Guala          5 0.00
## 16:    Badala          3 0.00
## 17:    Makuna          1 0.00
## 18:     Gania          1 0.00
## 19: Gbainkfay          0 0.00
## 20:    Petema          0 0.00


round(site.order$TotCapture / sum(site.order$TotCapture),2)

## It is clear that the captures are biased toward Bantou (37%), Tanganya(27%), and Bafodia(26%). 
## 90% of the presences comes from these 3 towns.


## Simple barplot of number of rodents captured
site.order$Site <- factor(site.order$Site, levels = unlist(site.order$Site, use.names = FALSE),
                    ordered = TRUE)
ggplot(site.order) + geom_bar(aes(x = Site, weight = TotCapture),
                        position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    scale_fill_manual("Mn captures",
                      values = c("1" = "blue", "0" = "orange")) + 
ggtitle('Mastomys natalensis captures') + ylab('Count')
ggsave(filename = 'Figures_Both/Barplot_Captures_Mn.png')

## Simple barplot of TS
ggplot(site.order) + geom_bar(aes(x = Site, weight = TS),
                        position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
ggtitle('Mastomys natalensis trap success') + ylab('TS')
ggsave(filename = 'Figures_Both/Barplot_TS_Mn.png')

## Look closer at individual environments
out = data.table::setDT(scaled.data)[,.(TotCapture = sum(Trap.weight[Mna==1]),
                                        TotTraps = sum(Trap.weight)),
                   by = list(Site, Type)]
out$Site <- factor(out$Site, levels = unlist(site.order$Site, use.names = FALSE),
                    ordered = TRUE)
out$TS <- with(out, TotCapture / TotTraps)
ggplot(out) + geom_bar(aes(x = Site, weight = TS, fill = Type),
                        position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
ggtitle('Mastomys natalensis trap success') + ylab('TS')
ggsave(filename = 'Figures_Both/Barplot_Type_TS_Mn.png')

## Naiawama and Talama have quite a few traps in House sites. 


library(tidyr)
temp = gather(scaled.data, variable, value,
           all_of(focal.variables), factor_key = TRUE)
temp = temp[,c('Site', 'Rodent', 'Trap.weight', 'Lassa', 'Date', 'Sp', 'Mna', 'variable', 'value')]
temp$type <- ifelse(temp$variable %in% full.sat.vars, 'sat',
             ifelse(temp$variable %in% modis.vars, 'modis', 'prec'))
f <- function(x){
    as.numeric(strsplit(as.character(x), '\\.')[[1]][[2]])
}
temp$radius <- NA
sat.mask <- temp$type=='sat'
temp$radius[sat.mask] <- sapply(temp$variable[sat.mask], FUN = f)
f <- function(x){
    strsplit(as.character(x), '\\.')[[1]][[1]]
}
temp$base.name <- NA
temp$base.name[sat.mask] <- sapply(temp$variable[sat.mask], FUN = f)

melt.scaled.data <- temp
rm(temp)

## --- TRAP SUCCESS PLOTS

## -Plot Mastomys natalensis only 
sat.dat <- subset(melt.scaled.data, type == 'sat')
sat.dat$base.name <- as.character(sat.dat$base.name)
sat.dat$base.name <- as.factor(sat.dat$base.name)
label.fxn <- as_labeller(function(y){paste0('Focal Radius: ', y)})
ggplot(data = sat.dat,
       aes(x = as.factor(base.name), y = value, fill = Mna, weight = Trap.weight)) +
    geom_boxplot(alpha = 0.65) + xlab('') + ylab('Scaled Value') +
    coord_cartesian(ylim = c(-5,10)) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    facet_wrap(radius~., ncol = 1, labeller = label.fxn,
               strip.position="top")
ggsave(filename = 'Figures_Both/Mn_Pres_Abs_Features.png',
       height = 7, width = 7)

## -Plot Mastomys natalensis in high-trap towns only
dir.create('Figures_Both/Mn_Pres_By_Site')
site.names <- c('Bantou', 'Tanganya', 'Bafodia', 'Naiawama', 'Talama')
##site.names <- unique(sat.dat$Site)
trunc.sat.dat <- sat.dat
val <- trunc.sat.dat$value
trunc.val <- 10
trunc.sat.dat$value[val > trunc.val] <- trunc.val
trunc.sat.dat$value[val < -trunc.val] <- -trunc.val
for(site in site.names){
    label.fxn <- as_labeller(function(y){paste0('Focal Radius: ', y)})
    ggplot(data = subset(trunc.sat.dat, Site==site),
           aes(x = base.name, y = value, fill = Mna, weight = Trap.weight)) +
        geom_boxplot(alpha = 0.65) + xlab('') + ylab('Scaled Value') +
            coord_cartesian(ylim = c(-trunc.val,trunc.val)) + 
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
        facet_wrap(radius~., ncol = 1, labeller = label.fxn)+
        ggtitle(paste0('Trap outcome in ', site)) 
    ggsave(filename = paste0('Figures_Both/Mn_Pres_By_Site/', site, '_Mn_Pres_Abs_Features.png'),
       height = 7, width = 7)
}



label.fxn <- as_labeller(function(y){paste0('Focal Radius: ', y)})
ggplot(data = subset(trunc.sat.dat,
                     Site%in% c('Bantou', 'Tanganya', 'Bafodia', 'Naiawama', 'Talama')),
       aes(x = base.name, y = value, fill = Mna, weight = Trap.weight)) +
    geom_boxplot(alpha = 0.65) + xlab('') + ylab('Scaled Value') +
        coord_cartesian(ylim = c(-trunc.val, trunc.val)) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    facet_wrap(radius~., ncol = 1, labeller = label.fxn)+
    ggtitle(paste0('Trap outcome in high Mn sites')) 
ggsave(filename = paste0('Figures_Both/Mn_Pres_By_Site/High_Mn_Pres_Abs_Features.png'),
       height = 7, width = 7)


dir.create('Figures_Both/Compare_Hist_Vars')
ss.melt <- subset(melt.scaled.data, Site %in% c('Bantou', 'Tanganya', 'Bafodia') &
                  radius %in% c(25, 50, 100,1000))

all.melt <- subset(melt.scaled.data, 
                  radius %in% c(25, 50, 100,1000))
all.melt$Site <- factor(paste(all.melt$Site),
                        levels = unlist(site.order$Site, use.names = FALSE),
                         ordered = TRUE)

for(var.name in sat.vars){
try({
    ggplot(subset(ss.melt, base.name==var.name),
           aes(x=value, y = ..density..,
               fill=as.factor(Mna), weight = Trap.weight)) +
        geom_histogram(alpha = 0.4, color = 'black',
                       position = 'identity', bins = 5) +
        guides(fill=guide_legend(title="Mn Present")) + 
    facet_grid(rows = vars(as.factor(Site)),
               cols = vars(as.factor(radius))) +
        theme(strip.text.y.right = element_text(angle = 0)) + 
        ggtitle(var.name)
    ggsave(paste0('Figures_Both/Compare_Hist_Vars/SS_',
                  paste(var.name), '_Hist.png'),
           height = 7, width = 7, units = 'in')

    ggplot(subset(all.melt, base.name==var.name),
           aes(x=value, y = ..density..,
               fill=as.factor(Mna), weight = Trap.weight)) +
        geom_histogram(alpha = 0.4, color = 'black',
                       position = 'identity', bins = 10) +
        guides(fill=guide_legend(title="Mn Present")) + 
    facet_grid(rows = vars(as.factor(Site)),
               cols = vars(as.factor(radius))) +
        theme(strip.text.y.right = element_text(angle = 0)) + 
        ggtitle(var.name)
    ggsave(paste0('Figures_Both/Compare_Hist_Vars/All_',
                  paste(var.name), '_Hist.png'),
           height = 14, width = 7, units = 'in')
    
})
}

wide.dataset <- rodent.data

## Function omit.names can be used to omit vector of strings
## (e.g. bare, tree) from glm.base.names
#glm.names <- omit.names(glm.base.names, 'tree')
glm.base.names <- c('Frac_bare', 'Frac_grass', 'Frac_tree', 'Frac_burn', 
                    'Frac_rice', 
                    'Density_Buildings')

glm.names <- glm.base.names
wide.dataset.norm <- create.norm.df(wide.dataset, glm.names)

wide.dataset.norm <- add.spatial.blocks(wide.dataset.norm)
wide.dataset.norm[,'block.i'] = factor(wide.dataset.norm[,'block.i'])

## --- Could adjust weights so that +/- are equal
## wi.pres <- wide.dataset.norm$Mna==1
## wi.abs <- wide.dataset.norm$Mna==0
## abs.adj <- sum(wi.pres) / sum(wi.abs)
## wide.dataset.norm$weight = 1
## wide.dataset.norm$weight[wi.abs] = abs.adj

## --- Perform simple ANOVA to evaluate within vs between site effects
simp.mod <- glm(Mna~ Site, data = wide.dataset.norm,
                family = 'binomial', weights = Trap.weight)
anova.mod <- anova(simp.mod)

simp.mod <- glm(Mna~ I(Site=='Bafodia'),
                data = subset(wide.dataset.norm, Source == 'PRE'),
                family = 'binomial', weights = Trap.weight)
anova.mod <- anova(simp.mod)
## 16%

simp.mod <- glm(Mna~ Site,
                data = subset(wide.dataset.norm, Source == 'PRE'),
                family = 'binomial', weights = Trap.weight)
anova.mod <- anova(simp.mod)
## 18% 

simp.mod <- glm(Mna~ I(Site%in%c('Bafodia', 'Tanganya', 'Bantou')),
                data = wide.dataset.norm,
                family = 'binomial', weights = Trap.weight)
anova.mod <- anova(simp.mod)
## 10%

## --- Rerun above GLM's with random effects

## RE for site
mod = glmer(paste('Mna ~ (1|Site)'),
            wide.dataset.norm, family = 'binomial',
            weights = Trap.weight)
summary(mod) ## AIC = 3687



## require('cAIC4')
## mod = glmer(paste('Mna ~ (1|Site)'),
##             wide.dataset.norm[1:20000,], family = 'binomial',
##             weights = Trap.weight)
## out = stepcAIC(mod, )


## Incorporate random effects for block and Site
mod = glmer(paste('Mna ~ (1|Site) + (1|block.i)'),
            wide.dataset.norm, family = 'binomial',
            weights = Trap.weight)
summary(mod) ## AIC = 3642



mod = glmer(paste(get.features(c(100), glm.names)$formula, '+ (1|Site) + (1|block.i)'),
    wide.dataset.norm, family = 'binomial', weights = Trap.weight)
summary(mod) ## AIC: 3592
## Random effects:
##  Groups  Name        Variance Std.Dev.
##  block.i (Intercept) 0.3542   0.5951  
##  Site    (Intercept) 1.4240   1.1933  
## Number of obs: 48135, groups:  block.i, 88; Site, 20

## Fixed effects:
##                       Estimate Std. Error z value Pr(>|z|)    
## (Intercept)           -5.11444    0.32002 -15.981  < 2e-16 ***
## Frac_bare.100          0.59786    1.79523   0.333    0.739    
## Frac_grass.100         0.10837    0.55469   0.195    0.845    
## Frac_tree.100          0.20535    1.64498   0.125    0.901    
## Frac_burn.100          0.05115    0.08749   0.585    0.559    
## Frac_rice.100          0.13591    0.23855   0.570    0.569    
## Density_Buildings.100  0.56219    0.14226   3.952 7.76e-05 ***

mod = glmer(paste(get.features(c(100), glm.names)$formula, '+ (1|Site)'),
    wide.dataset.norm, family = 'binomial', weights = Trap.weight)
summary(mod) ## AIC: 3566
## Random effects:
##  Groups Name        Variance Std.Dev.
##  Site   (Intercept) 1.566    1.251   
## Number of obs: 48135, groups:  Site, 20

## Fixed effects:
##                         Estimate Std. Error z value Pr(>|z|)    
## (Intercept)           -5.1092563  0.3182199 -16.056  < 2e-16 ***
## Frac_bare.100         -0.0842614  2.0147721  -0.042    0.967    
## Frac_grass.100        -0.3208447  0.6248629  -0.513    0.608    
## Frac_tree.100         -0.2957569  1.8611018  -0.159    0.874    
## Frac_burn.100          0.0456255  0.0729257   0.626    0.532    
## Frac_rice.100          0.0007181  0.2661323   0.003    0.998    
## Density_Buildings.100  0.4560159  0.1005820   4.534 5.79e-06 ***

mod = glmer(paste(get.features(c(200), glm.names)$formula, '+ (1|Site)'),
    wide.dataset.norm, family = 'binomial', weights = Trap.weight)
summary(mod) ## AIC: 3553
## Random effects:
##  Groups Name        Variance Std.Dev.
##  Site   (Intercept) 1.165    1.079   
## Number of obs: 48135, groups:  Site, 20

## Fixed effects:
##                       Estimate Std. Error z value Pr(>|z|)    
## (Intercept)           -4.92348    0.28092 -17.526  < 2e-16 ***
## Frac_bare.200         -0.02357    0.80095  -0.029    0.977    
## Frac_grass.200        -0.17379    0.25703  -0.676    0.499    
## Frac_tree.200         -0.14688    0.73598  -0.200    0.842    
## Frac_burn.200         -0.07020    0.13666  -0.514    0.607    
## Frac_rice.200          0.14255    0.11117   1.282    0.200    
## Density_Buildings.200  0.73157    0.18097   4.042 5.29e-05 ***

mod = glmer(paste(get.features(c(1000), glm.names)$formula, '+ (1|Site)'),
    wide.dataset.norm, family = 'binomial', weights = Trap.weight)
summary(mod) ## AIC: 3651
## Random effects:
##  Groups Name        Variance Std.Dev.
##  Site   (Intercept) 1.883    1.372   
## Number of obs: 48135, groups:  Site, 20
## Fixed effects:
##                        Estimate Std. Error z value Pr(>|z|)    
## (Intercept)             -6.0216     0.8025  -7.504 6.21e-14 ***
## Frac_bare.1000          -0.1849     0.7422  -0.249   0.8032    
## Frac_grass.1000         -1.3113     0.8424  -1.557   0.1195    
## Frac_tree.1000          -1.0280     0.9943  -1.034   0.3012    
## Frac_burn.1000          -3.8053     1.6919  -2.249   0.0245 *  
## Frac_rice.1000           0.3135     0.3118   1.005   0.3148    
## Density_Buildings.1000   0.9400     0.4546   2.068   0.0387 *  


modis.formula = paste('Mna~ ', paste(modis.vars, collapse = ' + '))
mod = glmer(paste(modis.formula, '+ (1|Site)'),
    wide.dataset.norm, family = 'binomial', weights = Trap.weight)
summary(mod) ## AIC: 




## Iterate to find a set of uncorrelated predictors. Threshold correlation
## is
thresh.corr = 0.7
radii <- c(1000,500,200,100,50,25)

## Only choose one set of EFC data
mask.data <- with(wide.dataset.norm, Source=='PRE' | (Source=='EFC' & Mna==1))
cor.data <- wide.dataset.norm[mask.data,]



all.preds <- full.sat.vars
tab.cor <- cor(cor.data[,all.preds])




## First remove predictors with no variation (correlation is NA)
mask.na <- is.na(tab.cor[,1])
foc.vars <- focal.variables[!mask.na]




## Load in all candidate variables
var.names <- get.features(c(25, 50,100,200,1000), sat.vars)$names

## Perform Wilcox test on each predictor in the candidate set
pvalues <- data.frame(var = var.names, p = NA)
for(i in 1:nrow(pvalues)){
    pvalues[i,'p'] = wilcox.test(get(paste(pvalues[i,'var']))~Mna,
                                 wide.dataset, weights = Trap.weight)$p.value
}

## Keep predictors that are significant at the p < 0.05 level
pvalues <- pvalues[order(pvalues$p, decreasing = FALSE),]




