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
preempt.data.path = '../Trap_Data/Data/Clean_PRE_Data_By_Trap.csv'
pre.dataset <- read.csv(preempt.data.path)
mask.presence <- pre.dataset$Sp%in%c('Mna')
pre.dataset$weight <- 1
pre.dataset$weight[!mask.presence] <- pre.dataset$pEmptyTrap[!mask.presence]
pre.dataset$Mna = 1*mask.presence

## --- Load EFC data
efc.dataset <- read.csv(efc.data.path)
efc.dataset <- subset(efc.dataset, Site %in% c('Bantou', 'Tanganya', 'Gania'))

## resample each trap location in efc dataset
efc.1 <- efc.dataset[efc.dataset$Mn==1,]
ws <- efc.1[,'weight']
efc.1$Mna = 1*(runif(length(ws)) < ws)


## --- Join EFC and PREEMPT data
rodent.data <- rbind.fill(efc.1, pre.dataset)


## Note that we omit cloud here
sat.vars <- c('Frac_bare', 'Frac_grass', 'Frac_tree', 
                     'Frac_burn', 'Frac_rice', 'Frac_water', 
                     'Frac_mound', 
                     'Density_Buildings', 'Density_Moderns',
                     'Density_Traditionals')
full.sat.vars <- get.features(rvec = c(50,100,200,1000), names = sat.vars)$names
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
for(r in c(50,100,200,1000)){
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
out = data.table::setDT(scaled.data)[,.(sum(Mna==1)),
                   by = list(Site)]
out[order(out$V1, decreasing = TRUE),]
site.order = out[order(out$V1, decreasing = TRUE),'Site']
##          Site  V1
##  1:    Bantou 336
##  2:   Bafodia 230
##  3:  Tanganya 209
##  4:     Largo  15
##  5:   Benduma  11
##  6:    Barlie   9
##  7:    Talama   8
##  8:   Kapethe   7
##  9:  Naiawama   7
## 10: Gbenikoro   6
## 11:     Guala   5
## 12:    Makump   5
## 13:    Badala   3
## 14:     Gania   2
## 15:  Njaguima   2
## 16:   Yekeyor   2
## 17:    Makuna   1
## 18:   Mokorie   1
## 19: Gbainkfay   0
## 20:    Petema   0


## Simple barplot of number of rodents captured
temp <- scaled.data#[scaled.data$Mna,]
temp$Site <- factor(temp$Site, levels = unlist(site.order, use.names = FALSE),
                    ordered = TRUE)
ggplot(temp) + geom_bar(aes(x = Site, fill = Mna),
                        position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    scale_fill_manual("Mn captures",
                      values = c("1" = "blue", "0" = "orange")) + 
ggtitle('Mastomys natalensis captures') + ylab('Count')
ggsave(filename = 'Figures_Both/Barplot_Captures_Mn.png')




library(tidyr)
temp = gather(scaled.data, variable, value,
           all_of(focal.variables), factor_key = TRUE)
temp = temp[,c('Site', 'Rodent', 'Lassa', 'Date', 'Sp', 'Mna', 'variable', 'value')]
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
       aes(x = as.factor(base.name), y = value, fill = Mna)) +
    geom_boxplot(alpha = 0.65) + xlab('') + ylab('Scaled Value') +
    coord_cartesian(ylim = c(-5,10)) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    facet_wrap(radius~., ncol = 1, labeller = label.fxn,
               strip.position="top")
ggsave(filename = 'Figures_Both/Mn_Pres_Abs_Features.png',
       height = 10, width = 14)


## -Plot Mastomys natalensis in high-trap towns only
dir.create('Figures_Both/Mn_Pres_By_Site')
site.names <- c('Bantou', 'Tanganya', 'Bafodia')
for(site in site.names){
    label.fxn <- as_labeller(function(y){paste0('Focal Radius: ', y)})
    ggplot(data = subset(sat.dat, Site==site),
           aes(x = base.name, y = value, fill = Mna)) +
        geom_boxplot(alpha = 0.65) + xlab('') + ylab('Scaled Value') +
        ##    coord_cartesian(ylim = c(-2.5,2.5)) + 
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
        facet_wrap(radius~., ncol = 1, labeller = label.fxn)+
        ggtitle(paste0('Trap outcome in ', site)) 
    ggsave(filename = paste0('Figures_Both/Mn_Pres_By_Site/', site, '_Mn_Pres_Abs_Features.png'),
       height = 10, width = 14)
}

label.fxn <- as_labeller(function(y){paste0('Focal Radius: ', y)})
ggplot(data = subset(sat.dat, Site%in% c('Bantou', 'Tanganya', 'Bafodia')),
       aes(x = base.name, y = value, fill = Mna)) +
    geom_boxplot(alpha = 0.65) + xlab('') + ylab('Scaled Value') +
    ##    coord_cartesian(ylim = c(-2.5,2.5)) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    facet_wrap(radius~., ncol = 1, labeller = label.fxn)+
    ggtitle(paste0('Trap outcome in ', site)) 
ggsave(filename = paste0('Figures_Both/Mn_Pres_By_Site/High_Mn_Pres_Abs_Features.png'),
       height = 10, width = 14)


dir.create('Figures_Both/Compare_Hist_Vars')
ss.melt <- subset(melt.scaled.data, Site %in% c('Bantou', 'Tanganya', 'Bafodia') &
                  radius %in% c(100,200,1000))

all.melt <- subset(melt.scaled.data, 
                  radius %in% c(100,200,1000))
all.melt$Site <- factor(paste(all.melt$Site), levels = unlist(site.order, use.names = FALSE),
                         ordered = TRUE)

for(var.name in sat.vars){
try({
    ggplot(subset(ss.melt, base.name==var.name),
           aes(x=value, y = ..density..,
               fill=as.factor(Mna))) +
        geom_histogram(alpha = 0.4, color = 'black',
                       position = 'identity', bins = 5) +
        guides(fill=guide_legend(title="Mn Present")) + 
    facet_grid(rows = vars(as.factor(Site)),
               cols = vars(as.factor(radius))) +
        theme(strip.text.y.right = element_text(angle = 0)) + 
        ggtitle(var.name)
    ggsave(paste0('Figures_Both/Compare_Hist_Vars/SS_',
                  paste(var.name), '_Hist.png'),
           height = 10, width = 10, units = 'in')

    ggplot(subset(all.melt, base.name==var.name),
           aes(x=value, y = ..density..,
               fill=as.factor(Mna))) +
        geom_histogram(alpha = 0.4, color = 'black',
                       position = 'identity', bins = 10) +
        guides(fill=guide_legend(title="Mn Present")) + 
    facet_grid(rows = vars(as.factor(Site)),
               cols = vars(as.factor(radius))) +
        theme(strip.text.y.right = element_text(angle = 0)) + 
        ggtitle(var.name)
    ggsave(paste0('Figures_Both/Compare_Hist_Vars/All_',
                  paste(var.name), '_Hist.png'),
           height = 10, width = 10, units = 'in')
    
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

wi.pres <- wide.dataset.norm$Mna==1
wi.abs <- wide.dataset.norm$Mna==0
abs.adj <- sum(wi.pres) / sum(wi.abs)
wide.dataset.norm$weight = 1
wide.dataset.norm$weight[wi.abs] = abs.adj

## --- Rerun above GLM's with random effects

## Incorporate random effects for block and Site
mod = glmer(paste('Mna ~ (1|Site) + (1|block.i)'),
    wide.dataset.norm, family = 'binomial')
summary(mod) ## 6418.4




mod = glmer(paste(get.features(c(100), glm.names)$formula, '+ (1|Site) + (1|block.i)'),
    wide.dataset.norm, family = 'binomial')
summary(mod) ## AIC: 5284
## Fixed effects:
##                       Estimate Std. Error z value Pr(>|z|)    
## (Intercept)           -5.15378    0.28509 -18.078  < 2e-16 ***
## Frac_bare.100          0.26923    1.86906   0.144    0.885    
## Frac_grass.100         0.01035    0.44509   0.023    0.981    
## Frac_tree.100         -0.16194    1.72032  -0.094    0.925    
## Frac_burn.100         -0.04333    0.14996  -0.289    0.773    
## Frac_rice.100          0.19337    0.29766   0.650    0.516    
## Density_Buildings.100  0.45437    0.07313   6.213 5.19e-10 ***


mod = glmer(paste(get.features(c(200), glm.names)$formula, '+ (1|Site) + (1|block.i)'),
    wide.dataset.norm, family = 'binomial')
summary(mod) ## AIC: 5280.4
## Fixed effects:
##                       Estimate Std. Error z value Pr(>|z|)    
## (Intercept)            -5.0232     0.2524 -19.904  < 2e-16 ***
## Frac_bare.200          -0.2910     0.7462  -0.390    0.697    
## Frac_grass.200         -0.1178     0.2030  -0.580    0.562    
## Frac_tree.200          -0.2638     0.6813  -0.387    0.699    
## Frac_burn.200          -0.1946     0.1567  -1.242    0.214    
## Frac_rice.200           0.1789     0.1234   1.450    0.147    
## Density_Buildings.200   0.8744     0.1234   7.087 1.37e-12 ***


mod = glmer(paste(get.features(c(1000), glm.names)$formula, '+ (1|Site) + (1|block.i)'),
    wide.dataset.norm, family = 'binomial')
summary(mod) ## AIC: 6384
## Fixed effects:
##                        Estimate Std. Error z value Pr(>|z|)    
## (Intercept)             -5.2835     0.4904 -10.773   <2e-16 ***
## Frac_bare.1000          -0.5178     1.0515  -0.492   0.6224    
## Frac_grass.1000         -1.7109     0.6954  -2.460   0.0139 *  
## Frac_tree.1000          -2.3261     1.1475  -2.027   0.0427 *  
## Frac_burn.1000          -2.4289     0.8484  -2.863   0.0042 ** 
## Frac_rice.1000           0.5201     0.3253   1.599   0.1099    
## Density_Buildings.1000   0.6505     0.2659   2.447   0.0144 *  



## Load in all candidate variables
var.names <- get.features(c(50,100,200,1000), sat.vars)$names

## Perform Wilcox test on each predictor in the candidate set
pvalues <- data.frame(var = var.names, p = NA)
for(i in 1:nrow(pvalues)){
    pvalues[i,'p'] = wilcox.test(get(paste(pvalues[i,'var']))~Mna, wide.dataset)$p.value
}

## Keep predictors that are significant at the p < 0.05 level
pvalues <- pvalues[order(pvalues$p, decreasing = FALSE),]




