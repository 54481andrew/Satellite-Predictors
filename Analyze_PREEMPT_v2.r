## v2 is built to read in wide data, convert to long (include MODIS),
## and make plots that show where Mna occur and do not.


## Figures and preliminary analyses of PREEMPT trap data with
## imagery-based predictors

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

## --- Load PREEMPT data with satellite predictors
preempt.data.path = '../Trap_Data/Data/Clean_PRE_Data_By_Trap.csv'
pre.dataset <- read.csv(preempt.data.path)
mask.presence <- pre.dataset$Sp%in%c('Mna')
pre.dataset$weight <- 1
pre.dataset$weight[!mask.presence] <- pre.dataset$pEmptyTrap[!mask.presence]
pre.dataset$Mna = 1*mask.presence

## --- Calculate some basic stats: number of traps total
dim(pre.dataset) ## 16623 ~ each row is one trap
## basic stats: number of traps in Bafodia
dim(subset(pre.dataset, Site=='Bafodia')) ## 2810 ~ each row is one trap

## Take into account the number of missing / CE traps by only adding up pEmptyTrap
temp <- pre.dataset[,c('Mna', 'pEmptyTrap', 'Rodent')] 
sum((1*(temp$Rodent==0)*temp$pEmptyTrap) + 1*(temp$Rodent==1))
## 15989

## Calculate number of Mna captures
sum(temp$Mna)
## 312

## --- Add columns, define important feature names



## Choose focal variables to normalize and (later) plot

## Note that we omit cloud here
sat.vars <- c('Frac_bare', 'Frac_grass', 'Frac_tree', 
                     'Frac_burn', 'Frac_rice', 'Frac_water', 
                     'Frac_mound', 
                     'Density_Buildings', 'Density_Moderns',
                     'Density_Traditionals')
lc.sat.vars = sat.vars[1:7]

full.sat.vars <- get.features(rvec = c(25,50,100,200,500,1000), names = sat.vars)$names
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

scaled.df <- pre.dataset
for(r in c(25,50,100,200,500,1000)){
    select.variables <- get.features(rvec = r, names = sat.vars)$names
    norm.out <- norm.preds(df = pre.dataset, foc.vars = select.variables)
    scaled.df[,select.variables] = norm.out[[1]]
}
nonsat.vars <- c(modis.vars, lc.vars, plag.vars)
norm.out <- norm.preds(df = pre.dataset, foc.vars = nonsat.vars)
scaled.df[,nonsat.vars] <- norm.out[[1]]


## Initialize new dataframe, scaled.df, with norm'd data and Site
scaled.data <- scaled.df
scaled.data$Mna <- factor(1*(scaled.data$Mna), levels = c('0','1'), ordered = TRUE)
scaled.data$ID <- 1:nrow(scaled.data)
scaled.data$Sp <- with(scaled.data, ifelse(Sp%in%c('Mna','Mer','Rra'), paste(Sp),
                                           ifelse(is.na(Sp), 'Absence', 'Other')))
scaled.data$Lassa <- factor(1*scaled.data$Lassa, levels = c('0','1'), ordered = TRUE)

## Which species are found where
table(scaled.data$Site, scaled.data$Sp)
##           Absence  Mer  Mna Other  Rra
## Badala       1889    2    3     7    5 
## Bafodia      2529   21  230    30    0
## Barlie       1083    1    9    26    6
## Benduma       939    0   11    12    2
## Gbainkfay     256    2    0     5    9
## Gbenikoro     355    0    6     9    0
## Guala        1663    0    5    32   15
## Kapethe      1039    2    7    20   10
## Largo        2052    1   15    16    2
## Makump        810    6    5    11    1
## Makuna        767    3    1    15    0
## Mokorie       269    2    1     3    9
## Naiawama      245    0    7     6    2
## Njaguima      384    3    2    14    1
## Petema        383    0    0     4    3
## Talama        662    3    8     7    0
## Yekeyor       648    1    2     4    5


## Badala: less buildings, more tree
## Bafodia: more buildings, less tree
## Barlie: weak
## Benduma: less tree, more building



## Simple barplot of number of rodents captured
not.murid <- c('Crocidura lamottei',
               'Crocidura olivieri', 'Crocidura buettikoferi',
               'Crocidura crossei', 'Crocidura therseae',
               'Crocidura grandiceps')
mask.shrew <- pre.dataset$Species_ID %in% not.murid
pre.dataset$Lassa[mask.shrew]
## Only 1 shrew had Lassa


not.murid <- c('Cla','Col','Cbu','Ccr','Cth', 'Cgr')
temp <- scaled.data[scaled.data$Rodent==1,]
temp$Sp <- with(temp, ifelse(Sp%in%c('Mna','Mer','Rra'), paste(Sp), 'Other'))
ggplot(temp) + geom_bar(aes(x = Site, fill = Sp),
                        position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    ggtitle('All captures') + ylab('Count')
ggsave(filename = 'Figures_PREEMPT/Barplot_Captures_All.png')

    
## Barplot of number of captures for total species rodents, shrews
data.table::setDT(pre.dataset)
temp = pre.dataset[,.(Total = length(unique(Sp[!is.na(Sp)])),
                      Rodents = length(unique(Sp[!is.na(Sp) & !(Sp%in%not.murid)])),
                      Shrews =  length(unique(Sp[!is.na(Sp) & (Sp%in%not.murid)]))),
                   by = list(Site)]
sp.site.order = temp$Site[order(temp$Total, decreasing = TRUE)]
temp = gather(data = temp, value, number, Total:Shrews, factor_key = TRUE)
temp$Site <- factor(temp$Site, levels = sp.site.order, ordered = TRUE)

ggplot(temp, aes(x = (Site), y = number, fill = value)) +
    geom_bar(stat = "identity",
             position = position_dodge())



## Simple barplot of number of rodents captured
temp <- scaled.data[scaled.data$Mna==1,]
ggplot(scaled.data) + geom_bar(aes(x = Site, fill = Sp==Mna),
                        position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
ggtitle('Mastomys natalensis captures') + ylab('Count')
ggsave(filename = 'Figures_PREEMPT/Barplot_Captures_Mn.png')

out = data.table::setDT(scaled.data)[,.(sum(Mna==1)),
                   by = list(Site)]
out[order(out$V1, decreasing = TRUE),]
##          Site  V1
##  1:   Bafodia 230
##  2:     Largo  15
##  3:   Benduma  11
##  4:    Barlie   9
##  5:    Talama   8
##  6:   Kapethe   7
##  7:  Naiawama   7
##  8: Gbenikoro   6
##  9:     Guala   5
## 10:    Makump   5
## 11:    Badala   3
## 12:  Njaguima   2
## 13:   Yekeyor   2
## 14:    Makuna   1
## 15:   Mokorie   1
## 16: Gbainkfay   0
## 17:    Petema   0
site.order = out[order(out$V1, decreasing = TRUE),'Site']

## Melt data so that it can be plotted in ggplot
##melt.cleaned.data <- melt(cleaned.data, id.vars = c('ID', 'Mna', 'Sp', 'Lassa', 'radius', 'Site'))
## --- Make a long version of the PRE dataset
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
    coord_cartesian(ylim = c(-2.5,2.5)) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    facet_wrap(radius~., ncol = 1, labeller = label.fxn,
               strip.position="top")
ggsave(filename = 'Figures_PREEMPT/Mn_Pres_Abs_Features.png',
       height = 10, width = 14)
## for(r in c(50,100,200,1000)){
##     ggplot(data = subset(sat.dat, radius==r),
##            aes(x = variable, y = value, fill = Mna)) +
##         geom_boxplot() + xlab('') + ylab('Scaled Value') +
##         theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
##         facet_wrap(radius~., ncol = 1, labeller = label.fxn) 
##     ggsave(filename = paste0('Figures_PREEMPT/Mn_Pres_Abs_Features_', r,'.png'),
##            height = 6, width = 14)
## }

ggplot(data = subset(sat.dat, Site!='Bafodia'),
       aes(x = as.factor(base.name), y = value, fill = Mna)) +
    geom_boxplot(alpha = 0.65) + xlab('') + ylab('Scaled Value') +
    coord_cartesian(ylim = c(-2.5,2.5)) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    facet_wrap(radius~., ncol = 1, labeller = label.fxn,
               strip.position="top")
ggsave(filename = 'Figures_PREEMPT/Mn_Pres_Abs_Features_NOT_Bafodia.png',
       height = 10, width = 14)


## -Plot Mastomys natalensis in Bafodia only
site.names <- unique(sat.dat$Site)
for(site in site.names){
    label.fxn <- as_labeller(function(y){paste0('Focal Radius: ', y)})
    ggplot(data = subset(sat.dat, Site==site),
           aes(x = base.name, y = value, fill = Mna)) +
        geom_boxplot(alpha = 0.65) + xlab('') + ylab('Scaled Value') +
            coord_cartesian(ylim = c(-2.5,2.5)) + 
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
        facet_wrap(radius~., ncol = 1, labeller = label.fxn)+
        ggtitle(paste0('Trap outcome in ', site)) 
    ggsave(filename = paste0('Figures_PREEMPT/Mn_Pres_By_Site/', site, '_Mn_Pres_Abs_Features.png'),
       height = 10, width = 14)
}

## create directory for Bafodia-centered analyses
dir.create('Figures_PREEMPT/BAF', showWarnings = FALSE)

## ggplot(data = subset(sat.dat, radius==1000 & base.name %in% c('Density_Buildings', 'Frac_rice', 'Frac_tree')),
##        aes(y = Site, x = value, fill = base.name)) +
##     geom_boxplot(alpha = 0.65) + xlab('') + ylab('Scaled Value') +
##     theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
##     ggtitle('Building density at 1000m') 
## ggsave(filename = 'Figures_PREEMPT/Site_Variation.png',
##        height = 10, width = 14)

## for(r in c(50,100,200,1000)){
##     ggplot(data = subset(sat.dat, radius==r & Site=='Bafodia'),
##            aes(x = base.name, y = value, fill = Mna)) +
##         geom_boxplot(alpha = 0.3) + xlab('') + ylab('Scaled Value') + 
##         theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
##         facet_wrap(radius~., ncol = 1, labeller = label.fxn) 
##     ggsave(filename = paste0('Figures_PREEMPT/BAF/BAF_Mn_Pres_Abs_Features_', r,'.png'),
##            height = 6, width = 14)
}

tbl[order(tbl, decreasing = TRUE)]
## Mna Rra Lsi Col Mer Cth Mmi Pro Mse Pda Cbu Lst Hsi Cgr Gsp Gke Uru Lfl Mmu 
## 312  70  50  47  47  27  27  22  12   8   7   5   4   3   3   2   2   1   1 



## -Plot across all species
label.fxn <- as_labeller(function(y){paste0('Focal Radius: ', y)})
ggplot(data = sat.dat,
       aes(x = base.name, y = value, fill = Sp)) +
    geom_boxplot(alpha = 0.65) + xlab('') + ylab('Scaled Value') +
    coord_cartesian(ylim = c(-2.5,2.5)) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    facet_wrap(radius~., ncol = 1, labeller = label.fxn)
ggsave(filename = 'Figures_PREEMPT/All_Species_Pres_Abs_Features.png',
       height = 10, width = 14)
## for(r in c(50,100,200,1000)){
##     ggplot(data = subset(sat.dat, radius==r),
##            aes(x = base.name, y = value, fill = Sp)) +
##         xlab('') + ylab('Scaled Value') + 
##         geom_boxplot() + 
##         theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
##         facet_wrap(radius~., ncol = 1, labeller = label.fxn)
##     ggsave(filename = paste0('Figures_PREEMPT/All_Species_Pres_Abs_Features_', r,'.png'),
##            height = 6, width = 14)
## }

## -Plot across species for Bafodia only
label.fxn <- as_labeller(function(y){paste0('Focal Radius: ', y)})
ggplot(data = subset(sat.dat, Site=='Bafodia'),
       aes(x = base.name, y = value, fill = Sp)) +
    geom_boxplot(alpha = 0.65) + xlab('') + ylab('Scaled Value') +
    coord_cartesian(ylim = c(-5,10)) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    facet_wrap(radius~., ncol = 1, labeller = label.fxn)
ggsave(filename = 'Figures_PREEMPT/BAF/BAF_All_Species_Pres_Abs_Features.png',
       height = 10, width = 14)
## for(r in c(50,100,200,1000)){
##     ggplot(data = subset(sat.dat, radius==r & Site=='Bafodia'),
##            aes(x = base.name, y = value, fill = Sp)) +
##         xlab('') + ylab('Scaled Value') + 
##         geom_boxplot(alpha = 0.65) + 
##         theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
##         facet_wrap(radius~., ncol = 1, labeller = label.fxn)
##     ggsave(filename = paste0('Figures_PREEMPT/BAF/BAF_All_Species_Pres_Abs_Features_', r,'.png'),
##            height = 6, width = 14)
## }


## --- LASSA PCR PLOTS
table(pre.dataset$Site[pre.dataset$Lassa==1], pre.dataset$Sp[pre.dataset$Lassa==1])
##           Col Lsi Mer Mmi Mna Rra
## Badala      0   0   0   2   0   0
## Bafodia     0   0   0   2  11   0
## Barlie      0   3   0   0   1   0
## Gbenikoro   0   0   0   1   0   0
## Guala       0   0   0   1   0   0
## Kapethe     0   0   0   0   2   1
## Largo       0   0   0   0   3   0
## Makump      1   0   0   0   1   0
## Talama      0   0   1   0   1   0


## -Mastomys natalensis only

temp <- scaled.data[scaled.data$Mna==1,]
ggplot(temp) + geom_bar(aes(x = Site, fill = Lassa),
                        position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
ggtitle('Mastomys natalensis Lassa +') + ylab('Count')
ggsave('Figures_PREEMPT/Lassa_Mn.png')
## 9 locations with Lassa

temp <- scaled.data[scaled.data$Rodent==1,]
ggplot(temp) + geom_bar(aes(x = Site, fill = Lassa),
                        position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
ggtitle('Rodent Lassa +') + ylab('Count')
ggsave('Figures_PREEMPT/Lassa_Rodent.png')
## 9 locations with Lassa

dir.create('Figures_PREEMPT/PCR', showWarnings = FALSE)
label.fxn <- as_labeller(function(y){paste0('Focal Radius: ', y)})
ggplot(data = subset(sat.dat, Mna==1 & !is.na(Lassa)),
       aes(x = base.name, y = value, fill = Lassa)) +
    geom_boxplot() + xlab('') + ylab('Scaled Value') +
    coord_cartesian(ylim = c(-2.5,2.5)) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    facet_wrap(radius~., ncol = 1, labeller = label.fxn) +
    ggtitle('Lassa in Mna')
ggsave(filename = 'Figures_PREEMPT/PCR/Mn_PCR_Features.png',
       height = 10, width = 14)
## for(r in c(50,100,200,1000)){
##     ggplot(data = subset(sat.dat, radius==r & Mna==1 & !is.na(Lassa)),
##            aes(x = base.name, y = value, fill = Lassa)) +
##         geom_boxplot() + xlab('') + ylab('Scaled Value') + 
##         theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
##         facet_wrap(radius~., ncol = 1, labeller = label.fxn)
##     ggsave(filename = paste0('Figures_PREEMPT/PCR/Mn_PCR_Features_', r,'.png'),
##            height = 6, width = 14)
## }

## -All species
label.fxn <- as_labeller(function(y){paste0('Focal Radius: ', y)})
ggplot(data = subset(sat.dat, !is.na(Lassa)),
       aes(x = base.name, y = value, fill = Lassa)) +
    geom_boxplot() + xlab('') + ylab('Scaled Value') +
    coord_cartesian(ylim = c(-2.5,2.5)) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    facet_wrap(radius~., ncol = 1, labeller = label.fxn) + 
    ggtitle('Lassa in Rodents')
ggsave(filename = 'Figures_PREEMPT/PCR/All_Species_PCR_Features.png',
       height = 10, width = 14)
## for(r in c(50,100,200,1000)){
##     ggplot(data = subset(sat.dat, radius==r & !is.na(Lassa)),
##            aes(x = base.name, y = value, fill = Lassa)) +
##         geom_boxplot() + xlab('') + ylab('Scaled Value') + 
##         theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
##         facet_wrap(radius~., ncol = 1, labeller = label.fxn)
##     ggsave(filename = paste0('Figures_PREEMPT/PCR/All_Species_PCR_Features_', r,'.png'),
##            height = 6, width = 14)
## }

## --- Make similar graphs for MODIS predictors here

## lc.vars

## plag.vars

## modis.vars




## --- Reshape data to make each radius x feature a unique predictor
wide.dataset <- pre.dataset
dir.create('Figures_PREEMPT/Corr', recursive = TRUE)
## Plot correlation matrices and predictor histograms
feature.names <- c()
for(r in c(50,100,200,1000)){
    feature.names <- get.features(rvec = r, sat.vars)$names
    ## Plot correlation matrices
    sd.df <- apply(wide.dataset[,feature.names], 2, sd, na.rm = TRUE)
    wide.df <- wide.dataset[,c(feature.names[!is.na(sd.df) & sd.df > 0])]
    cor.mat = round(cor(wide.df,
                        use = 'complete.obs'),1)
    p.mat <- cor_pmat(wide.df)
    ggcorrplot(cor.mat,
               p.mat = p.mat, colors = c('red', 'white', 'blue'),
               lab = TRUE) + ggtitle(paste0('Radius: ', r, ' m'))
    ggsave(paste0('Figures_PREEMPT/Corr/Correlation_Matrix_',r,'.png'))

}    

## The histograms indicate that, at the 1000m scale, the predictors
## that describe landcover types bare, grass, burn, rice, garden,
## and building predictors are all correlated with Mn presence/absence.
## We'll omit garden because this might not be predicted well. We'll only
## include a single set of building descriptors (type frac, for convenience). 
dir.create('Figures_PREEMPT/Compare_Hist_Vars')
melt.scaled.data$alt.Site = paste(ifelse(melt.scaled.data$Site=='Bafodia', 'Bafodia',
                                       'NotBafodia'))
temp <- melt.scaled.data
temp$alt.Site <- 'All'
melt.scaled.data.comb <- rbind(melt.scaled.data, temp)


for(var.name in sat.vars){
try({
    ggplot(subset(melt.scaled.data.comb, base.name==var.name),
           aes(x=value, y = ..density..,
               fill=as.factor(Mna))) +
        geom_histogram(alpha = 0.4, color = 'black',
                       position = 'identity', bins = 5) +
        guides(fill=guide_legend(title="Mn Present")) + 
    facet_grid(rows = vars(as.factor(Site)),
               cols = vars(as.factor(radius))) +
        theme(strip.text.y.right = element_text(angle = 0)) + 
        ggtitle(var.name)
    ggsave(paste0('Figures_PREEMPT/Compare_Hist_Vars/',
                  paste(var.name), '_Hist.png'),
           height = 10, width = 5, units = 'in')


    ggplot(subset(melt.scaled.data.comb, base.name==var.name),
           aes(x=value, y = ..density..,
               fill=as.factor(Mna))) +
        geom_histogram(alpha = 0.4, color = 'black',
                       position = 'identity', bins = 15) +
        guides(fill=guide_legend(title="Mn Present")) + 
        facet_grid(rows = vars(as.factor(alt.Site)),
                   cols = vars(as.factor(radius))) +
        theme(strip.text.y.right = element_text(angle = 0)) + 
        ggtitle(var.name)
    ggsave(paste0('Figures_PREEMPT/Compare_Hist_Vars/BAF_',
                  paste(var.name), '_Hist.png'))

})
}
        
writeLines('\n\n *** Stopping here; CHECK THIS CODE *** \n\n')
q()



## --- Fit generalized linear regression models

## Begin with a simple reasonable model: 1000m predictors are all that
## matter Restric analyses to only a few features that show an obvious
## relationship to Mn presence / absence in barplots.
glm.base.names <- c('Frac_bare', 'Frac_grass', 'Frac_tree', 'Frac_burn', 
                    'Frac_rice', 
                    'Density_Buildings')

## Function omit.names can be used to omit vector of strings
## (e.g. bare, tree) from glm.base.names
#glm.names <- omit.names(glm.base.names, 'tree')
glm.names <- glm.base.names
wide.dataset.norm <- create.norm.df(wide.dataset, glm.names)

## NULL model
glm('Mna~1', wide.dataset.norm, family = 'binomial') ## 3099

best.glm <- glm(get.features(1000, glm.names)$formula, wide.dataset.norm,
                family = 'binomial') 
summary(best.glm)
## Coefficients:
##                        Estimate Std. Error z value Pr(>|z|)    
## (Intercept)            -4.75640    0.10366 -45.883  < 2e-16 ***
## Frac_bare.1000         -0.52578    0.42612  -1.234 0.217255    
## Frac_grass.1000         0.40795    0.37851   1.078 0.281137    
## Frac_tree.1000          0.29811    0.65646   0.454 0.649738    
## Frac_burn.1000          0.57310    0.22093   2.594 0.009487 ** 
## Frac_rice.1000          0.54851    0.14665   3.740 0.000184 ***
## Density_Buildings.1000  0.65003    0.09282   7.003  2.5e-12 ***


## --- Use stepAIC to consider more complex models with finer spatial scales
## k = 2 denotes standard AIC
best.glm.step = stepAIC(best.glm, scope = list(lower = 'Mna ~ 1',
                                               upper = get.features(c(1000), glm.names)$formula),
                        k = 2)
summary(best.glm.step)
## Coefficients:
##                        Estimate Std. Error z value Pr(>|z|)    
## (Intercept)            -4.76229    0.10340 -46.057  < 2e-16 ***
## Frac_bare.1000         -0.71017    0.12363  -5.744 9.23e-09 ***
## Frac_grass.1000         0.24131    0.09560   2.524   0.0116 *  
## Frac_burn.1000          0.48250    0.09521   5.068 4.02e-07 ***
## Frac_rice.1000          0.48909    0.06394   7.650 2.02e-14 ***
## Density_Buildings.1000  0.64965    0.09302   6.984 2.87e-12 ***




## --- Incorporate random effects to wash out spatial autocorrelation

## Spatial autocorrelation biases the results above; 1000m predictors
## are essentially repeated 1000's of times per village; if villages
## A and B differ both in a confounding predictor like Bare_Area.1000,
## as well as a predictor that is consequential to Mna abundance,
## the regression will suggest that Bare_Area is heavily significant.

## Random effects at different spatial scales can wash out these
## confounding factors. The regression will correct for unknown
## spatial aspects by allowing each village's (or trap lines)
## intercept to come from a normal distribution.


## -Organize data into spatial blocks in each village


## Note on converting meters to lat / lon
## If your displacements
## aren't too great (less than a few kilometers) and you're not right at
## the poles, use the quick and dirty estimate that 111,111 meters
## (111.111 km) in the y direction is 1 degree (of latitude) and 111,111
## * cos(latitude) meters in the x direction is 1 degree (of longitude).


wide.dataset.norm <- add.spatial.blocks(wide.dataset.norm)
wide.dataset.norm[,'block.i'] = factor(wide.dataset.norm[,'block.i'])


## --- Rerun above GLM's with random effects


## Incorporate random effects for block and Site
mod = glmer(paste('Mna ~ (1|Site) + (1|block.i)'),
    wide.dataset.norm, family = 'binomial')
summary(mod) ## 1572


mod = glmer(paste(get.features(c(1000), glm.names)$formula, '+ (1|Site) + (1|block.i)'),
    wide.dataset.norm, family = 'binomial')
summary(mod) ## 1571
## Fixed effects:
##                        Estimate Std. Error z value Pr(>|z|)    
## (Intercept)            -5.28895    0.29848 -17.720   <2e-16 ***
## Frac_bare.1000         -0.69162    0.62987  -1.098    0.272    
## Frac_grass.1000        -0.58489    0.69340  -0.844    0.399    
## Frac_tree.1000         -0.99518    1.13036  -0.880    0.379    
## Frac_burn.1000         -0.38737    0.51683  -0.750    0.454    
## Frac_rice.1000          0.51965    0.29090   1.786    0.074 .  
## Density_Buildings.1000  0.06245    0.29740   0.210    0.834    


mod = glmer(paste(get.features(c(50,100,200,1000), glm.names)$formula, '+ (1|Site) + (1|block.i)'),
    wide.dataset.norm, family = 'binomial')
summary(mod) ## 1571
## Fixed effects:
##                        Estimate Std. Error z value Pr(>|z|)    
## (Intercept)            -5.20208    0.24663 -21.093   <2e-16 ***
## Frac_bare.50           -4.83326    4.06770  -1.188    0.235    
## Frac_grass.50          -0.99391    0.99578  -0.998    0.318    
## Frac_tree.50           -4.72594    3.71231  -1.273    0.203    
## Frac_burn.50           -0.30399    0.24776  -1.227    0.220    
## Frac_rice.50           -0.55710    0.42697  -1.305    0.192    
## Density_Buildings.50    0.30764    0.20748   1.483    0.138    
## Frac_bare.100           2.41695    5.93619   0.407    0.684    
## Frac_grass.100          0.44316    1.62428   0.273    0.785    
## Frac_tree.100           1.86435    5.36190   0.348    0.728    
## Frac_burn.100           0.20820    0.27234   0.764    0.445    
## Frac_rice.100           0.26749    0.70227   0.381    0.703    
## Density_Buildings.100  -0.49558    0.33631  -1.474    0.141    
## Frac_bare.200          -0.68173    0.84779  -0.804    0.421    
## Frac_grass.200         -0.02134    0.37488  -0.057    0.955    
## Frac_tree.200          -0.08975    0.71570  -0.125    0.900    
## Frac_burn.200          -0.10289    0.21435  -0.480    0.631    
## Frac_rice.200           0.20280    0.17918   1.132    0.258    
## Density_Buildings.200   0.70141    0.31501   2.227    0.026 *  
## Frac_bare.1000         -0.47556    0.61088  -0.778    0.436    
## Frac_grass.1000        -0.31440    0.64458  -0.488    0.626    
## Frac_tree.1000         -0.51007    1.06484  -0.479    0.632    
## Frac_burn.1000         -0.29886    0.64446  -0.464    0.643    
## Frac_rice.1000          0.25848    0.28605   0.904    0.366    
## Density_Buildings.1000  0.13461    0.25898   0.520    0.603    


mod = glmer(paste(get.features(c(200), glm.names)$formula, '+ (1|Site) + (1|block.i)'),
    wide.dataset.norm, family = 'binomial')
summary(mod) ## 2476
## Fixed effects:
##                       Estimate Std. Error z value Pr(>|z|)    
## (Intercept)           -5.12049    0.28784 -17.789  < 2e-16 ***
## Frac_bare.200         -0.07255    0.61612  -0.118  0.90627    
## Frac_grass.200        -0.02159    0.26789  -0.081  0.93578    
## Frac_tree.200         -0.11041    0.56190  -0.196  0.84422    
## Frac_burn.200         -0.14520    0.17335  -0.838  0.40223    
## Frac_rice.200          0.19039    0.14065   1.354  0.17584    
## Density_Buildings.200  0.47813    0.18430   2.594  0.00948 ** 



wi.pres <- wide.dataset.norm$Mna==1
wi.abs <- wide.dataset.norm$Mna==0
abs.adj <- sum(wi.pres) / sum(wi.abs)
wide.dataset.norm$weight = 1
wide.dataset.norm$weight[wi.abs] = abs.adj


## Just Bafodia
mod = glmer(paste(get.features(c(100), glm.names)$formula, ' + (1|block.i)'),
    subset(wide.dataset.norm, Site=='Bafodia'), family = 'binomial')
summary(mod) ## 1571
## Fixed effects:
##                       Estimate Std. Error z value Pr(>|z|)    
## (Intercept)            -3.3833     0.2070 -16.341  < 2e-16 ***
## Frac_bare.100          -7.0557    12.4161  -0.568    0.570    
## Frac_grass.100         -1.8681     3.3552  -0.557    0.578    
## Frac_tree.100          -5.5782    11.2935  -0.494    0.621    
## Frac_burn.100          -0.2218     0.5295  -0.419    0.675    
## Frac_rice.100          -0.6253     1.4719  -0.425    0.671    
## Density_Buildings.100   1.3802     0.2735   5.047 4.48e-07 ***


q() ## here

## ---- Spatial blocks with PCs
wide.dataset <- add.spatial.blocks(wide.dataset)
pc.out <- add.pcs(wide.dataset, pred.names = get.features(1000, glm.base.names)$names,
                  thresh = 0.95)
wide.dataset.pc <- pc.out[[1]]
mod = glmer(paste(pc.out$formula, '+ (1|Site) + (1|block.i)'),
    wide.dataset.pc, family = 'binomial')
summary(mod) ## 
## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  -5.2602     0.2960 -17.771   <2e-16 ***
## PC1           0.2573     0.1641   1.568   0.1169    
## PC2           0.4500     0.2016   2.232   0.0256 *  
## PC3          -0.2550     0.2406  -1.060   0.2893    
## PC4          -0.4733     0.2917  -1.622   0.1047    

## PC2 describes ares with high low trad, high burn, low bare,
## high modern buildings, low traditional buildings

dir.create('Figures_PREEMPT/PCs')
mat = pc.out$rot
st = 1.2
text.labels <- row.names(mat)
for(pci in 1:pc.out$num){
    jpeg(filename = paste0('Figures_PREEMPT/PCs/Pred_1000_PCs',pci,'.jpg'),
         width = 4, height = 4, units = 'in', res = 300)
    par(mfrow = c(1,1), mar = c(5,2,1,0))
    barplot(mat[,pci], yaxt = 'n', main = paste0('PC ',pci),
            horiz = T, names.arg = '',
            xlim = c(-.75, 0.75))
    text(y = (1:nrow(mat))*st -0.6,x = 0.001, labels = text.labels, xpd = NA,
         pos = 4, srt = 0, cex = 1)
    dev.off()
}


## --- Do PC's again but with all predictors!!

pc.out.all <- add.pcs(wide.dataset,
                      pred.names = get.features(c(25, 50, 100, 200, 500, 1000), glm.base.names)$names,
                  thresh = 0.95)
wide.dataset.pc <- pc.out.all[[1]]
mod.all = glmer(paste(pc.out.all$formula, '+ (1|Site) + (1|block.i)'),
    wide.dataset.pc, family = 'binomial')
summary(mod.all) ## 
## Fixed effects:
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -5.2001962  0.2419283 -21.495  < 2e-16 ***
## PC1         -0.1447844  0.0438808  -3.299 0.000969 ***
## PC2          0.0923872  0.0961774   0.961 0.336758    
## PC3          0.4141208  0.1540894   2.688 0.007198 ** 
## PC4         -0.0470893  0.1410880  -0.334 0.738562    
## PC5         -0.0007827  0.1031920  -0.008 0.993949    
## PC6         -0.2051790  0.1080706  -1.899 0.057622 .  
## PC7         -0.1506639  0.0982134  -1.534 0.125018    
## PC8         -0.0978294  0.1433311  -0.683 0.494897    
## PC9         -0.1896142  0.1229141  -1.543 0.122914    
## PC10        -0.0401318  0.1409965  -0.285 0.775928    
## PC11        -0.0798600  0.1276608  -0.626 0.531601    
## PC12         0.0206339  0.0973410   0.212 0.832126    
## PC13         0.0135435  0.1406623   0.096 0.923295    
## PC14         0.0235842  0.1248030   0.189 0.850115    
## PC15         0.0411900  0.1199239   0.343 0.731247    


best.mod.all = glmer(paste0('Mna~PC1 + PC3 + PC6', '+ (1|Site) + (1|block.i)'),
    wide.dataset.pc, family = 'binomial')
summary(best.mod.all) ## 
## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -5.13617    0.25927 -19.810  < 2e-16 ***
## PC1         -0.16121    0.02953  -5.460 4.76e-08 ***
## PC3          0.38360    0.08279   4.633 3.60e-06 ***
## PC6         -0.13614    0.07617  -1.787   0.0739 .  


st = 1.2
mat <- pc.out.all$rot
text.labels <- row.names(mat)
for(pci in c(1,3,6)){
    jpeg(filename = paste0('Figures_PREEMPT/PCs/Pred_all_PCs',pci,'.jpg'),
         width = 4, height = 5, units = 'in', res = 300)
    par(mfrow = c(1,1), mar = c(5,2,1,0))
    barplot(mat[,pci], yaxt = 'n', main = paste0('PC ',pci),
            horiz = T, names.arg = '',
            xlim = c(-.75, 1))
    text(y = (1:nrow(mat))*st -0.6,x = 0.001, labels = text.labels, xpd = NA,
         pos = 4, srt = 0, cex = 0.6)
    dev.off()
}

## -- hereherehere

## --- repeat above with just Bafodia
pc.out.baf <- add.pcs(wide.dataset,
                      pred.names = get.features(c(50, 100, 200), glm.base.names)$names,
                  thresh = 0.95)
wide.dataset.pc <- pc.out.baf[[1]]
mod = glmer(paste(pc.out.baf$formula, '+ (1|Site) + (1|block.i)'),
    wide.dataset.pc, family = 'binomial')
summary(mod) ## 
## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -5.25142    0.30055 -17.473  < 2e-16 ***
## PC1         -0.12337    0.07084  -1.741 0.081598 .  
## PC2         -0.17223    0.08919  -1.931 0.053471 .  
## PC3          0.07196    0.09505   0.757 0.448974    
## PC4          0.01283    0.16956   0.076 0.939677    
## PC5          0.11175    0.08717   1.282 0.199828    
## PC6         -0.25189    0.11266  -2.236 0.025362 *  
## PC7          0.14533    0.14860   0.978 0.328051    
## PC8         -0.03126    0.12911  -0.242 0.808667    
## PC9         -0.04741    0.15385  -0.308 0.757988    
## PC10        -0.47526    0.12243  -3.882 0.000104 ***
## PC11         0.13635    0.17662   0.772 0.440111    

mod.baf = glmer(paste('Mna~ PC1 + PC2 + PC6 + PC10 + (1|block.i)'),
    wide.dataset.pc, family = 'binomial')
summary(mod.baf) ## 
## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -5.19983    0.29488 -17.634  < 2e-16 ***
## PC1         -0.03036    0.06152  -0.494  0.62165    
## PC2         -0.10415    0.07195  -1.447  0.14778    
## PC6         -0.19544    0.08083  -2.418  0.01561 *  
## PC10        -0.32269    0.10357  -3.116  0.00184 ** 

st = 1.2
mat <- pc.out.baf$rot
text.labels <- row.names(mat)
for(pci in c(6,10)){
    jpeg(filename = paste0('Figures_PREEMPT/Pred_baf_PCs',pci,'.jpg'),
         width = 4, height = 5, units = 'in', res = 300)
    par(mfrow = c(1,1), mar = c(5,2,1,0))
    barplot(mat[,pci], yaxt = 'n', main = paste0('PC ',pci),
            horiz = T, names.arg = '',
            xlim = c(-.75, 1))
    text(y = (1:nrow(mat))*st -0.6,x = 0.001, labels = text.labels, xpd = NA,
         pos = 4, srt = 0, cex = 1)
    dev.off()
}

## PC6: (with coefficient < 0) more modern and more grass = more Mna
##200: low modern, high rice, high tree
##100: high trad low modern, high rice, high tree

## PC10 (with coefficient < 0) more modern, grass, rice = more Mna
## mixed on trad -- more Mna = less trad200/100, but more trad50


## --- Plot residuals

## Add in residuals to wide.dataset
summ.glmer.mod <- summary(best.mod.all)
summ.glmer.resid <- summ.glmer.mod$residuals
summ.glm.mod <- summary(best.glm.step)
summ.glm.resid <- summ.glm.mod$deviance.resid

wide.dataset$glm.resids = summ.glm.resid
wide.dataset$glmer.resids = summ.glmer.resid
wide.dataset$glm.fits = predict(best.glm.step, type = 'response')
wide.dataset$glmer.fits = predict(best.mod.all, type = 'response')


## Calculate residual for all traps in each block and village
agg.dev <- data.table::setDT(wide.dataset)[,.(glm.tot.resid = mean(glm.resids),
                                              glmer.tot.resid = mean(glmer.resids)),
                                           by = list(Site, block.i)]

cuts = c(-Inf, -0.025, -0.01, 0.01, 0.025, Inf)
cols = c('red', 'pink', 'black', 'lightblue', 'blue')

agg.focal.wide.data$glmer.cuts <- sapply(agg.focal.wide.data$glmer.tot.resids,
                                         function(x){max(which(x > cuts))})


## Fit all data with BRT

require(dismo)
require(raster)

wide.dataset$Mna = 1*(wide.dataset$Sp=='Mna')
raster.data.path <- '../../Storage/Raster_Data'
pred.stack <- raster::stack(x = paste0(raster.data.path, "/predictor_stack_05.grd"))
names.stack <- names(pred.stack)
keep.stack <- names.stack[!(names.stack %in% c('Water_Bodies', 'TS.House', 'TS.InTown'))]
pred.stack <- pred.stack[[keep.stack]]
hab.preds <- as.data.frame(extract(pred.stack, wide.dataset[,c('Longitude', 'Latitude')]))
## hab.sd <- apply(hab.preds, 2, sd) 
## hab.preds <- hab.preds[, !is.na(hab.sd) & hab.sd > 0]
hab.pred.names <- names(hab.preds)

comb.wide.dataset <- cbind(wide.dataset, hab.preds)

feature.names = get.features(c(50,100,200,1000), var.names)$names

x.names = c(hab.pred.names, feature.names)
imp.dataset <- comb.wide.dataset[,x.names]
standard.dev <- apply(imp.dataset, 2, sd)
nona.sd <- !is.na(standard.dev) & standard.dev > 0
nona.x.names <- x.names[nona.sd]
imp.dataset <- imp.dataset[,nona.x.names]
imp.names <- names(imp.dataset)

for(col in nona.x.names){
    x = imp.dataset[,col]
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    imp.dataset[,col] = x
}


## --- Implement smart fold (Just Bafodia!)
folded.data.baf <- add.cv.folds(train.data = subset(comb.wide.dataset, Site=='Bafodia'),
                            bloocv.r = 50)
#                            plot.path = 'Figures_PREEMPT/Map_Traps/CV_Folds/')

## ## Check between fold distances
## for(fold in unique(folded.data.baf$cvfold)){
##     ## Check distances
##     dist.mat <- distm(unique(folded.data.baf[folded.data.baf$cvfold == fold,c('Longitude', 'Latitude')]),
##                  unique(folded.data.baf[folded.data.baf$cvfold != fold,c('Longitude', 'Latitude')])
##                  )
##     print(min(dist.mat))
## }

## --- Fit Bafodia trap data with high resolution predictor of 50m
focal.variables <- c('Frac_bare', 'Frac_grass', 'Frac_tree', 
                     'Frac_burn', 'Frac_rice', 'Frac_water', 
                     'Frac_mound', 'Frac_garden', 
                     'Frac_modern_build', 'Frac_trad_build',
                     'Density_Buildings', 'Density_Moderns',
                     'Density_Traditionals')
brt.names <- get.features(c(1000, 200, 100,50), focal.variables)$names

## --- Set equal weights for Mna+ and Mna-
absent.weights <- subset(folded.data.baf, Mna==0 , weight)
present.weights <- subset(folded.data.baf, Mna==1 , weight)
adj.absent <- sum(present.weights) / sum(absent.weights)
folded.data.baf$adj.weight = ifelse(folded.data.baf$Mna==1, 1,
                                 adj.absent*folded.data.baf$weight)

gbm.mod.r <- NULL
gbm.mod.r <- mgbm.step(data = folded.data.baf, gbm.y = 'Mna',
                       gbm.x = brt.names,
                       site.weights = folded.data.baf$adj.weight,
                       learning.rate = 0.0001, max.trees = Inf,
                       fold.vector = folded.data.baf$cvfold,
                       n.folds = length(unique(folded.data.baf$cvfold)),
                       keep.fold.vector = TRUE, n.trees = 1, step.size = 1,
                             prev.stratify = FALSE, tree.complexity = 1,
                       bag.fraction = 0.5)



folded.data.all <- add.cv.folds(train.data = subset(comb.wide.dataset),
                            bloocv.r = 50)
## --- Set equal weights for Mna+ and Mna-
absent.weights <- subset(folded.data.all, Mna==0 , weight)
present.weights <- subset(folded.data.all, Mna==1 , weight)
adj.absent <- sum(present.weights) / sum(absent.weights)
folded.data.all$adj.weight = ifelse(folded.data.all$Mna==1, 1,
                                 adj.absent*folded.data.all$weight)

gbm.mod.r <- NULL
gbm.mod.r <- mgbm.step(data = folded.data.all, gbm.y = 'Mna',
                       gbm.x = brt.names,
                       site.weights = folded.data.all$adj.weight,
                       learning.rate = 0.0001, max.trees = Inf,
                       fold.vector = folded.data.all$cvfold,
                       n.folds = length(unique(folded.data.all$cvfold)),
                       keep.fold.vector = TRUE, n.trees = 1, step.size = 1,
                             prev.stratify = FALSE, tree.complexity = 1,
                       bag.fraction = 0.5)


## Neither of these model attempts dramatically reduces the mean deviance. 



## --- Plot villages with traps

library(png)
library(plotrix)

villages = unique(comb.wide.dataset$Site)

agg.wide.data <- data.table::setDT(comb.wide.dataset)[,
                                                 .(Mna = sum(Mna) > 0),
                                                 by = list(Longitude, Latitude,
                                                           Visit, Date, Site)]    

## Set parameters for image grid with units of meters
bs = 100 ## distance in meters
mlat <- mean(comb.wide.dataset$Latitude)
dx = bs / (111111*cos(pi/180*mlat))
dy = bs / 111111

for(village in villages){
    village.rast <- brick(paste0('../SkyLens_v3/Data/Rasters/Google_Raster_Images/',village,'.tif'))

    ## Build simple plot of where Mna captures occurred for each visit; 
    ## Aggregate over nights

    focal.wide.data <- subset(agg.wide.data, Site==village)
    
    visit.date <- unique(focal.wide.data[,c('Visit','Date')])
    visit.date <- visit.date[order(visit.date$Visit),]

    nvisits = max(focal.wide.data$Visit)
    for(vi in 1:nvisits){    
        png(file = paste0('Figures_PREEMPT/Map_Traps/Visits/',village,'_traps_v',vi,'.png'),
            width = 6, height = 6, units = 'in', res = 300)
        par(mai = c(0,0,0,0))
        plotRGB(village.rast)
        points(subset(focal.wide.data, Visit==vi , c(Longitude,Latitude)),
               pch = 1, cex = 0.25, 
               col = c('black', 'red')[focal.wide.data$Mna + 1])
        xseq = seq(village.rast@extent@xmin, village.rast@extent@xmax, by = dx)
        yseq = seq(village.rast@extent@ymin, village.rast@extent@ymax, by = dy)
        abline(v = xseq, lwd = 0.25)
        abline(h = yseq, lwd = 0.25)
        legend(x = 'topleft', legend = format(as.Date(visit.date$Date[vi]), format="%b-%Y"),
               adj = 0.15)
        dev.off()
    }

    ## --- Same, but aggregate across all visits
    agg.focal.wide.data <- data.table::setDT(focal.wide.data)[,
                                                   .(Mna = sum(Mna) > 0),
                                                   by = list(Longitude, Latitude)]    
    png(file = paste0('Figures_PREEMPT/Map_Traps/Aggregate/',village,'_traps.png'),
        width = 6, height = 6, units = 'in', res = 300)
    par(mai = c(0,0,1,0))
    plotRGB(village.rast)
    abline(v = xseq, lwd = 0.25)
    abline(h = yseq, lwd = 0.25)
    points(subset(agg.focal.wide.data, ,c(Longitude,Latitude)),
           pch = 1, cex = 0.25, 
           col = c('black', 'red')[agg.focal.wide.data$Mna + 1])
    dev.off()

    writeLines(paste0('Finished: ', village))
} ## End loop through villages



## --- Add in spatial block to pre.dataset, quantify lc diversity
sp.dataset <- add.spatial.blocks(pre.dataset)

for(rr in c(25,50,100,200, 500, 1000)){
    sp.dataset[,paste0('sh', rr)] <- NA
    features <- get.features(c(rr), lc.sat.vars)$names
    sp.dataset[,paste0('sh', rr)] <- -rowSums(sp.dataset[,..features, with = FALSE]*
                                 log(sp.dataset[,..features, with = FALSE]),
                                 na.rm = TRUE)
}



temp = sp.dataset[,.(Total = length(unique(Sp[!is.na(Sp)])),
                     Inv = length(unique(Sp[!is.na(Sp) & (Sp%in%c('Mmu','Rra'))])),
                     NInv = length(unique(Sp[!is.na(Sp) & !(Sp%in%c('Mmu','Rra'))])),
                     Mn = length(unique(Sp[!is.na(Sp) & (Sp%in%'Mna')])),
                     Rodents = length(unique(Sp[!is.na(Sp) & !(Sp%in%not.murid)])),
                     Shrews =  length(unique(Sp[!is.na(Sp) & (Sp%in%not.murid)])),
                     sh50 = mean(sh50),
                     sh100 = mean(sh100)
                     ),
           by = list(Site, block.i)]


plot(Total~sh100, temp)
plot(Rodents~sh100, temp)
plot(Shrews~sh100, temp)
plot(jitter(Inv)~sh100, temp)
plot(jitter(NInv)~sh100, temp)


temp = sp.dataset[,.(Total = length(unique(Sp[!is.na(Sp)])),
                     Inv = length(unique(Sp[!is.na(Sp) & (Sp%in%c('Mmu','Rra'))])),
                     NInv = length(unique(Sp[!is.na(Sp) & !(Sp%in%c('Mmu','Rra'))])),
                     Mn = length(unique(Sp[!is.na(Sp) & (Sp%in%'Mna')])),
                     Rodents = length(unique(Sp[!is.na(Sp) & !(Sp%in%not.murid)])),
                     Shrews =  length(unique(Sp[!is.na(Sp) & (Sp%in%not.murid)])),
                     sh25 = mean(sh25),
                     sh50 = mean(sh50),
                     sh100 = mean(sh100),
                     sh200 = mean(sh200),
                     sh500 = mean(sh500),
                     sh1000 = mean(sh1000)
                     ),
           by = list(Site)]

plot(Total~sh1000, temp)
plot(Rodents~sh1000, temp)
plot(Shrews~sh1000, temp)
plot(jitter(Inv)~sh1000, temp)
plot(jitter(NInv)~sh1000, temp)


temp <- sp.dataset[,.(Mna = 1*(sum(Sp %in% c('Mna')) > 0),
              Rra = 1*(sum(Sp %in% c('Rra')) > 0)),
           by = list(Site, block.i)]

mod = glm(Mna ~ Rra, temp, family = quasibinomial)
