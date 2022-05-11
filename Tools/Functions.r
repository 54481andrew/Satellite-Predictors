
## Function to calculate accuracy
calc.accuracy <- function(obs, pred, weights = rep(1, length(obs))){
    binary.pred = 1*(pred >= 0.5)
    return(1 - mean(binary.pred==obs))
}

## Function that computes Moran's I for a few spatial
## radii
compute.Moran <- function(y.name, block.name, df, dseq){
    y.mean <- mean(df[,y.name])
    df <- data.table::setDT(df)
    agg.df <- df[,.(Longitude = mean(Longitude),
                    Latitude = mean(Latitude),
                    y.mean = mean(get(y.name))
                    ),
                 by = list(get(block.name))]

    ## Compute distance matrix between blocks
    dist.mat <- distm(agg.df[,c('Longitude', 'Latitude')])
    moran.df <- data.frame(dist = dseq, i = NA, k = NA)
    nds <- length(dseq)
    
    ## Compute weights matrix
    for(di in 1:nds){
        d1 <- dseq[di]
        d2 <- dseq[di + 1]
        
        weights.mat <- dist.mat >= d1 & dist.mat < d2
        ntot <- nrow(weights.mat)
        diag(weights.mat) <- 0

    
        ## Convert binary weight matrix to neighbor list
        listw <- mat2listw(weights.mat)
        
        ## Compute Moran's I
        result <- moran(x = unlist(agg.df[,get]), listw = listw, n = nrow(weights.mat),
                        S0 = sum(weights.mat), zero.policy = TRUE)
        

        moran.df$i[di] <- result$i
        moran.df$k[di] <- result$k

    }## End loop through dseq
    return(result)
} ## End function




### *****************************

## Helper functions used in the various scrips. 

### *****************************

## Omit any names in pred.names that contain
## one of nix
omit.names <- function(pred.names, nix){
    new.names <- c()
    for(n in nix){
        keep.names <- which(!grepl(n, pred.names, ignore.case = TRUE))
        new.names <- c(new.names, pred.names[keep.names])
    }
    return(new.names)
}



get.radii <- function(names){
    return(
        as.numeric(sapply(strsplit(names, split = '[.]'), FUN = function(x){x[length(x)]}))
    )
}


## Find the bin of seq that point x is in
bin <- function(x, seq){
    xbin = max(which(x >= seq))    
    return(xbin)
}


## Add in principal components to dataframe
add.pcs <- function(df, pred.names, thresh = 0.9){

    prout <- prcomp(wide.dataset[,pred.names], retx=TRUE, center=TRUE, scale.=TRUE)
    
    ## Calculate explained variance
    cs <- cumsum(prout$sdev^2)
    scree.out <- cs / cs[length(cs)]

    include = sum(scree.out < thresh)

    df.new <- cbind(df, prout$x[,1:include])

    form <- paste('Mna~', paste(paste('PC', 1:include, sep = ''), collapse = ' + '))

    return(list(df.new, num = include, scree = scree.out, formula = form, rot = prout$rotation))
}


## Block the space of each village and assign folds. bs is block
## side length in meters.
add.spatial.blocks <- function(train.data, bs = 250){
    
    mlat <- mean(train.data$Latitude)
    dx = bs / (111111*cos(pi/180*mlat))
    dy = bs / 111111
    
    villages <- unique(train.data$Site)
    nvillages <- length(villages)
    fudge <- 1e-6 ## Prevent min / max points from falling on edge
    max.block = 0
    if('block.i' %in% names(train.data)){train.data <- subset(train.data, , -block.i)}
    train.data$block.i <- 0
    for(vi in 1:nvillages){
        mask <- train.data$Site==villages[vi]
        subdat <- train.data[mask,]
        min.lon <- min(subdat$Longitude) - fudge
        max.lon <- max(subdat$Longitude) + fudge
        min.lat <- min(subdat$Latitude) - fudge
        max.lat <- max(subdat$Latitude) + fudge
        
        lon.bins <- seq(min.lon, max.lon, by = dx)
        lat.bins <- seq(min.lat, max.lat, by = dy) 
        
        nlons = length(lon.bins) ## Number of bins in the x direction
        nlats = length(lat.bins) ## Number of bins in the y direction
        
        x.bins = sapply(subdat$Longitude, function(x){bin(x, lon.bins)})
        y.bins = sapply(subdat$Latitude, function(x){bin(x, lat.bins)})
        
        ## Assign folds
        blocks.i = x.bins + x.bins*(y.bins - 1) + max.block
        train.data[mask,'block.i'] = blocks.i
        
        ## Update block to keep block numbers between villages unique
        max.block <- max(blocks.i)
    }
    return(train.data)
} ## End function    

## Alternative to the function above. Only includes folds of
## traps that are bloocv.r distance (meters) away from any
## other fold. 
add.cv.folds <- function(train.data, bloocv.r = 50,
                         plot.path = NA){
    
    potential.traps <- train.data
    set.seed(42)
    fold.dataset <- c()
    ii = 1
    while(nrow(potential.traps) > 0){
        ## Choose trap
        focal.trap <- sample(1:nrow(potential.traps), size = 1)
        focal.trap.dat <- potential.traps[focal.trap,c('Longitude', 'Latitude')]
        
        ## Choose all traps within bloocv_r meters of the focal trap as the fold vector
        trap.distances <- distm(potential.traps[focal.trap,c('Longitude', 'Latitude')],
                                potential.traps[, c('Longitude', 'Latitude')])
        wi.traps <- which(trap.distances < bloocv.r)
        
        fold.traps <- potential.traps[wi.traps,]
        fold.traps[,'cvfold'] <- ii
        fold.dataset <- rbind(fold.dataset, fold.traps)
        
        ## remove those traps from potential traps    
        rem.traps <- potential.traps[wi.traps,]
        potential.traps <- potential.traps[-wi.traps,]
        
        ## Also remove any traps that overlap (are within 50m of the) folds traps
        ufold.traps <- unique(fold.traps[,c('Longitude', 'Latitude')])
        trap.distances <- distm(ufold.traps[,c('Longitude', 'Latitude')],
                                potential.traps[, c('Longitude', 'Latitude')])
        overlap.traps <- which(colSums(trap.distances < bloocv.r) > 0)
        
        overlapping.traps <- potential.traps[overlap.traps,]
        potential.traps <- potential.traps[-overlap.traps,]
        
        if(!is.na(plot.path)){
            dir.create(plot.path, recursive = TRUE, showWarnings = FALSE)
            png(file = paste0(plot.path, '/test_validation', ii, '.png'), width = 6, height = 6,
                units = 'in', res = 300)
            par(mai = c(0,0,0,0))
            plotRGB(village.rast)
            points(fold.dataset[,c('Longitude','Latitude')], pch = 1, cex = 0.5, 
                   col = 'black')
            points(fold.traps[,c('Longitude','Latitude')], pch = 1, cex = 0.5, 
                   col = 'orange')
            points(overlapping.traps[,c('Longitude','Latitude')], pch = 1, cex = 0.5,
                   col = 'green')
            points(focal.trap.dat[c('Longitude','Latitude')], pch = 1, cex = 1,
                   col = 'red')
            dev.off()
        }
        ii = ii + 1        
        print(ii)
    } ## end while loop

    fold.dataset
    return(fold.dataset)
} ## End function

## Create a weight matrix that weights presences and absences equally
## in each fold, for both train and validation subsets
## - adj.val: whether the validation
## - equal.within: 
create.weight.matrix <- function(train.data, fold.matrix, weight.col = 'weight',
                                 equal.within = TRUE){
                                 
    n.folds <- ncol(fold.matrix)
    weight.matrix <- matrix(train.data[,weight.col], nrow = nrow(fold.matrix), ncol = n.folds,
                            byrow = FALSE)
    
    if(equal.within){
        
        for(vi in 1:n.folds){
            ## Recall fold.matrix is encoded so that column vi describes the
            ## vi'th fold: validation set is masked by value vi, train is masked
            ## by vi + 1. Data not included is masked by NA. 
            mask.fold = fold.matrix[,vi] %in% c(vi, vi + 1)
            mask.foldpres = train.data$Mna==1 & mask.fold & !is.na(mask.fold)
            mask.foldabs = train.data$Mna==0 & mask.fold & !is.na(mask.fold)

            w.pres <- sum(train.data[mask.foldpres,weight.col])

            w.abs <- sum(train.data[mask.foldabs,weight.col])
            w.adj <- w.pres/w.abs
            weight.matrix[mask.foldabs,vi] <- w.adj*train.data[mask.foldabs,weight.col]

        } ## vi
    }## End if equal.weights

    return(weight.matrix)
}


## Alternative to the function above. Only includes folds of
## traps that are bloocv.r distance (meters) away from any
## other fold. 
create.fold.matrix <- function(train.data, n.folds = 10, r.min = 0, r.max = 50,
                               include.pres = TRUE, plot.path = NA,
                               fold.style = 'c', which.data = 'a',
                               verbose = FALSE){


    
    ## Buffered CV style fold
    if(fold.style=='c'){
        ## Each fold is constructed by choosing a point, then saving all points
        ## within an r.min [m] radius for evaluation, and any point outside of
        ## r.max [m] for training.
        
        ## include.pres: include at least 1 presence in each fold

        ## Store info on how the data are folded
        fold.matrix = matrix(NA, nrow = nrow(train.data), ncol = n.folds)

        train.data$trap.id <- 1:nrow(train.data)
        
        ## Track which traps are allowed to be fold focal traps (center of fold)
        if(include.pres){
            focal.traps <- which(train.data$Mna==1)
        }else{
            focal.traps <- 1:nrow(train.data)
        }
        
        ## Tracks which traps can be placed in a new CV fold. We'll make sure
        ## that a trap can only belong to one evaluation fold. 
        potential.traps <- 1:nrow(train.data)
        
        set.seed(42)
        fold.i = 0
        eval.traps = c()
        if(verbose){
            print(paste0('fold: ', NA,  ' rem:', length(potential.traps),
                         ' len:', NA))
        }
        while(fold.i < n.folds & length(potential.traps) > 0 & length(focal.traps) > 0){
            
            fold.i = fold.i + 1
            ## Choose trap
            ## Behavior of "sample" changes when focal.traps only has one
            ## element; handle that case separately. 
            if(length(focal.traps)>1){
                focal.trap <- sample(focal.traps, size = 1,
                                     prob = train.data$weight[focal.traps])
            }else{
                focal.trap = focal.traps
            }
            
            focal.trap.dat <- train.data[focal.trap, c('Longitude', 'Latitude')]
            
            ## Compute distances between traps and focal trap
            trap.distances <- distm(train.data[focal.trap, c('Longitude', 'Latitude')],
                                    train.data[, c('Longitude', 'Latitude')])
            
            ## Choose all traps within r.min meters of the focal trap as evaluation
            ## -- importantly, check that we are only choosing NEW validation traps
            wi.eval.traps <- which(trap.distances <= r.min &
                                   (train.data$trap.id %in% potential.traps))
            fold.matrix[wi.eval.traps, fold.i] <- fold.i
            
            ## traps that are outside the buffer distance get a non-NA fold number (!= fold.i)
            wi.train.traps <- which(trap.distances >= r.max)
            fold.matrix[wi.train.traps, fold.i] <- fold.i + 1

            ## Remove evaluation traps from further consideration
            focal.traps = focal.traps[-which(focal.traps %in% wi.eval.traps)]
            potential.traps = potential.traps[-which(potential.traps %in% wi.eval.traps)]

            if(verbose){
                print(paste0('fold: ', fold.i, ' rem:', length(potential.traps),
                             ' len:', length(wi.eval.traps)))
            }
            
            if(!is.na(plot.path)){
                
                village.name = unique(train.data[wi.eval.traps,'Site'])
                print(paste0('Focal trap village: ', village.name))
                village.rast <- brick(paste0('../SkyLens_v3/Data/Rasters/Google_Raster_Images/',
                                             village.name,'.tif'))
                dir.create(plot.path, recursive = TRUE, showWarnings = FALSE)
                png(file = paste0(plot.path, '/test_validation', fold.i, '.png'), width = 6, height = 6,
                    units = 'in', res = 300)
                par(mai = c(0,0,0,0))
                plotRGB(village.rast)
                points(train.data[fold.matrix[,fold.i]==fold.i, c('Longitude','Latitude')],
                       pch = 1, cex = 0.5, 
                       col = 'orange')
                points(train.data[fold.matrix[,fold.i]==(fold.i + 1), c('Longitude','Latitude')],
                       pch = 1, cex = 0.5, 
                       col = 'black')
                points(train.data[is.na(fold.matrix[,fold.i]), c('Longitude','Latitude')],
                       pch = 1, cex = 0.5,
                       col = 'green')
                points(train.data[focal.trap,c('Longitude','Latitude')],
                       pch = 1, cex = 0.5, 
                       col = ifelse(train.data[focal.trap,'Mna']==1, 'red', 'blue'))
                dev.off()
                
            }
            
        } ## end while loop
        
        ## Check if all folds were created
        if(fold.i < n.folds){
            n.folds = fold.i
            fold.matrix <- fold.matrix[,1:(n.folds)]
            print('WARNING: not all folds created')
        }        
    }## End fold.style=='c'

    ## Site fold    
    if(fold.style=='v'){
        ## Default is fold over all villages
        village.groups <- unique(train.data$Site)
        if(n.folds==4){
            village.groups <- list(c('Bafodia'), c('Guala', 'Makump', 'Badala', 'Talama'),
                                   c('Largo', 'Makuna', 'Njaguima'),
                                c('Barlie', 'Benduma'), 
                                c('Gbenikoro','Kapethe'))
        }
        if(n.folds==3){
            village.groups <- list(c('Bafodia'), c('Guala', 'Makump', 'Badala', 'Talama',
                                                   'Largo', 'Makuna', 'Njaguima'),
                                   c('Barlie', 'Benduma', 
                                     'Gbenikoro','Kapethe'))
        }
        if(n.folds==2){
            village.groups <- list(c('Bafodia'), c('Guala', 'Makump', 'Badala', 'Talama',
                                                   'Largo', 'Makuna', 'Njaguima',
                                                   'Barlie', 'Benduma', 
                                                     'Gbenikoro','Kapethe'))
        }
        
        nvillages <- length(village.groups)
        fold.matrix = matrix(NA, nrow = nrow(train.data), ncol = nvillages)
        for(vi in 1:nvillages){
            fold.matrix[,vi] = ifelse(train.data$Site %in% village.groups[vi][[1]],
                                      vi, vi+1)
        }
        n.folds <- ncol(fold.matrix)
    }

    
    ## Store simple stats on folds here
    stat.matrix = matrix(NA, nrow = 6, ncol = n.folds)

    ## Get simple stats on train and evaluation sets
    for(ii in 1:n.folds){
        col = fold.matrix[,ii]
        eval <- col==ii & !is.na(col)
        e.vals = train.data[eval,'Mna']
        train <- (col==ii+1) & !is.na(col)
        t.vals = train.data[train,'Mna']
        
        stat.matrix[,ii] <- c(sum(e.vals==0), sum(e.vals==1), sum(eval),
                              sum(t.vals==0), sum(t.vals==1), sum(train))
    }
    
    return(list(fold.matrix, stat.matrix))
} ## End function


## Function that bootstrap-samples indices from
## multiverse.EFC.dat. Indices of multiverse.EFC.dat
## are returned, with corresponding rows describing
## trap ID and (randomly sampled) Mna capture status. 
multiverse.subset <- function(single.EFC.copy, mask){

    ## Create list to store indices that describe the outcome of the sampled focal traps
    multiverse.index <- c()    

    ## Only proceed if any focal traps exist to sample
    if(sum(mask) > 0){
        ## Mask out values of interest
        mask.single.EFC <- single.EFC.copy[mask,]
        
        ## Create an array of masked indices for focal elements in single.EFC
        ind.array <- which(mask)
        
        ## Find total number of trap locations x nights to be sampled
        num.trap.nights <- mask.single.EFC$Nights ## Number of trap nights for each location
        tot.traps <- nrow(mask.single.EFC) ## Number of trap locations
        max.nights <- max(mask.single.EFC$Nights) 
        
        ## Sample trap outcome for all locations, for trap nights 1 to max.nights. Then,
        ## I'll set to NA trap locations that did not occur. 
        
        for(ni in 1:max.nights){
            ## For each trap location, use uniform number to determine trap outcome
            probs <- runif(tot.traps)
            Mna.status <- 1*(probs < mask.single.EFC$probMna)
            ## Set outcome to NA if this trap-night did not occur for this trap location
            Mna.status <- ifelse(num.trap.nights >= ni, Mna.status, NA)
            ## Save indices
            multiverse.index <- c(multiverse.index,
                                  ifelse(Mna.status==1, 2*ind.array - 1, 2*ind.array))
        }
        ## Remove away NA indices
        multiverse.index <- multiverse.index[!is.na(multiverse.index)]
    }## End if checking if focal traps exist in mask input

    return(multiverse.index)
}
    
## Analogous function of create.fold.matrix2 that is able to
## create village folds in EFC and PREEMPT datasets. Currently only
## supports fold across village name. 
create.fold.matrix2 <- function(train.data, single.EFC, fold.style = 'v'){
   
    pre.dat <- train.data[train.data$Source=='PRE',]
    efc.dat <- train.data[train.data$Source=='EFC',]
    
    ## Site fold
    if(fold.style=='v'){
        village.groups <- unique(train.data$Site)
        n.folds <- length(village.groups)
        fold.matrix = matrix(NA, nrow = nrow(train.data), ncol = n.folds)
        for(vi in 1:n.folds){

            ## Get index information for EFC folds (validation)
            mask.in <- single.EFC$Site %in% village.groups[vi][[1]]
            efc.dummy <- rep(NA, nrow(efc.dat))
            efc.in <- multiverse.subset(single.EFC, mask.in) 

            ## Get index information for EFC folds (train)
            efc.out <- multiverse.subset(single.EFC, !mask.in) 

            ## Get index information for PREEMPT folds (train and validation)
            fold.matrix[1:nrow(pre.dat),vi] = ifelse(pre.dat$Site %in% village.groups[vi][[1]],
                                      vi, vi+1)
            fold.matrix[nrow(pre.dat) + efc.in,vi] = vi
            fold.matrix[nrow(pre.dat) + efc.out,vi] = vi + 1
            
        }
        n.folds <- ncol(fold.matrix)
    }
    
    ## Store simple stats on folds here
    stat.matrix = matrix(NA, nrow = 6, ncol = n.folds)

    ## Get simple stats on train and evaluation sets
    for(ii in 1:n.folds){
        col = fold.matrix[,ii]
        eval <- col==ii & !is.na(col)
        e.vals = train.data[eval,'Mna']
        train <- (col==ii+1) & !is.na(col)
        t.vals = train.data[train,'Mna']
        
        stat.matrix[,ii] <- c(sum(e.vals==0), sum(e.vals==1), sum(eval),
                              sum(t.vals==0), sum(t.vals==1), sum(train))
    }
    
    return(list(fold.matrix, stat.matrix))
} ## End function


## Given a binary response vector yi,
## compute a weights vector so that 0's and 1's are equally weighted
correct.weights <- function(yi, weights){
    absent.weights <- weights[yi==0]
    present.weights <- weights[yi==1]
    adj.absent <- sum(present.weights) / sum(absent.weights)
    new.weights <- rep(1, length(weights))
    new.weights[yi==0] <- adj.absent*absent.weights
    return(new.weights)
}




## Return base feature names with radii pasted
## and also GLM formula
get.features <- function(rvec, names){
    feature.list <- c()
    for(r in rvec){
        feature.list <- c(feature.list, paste(names, r, sep = '.'))
    }
    form <- paste('Mna~', paste(feature.list, collapse = ' + '))
    return(list(formula = form, names = feature.list))
}


## Normalize the predictors in preparation for GLM
create.norm.df <- function(in.df, use.names){
    norm.df = in.df
    for(rr in c(50,100,200,1000)){
        temp.names <- get.features(rr, use.names)$names
        norm.df[,temp.names] = norm.preds(in.df,
                                          temp.names)[[1]]
    }
    return(norm.df)
}


## Build function to normalize over dataset
## - could do with "scale", but I want more control over NA values
norm.preds <- function(df, foc.vars){
    feature.mean <- colMeans(df[,foc.vars], na.rm = TRUE)
    feature.sd <- apply(df[,foc.vars], 2, sd, na.rm = TRUE)
    mean.matrix <- matrix(feature.mean,
                          nrow = nrow(df),
                          ncol = length(feature.mean),
                          byrow = TRUE)
    centered.df <- df[,foc.vars] - mean.matrix
    sd.matrix <- matrix(feature.sd,
                        nrow = nrow(df),
                        ncol = length(feature.mean),
                        byrow = TRUE)
    scaled.df <- centered.df/sd.matrix
    return(list(scaled.df, feature.mean, feature.sd))
}



## Binarize the data
binarize <- function(data, count, total, new.name, num = TRUE){
    classi.dat.bin <- c()
    for(ii in 1:nrow(data)){
        df <- data[ii,]
        if(!is.na(data[ii,total])){
            if(data[ii,total] > 0){
                tot <- df[,total]
                pos <- df[,count]
                neg <- tot - pos
                df.neg <- df[rep(seq_len(nrow(df)), each = neg),]
                df.neg[,new.name] <- rep('F',neg)
                df.pos <- df[rep(seq_len(nrow(df)), each = pos),]
                df.pos[,new.name] <- rep('S',pos)
                test.dat <- rbind(df.neg, df.pos)
                classi.dat.bin <- rbind(classi.dat.bin, test.dat)
            }
        }
    }
    if(!num){
        classi.dat.bin[,new.name] = as.factor(classi.dat.bin[,new.name])
    }else{
        classi.dat.bin[,new.name] = 1*(classi.dat.bin[,new.name]=='S')
    }
    return(classi.dat.bin)
} ## End function

## Use principal components info to calculate z-scores of dataset
Calc.ZScores <- function(dataset, pr.dat, npcs){
    hab.center.mat <- matrix(pr.dat$center, byrow = TRUE, ncol = length(pr.dat$center),
                             nrow = nrow(dataset))
    hab.scale.mat <- matrix(pr.dat$scale, byrow = TRUE, ncol = length(pr.dat$scale),
                            nrow = nrow(dataset))
    scaled.data <- (dataset[,names(pr.dat$center)] - hab.center.mat)/hab.scale.mat
    dataset[,paste0('PC',1:npcs)] <- (as.matrix(scaled.data) %*% hab.pc.mat)[,1:npcs]
    return(dataset)
}


## Apply elastic net regression model to predictor stack -- not used for BRT
glmnet.raster <- function(coef.dat, predictor.stack, rain.preds, verbose = FALSE){
    model.sum <- predictor.stack[[1]]*0
    names(model.sum) <- 'Out'
    terms <- subset(coef.dat, select = name, subset = val != 0)
    terms <- as.vector(unlist(terms))
    ## Fudge to make the custom.overlay function work correctly
    dummy.var <- rain.preds[which(rain.preds %in% coef.dat$name)[1]]
    for(term in terms){
        input.raster <- predictor.stack
        eval(parse(text = paste("custom.overlay <- function(",
                                paste(names(input.raster), collapse = ", "),
                                ") {",
                                term, '*(0*', dummy.var,' + 1)',
                                "}"
                                )))
        term.out <- overlay(predictor.stack, fun = custom.overlay, unstack = TRUE)
        model.sum <- model.sum +
            as.numeric(subset(coef.dat, subset = name==term, select = val))*term.out
    if(verbose){print(term); print(c(minValue(model.sum),maxValue(model.sum), model.sum[230,330]))}
    }
    if(nlayers(term.out) > 1){print("ERROR")}
    return(model.sum)
}

## Calculate deviance and likelihood of model fit
dev.fun <- function (data.prop, pred.prop, data.weight){
    LL.sat <- sum(dbinom(x = data.prop*data.weight,
                     size = data.weight, prob = data.prop,
                     log = TRUE))
    LL.mod <- sum(dbinom(x = data.prop*data.weight,
                     size = data.weight, prob = pred.prop,
                     log = TRUE))
    cv.dev <- 2*(LL.sat - LL.mod)
    return(c(cv.dev, LL.mod, LL.sat))
}


## legend html generator:
markerLegendHTML <- function(IconSet) {   
    ## container div:
    legendHtml <- "<div style='padding: 10px; padding-bottom: 10px;'><h4 style='padding-top:0; padding-bottom:10px; margin: 0;'> Marker Legend </h4>"
    n <- 1
    ## add each icon for font-awesome icons icons:
    for (Icon in IconSet) {
        if (Icon[["library"]] != "fa") {
            legendHtml<- paste0(legendHtml, "<div style='width: auto; height: 45px'>",
                                "<div style='position: relative; display: inline-block; width: 36px; height: 45px' class='awesome-marker-icon-", Icon[["markerColor"]][1], " awesome-marker'>",
                                "<i style='margin-left: 8px; margin-top: 11px; 'class= 'fa fa-", Icon[["icon"]], " fa-inverse'></i>",
                                "</div>",
                                "<p style='position: relative; top: 10px; display: inline-block; ' >", names(IconSet)[n] ,"</p>",
                                "</div>")    
        }
        n<- n + 1
    }
    paste0(legendHtml, "</div>")
}


## Return interpretable predictor names for plots
pretty.labels <- function(x){
    x.out <- c()
    for(xi in x){
        xi.name <- switch(xi,
                          TreeCover = 'Tree Cover',
                          ShrubCover = 'Shrub Cover',
                          Grasslands = 'Grassland',
                          Cropland = 'Cropland',
                          AquaticVeg = 'Aquatic Veg.',
                          SparseVeg = 'Sparse Veg.',
                          Bare = 'Bare',
                          BuiltUp = 'Built Up',
                          SnowIce = 'Snow and Ice',
                          OpenWater = 'Open Water',
                          Pmu = 'Mean Precipitation',
                          Tmu = 'Mean Temperature',
                          Nmu = 'NDVI',
                          Pmin = 'Minimum Precipitation',
                          Pmax = 'Maximum Precipitation',
                          Nmin = 'Minimum NDVI',
                          Nmax = 'Maximum NDVI',
                          Pcv = 'Precipitation Coeficient of Variation',
                          Ncv = 'NDVI Coef. of Variation',
                          Pc = 'Precipitation Constancy',
                          Pm = 'Precipitation Contingency',
                          Nc = 'NDVI Constancy',
                          Nm = 'NDVI Contingency',
                          Pdur = 'Dry Duration',
                          Ndur = 'Brown Duration',
                          Elev = 'Elevation',
                          Pop = 'Population',
                          TotCrop = 'Croplands',
                          Evergreen_Needleleaf_Forest = 'Evergreen Needle Forest',
                          Evergreen_Broadleaf_Forest = 'Evergreen Broadleaf Forest',
                          Deciduous_Needleleaf_Forest = 'Deciduous Needle Forest',
                          Deciduous_Broadleaf_Forest = 'Deciduous Broadleaf Forest',
                          Mixed_Forest= 'Mixed Forest',
                          Closed_Shrubland = 'Closed Shrubland',
                          Open_Shrubland = 'Open Shrubland',
                          Woody_Savanna = 'Woody Savanna',
                          Savannas = 'Savanna',
                          Grasslands = 'Grassland',
                          Permanent_Wetlands = 'Wetland',
                          Croplands = 'Cropland',
                          Urban_BuiltUp = 'Urban',
                          Cropland_Natural_Mosaic = 'Cropland Mosaic',
                          Permanent_Snow_Ice = 'Snow/Ice',
                          Barren = 'Barren',
                          Water_Bodies = 'Water',
                          Unclassified = 'NA')
        if(is.null(xi.name)){xi.name = xi}
        x.out <- c(x.out, xi.name)
    }
    return(x.out)
}

## Function that creates a confusion matrix from binary data
## rows 1 and 2 correspond to predictions of types (e.g., p1, p2),
## cols 1 and 2 correspond to true types (e.g., t1, t2)
confusion.matrix <- function(true_vec, pred_vec, types = c(0,1)){
    conf.mat <- matrix(NA, nrow = 2, ncol = 2)
    len.types <- length(types)
    for(i in 1:len.types){
        for(j in 1:len.types){
            conf.mat[i,j] <- sum(true_vec==types[j] & pred_vec==types[i])
        }
    }
    row.names(conf.mat) <- paste0('p',types)
    colnames(conf.mat) <- paste0('t',types)
    return(conf.mat)
}

## Returns a list of precision, recall, and F1 score
prf1 <- function(true_vec, pred_vec){
    conf.mat <- confusion.matrix(true_vec, pred_vec)
    tp <- conf.mat[2,2]
    fp <- conf.mat[2,1]
    tn <- conf.mat[1,1]
    fn <- conf.mat[1,2]
    if((tp + fp) == 0){
        prec <- 0
        rec <- 0
        f1 <- 0
    }else{
        prec <- tp/(tp + fp)
        rec <- tp/(tp + fn)
        f1 <- 2* (prec*rec) / (prec + rec)
    }
    return(list(tp = tp, fp = fp, tn = tn, fn = fn, accuracy = (tp + tn) / (length(true_vec)),
                precision = prec, recall = rec, f1 = f1))
}
