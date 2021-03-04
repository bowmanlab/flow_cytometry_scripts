setwd("~/bowman_lab/sandbox")

#### general parameters ####

stain = 'SG' # indicate AF or SG
output <- 'SCCOOS_SG'                   # identifier for output files
data.path <- './'                     # make sure this ends with "/"
f.list <- list.files(path = data.path,
                     pattern = '*fcs',
                     ignore.case = F)      # list of fcs files to analyze
					 
## what channels should be used?

if(stain == 'AF'){
  params <- c("FSC-HLin", "SSC-HLin", "BLU-V-HLin", "YEL-B-HLin", "RED-V-HLin", "RED-B-HLin")
}

if(stain == 'SG'){
  params <- c("FSC-HLin", "SSC-HLin", "BLU-V-HLin", "GRN-B-HLin")
}

#### functions ####

## Define a general plotting function for fcm data,
## takes as input plot title, dataframe, dataframe of bead events
## (can be NULL if not needed), x and y parameters.

plot.fcm <- function(name, fcm.dataframe, beads=NA, x='FSC-HLin', y='RED-V-HLin'){
  fcm.hex <- hexbin(log10(fcm.dataframe[,x]), log10(fcm.dataframe[,y]), 100)
  plot(fcm.hex@xcm,
       fcm.hex@ycm,
       col = BTC(100)[as.numeric(cut(fcm.hex@count, 100))],
       ylab = paste0('log10(', y, ')'),
       xlab = paste0('log10(', x, ')'),
       main = name,
       pch = 19,
       cex = 0.4)
  try({
  rect(c(min(log10(beads[,x])), min(log10(beads[,x]))),
       c(min(log10(beads[,y])), min(log10(beads[,y]))),
       c(max(log10(beads[,x])), max(log10(beads[,x]))),
       c(max(log10(beads[,y])), max(log10(beads[,y]))))
  }, silent = T)

}

#### QC parameters ####

## Lower limits for the key parameters that will be used for
## QC (assumes log10 scale).

FSC.llimit <- -1
SSC.llimit <- -1
RED.V.HLin.llimit <- -1 # only for AF
GRN.B.HLin.llimit <- 1.4 # only for SG

## Lower limits for the key parameters used to
## define beads (assumes log10 scale).

SSC.beads.llimit <- 3.3
FSC.beads.llimit <- 3.8
FL5.beads.llimit <- 3.3

#### aggregation and QC ####

library(hexbin)
library('flowCore')

if(stain == 'AF'){

  training.events <- data.frame(`FSC-HLin` = numeric(),
                                `SSC-HLin` = numeric(),
                                `BLU-V-HLin` = numeric(),
                                `YEL-B-HLin` = numeric(),
                                `RED-V-HLin` = numeric(),
                                `RED-B-HLin` = numeric()) # a dataframe to hold a selection of data for training the model
  
  colnames(training.events) <- c("FSC-HLin", "SSC-HLin", "BLU-V-HLin", "YEL-B-HLin", "RED-V-HLin", "RED-B-HLin")
}

if(stain == 'SG'){
  training.events <- data.frame(`FSC-HLin` = numeric(),
                                `SSC-HLin` = numeric(),
                                `BLU-V-HLin` = numeric(),
                                `GRN-B-HLin` = numeric()) # a dataframe to hold a selection of data for training the model
  
  colnames(training.events) <- c("FSC-HLin", "SSC-HLin", "BLU-V-HLin", "GRN-B-HLin")
}

## I've kept the sample size in here for now, but much better to select 1-2 representative
## samples for the model.  You need to avoid too many points that "fill in" the
## space between natural populations of cells.

## For future development new samples should be evaluated against the old model,
## and the model could be rebuilt if there are populations that are not included.

#sample.size <- 6000 # size to sample from each for training data
training.samples <- c('SG') # names of samples to add to training data.  Avoid adding too many.

## Iterate across all FCS files, performing QC, making plots,
## and taking a random selection of QC'd data for training.

#f.name <- f.list[25]
#grep('blank', f.list)

pdf(paste0(output, '_fcm_plots.pdf'),
    width = 5,
    height = 5)

for(fcs in c(f.list)){
  
  analyze.next <- T
  i <- 0

  while(analyze.next == T){
    i <- i + 1
    
    fcm <- read.FCS(paste0(data.path, fcs), emptyValue = F, dataset = i)
    f.name <- fcm@description$`GTI$SAMPLEID`
    
    f.name <- sub('SYBR', 'SG', f.name)
    
    ## If it's the last record in the file, stop loop after this iteration.
    
    if(as.integer(fcm@description$`$NEXTDATA`) == 0){
      analyze.next <- F}
    
    ## If analysis name has the stain in it, proceed.
    
    if(grepl(stain, f.name)){
      print(f.name)
      
      ## Convert data to dataframe.
      
      fcm <- as.data.frame(fcm@exprs)
      fcm <- fcm[,grep('HLin', colnames(fcm))]
	  
      ## Drop all events with a negative value in any channel.  Save these in case you need to look at.
      
      fcm.bad <- fcm[apply(fcm, 1, function(row) any(row < 0)),]
      
      for(param in params){
        fcm[fcm[param] <= 0,param] <- NA
      }
      
      fcm <- na.omit(fcm)
      
      ## Identify beads.  If you don't know where the beads are start with an empty beads dataframe and
      ## make some plots (SSC, FL5) to identify.
      
      fcm.beads <- data.frame()
      fcm.beads <- fcm[which(log10(fcm$`SSC-HLin`) > SSC.beads.llimit &
                               log10(fcm$`BLU-V-HLin`) > FL5.beads.llimit &
                               log10(fcm$`FSC-HLin`) > FSC.beads.llimit),]
      
      ## Make plots of all events.
      
      if(stain == 'AF'){
        plot.fcm(f.name, fcm, fcm.beads, x='FSC-HLin', y='RED-V-HLin')
        plot.fcm(f.name, fcm, fcm.beads, x = 'YEL-B-HLin', y = 'RED-V-HLin')
      }
      
      if(stain == 'SG'){
        plot.fcm(f.name, fcm, fcm.beads, x='FSC-HLin', y='GRN-B-HLin')
      }
    
      ## Remove events that are below limits (thresholds).
      
      fcm <- fcm[log10(fcm$`FSC-HLin`) > FSC.llimit,]
      fcm <- fcm[log10(fcm$`SSC-HLin`) > SSC.llimit,]
      
      if(stain == 'AF'){
        fcm <- fcm[log10(fcm$`RED-V-HLin`) > RED.V.HLin.llimit,]
      }
      
      if(stain == 'SG'){
        fcm <- fcm[log10(fcm$`GRN-B-HLin`) > GRN.B.HLin.llimit,]
      }
      
      ## Make plots of only those events remaining.
      
      if(stain == 'AF'){
        plot.fcm(paste(f.name, 'QC'), fcm, fcm.beads, x='FSC-HLin', y='RED-V-HLin')
        plot.fcm(paste(f.name, 'QC'), fcm, fcm.beads, x = 'YEL-B-HLin', y = 'RED-V-HLin')
      }
      
      if(stain == 'SG'){
        plot.fcm(paste(f.name, 'QC'), fcm, fcm.beads, x='FSC-HLin', y='GRN-B-HLin')
      }
      
      ## blanks and other very clean samples may not have enough points to donate
      ## to training dataset.  Currently turned off as random selection does not
      ## seem to be the best way to select the training set.
      
      # try({
      #   fcm.sample <- fcm[sample(1:length(fcm[,1]), sample.size),]
      #   training.events<- rbind(training.events, fcm.sample[,colnames(training.events)])
      #   }, silent = T)
      
      ## Check to see if sample should be added to training data.
      
      if(f.name %in% training.samples){
        training.events<- rbind(training.events, fcm[,colnames(training.events)])
      }
      
      write.csv(fcm, paste0(f.name, '.qc.csv'), quote = F)
    }
  }
}

dev.off()

write.csv(training.events, paste0(output, '.training_events.csv'), quote = F)

#### model training and selection ####

library(kohonen)
library(oce)
library(hexbin)

## Define a function to train a SOM for select parameters
## from training data selected during QC.

train.fcm <- function(event.file, params){

  events <- read.csv(event.file, stringsAsFactors = F)
  events <- events[sample.int(dim(events)[1], size = 5000),]

  for(param in params){
    plot.fcm('training', events, NULL, x = 'FSC.HLin', y = param)
  }

  sample.mat <- as.matrix(events[,params])
  colnames(sample.mat) <- params

  sample.mat <- log10(sample.mat)

  grid.size <- ceiling(dim(events)[1] ^ (1/2.5))

  som.grid <- somgrid(xdim = grid.size, ydim = grid.size, topo="hexagonal", toroidal = T)

  som.model <- som(sample.mat, 
                   grid = som.grid, 
                   rlen = 100,
                   alpha = c(0.05,0.01),
                   keep.data = TRUE)
  return(som.model)
}

## Execute the SOM function for select parameters.

event.file <- paste0(output, '.training_events.csv')

params <- gsub("-", ".", params)

pdf(paste0(output, '.som_model_training_plots.pdf'),
    width = 5,
    height = 5)

som.model <- train.fcm(event.file, params)

dev.off()

## Use kmeans and determine best number of clusters, following
## http://www.shanelynn.ie/self-organising-maps-for-customer-segmentation-using-r/.
## If you're really sure about the number of clusters you expect, skip this part and
## define as k.

## If you've saved a model earlier, and want to work on the clustering load the old model.

#load('sccoos_SG.som.Rdata')

## Guestimate the number of clusters based on within clusters sum of squares using k-means.

som.events <- som.model$codes[[1]]

wss <- rep(NA, 30)

for (i in 2:30) {
  wss[i] <- sum(kmeans(som.events, centers=i, iter.max = 20)$withinss)
}

plot(wss,
     pch = 19,
     ylab = 'Within-clusters sum of squares',
     xlab = 'K')

## Pick the elbow, this is the starting point for number of clusters.

k <- 4

## Now you need to manually evaluate some different clustering techniques,
## including k-means, model based, and hierarchical. In practice
## K-mean clustering has produced the most sensible clusters.

## Now use the classification model to identify populations in the
## training data.  You will want to evaluate a range of options for k
## (based on estimate provded by SoS analysis), and multiple clustering
## algorithms until the result "looks right".

## In practice kmeans clustering has consistently produced the most coherent
## clusters, however, you need to evaluate this for your own data!

library(vegan)
library(oce)

## Define a function to make a basic scatterplot, with events color-coded
## by cluster.

plot.clusters <- function(alg, paramx, paramy, som.model, som.cluster, flow.col, j){

  plot(som.model$data[[1]][,paramx], som.model$data[[1]][,paramy],
       type = 'n',
       xlab = paramx,
       ylab = paramy,
       main = paste(alg, ',', 'k =', j))
  
  for(p in 1:j){
    i <- which(som.cluster[som.model$unit.classif] == p)
    points(som.model$data[[1]][i,paramx], som.model$data[[1]][i,paramy],
           col = flow.col[p],
           pch = 19,
           cex = 0.4)
  }
  
  legend('topleft',
         legend = paste('Cluster', 1:j),
         pch = 19,
         col = flow.col,
         bg = 'white')
}

## Define a function to make SOM property plots.

som.property.plot <- function(som.model, som.cluster, property, title){
  plot(som.model, type = 'property', property = property, main = title)
  add.cluster.boundaries(som.model, som.cluster, lwd = 2)
}

## Generate a bunch of plots to guide selection of appropriate clustering
## algorithm.

cluster.tries <- list()

pdf(paste0(output, '.cluster_eval.pdf'), width = 5, height = 5)

for(j in (k-2):(k+2)){
  
  if(stain == 'AF'){
    plot.param.1 = 'YEL.B.HLin'
    plot.param.2 = 'RED.V.HLin'
  }
  
  if(stain == 'SG'){
    plot.param.1 = 'FSC.HLin'
    plot.param.2 = 'GRN.B.HLin'
  }

  som.cluster.k <- kmeans(som.events, centers = j, iter.max = 100, nstart = 10)$cluster # k-means
  som.dist <- vegdist(som.events) # hierarchical, step 1
  som.cluster.h <- cutree(hclust(som.dist), k = j) # hierarchical, step 2
  
  cluster.tries[[paste0('som.cluster.k.', j)]] <- som.cluster.k
  cluster.tries[[paste0('som.cluster.h.', j)]] <- som.cluster.h

  flow.col <- oce.colorsFreesurface(j)
  
  ## Plots for k-means.
  
  plot.clusters('kmeans', plot.param.1, plot.param.2, som.model, som.cluster.k, flow.col, j)
  som.property.plot(som.model, som.cluster.k, som.events[,1], paste0('log10(', plot.param.1, '),', 'kmeans, k =', j))
  som.property.plot(som.model, som.cluster.k, som.events[,3], paste0('log10(', plot.param.2, '),', 'kmeans, k =', j))
  
  plot(som.model,
       type = "mapping",
       property = som.cluster.k,
       main = paste0('Cluster locations,', 'kmeans, k =', j),
       bgcol = flow.col[som.cluster.k],
       col = NA)
  
  ## Plots for hierarchical.
  
  plot.clusters('hierarchical', plot.param.1, plot.param.2, som.model, som.cluster.h, flow.col, j)
  som.property.plot(som.model, som.cluster.h, som.events[,1], paste0('log10(', plot.param.1, '),', 'hclust, k =', j))
  som.property.plot(som.model, som.cluster.h, som.events[,3], paste0('log10(', plot.param.2, '),', 'hclust, k =', j))
  
  plot(som.model,
       type = "mapping",
       property = som.cluster.h,
       main = paste0('Cluster locations,', 'hclust, k =', j),
       bgcol = flow.col[som.cluster.h],
       col = NA)
}
  
dev.off()

#### describe selected model ####

## Select the clustering algorithm that you like best and final number of 
## clusters, and save the model.

k <- 4
make.time <- format(Sys.time(), "%Y%m%d_%H%M%S")
cluster.method <- 'k' #either k, pm, or h
som.cluster <- cluster.tries[[paste('som.cluster', cluster.method, k, sep = '.')]]
cluster.notes <- paste(cluster.method, 'k=', k)
save(list = c('som.model', 'som.cluster', 'k', 'cluster.notes'),
     file = paste0(output, '.', make.time, '.som.Rdata'))

## Refine cluster plots, if needed.

flow.col <- oce.colorsFreesurface(k)

pdf(paste0(output, '.', make.time, '.final_clusters.pdf'), width = 5, height = 5)

plot.clusters('kmeans', 'FSC.HLin', 'RED.V.HLin', som.model, som.cluster, flow.col, k)
plot.clusters('kmeans', 'SSC.HLin', 'RED.V.HLin', som.model, som.cluster, flow.col, k)
plot.clusters('kmeans', 'BLU.V.HLin', 'RED.V.HLin', som.model, som.cluster, flow.col, k)
plot.clusters('kmeans', 'YEL.B.HLin', 'RED.V.HLin', som.model, som.cluster, flow.col, k)
plot.clusters('kmeans', 'RED.B.HLin', 'RED.V.HLin', som.model, som.cluster, flow.col, k)
plot.clusters('kmeans', 'BLU.V.HLin', 'RED.V.HLin', som.model, som.cluster, flow.col, k)

plot.clusters('kmeans', 'FSC.HLin', 'GRN.B.HLin', som.model, som.cluster, flow.col, k)
plot.clusters('kmeans', 'FSC.HLin', 'SSC.HLin', som.model, som.cluster, flow.col, k)
plot.clusters('kmeans', 'BLU.V.HLin', 'GRN.B.HLin', som.model, som.cluster, flow.col, k)

dev.off()

## ignore

plot(fcm$`RED-V-HLin` ~ fcm$`YEL-B-HLin`, log = 'xy')