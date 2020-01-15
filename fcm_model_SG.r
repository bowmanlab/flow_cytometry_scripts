setwd("~/bowman_lab/bowman_lab_github/flow_cytometry_scripts")

#### general parameters ####

output <- 'test_SG'                      # identifier for output files
data.path <- './'                     # make sure this ends with "/"
f.list <- list.files(path = data.path,
                     pattern = '*fcs',
                     ignore.case = F)      # list of fcs files to analyze

## Define a general plotting function for fcm data,
## takes as input plot title, dataframe, dataframe of bead events
## (can be NULL if not needed), x and y parameters.

plot.fcm <- function(name, fcm.dataframe, beads=NA, x='SSC-HLin', y='GRN-B-HLin'){
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

FSC.llimit <- -0.9
SSC.llimit <- 0.1
FL1.llimit <- 0

## Lower limits for the key parameters used to
## define beads (assumes log10 scale).

SSC.beads.llimit <- 3.1
FSC.beads.llimit <- 4
FL5.beads.llimit <- 3.3

#### aggregation and QC ####

library(hexbin)
library('flowCore')

training.events <- data.frame(`FSC-HLin` = numeric(),
                              `SSC-HLin` = numeric(),
                              `GRN-B-HLin` = numeric(),
                              `BLU-V-HLin` = numeric()) # a dataframe to hold a selection of data for training the model

colnames(training.events) <- c("FSC-HLin", "SSC-HLin", "GRN-B-HLin", "BLU-V-HLin")

sample.size <- 1000 # size to sample from each for training data

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
    print(f.name)
    
    ## If it's the last record in the file, stop loop after this iteration.
    
    if(as.integer(fcm@description$`$NEXTDATA`) == 0){
      analyze.next <- F}
    
    ## Convert data to dataframe.
    
    fcm <- as.data.frame(fcm@exprs)
    fcm <- fcm[,grep('HLin', colnames(fcm))]
    fcm[fcm < 1] <- NA
    fcm <- na.omit(fcm)
    
    ## Identify beads.  If you don't know where the beads are start with an empty beads dataframe and
    ## make some plots (SSC, FL5) to identify.
    
    fcm.beads <- data.frame()
    fcm.beads <- fcm[which(log10(fcm$`SSC-HLin`) > SSC.beads.llimit &
                             log10(fcm$`BLU-V-HLin`) > FL5.beads.llimit &
                             log10(fcm$`FSC-HLin`) > FSC.beads.llimit),]
    
    ## Make plots of all events.
    
    plot.fcm(f.name, fcm, fcm.beads, y = "FSC-HLin")
    plot.fcm(f.name, fcm, fcm.beads)
  
    ## Remove events that are below limits (thresholds).
    
    fcm <- fcm[log10(fcm$`FSC-HLin`) > FSC.llimit,]
    fcm <- fcm[log10(fcm$`SSC-HLin`) > SSC.llimit,]
    fcm <- fcm[log10(fcm$`GRN-B-HLin`) > FL1.llimit,]
    
    ## Make plots of only those events remaining.
    
    plot.fcm(paste(f.name, 'QC'), fcm, fcm.beads, y = 'FSC-HLin')
    plot.fcm(paste(f.name, 'QC'), fcm, fcm.beads, y = 'GRN-B-HLin')
    plot.fcm(paste(f.name, 'QC'), fcm, fcm.beads, y = 'RED-V-HLin')
    
    ## blanks and other very clean samples may not have enough points to donate
    ## to training dataset
    
    try({
      fcm.sample <- fcm[sample(1:length(fcm[,1]), sample.size),]
      training.events<- rbind(training.events, fcm.sample[,colnames(training.events)])
      }, silent = T)
    
    write.csv(fcm, paste0(f.name, '.qc.csv'), quote = F)
    }
  
}

dev.off()

write.csv(training.events, paste0(output, '.training_events.csv'), quote = F)

#### model training and selection ####

library(kohonen)
library(oce)
library(hexbin)

training.events <- read.csv(paste0(output, '.training_events.csv'), row.names = 1)

## Define a function to train a SOM for select parameters
## from training data selected randomly during QC.

train.fcm <- function(event.file, params){

  events <- read.csv(event.file, stringsAsFactors = F)

  plot.fcm('training', events, NULL, x = 'SSC.HLin', y = 'FSC.HLin')
  plot.fcm('training', events, NULL, x = 'SSC.HLin', y = 'GRN.B.HLin')
  plot.fcm('training', events, NULL, x = 'SSC.HLin', y = 'BLU.V.HLin')

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
params <- c('FSC.HLin', 'SSC.HLin', 'GRN.B.HLin')

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

k <- 6

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
library(pmclust)
library(oce)

## Define a function to make a basic scatterplot, with events color-coded
## by cluster.

plot.clusters <- function(alg, param, som.model, som.cluster, flow.col, j){

  plot(som.model$data[[1]][,param], som.model$data[[1]][,'GRN.B.HLin'],
       type = 'n',
       xlab = param,
       ylab = 'GRN.B.HLin',
       main = paste(alg, ',', 'k =', j))
  
  for(p in 1:j){
    i <- which(som.cluster[som.model$unit.classif] == p)
    points(som.model$data[[1]][i,param], som.model$data[[1]][i,'GRN.B.HLin'],
           col = flow.col[p],
           pch = 19,
           cex = 0.4)
  }
  
  legend('topleft',
         legend = paste('Cluster', 1:j),
         pch = 19,
         col = flow.col)
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

  som.cluster.pm <- pmclust(som.events, K = j, algorithm = 'apecma')$class # model based
  som.cluster.k <- kmeans(som.events, centers = j, iter.max = 100, nstart = 10)$cluster # k-means
  som.dist <- vegdist(som.events) # hierarchical, step 1
  som.cluster.h <- cutree(hclust(som.dist), k = j) # hierarchical, step 2
  
  cluster.tries[[paste0('som.cluster.pm.', j)]] <- som.cluster.pm
  cluster.tries[[paste0('som.cluster.k.', j)]] <- som.cluster.k
  cluster.tries[[paste0('som.cluster.h.', j)]] <- som.cluster.h

  flow.col <- oce.colorsFreesurface(j)
  
  ## Plots for pm.
  
  plot.clusters('pmclust', 'SSC.HLin', som.model, som.cluster.pm, flow.col, j)
  som.property.plot(som.model, som.cluster.pm, som.events[,1], paste0('log10(SSC.HLin),', 'pmclust, k =', j))
  som.property.plot(som.model, som.cluster.pm, som.events[,3], paste0('log10(GRN.B.HLin),', 'pmclust, k =', j))
  
  plot(som.model,
       type = "mapping",
       property = som.cluster.pm,
       main = paste0('Cluster locations,', 'pmclust, k =', j),
       bgcol = flow.col[som.cluster.pm],
       col = NA)
  
  ## Plots for k-means.
  
  plot.clusters('kmeans', 'SSC.HLin', som.model, som.cluster.k, flow.col, j)
  som.property.plot(som.model, som.cluster.k, som.events[,1], paste0('log10(SSC.HLin),', 'kmeans, k =', j))
  som.property.plot(som.model, som.cluster.k, som.events[,3], paste0('log10(GRN.B.HLin),', 'kmeans, k =', j))
  
  plot(som.model,
       type = "mapping",
       property = som.cluster.k,
       main = paste0('Cluster locations,', 'kmeans, k =', j),
       bgcol = flow.col[som.cluster.k],
       col = NA)
  
  ## Plots for hierarchical.
  
  plot.clusters('hierarchical', 'SSC.HLin', som.model, som.cluster.h, flow.col, j)
  som.property.plot(som.model, som.cluster.h, som.events[,1], paste0('log10(SSC.HLin),', 'hclust, k =', j))
  som.property.plot(som.model, som.cluster.h, som.events[,3], paste0('log10(GRN.B.HLin),', 'hclust, k =', j))
  
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
cluster.method <- 'k' #either k, pm, or h
som.cluster <- cluster.tries[[paste('som.cluster', cluster.method, k, sep = '.')]]
cluster.notes <- paste(cluster.method, 'k=', k)
save(list = c('som.model', 'som.cluster', 'k', 'cluster.notes'), file = paste0(output, '.som.Rdata'))

## Refine cluster plots, if needed.

flow.col <- oce.colorsFreesurface(k)

plot.clusters('kmeans', 'FSC.HLin', som.model, som.cluster, flow.col, k)
plot.clusters('kmeans', 'SSC.HLin', som.model, som.cluster, flow.col, k)