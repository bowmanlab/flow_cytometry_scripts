setwd("~/bowman_lab/mosaic/utqiagvik/fcm")

#### general parameters ####

output <- 'test_AF'                      # identifier for output files
data.path <- './'                     # make sure this ends with "/"
f.list <- list.files(path = data.path,
                     pattern = '*AF.*fcs',
                     ignore.case = T)      # list of fcs files to analyze

## Define a general plotting function for fcm data,
## takes as input plot title, dataframe, dataframe of bead events
## (can be NULL if not needed), x and y parameters.

plot.fcm <- function(name, fcm.dataframe, beads, x, y){
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
SSC.llimit <- -0.7
FL3.llimit <- 0
FL6.llimit <- 0

## Lower limits for the key parameters used to
## define beads (assumes log10 scale).

SSC.beads.llimit <- 2
FSC.beads.llimit <- 2
FL5.beads.llimit <- 1.4

#### aggregation and QC ####

library(hexbin)
library('flowCore')

training.events <- data.frame(FSC = numeric(),
                              SSC = numeric(),
                              FL1 = numeric(),
                              FL2 = numeric(),
                              FL3 = numeric(),
                              FL4 = numeric(),
                              FL5 = numeric(),
                              FL6 = numeric()) # a dataframe to hold a selection of data for training the model

sample.size <- 200 # size to sample from each for training data

## Iterate across all FCS files, performing QC, making plots,
## and taking a random selection of QC'd data for training.

#f.name <- f.list[1]
#grep('blank', f.list)

pdf(paste0(output, '_fcm_plots.pdf'),
    width = 5,
    height = 5)

for(f.name in f.list){
  print(f.name)
  
  fcm <- read.FCS(paste0(data.path, f.name))
  fcm <- as.data.frame((exprs(fcm)))
  
  fcm[fcm == 0] <- NA
  fcm <- na.omit(fcm)
  
  ## Identify beads.  If you don't know where the beads are start with an empty beads dataframe and
  ## make some plots (SSC, FL5) to identify.
  
  fcm.beads <- data.frame()
  fcm.beads <- fcm[which(log10(fcm$SSC) > SSC.beads.llimit &
                           log10(fcm$FL5) > FL5.beads.llimit &
                           log10(fcm$FSC) > FSC.beads.llimit),]
  
  ## Make plots of all events.
  
  plot.fcm(f.name, fcm, fcm.beads, 'FSC', 'SSC')
  plot.fcm(f.name, fcm, fcm.beads, 'FSC', 'FL3')
  plot.fcm(f.name, fcm, fcm.beads, 'FSC', 'FL6')
  plot.fcm(f.name, fcm, fcm.beads, 'FL6', 'FL2')
  plot.fcm(f.name, fcm, fcm.beads, 'FL6', 'FL3')
  plot.fcm(f.name, fcm, fcm.beads, 'FL6', 'FL4')

  ## Remove events that are below limits (thresholds).
  
  fcm <- fcm[log10(fcm$FSC) > FSC.llimit,]
  fcm <- fcm[log10(fcm$SSC) > SSC.llimit,]
  fcm <- fcm[log10(fcm$FL3) > FL3.llimit,]
  fcm <- fcm[log10(fcm$FL6) > FL6.llimit,]
  
  ## Make plots of only those events remaining.
  
  plot.fcm(paste(f.name, 'QC'), fcm, fcm.beads, 'FSC', 'SSC')
  plot.fcm(paste(f.name, 'QC'), fcm, fcm.beads, 'FSC', 'FL3')
  plot.fcm(paste(f.name, 'QC'), fcm, fcm.beads, 'FSC', 'FL5')
  plot.fcm(paste(f.name, 'QC'), fcm, fcm.beads, 'FSC', 'FL6')

  
  ## blanks and other very clean samples may not have enough points to donate
  ## to training dataset
  
  try({fcm.sample <- fcm[sample(1:length(fcm[,1]), sample.size),]}, silent = T)
  
  training.events<- rbind(training.events, fcm.sample[,colnames(training.events)])
  
  write.csv(fcm, sub('FCS', 'qc.csv', f.name), quote = F)
}

dev.off()

write.csv(training.events, paste0(output, '.training_events.csv'), quote = F)

#### model training and selection ####

library(kohonen)
library(oce)
library(hexbin)

## Define a function to train a SOM for select parameters
## from training data selected randomly during QC.

train.fcm <- function(event.file, params){

  events <- read.csv(event.file, stringsAsFactors = F)

  plot.fcm('training', events, NULL, 'FSC', 'SSC')
  plot.fcm('training', events, NULL, 'FSC', 'FL6')
  plot.fcm('training', events, NULL, 'FSC', 'FL5')
  plot.fcm('training', events, NULL, 'FL6', 'FL2')
  plot.fcm('training', events, NULL, 'FL6', 'FL3')
  plot.fcm('training', events, NULL, 'FL6', 'FL4')

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
params <- c('FSC', 'SSC', 'FL2', 'FL3', 'FL4', 'FL5', 'FL6')

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

#load('test_AF.som.Rdata')

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

k <- 5

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

  plot(som.model$data[[1]][,param], som.model$data[[1]][,'FL6'],
       type = 'n',
       xlab = param,
       ylab = 'FL6',
       main = paste0(alg, ', ', 'k = ', j))
  
  for(p in 1:j){
    i <- which(som.cluster[som.model$unit.classif] == p)
    points(som.model$data[[1]][i,param], som.model$data[[1]][i,'FL6'],
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
  
  for(param in params){
    plot.clusters('pmclust', param, som.model, som.cluster.pm, flow.col, j)
  }
  
  som.property.plot(som.model, som.cluster.pm, som.events[,'FSC'], paste0('log10(FSC),', 'pmclust, k =', j))
  som.property.plot(som.model, som.cluster.pm, som.events[,'FL6'], paste0('log10(FL6),', 'pmclust, k =', j))
  
  plot(som.model,
       type = "mapping",
       property = som.cluster.pm,
       main = paste0('Cluster locations,', 'pmclust, k =', j),
       bgcol = flow.col[som.cluster.pm],
       col = NA)
  
  ## Plots for k-means.
  
  for(param in params){
    plot.clusters('kmeans', param, som.model, som.cluster.k, flow.col, j)
  }
  
  som.property.plot(som.model, som.cluster.k, som.events[,'FSC'], paste0('log10(FSC),', 'kmeans, k =', j))
  som.property.plot(som.model, som.cluster.k, som.events[,'FL6'], paste0('log10(FL6),', 'kmeans, k =', j))
  
  plot(som.model,
       type = "mapping",
       property = som.cluster.k,
       main = paste0('Cluster locations,', 'kmeans, k =', j),
       bgcol = flow.col[som.cluster.k],
       col = NA)
  
  ## Plots for hierarchical.
  
  for(param in params){
    plot.clusters('hierarchical', param, som.model, som.cluster.h, flow.col, j)
  }
  
  som.property.plot(som.model, som.cluster.h, som.events[,'FSC'], paste0('log10(FSC),', 'hclust, k =', j))
  som.property.plot(som.model, som.cluster.h, som.events[,'FL6'], paste0('log10(FL6),', 'hclust, k =', j))
  
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

k <- 6
cluster.method <- 'k'
som.cluster <- cluster.tries[[paste('som.cluster', cluster.method, k, sep = '.')]]
cluster.notes <- paste(cluster.method, 'k=', k)
save(list = c('som.model', 'som.cluster', 'k', 'cluster.notes'), file = paste0(output, '.som.Rdata'))

## Refine cluster plots, if needed.

flow.col <- oce.colorsFreesurface(k)

plot.clusters('kmeans', 'FSC', som.model, som.cluster, flow.col, k)
plot.clusters('kmeans', 'SSC', som.model, som.cluster, flow.col, k)

## pca on nodes ##

node.pca <- prcomp(som.model$codes[[1]])

pdf(paste0(output, '.node_pca.pdf'),
    width = 5,
    height = 5)

plot(node.pca$x[,1], node.pca$x[,2], type = 'n',
     xlab = 'PC1',
     ylab = 'PC2')

for(p in 1:k){
  points(node.pca$x[which(som.cluster == p),1], node.pca$x[which(som.cluster == p),2],
         pch = 19,
         col = flow.col[p],
         cex = 0.6)
}

scaling.factor = 3

arrows(0,
       0,
       node.pca$rotation[,1] * scaling.factor,
       node.pca$rotation[,2] * scaling.factor,
       lwd = 2)

text(node.pca$rotation[,1] * scaling.factor,
     node.pca$rotation[,2] * scaling.factor,
     labels = rownames(node.pca$rotation),
     adj = c(1,1),
     font = 2)

dev.off()