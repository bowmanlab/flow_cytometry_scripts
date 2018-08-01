setwd('~/bowman_lab/cyflow_space/ecuador_2018')

#### general parameters ####

output <- 'ecuador_2018_AF'
plot.fcm <- function(name, fcm.dataframe, x, y){
  hist2d(log10(fcm.dataframe[,x]),
         log10(fcm.dataframe[,y]),
         col = c('grey', colorRampPalette(c('white', 'lightgoldenrod1', 'darkgreen'))(100)),
         nbins = 100,
         bg = 'grey',
         ylab = paste0('log10(', y, ')'),
         xlab = paste0('log10(', x, ')'),
         main = name,
         cex.main = 0.8)
  
  box()
}

#### QC parameters ####

f.list <- list.files(path = '.', pattern = '*AF.fcs', ignore.case = T)
FSC.limit <- 5
SSC.limit <- 1
FL6.limit <- 1

#### aggregation and QC ####

library(gplots)

training.events <- data.frame(FSC = numeric(),
                                 SSC = numeric(),
                                 FL2 = numeric(),
                                 FL3 = numeric(),
                                 FL4 = numeric(),
                                 FL6 = numeric())

library('flowCore')

#f.name <- f.list[7]
#grep('blank', f.list)

sample.size <- 1000 # size to sample from each for training data

pdf(paste0(output, '_fcm_plots.pdf'),
    width = 5,
    height = 5)

for(f.name in f.list){
  print(f.name)
  
  fcm <- read.FCS(f.name)
  fcm <- as.data.frame((exprs(fcm)))
  fcm$FL5 <- NULL # FL5 channel does not contain biologically useful information
  
  fcm[fcm == 0] <- NA
  fcm <- na.omit(fcm)
  
  plot.fcm(f.name, fcm, 'SSC', 'FSC')
  plot.fcm(f.name, fcm, 'SSC', 'FL6')
  plot.fcm(f.name, fcm, 'FL2', 'FL6')
  plot.fcm(f.name, fcm, 'FL4', 'FL6')

  fcm <- fcm[fcm$FSC > FSC.limit,]
  fcm <- fcm[fcm$SSC > SSC.limit,]
  fcm <- fcm[fcm$FL6 > FL6.limit,]

  plot.fcm(paste(f.name, 'QC'), fcm, 'SSC', 'FL6')

  ## blanks and other very clean samples may not have enough points to donate
  ## to training dataset
  
  try({fcm.sample <- fcm[sample(1:length(fcm[,1]), sample.size),]}, silent = T)
  
  training.events<- rbind(training.events, fcm.sample)
  
  write.csv(fcm, sub('FCS', 'qc.csv', f.name), quote = F)
}

dev.off()

write.csv(training.events, paste0(output, '.training_events.csv'), quote = F)

#### training ####

library(kohonen)
library(oce)
library(gplots)

train.fcm <- function(event.file, params){

  events <- read.csv(event.file, stringsAsFactors = F)

  plot.fcm('training', events, 'SSC', 'FSC')
  plot.fcm('training', events, 'SSC', 'FL6')
  plot.fcm('training', events, 'FL2', 'FL6')
  plot.fcm('training', events, 'FL4', 'FL6')

  sample.mat <- as.matrix(events[,params])
  colnames(sample.mat) <- params

  sample.mat <- log10(sample.mat)

  grid.size <- ceiling(dim(events)[1] ^ (1/2.5))

  som.grid <- somgrid(xdim = grid.size, ydim = grid.size, topo="hexagonal")

  som.model <- som(sample.mat, 
                   grid = som.grid, 
                   rlen = 100,
                   alpha = c(0.05,0.01),
                   keep.data = TRUE)
  return(som.model)
}

event.file <- paste0(output, '.training_events.csv')
params <- c('SSC', 'FSC', 'FL2', 'FL3', 'FL4', 'FL6')

pdf(paste0(output, '.som_mode_training_plots.pdf'),
    width = 5,
    height = 5)

som.model <- train.fcm(event.file, params)

dev.off()

## use kmeans and determine best number of clusters
## http://www.shanelynn.ie/self-organising-maps-for-customer-segmentation-using-r/

## if you want to work with previous models, load
#load('ecuador_2018.af.som.Rdata')

som.events <- som.model$codes[[1]]

wss <- rep(NA, 30)

for (i in 2:30) {
  wss[i] <- sum(kmeans(som.events, centers=i, iter.max = 20)$withinss)
}

plot(wss,
     pch = 19,
     ylab = 'Within-clusters sum of squares',
     xlab = 'K')

## pick the elbow

library(vegan)
library(pmclust)

## K-mean clustering currently providing much more sensible clusters

k <- 10
#som.cluster.sg <- kmeans(som.events.sg, centers = k.sg)$cluster
som.cluster <- pmclust(som.events, K = k)$class
#som.dist <- vegdist(som.events)
#som.cluster <- cutree(hclust(som.dist), k = k)

### populations in the training data ###

library(oce)

flow.col <- oce.colorsFreesurface(k)

pdf(paste0(output, '.param_plots.pdf'), width = 5, height = 5)

for(param in colnames(som.model$data[[1]])){
  
  plot(som.model$data[[1]][,param], som.model$data[[1]][,'FL6'],
       type = 'n',
       xlab = param,
       ylab = 'FL6')
  
  for(p in 1:k){
    i <- which(som.cluster[som.model$unit.classif] == p)
    points(som.model$data[[1]][i,param], som.model$data[[1]][i,'FL6'],
           col = flow.col[p],
           pch = 19,
           cex = 0.4)
  }
}

dev.off()

### save the models ###

save(list = c('som.model', 'som.cluster', 'k'), file = paste0(output, '.som.Rdata'))

### som property plots ###

som.property.plot <- function(som.model, som.cluster, property, title){
  plot(som.model, type = 'property', property = property, main = title)
  add.cluster.boundaries(som.model, som.cluster, lwd = 2)
}

pdf(paste0(output, '.properties.pdf'),
    width = 6,
    height = 5)

plot(som.model, type = 'counts')
add.cluster.boundaries(som.model, som.cluster)

som.property.plot(som.model, som.cluster, som.events[,1], 'log10(FSC)')
som.property.plot(som.model, som.cluster, som.events[,2], 'log10(SSC)')
som.property.plot(som.model, som.cluster, som.events[,3], 'log10(FL2)')
som.property.plot(som.model, som.cluster, som.events[,4], 'log10(FL3)')
som.property.plot(som.model, som.cluster, som.events[,4], 'log10(FL4)')
som.property.plot(som.model, som.cluster, som.events[,4], 'log10(FL6)')

plot(som.model,
     type = "mapping",
     property = som.cluster,
     main = 'Cluster locations',
     #main = NULL,
     bgcol = flow.col[som.cluster],
     col = NA)

dev.off()

### pca on nodes ###

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