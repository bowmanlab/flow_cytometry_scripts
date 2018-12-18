setwd('~/sccoos_fcm')

#### parameters ####

f.list <- list.files(path = '.', pattern = '_SG.qc.csv', ignore.case = T)
output <- 'test.sg'
input <- 'test_SG.som.Rdata'
label <- 'sg'
paramx <- 'FSC' # SSC best indicator for size, pg. 422 "In Living Color"?
paramy <- 'FL1'

#### classify ####

library(kohonen)
library(oce)
library(plotrix)

load(input)

cluster.tally <- matrix(nrow = length(f.list), ncol = k)
colnames(cluster.tally) <- 1:k
row.names(cluster.tally) <- f.list
flow.col <- oce.colorsFreesurface(k)

## Set variables for testing only.

paramx <- 'FSC'
paramy <- 'FL1'
label <- 'sg'
sample <- f.list[1]
som.model <- som.model
cluster.tally.df <- cluster.tally
cluster.vector <- som.cluster

## Classify events from all samples, making a single compiled pdf
## showing size and layout of clusters for each sample.

classify.fcm <- function(sample, som.model, cluster.vector, paramx, paramy, label, flow.col, k){
  
  print(sample)
  
  params <- colnames(som.model$data[[1]])
  
  sample.df <- read.csv(sample, row.names = 1, header = T)
  sample.mat <- as.matrix(sample.df[,params])
  sample.mat <- log10(sample.mat)
  
  sample.predict <- predict(som.model, sample.mat)
  sample.df[paste0('cluster.', label)] <- cluster.vector[sample.predict$unit.classif]
  
  out <- vector(length = k)
  
  ## regular plot
  
  plot(sample.mat[,paramy] ~ sample.mat[,paramx],
       type = 'n',
       main = sample,
       xlab = paramx,
       ylab = paramy)
  
  for(cluster in 1:k){
    r = which(sample.df[paste0('cluster.', label)] == cluster)
    out[cluster] <- length(r)

    temp.sd.x <- sd(sample.mat[r, paramx])
    temp.sd.y <- sd(sample.mat[r, paramy])
    
    temp.mean.x <- mean(sample.mat[r, paramx])
    temp.mean.y <- mean(sample.mat[r, paramy])
    
    draw.ellipse(temp.mean.x, temp.mean.y, a = temp.sd.x, b = temp.sd.y, border = flow.col[cluster])
    
    text(temp.mean.x, temp.mean.y, length(r))
    
  }
  
  legend('topleft',
         legend = paste('Cluster', 1:k),
         pch = 1,
         col = flow.col)
  
  write.csv(sample.df, sample, quote = F)
  return(out)
}

pdf(paste0(output, '.clusters.pdf'))

for(sample in f.list){
  cluster.tally[sample,] <- classify.fcm(sample, som.model, som.cluster, paramx, paramy, label, flow.col, k)
}

dev.off()

write.csv(cluster.tally, paste0(output, '.cluster.tally.csv'), quote = F)

## Define a function to make a detailed plot of all events,
## color-coded by cluster, for a given sample.

event.plot <- function(sample, paramx, paramy, k, label, flow.col){
  
  sample.df <- read.csv(sample, row.names = 1, header = T)
  
  plot(sample.df[,paramy] ~ sample.df[,paramx],
       type = 'n',
       main = sample,
       xlab = paramx,
       ylab = paramy,
       log = 'xy')
  
  for(cluster in 1:k){
    r = which(sample.df[paste0('cluster.', label)] == cluster)
    out[cluster] <- length(r)
    
    points(sample.df[r, paramy] ~ sample.df[r, paramx],
           pch = 19,
           cex = 0.3,
           col = flow.col[cluster])
  }
  
  legend('topleft',
         legend = paste('Cluster', 1:k),
         pch = 19,
         col = flow.col)
}

for(sample in f.list){
  pdf(paste0(sample, '.pdf'))
  event.plot(sample, 'FSC', 'FL1', k, label, flow.col)
  dev.off()
}
  