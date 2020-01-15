setwd("~/bowman_lab/bowman_lab_github/flow_cytometry_scripts")

#### parameters ####

f.list <- list.files(path = '.', pattern = '.qc.csv', ignore.case = T)
output <- 'test.sg'
input <- 'test_SG.som.Rdata'
label <- 'sg'
paramx <- 'SSC.HLin' # SSC best indicator for size, pg. 422 "In Living Color"?
paramy <- 'GRN.B.HLin'
beads.added <- 0.5 * 10^4 # 50 ul of 10^5 per ml
SSC.beads.llimit <- 3.1
FSC.beads.llimit <- 4
FL5.beads.llimit <- 3.3

#### classify ####

library(kohonen)
library(oce)
library(plotrix)

load(input)

cluster.tally <- matrix(nrow = length(f.list), ncol = k + 2)
colnames(cluster.tally) <- c(1:k, 'beads', 'correction_factor_applied')
row.names(cluster.tally) <- f.list
flow.col <- oce.colorsFreesurface(k)

## Set variables for testing only.

paramx <- 'SSC.HLin'
paramy <- 'GRN.B.HLin'
label <- 'sg'
sample <- f.list[1]
sample <- 'NC1.qc.csv'
som.model <- som.model
cluster.tally.df <- cluster.tally
cluster.vector <- som.cluster

## Classify events from all samples, making a single compiled pdf
## showing size and layout of clusters for each sample.

classify.fcm <- function(sample,
                         som.model,
                         cluster.vector,
                         paramx,
                         paramy,
                         label,
                         flow.col,
                         k,
                         SSC.beads.llimit,
                         FSC.beads.llimit,
                         FL5.beads.llimit){
  
  print(sample)
  
  params <- colnames(som.model$data[[1]])
  
  sample.df <- read.csv(sample, row.names = 1, header = T)
  sample.mat <- as.matrix(sample.df[,params])
  sample.mat <- log10(sample.mat)
  
  sample.predict <- predict(som.model, sample.mat)
  sample.df[paste0('cluster.', label)] <- cluster.vector[sample.predict$unit.classif]
  
  sample.beads <- which(log10(sample.df$SSC.HLin) > SSC.beads.llimit &
                           log10(sample.df$BLU.V.HLin) > FL5.beads.llimit &
                           log10(sample.df$FSC.HLin) > FSC.beads.llimit)
  
  ## beads are assigned cluster 0
  
  sample.df[sample.beads, paste0('cluster.', label)] <- 0
  
  out <- vector(length = k + 2)
  
  ## regular plot
  
  plot(sample.mat[,paramy] ~ sample.mat[,paramx],
       type = 'n',
       main = sample,
       xlab = paramx,
       ylab = paramy)
  
  for(cluster in c(0:k)){
    
    r = which(sample.df[paste0('cluster.', label)] == cluster)
    out[cluster] <- length(r)

    temp.sd.x <- sd(sample.mat[r, paramx])
    temp.sd.y <- sd(sample.mat[r, paramy])
    
    temp.mean.x <- mean(sample.mat[r, paramx])
    temp.mean.y <- mean(sample.mat[r, paramy])
    
    if(cluster == 0){
      draw.ellipse(temp.mean.x, temp.mean.y, a = temp.sd.x, b = temp.sd.y)
    }else(draw.ellipse(temp.mean.x, temp.mean.y, a = temp.sd.x, b = temp.sd.y, border = flow.col[cluster]))

    text(temp.mean.x, temp.mean.y, length(r))
  }
  
  legend('topleft',
         legend = c(paste('Cluster', 1:k), 'Beads'),
         pch = 1,
         col = c(flow.col, 'black'))
  
  ## Convert to events ml^-1.  This script was originally setup to assume the beads formed a coherent and exclusive cluster, however,
  ## there are too often other events in this cluster, resulting in an overcount.  Now selecting based on user-defined limits, which
  ## should be the same as those selected for fcm_model_SG 
  
  beads.counted <- length(sample.beads)
  count.cf <- beads.added / beads.counted
  out[k + 1] <- beads.counted
  out <- out * count.cf
  out[k + 2] <- count.cf
  
  write.csv(sample.df, sample, quote = F)
  return(out)
}

pdf(paste0(output, '.clusters.pdf'))

for(sample in f.list){
  cluster.tally[sample,] <- classify.fcm(sample,
                                         som.model,
                                         som.cluster,
                                         paramx,
                                         paramy,
                                         label,
                                         flow.col,
                                         k,
                                         SSC.beads.llimit,
                                         FSC.beads.llimit,
                                         FL5.beads.llimit)
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
  