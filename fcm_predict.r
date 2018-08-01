setwd('~/bowman_lab/cyflow_space/ecuador_2018')

#### parameters ####

f.list <- list.files(path = '.', pattern = '_SG.qc.csv', ignore.case = T)
output <- 'ecuador_2018.sg'
input <- 'ecuador_2018_SG.som.Rdata'
label <- 'sg'
paramx <- 'SSC' #SSC best indicator for size, pg. 422 "In Living Color"
paramy <- 'FL1'

### classify ###

library(kohonen)
library(oce)

load(input)
k <- max(som.cluster)

cluster.tally <- matrix(nrow = length(f.list), ncol = k)
colnames(cluster.tally) <- 1:k
row.names(cluster.tally) <- f.list
flow.col <- oce.colorsFreesurface(k)

## variables for testing 

# paramx <- 'SSC'
# paramy <- 'FL1'
# label <- 'sg'
# sample <- f.list[1]
# som.model <- som.model
# cluster.tally.df <- cluster.tally
# cluster.vector <- som.cluster

## end

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
    points(sample.mat[r, paramy] ~ sample.mat[r, paramx],
           pch = 19,
           cex = 0.3,
           col = flow.col[cluster])
  }
  
  write.csv(sample.df, sample, quote = F)
  return(out)
}

pdf(paste0(output, '.clusters.pdf'))

for(sample in f.list){
  cluster.tally[sample,] <- classify.fcm(sample, som.model, som.cluster, paramx, paramy, label, flow.col, k)
}

dev.off()

## Keep commented unless you want to reassign clusters
write.csv(cluster.tally, paste0(output, '.cluster.tally.csv'), quote = F)