# working script for running BIGCAAT locally 
# v 0.9
require(tidyverse)
require(doParallel)
require(parallel)
require(foreach)
require(tibble)
require(HLAtools)
require(SSHAARP)
require(dplyr)
require(gtools)
require(BIGDAWG)
require(stringr)
require(data.table)

scripts<-list.files(paste(getwd(), 'scripts', sep ='/'), full.names=TRUE)

for(i in 1:length(scripts)){
  source(scripts[i])
}

# usage
dataset <- ""
loci <- c()

BIDS_results <- BIDS(loci, dataset)


