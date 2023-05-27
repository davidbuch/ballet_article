library(tidyverse)
library(PTT)
source("R/scale_box.R")
source("R/fn_from_apt.R")


## we seem to have some memory overflow issues when fitting with 
## resolutionsgreater than 12
galaxy_data <- read.csv("data/intermediate_data/primary/edsgc.csv")
x <- scale_box(as.matrix(galaxy_data))

nsims <- 1000
max.resol = 12
start <- Sys.time()
res <- apt(X = x, Xpred = x, n.post.samples = nsims, max.resol = max.resol)
end <- Sys.time()
print(end - start)

# currently extracting samples 990-1000... takes several minutes
fsamps <- fn_from_apt(res, x, max.resol)

