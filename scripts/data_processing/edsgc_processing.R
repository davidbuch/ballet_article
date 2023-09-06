# The raw data were provided by Woncheol Jang
# He used the dataset in his 2004 article:
# "Nonparametric density estimation and clustering in astronomical sky surveys"
library(tidyverse)
dir.create('data/intermediate_data', recursive = TRUE, showWarnings = FALSE)

galaxy_data <- read.table("data/raw_data/edsgc.dat", skip = 1) %>% 
  dplyr::rename(x = V3, y = V4) %>% 
  select(x,y)
abell_clusters <- read.table("data/raw_data/abell.dat", skip = 9) %>% 
  dplyr::rename(x = V2, y = V3) %>%
  select(x,y)
edcci_clusters <- read.table("data/raw_data/edcci.dat") %>%
  dplyr::rename(x = V2, y = V3) %>%
  select(x,y)

write.csv(galaxy_data, 
          "data/intermediate_data/edsgc.csv",
          row.names = FALSE)
write.csv(abell_clusters, 
          "data/intermediate_data/abell_clusters.csv",
          row.names = FALSE)
write.csv(edcci_clusters, 
          "data/intermediate_data/edcci_clusters.csv",
          row.names = FALSE)
