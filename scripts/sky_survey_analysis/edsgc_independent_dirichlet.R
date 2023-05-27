library(tidyverse)

galaxy_data <- read.csv("data/intermediate_data/primary/edsgc.csv")
galaxy_clusters <- read.csv("data/intermediate_data/primary")

N <- nrow(galaxy_data)
S <- 1000
bin_resolution <- 100

splitrange <- function(v, n = 10) {rv <- range(v); rv[1] + ((rv[2] - rv[1])/n)*c(-1e-6,1:n)}

galaxy_data <- galaxy_data %>%
  mutate(x = RA,
         y = -DEC,
         xbin = cut(x, splitrange(x, n = bin_resolution)),
         ybin = cut(y, splitrange(y, n = bin_resolution)),
         xybin = factor(paste(xbin,ybin, sep = ", ")))


bin_counts <- galaxy_data %>% count(xybin) %>% rename(bin_count = n)
galaxy_data <- galaxy_data %>% left_join(bin_counts, by = 'xybin')


x <- galaxy_data %>% select(x,y) %>% as.matrix
bin_count_vec <- galaxy_data %>% pull(bin_count)

density_samples <- matrix(nrow = S, ncol = N)
for(s in 1:S){
  density_samples[s,] <- rgamma(N, shape = bin_count_vec)
}
density_samples <- density_samples / rowSums(density_samples)

