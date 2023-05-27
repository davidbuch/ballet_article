
library(salso)
min_vi_part <- salso(clustering_samps)

ggplot() + 
  geom_point(aes(x = Y[,1], y = Y[,2], color = factor(min_vi_part)))

for(i in sample(1:S, 5)){
  p <- ggplot() + 
    geom_point(aes(x = Y[,1], y = Y[,2], color = factor(clustering_samps[i,]))) + 
    guides(color = "none")
  print(p)
}
