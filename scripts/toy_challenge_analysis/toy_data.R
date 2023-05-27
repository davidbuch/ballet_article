library(tidyverse)
dir.create('output/toy_challenge', recursive = TRUE, showWarnings = FALSE)

# Load and Plot the Three Toy Challenge Datasets
two_moons <- readRDS("data/clean_data/two_moons.rds")
circles <- readRDS("data/clean_data/noisy_circles.rds")
tsne <- readRDS("data/clean_data/tsne.rds")

two_moons <- two_moons %>% 
  mutate(dataset = 'two_moons')
circles <- circles %>%
  mutate(dataset = 'circles')
tsne <- tsne %>%
  mutate(dataset = 'tsne')

plot_data <- rbind(two_moons,
                   tsne,
                   circles)
plot_data <- plot_data %>%
  mutate(dataset = factor(dataset,
                          levels = c('two_moons', 'circles', 'tsne')))

png("output/toy_challenge/toy_challenge_data.png", 
    width = 12, height = 6, units = 'in', res = 300)
ggplot(plot_data) + 
  geom_point(aes(x = x, y = y), size = 0.25) + 
  facet_wrap(~ dataset, nrow = 1) + 
  labs(title = "Toy Challenge Datasets")
dev.off()
