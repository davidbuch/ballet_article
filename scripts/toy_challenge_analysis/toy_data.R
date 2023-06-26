library(tidyverse)
dir.create('output/toy_challenge', recursive = TRUE, showWarnings = FALSE)

# Load and Plot the Three Toy Challenge Datasets
two_moons <- readRDS("data/clean_data/two_moons.rds")
circles <- readRDS("data/clean_data/noisy_circles.rds")
tsne <- readRDS("data/clean_data/tsne.rds")

two_moons <- two_moons %>% 
  mutate(dataset = 'Two Moons')
circles <- circles %>%
  mutate(dataset = 'Noisy Circles')
tsne <- tsne %>%
  mutate(dataset = 't-SNE')

plot_data <- rbind(two_moons,
                   tsne,
                   circles)
plot_data <- plot_data %>%
  mutate(dataset = factor(dataset,
                          levels = c('Two Moons', 'Noisy Circles', 't-SNE')))

png("output/toy_challenge/toy_challenge_data.png", 
    width = 12, height = 4, units = 'in', res = 300)
print(
  ggplot(plot_data) + 
    geom_point(aes(x = x, y = y), size = 0.25) + 
    facet_wrap(~ dataset, nrow = 1) + 
    labs(title = "Toy Challenge Datasets")
)
dev.off()
