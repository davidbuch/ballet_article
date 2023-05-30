library(dplyr)
library(xtable)

dir_name <- "output/sky_survey_analysis/sim_study"
filename_list <- list.files(path="output/sky_survey_analysis/sim_study")
filename_list <- filename_list[str_detect(filename_list, '.rds')]

files <- list()
for(i in 1:length(filename_list)){
  files[[i]] <- readRDS(filename_list[[i]])
}

test_results <- do.call(rbind, files)

test_results %>% 
  group_by(type) %>%
  summarise(sensitivity = mean(sensitivity),
            specificity = mean(specificity)) %>%
  xtable()
