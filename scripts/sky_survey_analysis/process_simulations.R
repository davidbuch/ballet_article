library(tidyverse)
library(xtable)

dir_name <- "output/sky_survey_analysis/sim_study"
filename_list <- list.files(path="output/sky_survey_analysis/sim_study")
filename_list <- filename_list[str_detect(filename_list, 'accuracy_\\d+.rds')]

files <- list()
for(i in 1:length(filename_list)){
  filename <- paste0(dir_name, '/',filename_list[[i]])
  cat(paste("Loading", filename, "\n"))
  files[[i]] <- readRDS(filename)
}

test_results <- do.call(rbind, files)

result_summary <- test_results %>% 
  group_by(type) %>%
  summarise(sensitivity = mean(sensitivity),
            specificity = mean(specificity),
            exact_match_frac = mean(exact_match_frac)) %>%
  xtable()

print(result_summary, 
      file = "output/sky_survey_analysis/sim_study_test_results.txt")

print(result_summary)
