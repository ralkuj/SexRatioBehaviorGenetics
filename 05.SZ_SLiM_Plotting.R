#### Packages and setup ####

library(dplyr)
library(tidyr) 
library(ggplot2)
library(data.table)

options(scipen = 999) #to avoid inconsistencies with scientific notation

#Set working directory containing combined results file
#setwd("")

#### Upload and organise data ####

# #Code used to originally combine all files into master output
# #Grab the set wd
# stated_wd <- getwd()
# #Set working directory that has all of the simulation results files
# wd <- paste0(stated_wd, "/SongZhang-original_Simulation_SLiM_human_directional/")
# setwd(wd)
# # Get a list of files matching the pattern, ? means optional, * means 0+ digits
# folder_list <- Filter(function(x) grepl("Human_mr1e-\\d+_ms0\\.\\d+$", x),
#                       list.dirs(path = ".", full.names = TRUE, recursive = TRUE))
# # Initialize an empty data frame to store the combined data
# simulation_data <- data.frame()
# # Set up list of files, AF_O.txt to AF_29.txt
# file_list <- sprintf("Summary_all_%d.txt", 0:29)
# # Loop through each folder and get all files
# for(folder_name in folder_list){
#   # Extract the mutation rate and size from the file name using regular expression
#   mr <- as.numeric(sub("./Human_mr1e-(\\d+)_ms0\\.\\d+", "\\1", folder_name))
#   ms <- as.numeric(sub("./Human_mr1e-\\d+_ms0\\.(\\d+)", "\\1", folder_name))
#   # Change working directory to that folder
#   wd <- paste0(stated_wd, "/SongZhang-original_Simulation_SLiM_human_directional", folder_name, "/")
#   setwd(wd)
#   # Loop through each file and read the data
#   for (file_name in file_list) {
#     temp_data <- fread(file_name) # Read the CSV file
#     # Only retain final generation data
#     temp_data <- temp_data[which.max(temp_data$N_gen), , drop = FALSE]
#     # Add mutation rate and mutation size info
#     temp_data$mr <- mr
#     temp_data$ms <- ms
#     simulation_data <- rbind(simulation_data, temp_data)  # Combine the data
#   }
# }
# # Save combined data file
# write.csv(simulation_data, file = "SongZhang-original_Simulation_SLiM_human_directional_combined.csv", row.names = FALSE)

# Read combined data file
simulation_data <- read.csv("SongZhang-original_Simulation_SLiM_human_directional_combined.csv")

#Fix the ms and mr values
simulation_data <- simulation_data %>%
  mutate(ms = case_when(
    ms == 1 ~ 0.01,
    ms == 2 ~ 0.02,
    ms == 4 ~ 0.04,
    ms == 8 ~ 0.08,
    ms == 16 ~ 0.16,
    ms == 5 ~ 0.005,
    ms == 25 ~ 0.0025,
    ms == 125 ~ 0.00125,
    TRUE ~ ms  # Preserve original value if not matched
  ))
simulation_data <- simulation_data %>%
  mutate(mr = case_when(
    mr == 2 ~ 1e-02,
    mr == 3 ~ 1e-03,
    mr == 4 ~ 1e-04,
    mr == 5 ~ 1e-05,
    mr == 6 ~ 1e-06,
    TRUE ~ mr  # Preserve original value if not matched
  ))

# Average for each parameter combo
  simulation_data <- simulation_data %>%
  group_by(mr, ms) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

#### Heat map ####

# Convert mr and ms to factors to treat them as categorical variables
simulation_data$mr <- as.factor(simulation_data$mr)
simulation_data$ms <- as.factor(simulation_data$ms)
  
# Set up steps showing simulations where h2 < .00038 (upper 95CI)
# + 1 to each parameter value since the line will need to be shifted down in ggplot
steps <- data.frame(mr = factor(c(1e-06, 1e-05, 1e-04, 1e-03, 1e-02)), ms = factor(c(.02, .01, .005, .0025, .00125)))

  # Create the plot
  ggplot(simulation_data, aes(x = as.factor(mr), y = as.factor(ms), fill = Sex_ratio_European)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("red", "white", "darkgreen"),
      #Sets the midway point according to Song and Zhang's "optimal value" of 0.524
      values = c(0, 0.4615, 1), 
      na.value = NA,
      limits = range(simulation_data$Sex_ratio_European, na.rm = TRUE)
    ) +
    # Add black step
    geom_step(data = steps, mapping = aes(x=mr, y=ms), group = 1, 
              inherit.aes = F, position = position_nudge(x = -0.5, y = -0.5)) +
    labs(
      x = "Mutation Rate",
      y = "Mutation Size",
      fill = "Resulting\nSex Ratio"
    ) +
    # Add *s
    annotate("text", x = as.factor(0.00001), y = as.factor(0.005), 
             label = "*", color = "black", 
             size = 6) +
    annotate("text", x = as.factor(0.0001), y = as.factor(0.0025), 
             label = "*", color = "black", 
             size = 6) +
    annotate("text", x = as.factor(0.001), y = as.factor(0.00125), 
             label = "*", color = "black", 
             size = 6) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White panel background
      panel.grid = element_blank(),  # Remove grid lines
      panel.border = element_blank(),  # Remove any border
      axis.ticks = element_blank(),  # Remove axis ticks
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels if needed
    )
  
#### Checking the h2 values against the 95%CI ####  
  
# Summarise European h2 values for each mr and ms combination
h2_summary <- simulation_data %>%
  select(mr, ms, Herit_European) %>%
  arrange(mr, ms)

# Create a wide format table with mr as rows and ms as columns
h2_wide <- simulation_data %>%
  select(mr, ms, Herit_European) %>%
  pivot_wider(names_from = ms, values_from = Herit_European, names_prefix = "ms_")

# Add flag column to identify h2 values < 0.00038 (upper 95CI threshold)
h2_summary_flagged <- simulation_data %>%
  select(mr, ms, Herit_European) %>%
  mutate(below_threshold = ifelse(Herit_European < 0.00038, "*", "")) %>%
  arrange(mr, ms)

# Create wide format with asterisks marking values < 0.00038
h2_wide_flagged <- simulation_data %>%
  select(mr, ms, Herit_European) %>%
  mutate(h2_display = ifelse(Herit_European < 0.00038, 
                              paste0(format(Herit_European, scientific = FALSE), "*"),
                              format(Herit_European, scientific = FALSE))) %>%
  select(mr, ms, h2_display) %>%
  pivot_wider(names_from = ms, values_from = h2_display, names_prefix = "ms_")

print(h2_wide_flagged)

  
