library(ggplot2)
library(gganimate)
library(tidyverse)

data <- read.csv("fdtd_output.csv")

# Filters data to the center of the y-axis
data_filtered <- data %>% filter(j == 20)

p <- ggplot(data_filtered, aes(x = i, y = k, fill = Hz)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = 'Time: {frame_time}', x = 'X', y = 'Z') +
  theme_minimal()

anim <- p + transition_time(T) +
  labs(title = 'Hz field at Time: {frame_time}')

anim_save("Hz_animation.gif", animation = anim)