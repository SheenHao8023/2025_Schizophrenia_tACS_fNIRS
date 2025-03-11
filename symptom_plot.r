pkg_list <- c("ggplot2", "dplyr", "tidyr", "readxl")  
sapply(pkg_list, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
})

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

file_path <- "C:/Users/XinHao/Desktop/tES_SZ_Symptom/symptom_data.xlsx"
data <- read_excel(file_path)
long_data <- data %>% 
  #  filter(Group == 1) %>% # Experimental
  #  filter(Group == 2) %>% # ControlActive
  #  filter(Group == 3) %>% # ControlSham
  #  filter(Group == 4) %>% # ControlResting
  #filter(ID == 'HZ001' & Type == 'SZ') %>% # Individual
  filter(ID %in% c('HZ001', 'HZ002') & Type == 'SZ') %>%
  pivot_longer(cols = 9:ncol(data), names_to = "Measure", values_to = "Score")
long_data$Session <- factor(long_data$Session, levels = c("Pre", "Post"))
theme_custom <- theme(
  panel.grid = element_blank(),
  axis.line = element_line(color = 'black', size = 0.75),
  axis.ticks = element_line(color = 'black', size = 0.75),
  axis.text = element_text(size = 12, color = "black"),
  axis.title = element_text(size = 18),
  legend.position = "none"
)
color_map <- c("Pre" = "#f47c7c", "Post" = "#80d6ff")
measure_names <- unique(long_data$Measure) 

for (measure in measure_names) {
  p <- ggplot(long_data %>% filter(Measure == measure), aes(x = Session, y = Score, fill = Session)) +
    geom_bar(stat = "summary", fun = "mean", width = 0.6, colour = "black", size = 0.5) +
    geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, size = 0.8) +
    #geom_point(stat = "summary", fun = "mean", size = 2, shape = 16) +
    #geom_line(stat = "summary", fun = "mean", color = "gray70", size = 0.75, linetype = "solid") +
    scale_fill_manual(values = color_map) +
    labs(x = "Session", y = "Score", title = measure) +
    theme_custom
  output_path <- paste0(dirname(file_path), "/", measure, "_bar_plot.png")
  ggsave(output_path, plot = p, width = 4, height = 6, dpi = 300)
}