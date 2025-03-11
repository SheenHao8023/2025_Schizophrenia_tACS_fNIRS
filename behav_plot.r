## matlab跑完得到的excel文件用于此段代码
pkg_list <- c("ggplot2", "readxl", "tidyr", "dplyr", "ggprism")  # 你的包列表
sapply(pkg_list, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
})
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggprism)

file_path <- "C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx"
sheets <- c("IC", "WS", "RD", "TV")  
conditions <- c("Hearing each other", "Hearing A", "Hearing B")
output_dir <- "C:/Users/XinHao/Desktop/tES_SZ_Behav/"
color_map <- list("IC" = rep('#f47c7c', 5),"WS" = rep('#80d6ff', 5),"RD" = rep('#ffa07a', 5),"TV" = rep('#9370db', 5))
x_labels <- c("B1" = "Pre-test", "B2" = "Online1", "B3" = "Online2", "B4" = "Offline", "B5" = "Post-test")
y_limits <- list("IC" = c(0, 1), "WS" = c(0, 0.2), "RD" = c(0, 200), "TV" = c(0, 1000))

theme_custom <- theme_prism(axis_text_angle = 0) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = 'black', size = 0.75),
    axis.ticks = element_line(color = 'black', size = 0.75),
    axis.text = element_text(size = 12, color = "black", face = 'plain'),
    axis.title = element_text(size = 14, face = 'plain'),
    legend.position = "none"  
  )

for (sheet in sheets) {
  for (condition in conditions) {
    data <- read_excel(file_path, col_names = TRUE, sheet = sheet)
    # data <- data %>% filter(Group == 1) # Experimental
    # data <- data %>% filter(Group == 2) # ControlActive
    # data <- data %>% filter(Group == 3) # ControlSham
    # data <- data %>% filter(Group == 4) # ControlResting
    # data <- data %>% filter(ID == 1181) # Individual
    col_suffix <- ifelse(condition == "Hearing each other", "C1", 
                         ifelse(condition == "Hearing A", "C2", "C3"))
    selected_cols <- c("ID", paste0("B", 1:5, col_suffix))
    data_condition <- data %>% select(all_of(selected_cols))
    colnames(data_condition)[1] <- "ID"
    data_long <- data_condition %>% 
      pivot_longer(cols = -ID, names_to = "Condition", values_to = "value") %>%
      mutate(Condition = factor(gsub("B(\\d)C[1-3]", "B\\1", Condition), levels = names(x_labels), labels = x_labels))
    summary_data <- data_long %>%
      group_by(Condition) %>%
      summarise(mean_value = mean(value, na.rm = TRUE),sd_value = sd(value, na.rm = TRUE))

    p <- ggplot() +
      geom_bar(data = summary_data, aes(x = Condition, y = mean_value, fill = Condition), 
               stat = "identity", width = 0.7, colour = "black", size = 0.5, show.legend = FALSE) +
      geom_errorbar(data = summary_data, 
                    aes(x = Condition, ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                    width = 0.4, size = 0.8) +
      geom_point(data = data_long, aes(x = Condition, y = value), position = position_nudge(x = 0), size = 2, color = "black") +
      geom_line(data = data_long, aes(x = Condition, y = value, group = ID), color = "gray70", size = 0.5, na.rm = TRUE) +
      theme_custom +
      coord_cartesian(ylim = y_limits[[sheet]]) +
      labs(x = condition, y = sheet) +
      scale_fill_manual(values = color_map[[sheet]])
    ggsave(filename = paste0(output_dir, sheet, "_", gsub(" ", "_", condition), "_plot.png"), plot = p, width = 5, height = 5, dpi = 300)
  }
}