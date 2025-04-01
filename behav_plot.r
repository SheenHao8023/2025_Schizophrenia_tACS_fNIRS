# matlab跑通后导出的excel直接用于此脚本
# ----------------- 1. 安装 & 加载所需 R 包 -----------------
pkg_list <- c("ggplot2", "readxl", "tidyr", "dplyr", "ggprism", "patchwork")  
sapply(pkg_list, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
})

# ----------------- 2. 读取 Excel 数据 -----------------
file_path <- "C:/Users/XinHao/Desktop/tES_SZ_Behav/behavior_data.xlsx"
output_dir <- "C:/Users/XinHao/Desktop/tES_SZ_Behav/"

sheets <- c("IC", "WS", "RD", "TV")  # 4 种度量
conditions <- c("Hearing each other", "Hearing A", "Hearing B")  # 3 个实验条件

# 颜色与组名映射表
color_map <- list("IC" = rep('#f47c7c', 5), 
                  "WS" = rep('#80d6ff', 5), 
                  "RD" = rep('#ffa07a', 5), 
                  "TV" = rep('#9370db', 5))
group_names <- c("1" = "Experimental", 
                 "2" = "Control active", 
                 "3" = "Control sham", 
                 "4" = "Control resting")

# 横坐标标签映射
x_labels <- c("B1" = "Pre-test", "B2" = "Online1", "B3" = "Online2", "B4" = "Offline", "B5" = "Post-test")

# y 轴范围设置
y_limits <- list("IC" = c(0, 1), "WS" = c(0, 0.2), "RD" = c(0, 200), "TV" = c(0, 1000))

# 自定义主题
theme_custom <- theme_prism(axis_text_angle = 0) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = 'black', size = 0.75),
    axis.ticks = element_line(color = 'black', size = 0.75),
    axis.text = element_text(size = 12, color = "black", face = 'plain'),
    axis.title = element_text(size = 14, face = 'plain'),
    legend.position = "none"  
  )

# ----------------- 3. 读取数据并获取所有 Group -----------------
all_data <- read_excel(file_path, col_names = TRUE, sheet = sheets[1])  # 读取其中一个 sheet 确定 Group
unique_groups <- unique(all_data$Group)  # 获取所有独立 Group 值

# ----------------- 4. 生成每个 Group 的大图 -----------------
for (group in unique_groups) {
  plot_list <- list()  # 存储当前 Group 的所有子图
  
  for (sheet in sheets) {
    for (condition in conditions) {
      data <- read_excel(file_path, col_names = TRUE, sheet = sheet)  # 读取数据
      data_group <- data %>% filter(Group == group)  # 筛选当前组的数据
      
      # 选取不同实验条件对应的列
      col_suffix <- ifelse(condition == "Hearing each other", "C1", 
                           ifelse(condition == "Hearing A", "C2", "C3"))
      selected_cols <- c("ID", paste0("B", 1:5, col_suffix))
      data_condition <- data_group %>% select(all_of(selected_cols))
      colnames(data_condition)[1] <- "ID"  # 重新命名 ID 列
      
      # 数据变长格式
      data_long <- data_condition %>% 
        pivot_longer(cols = -ID, names_to = "Condition", values_to = "value") %>%
        mutate(Condition = factor(gsub("B(\\d)C[1-3]", "B\\1", Condition), 
                                  levels = names(x_labels), labels = x_labels))
      
      # 计算均值 & 标准差
      summary_data <- data_long %>%
        group_by(Condition) %>%
        summarise(mean_value = mean(value, na.rm = TRUE), sd_value = sd(value, na.rm = TRUE))
      
      # 绘制 ggplot 图
      p <- ggplot() +
        geom_bar(data = summary_data, aes(x = Condition, y = mean_value, fill = Condition), 
                 stat = "identity", width = 0.7, colour = "black", size = 0.5, show.legend = FALSE) +
        geom_errorbar(data = summary_data, 
                      aes(x = Condition, ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                      width = 0.4, size = 0.8) +
        geom_point(data = data_long, aes(x = Condition, y = value), 
                   position = position_nudge(x = 0), size = 2, color = "black") +
        geom_line(data = data_long, aes(x = Condition, y = value, group = ID), 
                  color = "gray70", size = 0.5, na.rm = TRUE) +
        theme_custom +
        coord_cartesian(ylim = y_limits[[sheet]]) +
        labs(x = condition, y = sheet, title = paste(sheet, "-", condition)) +
        scale_fill_manual(values = color_map[[sheet]])
      
      # 将子图存入列表
      plot_list[[length(plot_list) + 1]] <- p  
    }
  }
  
  # **将 12 张图排列成 4×3**
  combined_plot <- wrap_plots(plot_list, ncol = 3) + 
    plot_annotation(title = paste("Group:", group_names[as.character(group)]))  # 添加大图标题
  
  # 保存大图
  ggsave(filename = paste0(output_dir, "Group", group_names[as.character(group)], ".png"), 
         plot = combined_plot, width = 12, height = 16, dpi = 300)
}
