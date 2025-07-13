install_if_missing <- function(packages) {
  to_install <- packages[!packages %in% installed.packages()[, "Package"]]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
  invisible(lapply(packages, library, character.only = TRUE))
}

install_if_missing(c("ggplot2", "readxl", "viridis", "reshape2", "patchwork"))

file_path <- "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/4.glm/activation.xlsx"

# Step 1: 读取原始数据
raw_data <- suppressMessages(read_excel(file_path, col_names = FALSE))

# Step 2: 提取组号（第2列）
group_numbers <- as.integer(unlist(raw_data[, 2]))

# Step 3: 构建自动分组
group_map <- split(1:nrow(raw_data), group_numbers)
group_names <- c(
  "1" = "Experimental",
  "2" = "ShamControl",
  "3" = "RestingControl",
  "4" = "ActiveControl",
  "5" = "BehaviorControl"
)

# Step 4: 提取数值矩阵（第3列开始为数据）
data_matrix <- as.matrix(raw_data[, 3:ncol(raw_data)])
data_matrix <- apply(data_matrix, 2, as.numeric)

# 检查行数一致
if (nrow(data_matrix) != nrow(raw_data)) stop("数据矩阵行数和原始数据不一致")

Channels <- c("rFG", "rSTG", "rTPJ", "rPMC", "rSFG", "rDLPFC", "rFPC", "DMPFC", 
              "lFPC", "lDLPFC", "lSFG", "lPMC", "lTPJ", "lSTG", "lFG")
Conditions <- c("R1", "HA1", "HB1", "HEO1", "HA2", "HB2", "HEO2", "HA3", "HB3", 
                "HEO3", "HA4", "HB4", "HEO4", "R2", "HA5", "HB5", "HEO5")
custom_colors <- colorRampPalette(c('#80d6ff', 'white', '#f47c7c'))(6)

for (group_name in names(all_groups)) {
  group_indices <- all_groups[[group_name]]
  
  if (length(group_indices) == 0) next  # 跳过空组
  
  # 获取该组的平均激活数据
  data_matrix_avg_group <- colMeans(data_matrix[group_indices, , drop = FALSE], na.rm = TRUE)
  
  # 重构成 17 条件 × 15 ROI 的矩阵
  data_matrix_avg_group <- matrix(data_matrix_avg_group, nrow = 17, ncol = 15, byrow = FALSE)
  
  # ✅ 画图 1：Online activation（第 2~13 条件）
  heatmap_data1 <- melt(data_matrix_avg_group[2:13, ])
  max_abs_value <- max(abs(range(heatmap_data1$value, na.rm = TRUE)))
  heatmap_plot1 <- ggplot(heatmap_data1, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = custom_colors, limits = c(-max_abs_value, max_abs_value)) +
    theme_minimal() +
    labs(title = "Activation (Online)", x = "Condition", y = "ROI") +
    theme(panel.grid = element_blank(),  
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12), 
          axis.text.y = element_text(angle = 0, hjust = 1, color = "black", size = 12),
          plot.background = element_rect(fill = "white"), 
          panel.background = element_rect(fill = "white")) +
    scale_x_continuous(breaks = 1:12, labels = Conditions[2:13]) +
    scale_y_continuous(breaks = 1:15, labels = Channels) +
    coord_fixed(ratio = 1)
  
  # ✅ 画图 2：RCT（第 1-4 和 14-17 条件）
  heatmap_data2 <- melt(data_matrix_avg_group[c(1:4, 14:17), ])
  max_abs_value <- max(abs(range(heatmap_data2$value, na.rm = TRUE)))
  heatmap_plot2 <- ggplot(heatmap_data2, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = custom_colors, limits = c(-max_abs_value, max_abs_value)) +
    theme_minimal() +
    labs(title = "Activation (RCT)", x = "Condition", y = "ROI") +
    theme(panel.grid = element_blank(),  
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12), 
          axis.text.y = element_text(angle = 0, hjust = 1, color = "black", size = 12),
          plot.background = element_rect(fill = "white"), 
          panel.background = element_rect(fill = "white")) +
    scale_x_continuous(breaks = 1:8, labels = Conditions[c(1:4, 14:17)]) +
    scale_y_continuous(breaks = 1:15, labels = Channels) +
    coord_fixed(ratio = 1)
  
  combined_plot <- heatmap_plot1 | heatmap_plot2
  
  ggsave(paste0("C:/Users/XinHao/Desktop/tES_SZ_fNIRS/4.glm/activation_heatmap_", group_name, ".png"),
         plot = combined_plot, width = 12, height = 8, dpi = 300)
}
