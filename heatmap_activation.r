install_if_missing <- function(packages) {
  to_install <- packages[!packages %in% installed.packages()[, "Package"]]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
  invisible(lapply(packages, library, character.only = TRUE))
}

install_if_missing(c("ggplot2", "readxl", "viridis", "reshape2", "patchwork"))

file_path <- "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/5.glm/activation.xlsx"
all_groups <- list(
  Experimental = c(1,2,4,5), 
  RestingControl = c(3), 
  ShamControl = c(), 
  ActiveControl = c())

data <- suppressMessages(read_excel(file_path, col_names = FALSE))
data <- as.matrix(data)
data <- apply(data[1:5, 3:ncol(data)], 2, as.numeric)

Channels <- c("rFG", "rSTG", "rTPJ", "rPMC", "rSFG", "rDLPFC", "rFPC", "DMPFC", 
              "lFPC", "lDLPFC", "lSFG", "lPMC", "lTPJ", "lSTG", "lFG")
Conditions <- c("R1", "HA1", "HB1", "HEO1", "HA2", "HB2", "HEO2", "HA3", "HB3", 
                "HEO3", "HA4", "HB4", "HEO4", "R2", "HA5", "HB5", "HEO5")
custom_colors <- colorRampPalette(c('#80d6ff', 'white', '#f47c7c'))(6)

for (group_name in names(all_groups)) {
  group_indices <- all_groups[[group_name]]
  data_matrix_avg_group <- colMeans(data[group_indices, 3:ncol(data), drop = FALSE], na.rm = TRUE)
  data_matrix_avg_group <- matrix(data_matrix_avg_group, nrow = 17, ncol = 15, byrow = FALSE)
  
  # Baseline-online-offline
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
  
  # RCT
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
  
  # 组合两个图
  combined_plot <- heatmap_plot1 | heatmap_plot2  # 使用 patchwork 进行组合
  
  # 保存图像
  ggsave(paste0("C:/Users/XinHao/Desktop/tES_SZ_fNIRS/5.glm/activation_heatmap_", group_name, ".png"),
         plot = combined_plot, width = 12, height = 8, dpi = 300)
}
