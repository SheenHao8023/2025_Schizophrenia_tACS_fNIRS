install_if_missing <- function(packages) {
  to_install <- packages[!packages %in% installed.packages()[, "Package"]]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
  invisible(lapply(packages, library, character.only = TRUE))
}

install_if_missing(c("ggplot2", "readxl", "viridis"))

file_path <- "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/5.glm/activation.xlsx"
data <- suppressMessages(read_excel(file_path, col_names = FALSE))
data <- as.matrix(data)

data_matrix_avg_group <- colMeans(data)
data_matrix_avg_group <- matrix(data_matrix_avg_group, nrow = 17, ncol = 15, byrow = FALSE)

Channels <- c("rFG", "rSTG", "rTPJ", "rPMC", "rSFG", "rDLPFC", "rFPC", "DMPFC", 
              "lFPC", "lDLPFC", "lSFG", "lPMC", "lTPJ", "lSTG", "lFG")
Conditions <- c("R1", "HA1", "HB1", "HEO1", "HA2", "HB2", "HEO2", "HA3", "HB3", 
                "HEO3", "HA4", "HB4", "HEO4", "R2", "HA5", "HB5", "HEO5")
custom_colors <- colorRampPalette(c('#80d6ff', 'white', '#f47c7c'))(6)

# Baseline-online-offline
heatmap_data <- melt(data_matrix_avg_group[2:13, ])
heatmap_plot <- ggplot(heatmap_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = custom_colors, limits = c(-0.1, 0.1)) +
  theme_minimal() +
  labs(title = "Activation", x = "Condition", y = "ROI") +
  theme(panel.grid = element_blank(),  
      axis.text.x = element_text(angle = 0, hjust = 1, color = "black", size = 14), 
      axis.text.y = element_text(angle = 0, hjust = 1, color = "black", size = 14),
      plot.background = element_rect(fill = "white"), 
      panel.background = element_rect(fill = "white")) +
      scale_x_continuous(breaks = 1:12, labels = Conditions[2:13]) +
      scale_y_continuous(breaks = 1:15, labels = Channels) +
      coord_fixed(ratio = 1)
ggsave("C:/Users/XinHao/Desktop/tES_SZ_fNIRS/5.glm/activation_heatmap_onetime.png",
       plot = heatmap_plot, width = 12, height = 12, dpi = 300)

# RCT
heatmap_data <- melt(data_matrix_avg_group[c(1:4, 14:17), ])
heatmap_plot <- ggplot(heatmap_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = custom_colors, limits = c(-0.1, 0.1)) +
  theme_minimal() +
  labs(title = "Activation", x = "Condition", y = "ROI") +
  theme(panel.grid = element_blank(),  
      axis.text.x = element_text(angle = 0, hjust = 1, color = "black", size = 14), 
      axis.text.y = element_text(angle = 0, hjust = 1, color = "black", size = 14),
      plot.background = element_rect(fill = "white"), 
      panel.background = element_rect(fill = "white")) +
      scale_x_continuous(breaks = 1:8, labels = Conditions[c(1:4, 14:17)]) +
      scale_y_continuous(breaks = 1:15, labels = Channels) +
      coord_fixed(ratio = 1)
ggsave("C:/Users/XinHao/Desktop/tES_SZ_fNIRS/5.glm/activation_heatmap_tentime.png",
       plot = heatmap_plot, width = 12, height = 12, dpi = 300)