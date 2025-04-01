install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}
install_if_missing(c("R.matlab", "ggplot2", "reshape2", "dplyr", "purrr"))

folder_paths <- c(
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/6.coh_pair",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/6.coh_sz",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/7.gc_pairab",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/7.gc_pairba",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/7.gc_sz"
)
connectivity_name <- c('SZ--HC', 'SZ--SZ', 'SZ-->HC', 'HC-->SZ', 'SZ<-->SZ')
Y_labels <- c('ROI', 'ROI', 'SZ', 'HC', 'ROI')
limits_list <- list(c(0.15,0.5), c(0.3,0.9), c(0,0.035), c(0,0.025), c(0,0.125))
Conditions <- c("R1", "HA1", "HB1", "HEO1", "HA2", "HB2", "HEO2", "HA3", "HB3", 
                "HEO3", "HA4", "HB4", "HEO4", "R2", "HA5", "HB5", "HEO5")
roi_names <- c("rFG", "rSTG", "rTPJ", "rPMC", "rSFG", "rDLPFC", "rFPC", "DMPFC", 
               "lFPC", "lDLPFC", "lSFG", "lPMC", "lTPJ", "lSTG", "lFG")

custom_colors <- colorRampPalette(c('#f7fcb9', '#addd8e', '#31a354'))(5)

for (folder_path in folder_paths) {
  mat_files <- list.files(folder_path, pattern = "*.mat", full.names = TRUE)
  group_matrices <- list()
  
  for (condition in Conditions) {
    condition_files <- mat_files[grepl(paste0("_", condition, "\\.mat$"), mat_files)]
    if (length(condition_files) > 0) {
      matrices <- lapply(condition_files, function(file) {
        data <- readMat(file)
        matrix_name <- names(data)[1]  
        as.matrix(data[[matrix_name]])
      })
      
      group_matrices[[condition]] <- Reduce("+", matrices) / length(matrices)
    }
  }
  
  for (condition in Conditions) {
    if (!is.null(group_matrices[[condition]])) {
      avg_matrix <- group_matrices[[condition]]
      melted_data <- melt(avg_matrix)
      p <- ggplot(melted_data, aes(Var1, Var2, fill = value)) +
        geom_tile(color = "lightgray") +
        scale_fill_gradientn(colors = custom_colors, limits = limits_list[[which(folder_paths == folder_path)]]) +
        labs(title = paste(Conditions[which(Conditions == condition)]), 
             x = paste(connectivity_name[which(folder_paths == folder_path)]), 
             y = paste(Y_labels[which(folder_paths == folder_path)])) +
        theme_minimal() +
        theme(panel.grid = element_blank(),  
              axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
              axis.text.y = element_text(angle = 0, hjust = 1, color = "black"),
              plot.background = element_rect(fill = "white"), 
              panel.background = element_rect(fill = "white")) +
        scale_x_continuous(breaks = 1:15, labels = roi_names) +
        scale_y_continuous(breaks = 1:15, labels = roi_names)
      
      ggsave(filename = file.path(folder_path, paste0("heatmap_", condition, ".png")), plot = p, width = 6, height = 6)
    }
  }
}
