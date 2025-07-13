# 自动安装并加载包
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

install_if_missing(c("R.matlab", "ggplot2", "reshape2", "dplyr", "purrr", "patchwork", "tools"))

# 🧠 分组列表，按组写上被试编号前缀即可
group_subjects_list <- list(
  Experimental     = c("SX001", "SX007", "SX008"),
  ShamControl      = c(),
  RestingControl   = c("SX002"),
  ActiveControl    = c("SX006"),
  BehaviorControl  = c("SX003", "SX004", "SX005")
)

group_names <- names(group_subjects_list)

# 🔧 配置参数
folder_paths <- c(
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/5.coh_pair",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/5.coh_sz",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/5.gc_pairab",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/5.gc_pairba",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/5.gc_sz",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/6.mi_pair",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/6.mi_sz",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/6.te_pairab",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/6.te_pairba",
  "C:/Users/XinHao/Desktop/tES_SZ_fNIRS/6.te_sz"
)

connectivity_name <- c('SZ--HC', 'SZ--SZ', 'SZ-->HC', 'HC-->SZ', 'SZ<-->SZ')
Y_labels <- c('ROI', 'ROI', 'SZ', 'HC', 'ROI')

Conditions <- c("R1", "HA1", "HB1", "HEO1", "HA2", "HB2", "HEO2", "HA3", "HB3", 
                "HEO3", "HA4", "HB4", "HEO4", "R2", "HA5", "HB5", "HEO5")
roi_names <- c("rFG", "rSTG", "rTPJ", "rPMC", "rSFG", "rDLPFC", "rFPC", "DMPFC", 
               "lFPC", "lDLPFC", "lSFG", "lPMC", "lTPJ", "lSTG", "lFG")

custom_colors <- colorRampPalette(c('#f7fcb9', '#addd8e', '#31a354'))(5)

plot_order <- c(
  "R1",     NA,     NA,     NA,     "R2",
  "HA1",   "HA2",  "HA3",  "HA4",  "HA5",
  "HB1",   "HB2",  "HB3",  "HB4",  "HB5",
  "HEO1",  "HEO2", "HEO3", "HEO4", "HEO5"
)

# 主循环：每个文件夹
for (folder_path in folder_paths) {
  mat_files <- list.files(folder_path, pattern = "*.mat", full.names = TRUE)
  folder_index <- ((which(folder_paths == folder_path) - 1) %% length(connectivity_name)) + 1
  
  # 每组处理
  for (group in group_names) {
    subject_ids <- group_subjects_list[[group]]
    group_condition_matrices <- list()
    
    # 遍历条件
    for (condition in Conditions) {
      condition_matrices <- list()
      
      for (subject in subject_ids) {
        matching_files <- mat_files[grepl(paste0(subject, ".*_", condition, "\\.mat$"), basename(mat_files))]
        if (length(matching_files) > 0) {
          for (file in matching_files) {
            data <- readMat(file)
            matrix_name <- names(data)[1]
            mat <- as.matrix(data[[matrix_name]])
            condition_matrices[[length(condition_matrices) + 1]] <- mat
          }
        }
      }
      
      if (length(condition_matrices) > 0) {
        avg_mat <- Reduce("+", condition_matrices) / length(condition_matrices)
        group_condition_matrices[[condition]] <- avg_mat
      }
    }
    
    # 🧠 统一色阶范围（SZ组所有 condition 的 value除对角线）
    is_sz_folder <- grepl("(_|\\.)sz$", basename(folder_path), ignore.case = TRUE)
    
    all_values <- unlist(
      lapply(group_condition_matrices, function(mat) {
        if (is_sz_folder) {
          mat[row(mat) + col(mat) == nrow(mat) + 1] <- NA  # 忽略副对角线
        }
        as.vector(mat)
      })
    )
    
    global_min <- min(all_values, na.rm = TRUE)
    global_max <- max(all_values, na.rm = TRUE)
    
    
    # 生成拼图
    plot_list_filled <- lapply(plot_order, function(condition) {
      if (is.na(condition)) {
        ggplot() + theme_void()
      } else if (!is.null(group_condition_matrices[[condition]])) {
        avg_matrix <- group_condition_matrices[[condition]]
        melted_data <- melt(avg_matrix)
        
        ggplot(melted_data, aes(Var1, Var2, fill = value)) +
          geom_tile(color = "lightgray") +
          scale_fill_gradientn(
            colors = custom_colors,
            limits = c(global_min, global_max)
          ) +
          labs(title = condition,
               x = connectivity_name[[folder_index]],
               y = Y_labels[[folder_index]]) +
          theme_minimal() +
          theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(angle = 0, hjust = 1, color = "black"),
            plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white")
          ) +
          scale_x_continuous(breaks = 1:15, labels = roi_names) +
          scale_y_continuous(breaks = 1:15, labels = roi_names)
      } else {
        ggplot() + theme_void() + labs(title = paste(condition, "(Missing)"))
      }
    })
    
    final_plot <- wrap_plots(plot_list_filled, ncol = 5)
    group_file <- paste0("combined_heatmap_", group, ".png")
    
    ggsave(
      filename = file.path(folder_path, group_file),
      plot = final_plot,
      width = 20, height = 16
    )
  }
}
