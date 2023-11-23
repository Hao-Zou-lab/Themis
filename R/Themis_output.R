#' @title Bubble chart visualization
#'
#' @param directory Directory containing files that summarise the p-values for each method under different quantitative criteria in the same data
#' @export
generate_Bubble_chart <- function(directory) {
  required_packages <- c("ggplot2", "reshape2", "scales", "ggtree", "dplyr", "tidyverse")
  for (package_name in required_packages) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      install.packages(package_name)
      library(package_name, character.only = TRUE)
    } else {
      library(package_name, character.only = TRUE)
    }
  }
  # Get the names of all files in the directory, suffixed by _P_sum.txt
  file_names <- list.files(directory, pattern = "_P_sum\\.txt$", full.names = TRUE)

  # Iterate through each file
  for (file_name in file_names) {
    # Read data frame A
    A <- read.table(file_name, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

    # Replace values greater than 0.05 with NA
    A[A > 0.05] <- NA

    # Calculation B
    B <- -log10(A)

    # Convert to dataframe C and add row name columns
    C <- as.data.frame(B)
    C <- C %>% rownames_to_column(var = "name")

    # Melt operation
    dt2 <- melt(C)

    # Sort by name
    dt2$name <- factor(dt2$name, levels = rev(unique(dt2$name)), ordered = TRUE)

    # Draw and save the graph
    p <- ggplot(dt2, aes(variable, name)) +
      scale_color_gradientn(values = seq(0, 1, 0.2), colors = c('#FFD2D2', '#F50504', '#910504')) +
      theme_bw() +
      geom_point(aes(size = 10, color = `value`)) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
      xlab(NULL) + ylab(NULL) +
      theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
      geom_point(shape = 1, fill = dt2$value, stroke = 0.5, color = "black", size = 4.5)

    output_file_name <- paste0(dirname(file_name), "/", tools::file_path_sans_ext(basename(file_name)), "_bubble.tiff")
    ggsave(output_file_name, p, width = 4, height = 3, dpi = 600)
  }
}

#' @title Radargram quantification
#'
#' @param directory Directory containing files that summarise the p-values for each method under different quantitative criteria in the same data
#' @export
generate_radar_chart <- function(directory) {
  required_packages <- c("ggradar", "dplyr", "scales", "tibble", "fmsb", "ggplot2")
  for (package_name in required_packages) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      install.packages(package_name)
      library(package_name, character.only = TRUE)
    } else {
      library(package_name, character.only = TRUE)
    }
  }
  rank_columns <- function(data) {
    for (col in colnames(data)) {
      # Check for missing values and replace them with NA in the original data
      data[[col]][is.na(data[[col]])] <- NA

      # Rank the column while ignoring the missing values
      rank_col <- rank(-data[[col]], ties.method = "max", na.last = "keep")

      # Update the column with the ranks
      data[[col]] <- rank_col
    }
    return(data)
  }

  # Read input data
  file_names <- list.files(pattern = "_P_sum.txt", full.names = TRUE)

  for (input_file in file_names) {
    data <- read.table(input_file, sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
    df <- rank_columns(data)
    df[is.na(df)] <- 0

    # Create radar chart and set legend size and position
    df <- data.frame(Method = rownames(df), df)

    p1 <- ggradar(
      df %>%
        tail(nrow(df)) %>%
        select(1:ncol(df)),
      values.radar = c("0", "4", "8"),
      grid.min = 0, grid.mid = 4, grid.max = 9,
      axis.label.size = 4,
      group.line.width = 0.8,
      group.point.size = 3,
      group.colours = NULL,
      background.circle.colour = "#D7D6D1",
      gridline.mid.colour = "grey",
      plot.legend = if (nrow(df) > 1) TRUE else FALSE,
      legend.title = "Method",
      legend.position = "right",
    ) +
      ggtitle(gsub("_P_sum.txt", "", basename(input_file))) +
      theme(plot.title = element_text(size = 14, hjust = 0.5),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.position = c(1.1, 0.5))

    # Save radar chart as a TIFF file
    output_file <- paste0(gsub("_P_sum.txt", "", input_file), ".tiff")
    ggsave(output_file, p1, width = 15, height = 8, dpi = 600)

    df <- df[-1]

    # Calculate radar chart areas
    radar_area <- function(data) {
      library(sf)
      library(tidyverse)

      areas <- c()
      row_names <- rownames(data)
      for (i in 1:nrow(data)) {
        d1 <- data[i, ] %>% as.numeric()
        center <- c(x = 2.1, y = 2.1)
        half <- seq(0, pi, length.out = 180)
        dfx <- data.frame(x = NA, y = NA)

        for (j in 1:length(d1)) {
          x = d1[j] * sin(j * (2 * pi) / length(d1)) + center["x"]
          y = d1[j] * cos(j * (2 * pi) / length(d1)) + center["y"]
          dfa = data.frame(x = x, y = y)
          dfx = dfa %>% bind_rows(dfx)
        }

        dfx[length(d1) + 1, ] = dfx[1, ]
        pol = st_polygon(list(as.matrix(dfx)))
        areas[i] <- st_area(pol)
      }

      result_df <- data.frame(RowName = row_names, RadarArea = areas)
      return(result_df)
    }

    # Calculate radar chart areas and save to a CSV file
    result_df <- radar_area(df)
    output_file <- paste0(gsub("_P_sum.txt", "", input_file), ".csv")
    write.csv(result_df, file = output_file, row.names = FALSE)

    # Create bar chart using ggplot2
    library(ggplot2)
    summary_data <- read.csv(output_file)

    # Custom color scheme
    custom_colors <- c("#FF5A5F", "#FFB400", "#007A87", "#8CE085", "#7B0051", "#00D1C1", "#FFA98F", "#B4A76C", "#494646")

    # Create bar chart
    p <- ggplot(summary_data, aes(x = RowName, y = RadarArea, fill = RowName)) +
      geom_col(width = 0.5) +
      scale_fill_manual(values = custom_colors) +
      labs(title = paste0(gsub("_P_sum.txt", "", basename(input_file)), " Radar Chart Areas"),
           x = "Method",
           y = "Radar Area") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            panel.background = element_rect(fill = "white"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"))

    # Save bar chart as a TIFF file
    output_file <- paste0(gsub("_P_sum.txt", "", input_file), "_bar.tiff")
    ggsave(output_file, p, width = 6, height = 4, dpi = 600)
  }
}

#' Themis output
#'
#' @param expFile a m*n matrix with m genes and n samples
#' @param cliFile Results of various analyses on the sample data, including Stemness Score (mRNAsi), immune infiltration scores (ESTIMATEScore, ImmuneScore, and StromalScore), immune therapy response prediction (TIDE), chemotherapy drug sensitivity analysis, as well as clinical outcome information (OS/PFS) and merged clinical risk factors.
#' @param cluster_file_path Result file of tumor molecular typing analysis for each sample in the data.
#' @param directory Working path for storing files
#' @param custom.data A logical value indicating whether custom data files are used.
#' @param custom.Method A logical value indicating whether custom methods are used.
#' @param covariates List of covariates required for univariate analysis
#' @export
Themis_output <- function(directory) {
  directory = directory
  setwd(directory)
  message(paste0("\n",">>> Running ", "Bubble chart"))
  generate_Bubble_chart(directory)
  message(paste0("\n",">>> Running ", "radar chart"))
  generate_radar_chart(directory)
  message("Analysis Ending!")
  }

