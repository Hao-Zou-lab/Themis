#' @title Kaplan Meier analysis
#'
#' @param directory Default working directory, where all "_MetaMerge.txt" files are located
#'
#' @return
#' @export
#'
#' @examples
analyzeBatchSurvival <- function(directory) {
  # Check and install packages
  required_packages <- c("survival", "survminer")
  for (package_name in required_packages) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      install.packages(package_name)
      library(package_name, character.only = TRUE)
    } else {
      library(package_name, character.only = TRUE)
    }
  }

  # Get all "_MetaMerge.txt" files in the directory
  clusterFiles <- list.files(directory, pattern = "_MetaMerge.txt$", full.names = TRUE)

  # Create an empty result data frame
  result_df <- data.frame(Method = character(0), OS_pValue = numeric(0), PFS_pValue = numeric(0))

  # Iterate through each file and perform survival analysis
  for (file in clusterFiles) {
    # Extract filename prefix
    file_prefix <- sub("_MetaMerge.txt$", "", basename(file))
    file_prefix <- gsub(".txt$", "", file_prefix)

    # Read data file
    rt <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

    # Processing OSfutime and PFSfutime
    rt$OSfutime <- rt$OSfutime / 365
    rt$PFSfutime <- rt$PFSfutime / 365

    # OS Analysis
    length <- length(levels(factor(rt$Cluster)))
    diff <- survdiff(Surv(OSfutime, OSfustat) ~ Cluster, data = rt)
    OS_pValue <- 1 - pchisq(diff$chisq, df = length - 1)

    # Filter out samples with missing values for PFS time or PFS state
    rt <- rt[!is.na(rt$PFSfutime) & !is.na(rt$PFSfustat), ]

    # Check if there are still enough samples to analyse
    if (nrow(rt) > 0) {
      # Perform PFS analysis
      diff <- survdiff(Surv(PFSfutime, PFSfustat) ~ Cluster, data = rt)
      PFS_pValue <- 1 - pchisq(diff$chisq, df = length - 1)
    } else {
      # Not enough samples for PFS analysis, set PFS_pValue to NA or other appropriate value
      PFS_pValue <- NA
    }

    # Store the results in the result data frame
    result <- data.frame(Method = file_prefix, OS_pValue = OS_pValue, PFS_pValue = PFS_pValue)
    result_df <- rbind(result_df, result)
    result_file <- file.path(directory, "batch_KManalysis.txt")
    write.table(result_df, file = result_file, sep = "\t", quote = FALSE, row.names = F)
  }

  # Returns the result data frame of batch processing
  return(result_df)
}

#' @title Univariate cox analysis based on OS
#'
#' @param directory Default working directory, where all "_MetaMerge.txt" files are located
#'
#' @return Output the results of the univariate survival analysis for each file and save them in a file ending with "_OSunivariate.txt".Extract the p-values of the Cluster variables of all univariate analysis result files and save them in a file named "batch_OSunivariate.txt".
#' @export
#'
#' @examples
perform_OS_cox <- function(directory) {
  required_packages <- c("survival", "survminer", "tidyverse")
  for (package_name in required_packages) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      install.packages(package_name)
      library(package_name, character.only = TRUE)
    } else {
      library(package_name, character.only = TRUE)
    }
  }
  covariates <- c("Age", "gender", "stage", "tumor_event_type", "treatment_outcome", "Residual_tumo", "Cluster")
  # Perform Cox analysis
  perform_OScox_analysis <- function(covariates) {
    file_list <- list.files(directory, pattern = "_MetaMerge.txt", full.names = TRUE)

    for (file in file_list) {
      # Read data from the file
      data <- read.table(file, header = TRUE, sep = "\t")
      data <- data[, colSums(is.na(data)) != nrow(data)]
      # Create formulas for Cox analysis
      covariates <- c("Age", "gender", "stage", "tumor_event_type", "treatment_outcome", "Residual_tumo", "Cluster")
      covariates <- intersect(covariates, colnames(data))
      univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OSfutime, OSfustat)~', x)))

      # Perform Cox analysis for each covariate
      univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data)})

      # Extract and format results
      univ_results <- data.frame()
      for(i in 1:length(univ_models)){
        if(length(univ_models[[i]]$coefficients)== 1){
          x <- summary(univ_models[[i]])
          p.value <- signif(x$wald["pvalue"], 2)
          HR <- signif(x$conf.int[, "exp(coef)"], 2)
          HR.confint.lower <- signif(x$conf.int[, "lower .95"], 2)
          HR.confint.upper <- signif(x$conf.int[, "upper .95"], 2)
          HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
          res <- data.frame(p.value, HR)
          colnames(res) <- c("p.value", "HR(95% CI for HR)")
          rownames(res) <- names(univ_models[i])
        }else{
          x <- summary(univ_models[[i]])
          p.value <- signif(x[["coefficients"]][,"Pr(>|z|)"], 2)
          HR <- signif(x$conf.int[, "exp(coef)"], 2)
          HR.confint.lower <- signif(x$conf.int[, "lower .95"], 2)
          HR.confint.upper <- signif(x$conf.int[, "upper .95"], 2)
          HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
          names(HR) <- names(HR.confint.lower)
          res <- data.frame(p.value, HR)
          colnames(res) <- c("p.value", "HR(95% CI for HR)")
        }
        univ_results <- rbind(univ_results,res)
      }

      # Combine results into a data frame
      a <- rownames_to_column(univ_results,var="factor")

      # Save results to file
      file_prefix <- sub("_MetaMerge.txt", "", basename(file))
      output_file <- paste0(file_prefix, "_OSunivariate.txt")
      write.table(a, file = output_file, sep = "\t", row.names = F, quote = FALSE)
    }
  }
  # Perform Cox analysis
  perform_OScox_analysis(covariates)
  file_suffix <- "_OSunivariate.txt"
  # Extract and save results
  result_df <- data.frame()
  for (file_name in list.files(directory, pattern = file_suffix)) {
    file_prefix <- str_remove(file_name, file_suffix)
    file_data <- read.table(paste0(directory, "/", file_name), header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    p_value <- tail(file_data[["p.value"]], 1)
    result_df <- rbind(result_df, p_value)
  }

  row_names <- str_replace(list.files(directory, pattern = file_suffix), file_suffix, "")
  rownames(result_df) <- row_names

  col_names <- "OSunivariate cox"
  colnames(result_df) <- col_names

  result_file <- file.path(directory, "batch_OSunivariate.txt")
  write.table(result_df, file = result_file, sep = "\t", quote = FALSE)
}

#' @title Univariate cox analysis based on PFS
#'
#' @param directory Default working directory, where all "_MetaMerge.txt" files are located
#'
#' @return Output the results of the univariate survival analysis for each file and save them in a file ending with "_PFSunivariate.txt".Extract the p-values of the Cluster variables of all univariate analysis result files and save them in a file named "batch_PFSunivariate.txt".
#' @export
#'
#' @examples
perform_PFS_cox <- function(directory) {
  required_packages <- c("survival", "survminer", "tidyverse")
  for (package_name in required_packages) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      install.packages(package_name)
      library(package_name, character.only = TRUE)
    } else {
      library(package_name, character.only = TRUE)
    }
  }
  covariates <- c("Age", "gender", "stage", "tumor_event_type", "treatment_outcome", "Residual_tumo", "Cluster")
  # Perform Cox analysis
  perform_PFScox_analysis <- function(covariates) {
    file_list <- list.files(directory, pattern = "_MetaMerge.txt", full.names = TRUE)

    for (file in file_list) {
      # Read data from the file
      data <- read.table(file, header = TRUE, sep = "\t")
      data <- data[, colSums(is.na(data)) != nrow(data)]
      # Create formulas for Cox analysis
      covariates <- c("Age", "gender", "stage", "tumor_event_type", "treatment_outcome", "Residual_tumo", "Cluster")
      covariates <- intersect(covariates, colnames(data))
      univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(PFSfutime, PFSfustat)~', x)))

      # Perform Cox analysis for each covariate
      univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data)})

      # Extract and format results
      univ_results <- data.frame()
      for(i in 1:length(univ_models)){
        if(length(univ_models[[i]]$coefficients)== 1){
          x <- summary(univ_models[[i]])
          p.value <- signif(x$wald["pvalue"], 2)
          HR <- signif(x$conf.int[, "exp(coef)"], 2)
          HR.confint.lower <- signif(x$conf.int[, "lower .95"], 2)
          HR.confint.upper <- signif(x$conf.int[, "upper .95"], 2)
          HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
          res <- data.frame(p.value, HR)
          colnames(res) <- c("p.value", "HR(95% CI for HR)")
          rownames(res) <- names(univ_models[i])
        }else{
          x <- summary(univ_models[[i]])
          p.value <- signif(x[["coefficients"]][,"Pr(>|z|)"], 2)
          HR <- signif(x$conf.int[, "exp(coef)"], 2)
          HR.confint.lower <- signif(x$conf.int[, "lower .95"], 2)
          HR.confint.upper <- signif(x$conf.int[, "upper .95"], 2)
          HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
          names(HR) <- names(HR.confint.lower)
          res <- data.frame(p.value, HR)
          colnames(res) <- c("p.value", "HR(95% CI for HR)")
        }
        univ_results <- rbind(univ_results,res)
      }

      # Combine results into a data frame
      a <- rownames_to_column(univ_results,var="factor")

      # Save results to file
      file_prefix <- sub("_MetaMerge.txt", "", basename(file))
      output_file <- paste0(file_prefix, "_PFSunivariate.txt")
      write.table(a, file = output_file, sep = "\t", row.names = F, quote = FALSE)
    }
  }
  # Perform Cox analysis
  perform_PFScox_analysis(covariates)
  file_suffix <- "_PFSunivariate.txt"
  # Extract and save results
  result_df <- data.frame()
  for (file_name in list.files(directory, pattern = file_suffix)) {
    file_prefix <- str_remove(file_name, file_suffix)
    file_data <- read.table(paste0(directory, "/", file_name), header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    p_value <- tail(file_data[["p.value"]], 1)
    result_df <- rbind(result_df, p_value)
  }

  row_names <- str_replace(list.files(directory, pattern = file_suffix), file_suffix, "")
  rownames(result_df) <- row_names

  col_names <- "PFSunivariate cox"
  colnames(result_df) <- col_names

  result_file <- file.path(directory, "batch_PFSunivariate.txt")
  write.table(result_df, file = result_file, sep = "\t", quote = FALSE)
}

#' @title Comparative analysis of TIDE score, stemness score and estimate score
#'
#' @param directory Default working directory, where all "_MetaMerge.txt" files are located
#'
#' @return
#' @export
#'
#' @examples
process_immune_scores <- function(directory) {
  files <- list.files(directory, pattern = "_MetaMerge.txt$", full.names = TRUE)
  p_value_results <- data.frame()
  required_packages <- c("dplyr","ggplot2", "ggpubr")
  for (package_name in required_packages) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      install.packages(package_name)
      library(package_name, character.only = TRUE)
    } else {
      library(package_name, character.only = TRUE)
    }
  }
  processData <- function(df, col_names) {
    p_values <- c()

    for (col_name in col_names) {
      plot_data <- df %>%
        select(Cluster, all_of(col_name))

      p <- ggplot(plot_data, aes(x = Cluster, y = df[[col_name]], fill = Cluster)) +
        geom_boxplot() +
        geom_jitter(position = position_jitter(0.2), color = "gray27", size = 1) +
        ggtitle(paste0(col_name)) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x = element_text(size = 12, colour = "black"),
              axis.text.y = element_text(size = 12, colour = "black"),
              axis.title.x = element_text(size = 14, face = "bold"),
              axis.title.y = element_text(size = 14, face = "bold")) +
        stat_compare_means(method = "anova", aes(group = Cluster)) +
        scale_fill_manual(values = c("#DB9B34", "#108B96", "#B72230", "#317CB7")) +
        xlab("Cluster") + ylab(col_name)

      anova_result <- anova(lm(plot_data[[col_name]] ~ Cluster, data = plot_data))
      p_value <- anova_result$"Pr(>F)"[1]

      p_values <- c(p_values, p_value)

      file_prefix <- sub("_cluster.txt$", "", basename(file))
      file_prefix <- gsub(".txt$", "", file_prefix)
      file_name <- paste0(file_prefix, "_", col_name, ".tiff")
      ggsave(file_name, p, width = 5, height = 5, dpi = 600)
    }

    return(p_values)
  }
  for (file in files) {
    df <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    file_prefix <- sub("_MetaMerge.txt$", "", basename(file))
    file_prefix <- gsub(".txt$", "", file_prefix)

    col_names <- c("TIDE", "StromalScore", "ImmuneScore", "ESTIMATEScore", "mRNAsi_score")
    p_values <- processData(df, col_names)

    p_value_results <- rbind(p_value_results, data.frame(Method = file_prefix,
                                                         TIDE = p_values[1],
                                                         StromalScore = p_values[2],
                                                         ImmuneScore = p_values[3],
                                                         ESTIMATEScore = p_values[4],
                                                         mRNAsi_score = p_values[5]))
  }

  write.table(p_value_results, file = paste0(directory, "/batch_TMEscore.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

#' @title P value's summary
#'
#' @param directory Path to the file prefixed with "batch_" that contains the p-value associated with the clinical ration.
#'
#' @return Summary file containing p-value results for different methods in different clinical quantitative criteria
#' @export
#'
#' @examples
merge_all_p_values <- function(directory) {
  #Batch read file data into data frame
  file_list <- list.files(directory, pattern = "^batch_", full.names = TRUE)
  data_list <- lapply(file_list, function(file) {
    read.table(file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  })
  merged_data <- do.call(cbind, data_list)
  # Load all_pvalue_sum.Rdata
  load(file ="./sum data/all_pvalue_sum.RData")

  # Merge p_sum.txt with all_pvalue_sum.Rdata by column names
  merged_summary <- rbind(merged_data, data)

  # Extract common suffix from row names and combine them into separate data frames
  suffixes <- sub("^.*_", "", row.names(merged_summary))
  unique_suffixes <- unique(suffixes)
  for (suffix in unique_suffixes) {
    rows <- grepl(paste0("_", suffix, "$"), row.names(merged_summary), perl = TRUE)
    sub_df <- merged_summary[rows, ]

    # Output sub_df to file "_P_sum.txt"
    out_file <- file.path(directory, paste0(suffix, "_P_sum.txt"))
    write.table(sub_df, file = out_file, sep = "\t", quote = FALSE)
  }

  # Output merge result to file All_P_summary.txt
  out_file <- file.path(directory, "All_P_summary.txt")
  write.table(merged_summary, file = out_file, sep = "\t", quote = FALSE)
}

#' Themis'process
#'
#' @param directory Working path for storing files
#' @export
Themis_process1 <- function(directory){
  directory = directory
  setwd(directory)
  message(paste0("\n",">>> Running ", "K-M analysis(OS & PFS)"))
  analyzeBatchSurvival(directory)
  message(paste0("\n",">>> Running ", "univariate Cox analysis(OS)"))
  perform_OS_cox(directory)
  message(paste0("\n",">>> Running ", "univariate Cox analysis(PFS)"))
  perform_PFS_cox(directory)
  message(paste0("\n",">>> Running ", "stemness/immune/TIDE scores"))
  process_immune_scores(directory)
  message(paste0("\n",">>> Running ", "p values merge"))
  merge_all_p_values(directory)
}

