#' Merge subtyping results and metadata files
#'
#' @param directory Path to store the _cluster.txt file
#' @param path Path for storing 8 test datasets and 30 TCGA tumor metadata data
#'
#' @return Outputs a merged metadata file and cluster file for each piece of data
#' @export
#'
#' @examples
Themis_custom_method <- function(directory, path) {
  # Determine if the file name contains "Data".
  if (length(list.files(directory, pattern = "_cluster.txt$")) > 0 && grepl("Data", list.files(directory, pattern = "_cluster.txt$"))) {
    # Load the dataset data_meta_sum.RData
    load(file.path(path, "data_meta_sum.RData"))
    data_sum <- data
  } else {
    # Load the dataset tcga_meta_sum.RData
    load(file.path(path, "tcga_meta_sum.RData"))
    data_sum <- data
  }

  merge_files <- function(data_sum) {
    # Get all txt files under the specified path
    files <- list.files(directory, pattern = "_cluster.txt$")

    for(file in files){
      # Read txt file
      df <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

      # Merging of data
      merge_data <- merge(df, data_sum, by = "row.names")
      rownames(merge_data) <- merge_data$row.names
      merge_data$row.names <- NULL
      # Output results file
      output_file_name <- paste0(sub("_cluster.txt", "", file), "_MetaMerge.txt")
      write.table(merge_data, output_file_name, quote = FALSE, sep = "\t", row.names = F, col.names = TRUE)
    }
  }

  # Calling internal functions
  merge_files(data_sum)
}
