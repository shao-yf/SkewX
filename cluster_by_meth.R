#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(NanoMethViz))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

# Modified cluster_reads function with binary thresholding
cluster_reads_with_threshold <- function(x, chr, start, end, min_pts = 5, threshold = 0.5) {
    assertthat::assert_that(
        is(x, "ModBamResult"),
        assertthat::is.string(chr) || (is.factor(chr) && assertthat::is.scalar(chr)),
        assertthat::is.number(start) && assertthat::is.number(end),
        assertthat::is.number(min_pts) && min_pts >= 1,
        assertthat::is.number(threshold) && threshold >= 0 && threshold <= 1
    )

    # query data
    methy_data <- query_methy(x, chr, start, end)

    if (nrow(methy_data) == 0) {
        stop(glue::glue(
            "No reads containing methylation data found in region {chr}:{start}-{end}.\n",
            "This could be due to:\n",
            "  - No data in this genomic region\n",
            "  - Incorrect chromosome naming\n",
            "  - Region outside of data coverage\n",
            "Please check your genomic coordinates and data."
        ))
    }

    methy_data <- methy_data %>%
        dplyr::filter(.data$pos >= start & .data$pos < end)

    # Apply binary thresholding to mod_prob BEFORE any other processing
    methy_data <- methy_data %>%
        dplyr::mutate(mod_prob = ifelse(.data$mod_prob >= threshold, 1, 0))

    read_stats <- get_read_stats(methy_data)

    # identify the read names whose span is at least 90% the length of maximum span
    # filter methylation data for only those reads that meet the above condition of span
    max_span <- max(read_stats$span)
    keep_reads <- read_stats$read_name[read_stats$span > 0.9 * max_span]
    methy_data <- methy_data %>%
        dplyr::filter(.data$read_name %in% keep_reads)

    # convert methylation data into a matrix with one row for each read name
    mod_mat <- methy_data %>%
        dplyr::select("read_name", "pos", "mod_prob") %>%
        dplyr::arrange(.data$pos) %>%
        tidyr::pivot_wider(names_from = "pos", values_from = "mod_prob") %>%
        NanoMethViz:::df_to_matrix()

    # pre-check before filtering
    if (nrow(mod_mat) < min_pts) {
        stop(glue::glue(
            "Insufficient reads for clustering: found {nrow(mod_mat)} reads but need at least {min_pts}.\n",
            "Try reducing 'min_pts' parameter or expanding the genomic region."
        ))
    }

    # remove positions with high missingness (>60%) then reads with high missingness (>30%)
    mod_mat_filled <- mod_mat[order(rownames(mod_mat)), ]
    col_missingness <- NanoMethViz:::mat_col_map(mod_mat_filled, NanoMethViz:::missingness)
    mod_mat_filled <- mod_mat_filled[, col_missingness < 0.6]
    row_missingness <- NanoMethViz:::mat_row_map(mod_mat_filled, NanoMethViz:::missingness)
    mod_mat_filled <- mod_mat_filled[row_missingness < 0.3, ]

    # For binary data, fill missing values with the mode (most common value) of that read
    # or if tied, use the overall mode across the region
    for (i in seq_len(nrow(mod_mat_filled))) {
        missing_indices <- is.na(mod_mat_filled[i, ])
        if (any(missing_indices)) {
            row_values <- mod_mat_filled[i, !missing_indices]
            if (length(row_values) > 0) {
                # Calculate mode for this read
                read_mode <- as.numeric(names(sort(table(row_values), decreasing = TRUE))[1])
                mod_mat_filled[i, missing_indices] <- read_mode
            } else {
                # If all values are missing for this read, use 0 as default
                mod_mat_filled[i, missing_indices] <- 0
            }
        }
    }

    # post-check before filtering
    if (nrow(mod_mat_filled) < min_pts) {
        stop(glue::glue(
            "Insufficient reads after filtering: {nrow(mod_mat_filled)} reads remaining but need at least {min_pts}.\n",
            "Try reducing 'min_pts' parameter or adjusting filtering criteria."
        ))
    }

    # cluster reads using HDBSCAN algorithm with specified minimum number of points
    dbsc <- dbscan::hdbscan(mod_mat_filled, minPts = min_pts)
    clust_df <- data.frame(read_name = rownames(mod_mat_filled), cluster_id = dbsc$cluster)

    # merge and process results of cluster analysis and read statistics
    clust_df %>%
        dplyr::inner_join(read_stats, by = "read_name") %>%
        dplyr::arrange(.data$cluster_id) %>%
        dplyr::mutate(
            cluster_id = as.factor(.data$cluster_id),
            start = as.integer(.data$start),
            end = as.integer(.data$end),
            span = as.integer(.data$span)
        )
}

# Use the original get_read_stats function from NanoMethViz
get_read_stats <- function(methy_data) {
    methy_data %>%
        dplyr::group_by(.data$read_name) %>%
        dplyr::summarise(
            start = min(.data$pos),
            end = max(.data$pos),
            mean = mean(.data$mod_prob, na.rm = TRUE),
            span = .data$end - .data$start,
            strand = unique(.data$strand),
            .groups = 'drop'
        ) %>%
        dplyr::arrange(.data$strand)
}

apply_cluster_reads_parallel <- function(mbr, bed, min_pts, threshold = 0.5, num_cores = 4) {
  # Initialize a parallel backend
  #cl <- makeCluster(detectCores())
  #try with 4 cores first
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  # Simple progress message instead of deprecated progress_estimated
  message("Starting parallel clustering...")

  # Define a function to process each row
  process_row <- function(row) {
    # Get the current row
    current_row <- bed[row, ]

    # Apply the cluster_reads function with binary thresholding to the current row
    row_cluster <- tryCatch({
      cluster_reads_with_threshold(mbr, current_row$chr, current_row$start, current_row$end, min_pts = min_pts, threshold = threshold)}, 
      error = function(err){
        # Handle the error (e.g., print a message)
        message(paste("Error in row", row, ":", err$message))
        # Skip to the next row
        return(NA)
      })

    if (all(is.na(row_cluster))) {
      return(NULL)
    }

    # Add the CGI_id and chr, start, end
    row_cluster <- row_cluster %>% mutate(CGI_id = paste0(current_row$chr, ":", current_row$start, "-", current_row$end), chr = current_row$chr, start = current_row$start, end = current_row$end)

    # Calculate the average methylation by cluster_id and add it as a new column
    row_cluster <- row_cluster %>% group_by(cluster_id) %>% mutate(avg_cluster_methylation = mean(mean))

    # If there are exactly 2 clusters, assign the cluster with the lowest average methylation to Xa and the other to Xi
    if (nlevels(row_cluster$cluster_id) == 2) { # Watch out for NA clusters...
      low_mC <- row_cluster %>% filter(cluster_id %in% c("1", "2")) %>% pull(avg_cluster_methylation) %>% min()
      high_mC <- row_cluster %>% filter(cluster_id %in% c("1", "2")) %>% pull(avg_cluster_methylation) %>% max()
      row_cluster <- row_cluster %>% mutate(assigned_X = case_when(avg_cluster_methylation == low_mC ~ "Xa", avg_cluster_methylation == high_mC ~ "Xi", TRUE ~ "NA"))
    } else {
      row_cluster$assigned_X <- NA
    }

    return(row_cluster)
  }

  # Apply the function to each row in parallel
  res <- foreach(row = 1:nrow(bed), .combine = bind_rows, .packages = c("tidyverse", "NanoMethViz")) %dopar% {
    # Optionally print progress every 100 rows
    if (row %% 100 == 0) message(paste("Processed", row, "rows"))
    process_row(row)
  }

  # Stop the parallel backend
  stopCluster(cl)

  # Combine the results
  res <- bind_rows(res)

  return(res)
}

calculate_skew_by_block <- function(clustered_reads, haplotyped_reads){
  #remove uninformative reads that don’t clusters or reads from CGIs that don’t have exactly 2 clusters
  clustered_reads <- clustered_reads %>% filter(assigned_X %in% c("Xa","Xi"))
  #remove reads that appear multiple times because they span several CGIs
  clustered_reads <- clustered_reads %>% distinct(read_name, .keep_all = TRUE)
  #merge methylation cluster information with haplotype and phase set information
  df2 <- left_join(clustered_reads,haplotyped_reads)
  #remove the reads that couldn’t be haplotyped
  df2 <- df2 %>% filter(!is.na(HP))

  #count by haplotype blocks
  counts_by_block <- df2 %>% group_by(PS, assigned_X, HP) %>% summarise(counts = n())

  skew_by_block <- counts_by_block %>%
    unite(combi, assigned_X, HP) %>%
    mutate(combi = recode(combi, "Xa_1" = "H1_Xa", "Xa_2" = "H2_Xa", "Xi_1" = "H1_Xi", "Xi_2" = "H2_Xi")) %>%
    pivot_wider(id_cols = PS, names_from = combi, values_from = counts, values_fill = 0) %>%
    mutate(H1_Xa_skew = (H1_Xa + H2_Xi) / (H1_Xa + H1_Xi + H2_Xa + H2_Xi))

  return(skew_by_block)
}
#This updated version uses the foreach function to iterate over the rows in parallel, and the pb$tick() line updates the progress bar for each iteration. The results are combined using the .combine = bind_rows argument, and the required packages are specified using the .packages argument.


# Retrieve the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Access the argument(s)
lib <- args[1]
bam <- args[2]
BED <- read_tsv(args[3], col_names = c("chr", "start", "end"))
haplotyped_reads <- read_tsv(args[4], col_names = c("read_name", "HP", "PS"))
ncpus <- strtoi(args[5])

#create the ModBamResult object
mbr <- ModBamResult(
    methy = ModBamFiles(
        samples = lib,
        paths = bam
    ),
    samples = data.frame(
        sample = lib,
        group = 1
    )
)

# Apply clustering with binary thresholding
clustered_reads <- apply_cluster_reads_parallel(mbr, BED, min_pts = 5, num_cores = ncpus)
write_tsv(clustered_reads, paste0(lib, "_CGIX_clustered_reads.tsv.gz"))

skew <- calculate_skew_by_block(clustered_reads, haplotyped_reads)

write_tsv(skew, paste0(lib,"_CGIX_skew.tsv.gz"))
