# Since the prior script runs 20 instances in parallel, outputs are written
# to a single .rds file for each run. This script collates the outputs of
# every one of those runs into one big tsv

library(dplyr)
library(purrr)
library(stringr)


OUTPUT_DIR <- "bootstrap_markers_ttmt/"
all_files <- list.files(OUTPUT_DIR, pattern = "\\.rds$", full.names = TRUE)

# Read all files, combine, and add run/cluster information
final_markers_table <- purrr::map_dfr(all_files, function(file_path) {

    # Extract the sequential run number (i) which is the last number before .rds
    # This is the most robust way to handle the file naming from the previous steps.
    run_number <- stringr::str_extract(basename(file_path), "(?<=_s_)\\d+")

    # Read the list of marker data frames
    run_list <- readRDS(file_path)

    # Combine the list elements (clusters) into one data frame, adding cluster and run info
    purrr::imap_dfr(run_list, ~ .x %>%
                        # Use tibble::rownames_to_column for a safer operation
                        tibble::rownames_to_column(var = "gene") %>%
                        dplyr::mutate(
                            cluster = .y,
                            # Convert run_number to integer for the final table
                            run_id = as.integer(run_number)
                        )
                    )
})

# Save the final combined table as a language-agnostic TSV file
write.table(
    final_markers_table,
    file = "thc_1000_bootstrap_ttmt_aggregated.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
