library(dplyr)
library(purrr)
library(stringr)


OUTPUT_DIR <- "bootstrp_markers_geno/"
all_files <- list.files(OUTPUT_DIR, pattern = "\\.rds$", full.names = TRUE)
# only extract two genes
genes_to_save <- c("Reln","Cadm2")

# Read all files, combine, and add run/cluster informatio
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

filtered_markers_table <- final_markers_table %>%
    dplyr::filter(gene %in% genes_to_save)


# Save the final combined table as a language-agnostic TSV file
write.table(
    filtered_markers_table,
    file = "thc_1000_bootstrap_geno_aggregated.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
