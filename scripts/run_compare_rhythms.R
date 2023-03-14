# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.15")   # This installs bioconductor
#BiocManager::install(c("SummarizedExperiment", "DESeq2", "edgeR", "limma", "rain")) # Packages needed by compareRhythms
#install.packages("devtools")    # if it is not already installed
#devtools::install_github("cran/npsm")   # Package archived by CRAN
#devtools::install_github("cran/DODR")   # Package archived by CRAN
#devtools::install_github("bharathananth/compareRhythms", build_vignettes = TRUE, dependencies = TRUE)

library(tibble)
library(dplyr)
library(readr)

input_file1 <- snakemake@input[[1]]
input_file2 <- snakemake@input[[2]]

outfile <- snakemake@output[[1]]

timepoints1 <- snakemake@params[['timepoints1']]
timepoints2 <- snakemake@params[['timepoints2']]

print("Preparing compareRhythms")
print(input_file1)
print(input_file2)
print(snakemake@input[['outliers1']])
print(snakemake@input[['outliers2']])
print(outfile)
print(timepoints1)
print(timepoints2)

library(compareRhythms)

# Load the two datasets and combine them
print("Loading data")
data1 <-  read_tsv(input_file1)
data2 <-  read_tsv(input_file2)

# Select just the non-outlier samples
print("Loading outliers1")
outliers1 <- read_tsv(snakemake@input[["outliers1"]], col_names = c("sample_id"))
print("Loading outliers2")
outliers2 <- read_tsv(snakemake@input[["outliers2"]], col_names = c("sample_id"))
data1_outliers <- data1 %>% select(starts_with("GSM")) %>% colnames %in% outliers1$sample_id
data2_outliers <- data2 %>% select(starts_with("GSM")) %>% colnames %in% outliers2$sample_id
is_outlier <- c(data1_outliers, data2_outliers)

# Join the two datasets together
print("Joining data")
joined_data <- full_join(
    data1,
    data2,
    by = "Name",
    suffix = c("_A", "_B"),
) %>% column_to_rownames("Name")
print(head(as.tibble(joined_data)))

# Make the design matrix of time and groupings
design <- dplyr::bind_rows(
   tibble(
        time = timepoints1,
        group = factor("A", c("A", "B"))
    ),
    tibble(
        time = timepoints2,
        group = factor("B", c("A", "B"))
    )
)


# Find genes with a minimum expression amount
print("Selecting genes")
MIN_EXPRESSION = 20 # read counts
total_counts <- rowSums(joined_data %>% select(starts_with("GSM")))
selected_genes <- total_counts > MIN_EXPRESSION

# Get just the parts of the data that we care about
selected_data <- as.matrix(joined_data[selected_genes, !is_outlier])
selected_design <- design # We already have dropped outliers
print(head(selected_data))
print(selected_design)

# Run with default period range
print("Running compareRhythms")
res <- compareRhythms(
    selected_data,
    selected_design,
    method = "voom",
    period = 24,
)
write_tsv(res, outfile)
