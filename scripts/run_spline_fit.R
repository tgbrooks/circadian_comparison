library(readr)
library(dplyr)

source("scripts/spline_fit.R")

## Load the raw data
# Load with read.table not read_tsv: for some reason read_tsv inserts a extra, unwanted column
data_table <- as_tibble(read.table(snakemake@input[['tpm']], sep="\t", header=TRUE))

sample_table <- read_tsv(snakeamek@input[['sample_info']]) %>% 
    rename_with(function(x) { "sample" }, 1)

# We don't use all studies
# Some have a very different seuqencing methodoloy or bad alignment stats
drop_studies <- c("Greenwell19_AdLIb", "Greenwell19_NightFeed", "Janich15", "Manella21_Liver")
selected_samples <- (sample_table %>% filter(!(study %in% drop_studies)))

# Find genes passing an expression cutoff
# require at least 50 samples to have non-zero measurements
non_zero <- rowSums(data_table %>% select(selected_samples$sample) != 0)
passing_gene <- rs > 50


times <- selected_samples$time
study <- selected_samples$study
summary <- list()
curves <- list()
re_structure <- list()
re <- list()
for(gene in snakemake@params['genes']) {
    # Gather data for one gene
    gene_values <- data_table %>%
        filter(Name == gene) %>%
        select(selected_samples$sample)

    # Fit the spline to that gene
    res <- fit_splines(t(gene_values), study, times)
    summary[gene] <- res$summary
    curves[gene] <- res$curves
    re <- res$random_effects
    re_structure <- res$random_effects_structure
}
