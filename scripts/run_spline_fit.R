library(readr)
library(dplyr)
library(tibble)
library(tidyr)

source("scripts/spline_fit.R")

# Small value to add to gene before log transform
PSEUDOCOUNT <- 0.01

input <- snakemake@input
output <- snakemake@output
num_batches <- snakemake@params[['num_batches']]
batch <- as.integer(snakemake@wildcards[['batch']])

#input <- list(
#    tpm = "results/Liver/tpm_all_samples.txt",
#    sample_info = "results/Liver/all_samples_info.txt"
#)
#output <- list(
#    summary = "temp.summary.txt",
#    curves_fit = "temp.curves.fit.txt",
#    curves_pstd = "temp.curves.pstd.txt",
#    re = "temp.re.txt",
#    re_structure = "temp.re_structure.txt"
#)
#num_batches <- 3000
#batch <- 1

## Load the raw data
# Load with read.table not read_tsv: for some reason read_tsv inserts a extra, unwanted column
data_table <- as_tibble(read.table(input[['tpm']], sep="\t", header=TRUE))

sample_table <- read_tsv(input[['sample_info']], show_col_types=FALSE) %>% 
    rename_with(function(x) { "sample" }, 1)

# We don't use all studies
# Some have a very different seuqencing methodoloy or bad alignment stats
drop_studies <- c("Greenwell19_AdLib", "Greenwell19_NightFeed", "Janich15", "Manella21_Liver")
# and some samples are outliers to be dropped
drop_samples <- read.csv2(input[['outliers']], header=FALSE)$V1
selected_samples <- (sample_table %>% filter(!(study %in% drop_studies)) %>% filter(!(sample %in% drop_samples)))

# Find genes passing an expression cutoff
# require at least 33% of samples to have non-zero measurements
data_table$non_zero <- rowSums(data_table %>% select(selected_samples$sample) != 0)
data_table$passing_gene <- data_table$non_zero > (nrow(selected_samples)%/%3)
selected_data_table <- data_table %>% filter(passing_gene)


# Find the genes in our subset
num_genes_total <- sum(data_table$passing_gene)
batch_size <- (num_genes_total %/% num_batches) + 1
fst <- 1 + batch_size * batch
lst <- batch_size * (batch+1)
selected_genes <-  (selected_data_table %>% slice(fst:lst))$Name
#selected_genes <- c("ENSMUSG00000055866")# XXX TODO USE BATCHES


times <- selected_samples$time
study <- selected_samples$study
summary <- NULL
curves <- NULL
re_structure <- NULL
re <- NULL
for(gene in selected_genes) {
    # Gather data for one gene
    gene_values <- data_table %>%
        filter(Name == gene) %>%
        select(selected_samples$sample)

    # Fit the spline to that gene
    tryCatch({
        res <- fit_splines(log(t(gene_values)+PSEUDOCOUNT), study, times)

        if (is.null(summary)) {
            # If first iteration, must generate the tables we store these in
            summary <- as_tibble(res$summary) %>%
                        mutate(gene=gene, .before=)
            curves <- as_tibble(res$curves) %>%
                        mutate(gene=gene)
            re <- as_tibble(res$random_effects) %>%
                        mutate(study=rownames(res$random_effects),
                               gene=gene)
            re_structure <- as_tibble(res$random_effects_structure) %>%
                                mutate(var=rownames(res$random_effects_structure),
                                       gene=gene)
        } else {
            # Store the results
            summary <- bind_rows(
                summary,
                as_tibble(res$summary) %>%
                    mutate(gene=gene))
            curves <- bind_rows(
                curves,
                as_tibble(res$curves) %>%
                    mutate(gene=gene))
            re <- bind_rows(
                re,
                as_tibble(res$random_effects) %>%
                    mutate(study=rownames(res$random_effects),
                           gene=gene))
            re_structure <- bind_rows(
                re_structure,
                as_tibble(res$random_effects_structure) %>%
                    mutate(var=rownames(res$random_effects_structure),
                           gene=gene))
        }
    },
    error=function(e) {
        print("Failed on gene:")
        print(gene)
    })
    print("Done with gene:")
    print(gene)
}

write_tsv(summary %>% relocate(gene), output[['summary']])
write_tsv(curves %>% select(gene, u, fit) %>% spread(u,fit), output[['curves_fit']])
write_tsv(curves %>% select(gene, u, pstd) %>% spread(u,pstd), output[['curves_pstd']])
write_tsv(re %>% relocate(gene, study), output[['re']])
write_tsv(re_structure %>% relocate(gene, var), output[['re_structure']])
