library(dplyr)
library(readr)

input_file <- snakemake@input[[1]]
output_dir <- snakemake@params[['out_dir']]

timepoints <- snakemake@params[['timepoints']]
period <- snakemake@params[['period']]

print("Running MetaCycle JTK with following parameters:")
print(input_file)
print(output_dir)
print(timepoints)
print(period)

library(MetaCycle)
print("Loaded MetaCycle")

if (period == "default") {
    # Run with default period range
    meta2d(input_file, output_dir, filestyle="txt", timepoints, cycMethod=c("JTK"))
} else {
    period <- as.integer(period)
    num_timepoints <- length(na.exclude(unique(timepoints %% 24)))
    if ((period == 8) & (num_timepoints < 6)) {
        # Some studies have too few timepoints to run period 8
        # we output all non-significant for them
        output_file <- paste(output_dir, "/JTKresult_expression.tpm.txt", sep='');
        input <- read_tsv(input_file);
        output <- tibble(CycID = input$Name, BH.Q = 1, ADJ.P = 1, PER = 0, LAG = 0, AMP = 0);
        print(output)
        print(output_file)
        write_tsv(output, output_file)
    }  else {
        # Run study through meta2d to get JTK results
        meta2d(input_file, output_dir, filestyle="txt", timepoints, cycMethod=c("JTK"), minper=period, maxper=period)
    }
}
