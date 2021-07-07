input_file <- snakemake@input[[1]]
output_dir <- snakemake@params[['out_dir']]

timepoints <- snakemake@params[['timepoints']]

print("Running MetaCycle JTK with following parameters:")
print(input_file)
print(output_dir)
print(timepoints)

library(MetaCycle)
print("Loaded MetaCycle")

meta2d(input_file, output_dir, filestyle="txt", timepoints, cycMethod=c("JTK"))
