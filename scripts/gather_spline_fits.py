import pathlib
import pandas

batch_dir = pathlib.Path("results/Liver/spline_fit/batches")
file_types = ["summary.txt", "curves_fit.txt", "curves_pstd.txt", "re.txt", "re_structure.txt",]
for file_type in file_types:
    print(f"Processig {file_type}")
    data_list = []
    for path in batch_dir.glob(f"*.{file_type}"):
        data_list.append(pandas.read_csv(path, sep="\t"))
    df = pandas.concat(data_list, axis=0)
    print(f"Found {len(data_list)} files, with total of {len(df)} lines")
    dupe_columns = df.columns.intersection(["gene", "study", "var"])
    df.drop_duplicates(dupe_columns, inplace=True)
    df.sort_values(by="gene", inplace=True)
    df.to_csv(f"results/Liver/spline_fit/{file_type}", sep="\t", index=False)
