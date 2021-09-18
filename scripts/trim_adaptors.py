import subprocess
import pathlib

fastq_dir = pathlib.Path(snakemake.input[0])
trimmed_dir = pathlib.Path(snakemake.output[0])
trimmed_dir.mkdir(exist_ok=True)

for input_fastq in fastq_dir.glob("*.fastq"):
    output_fastq = trimmed_dir / (input_fastq.name)
    cmd = f"cutadapt {snakemake.params.args} -o {output_fastq} {input_fastq}"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
