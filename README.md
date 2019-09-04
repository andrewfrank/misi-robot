# misi-robot

A short python script to automatically submit jobs to the JGI MiSI tool.

## Usage

Query and reference fasta files should contain the predicted genes (DNA) of 
specified genomes. No tRNA or rRNA features should be included. No fasta files 
greater than 500MB in size will be accepted.

If your output directory contains existing misi-robot results, new results will
be appended to the existing results file.

```{python}
misi-robot.py [-h] [--start NUMBER] qry_dir_path ref_dir_path out_dir_path

positional arguments:
  qry_dir_path    The full path to a directory containing query genomes.
  ref_dir_path    The full path to a directory containing reference genomes.
  out_dir_path    The full path to a directory for the output.

optional arguments:
  -h, --help      show this help message and exit
  --start NUMBER  Start submitting jobs at this job number. Useful for
                  restarting this script. Defaults to 1.
```

### Example Command

```{bash}
python "D:\scripts\misi-robot.py" "E:\experiments\query_genomes" "E:\experiments\reference_genomes" "E:\experiments\output"
```

### Example for Restarting a Submission Job

```{bash}
python "D:\scripts\misi-robot.py" --start 145 "E:\experiments\query_genomes" "E:\experiments\reference_genomes" "E:\experiments\output"
```
