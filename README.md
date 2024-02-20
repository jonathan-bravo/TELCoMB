# TELCoMB
A workflow for contextualization of antibiotic resistance in microbiomes.

### Requirements 
We manage our dependencies trough [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Install Snakemake and Clone the Repository

Create the environment for `telcomb` and install snakemake.

```bash
conda create -c conda-forge -c bioconda -n telcomb snakemake git
```

Clone the repository.

```bash
conda activate telcomb
git clone https://github.com/jonathan-bravo/TELCoMB
```

### Usage on local desktop

The `telcomb` workflows assumes that the `fastq` files will be stored in a directory called `samples` in the working directory. Here we show the usage and the directories structure that can be used with the default `config.json` file.

If all the databases are already available it is possible to avoid re-downloading them by specifying the directory in the config file. The names of the database have to be as in the following table. There is no need to copy them, a soft link is sufficient (`ln -s`).

| Database               | File name in `DATABASES_DIR`         |
|------------------------|--------------------------------------|
| MGEs combined database | `mges_combined.fasta`                |
| MEGARes Database       | `megares_full_database.fasta`        |
| MEGARes Ontology       | `megares_full_annotations.csv`       |



```bash
cd TELCoMB
mamba activate telcomb

# Create the directories structure
mkdir -p work_dir/samples work_dir/logs 

# Move the fastq files in the samples directory
mv <your_data>.fastq work_dir/samples

# Run the workflow
snakemake -c <number of threads available> --use-conda --conda-frontend conda
```

### Usage on slurm cluster

Edit the `cluster.json` file in order to fit your resources.

```bash
cd TELCoMB

# Create the directories structure
mkdir -p work_dir/samples work_dir/logs 

# Move the fastq files in the samples directory
mv <your_data>.fastq work_dir/samples

# Run the workflow
mkdir -p logs
sbatch run.sh
```