# TELSeq
A workflow for contextualization of antibiotic resistance in microbiomes.

### Requirements 
We manage our dependencies trough [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Install Snakemake and Clone the Repository

Create the environment for `telseq` and install snakemake.

```bash
conda create -c conda-forge -n mamba_base mamba
conda activate mamba_base
mamba create -c conda-forge -c bioconda -n telseq snakemake git
mamba activate telseq
```

Clone the repository.

```bash
mamba activate telseq
git clone https://github.com/marco-oliva/telseq.git
```

### Usage on local desktop

The `telseq` workflows assumes that the `fastq` files will be stored in a directory called `samples` in the working directory. Here we show the usage and the directories structure that can be used with the default `config.json` file.

If al the databases are already available it is possible to avoid re-downloading them by specifying the directory in the config file. The names of the database have to be as in the following table. There is no need to copy them, a soft link is sufficient (`ln -s`).

| Database               | File name in `DATABASES_DIR`         |
|------------------------|--------------------------------------|
| KEGG genes             | `kegg_genes.fasta`                   |
| MGEs combined database | `mges_combined.fasta`                |
| MEGARes Database       | `megares_full_database_v2.00.fasta`  |
| MEGARes Ontology       | `megares_full_annotations_v2.00.csv` |



```bash
cd telseq
mamba activate telseq

# Create the directories structure
mkdir -p work_dir/samples work_dir/logs 

# Move the fastq files in the samples directory
mv <your_data>.fastq work_dir/samples

# Run the workflow
snakemake -j <number of threads available> --use-conda --conda-frontend mamba
```

#### Plots
```bash
cd telseq

# Spawn the jupyter notebook for the violin plot
snakemake -j1 --use-conda --conda-frontend mamba --edit-notebook violin_plot_all_samples.pdf

# Spawn the jupyter notebook for the resitome heatmap
snakemake -j1 --use-conda --conda-frontend mamba --edit-notebook heatmap_all_samples.pdf

# Spawn the jupyter notebook for the colocalizations plot for a specific sample
snakemake -j1 --use-conda --conda-frontend mamba --edit-notebook <your sample name>_colocalizations_plot.pdf"
```

### Usage on slurm cluster

Edit the `cluster.json` file in order to fit your resources.

```bash
cd telseq

# Create the directories structure
mkdir -p work_dir/samples work_dir/logs 

# Move the fastq files in the samples directory
mv <your_data>.fastq work_dir/samples

# Run the workflow
mkdir -p logs
sbatch run.sh
```

#### Plots
This will run on a node in the cluster so you will need an ssh tunnel to the node to be able to connect to the jupyter notebook.

```bash
cd telseq

# Spawn the jupyter notebook for the violin plot
snakemake -j1 --use-conda --conda-frontend mamba --edit-notebook violin_plot_all_samples.pdf

# Spawn the jupyter notebook for the resitome heatmap
snakemake -j1 --use-conda --conda-frontend mamba --edit-notebook heatmap_all_samples.pdf

# Spawn the jupyter notebook for the colocalizations plot for a specific sample
snakemake -j1 --use-conda --conda-frontend mamba --edit-notebook <your sample name>_colocalizations_plot.pdf
```

