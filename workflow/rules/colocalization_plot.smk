# Going to have to change code here as well
rule colocalization_visualizations_notebook:
    input:
#        megares_db=f"{DATABASES}/megares_database_v3.00.fasta",
#        megares_annotation=f"{DATABASES}/megares_annotations_v3.00.csv",
#        mges_db=f"{DATABASES}/mges_combined.fasta", 
        reads_length=f"{OUTDIR}/{{sample}}_non_host_reads_lengths.json",
        colocalizations=f"{OUTDIR}/{{sample}}_deduplicated_colocalizations.csv",
#        config_file=f"config/config.ini",
    output:
       touch(f"{OUTDIR}/{{sample}}_colocalizations_plots.pdf"),
    conda:
        "../envs/plots.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_colocalizations_plot.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_colocalizations_plot_snakemake.log",
    shell:
        "python workflow/scripts/colocalizations_notebook_modified.py "
#        "--config_file {input.config_file} "
        "--read_lengths {input.reads_length} "
        "--colocalizations {input.colocalizations} "
        "--output_plot {output}"
    # notebook:
    #     "workflow/notebooks/colocalizations_notebook.py.ipynb"
