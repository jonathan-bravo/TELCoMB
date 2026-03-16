rule heatmap_notebook:
    input:
        megares_db=f"{DATABASES}/megares_database_v3.00.fasta",
        megares_annotation=f"{DATABASES}/megares_annotations_v3.00.csv",
        config_file=f"config/config.ini",
        data=expand("{output}" + "/" + "{sample_name}" + "_deduplicated_amr_features.csv", sample_name=SAMPLES,output=OUTDIR),
    output:
        touch(f"{OUTDIR}/heatmap_all_samples.svg"),
    params:
        samples_list=SAMPLES,
    conda:
        "../envs/plots.yaml"
    benchmark:
        f"{BENCHDIR}/heatmap_plot.benchmark"
    log:
        f"{LOGDIR}/heatmap_plot_snakemake.log",
    shell:
        "python workflow/scripts/heatmap_notebook.py "
        "--dedup_string _deduplicated "
        "--config_file {input.config_file} "
        "--output_plot {output} "
        "--samples_list {params.samples_list} "
    # notebook:
    #     "workflow/notebooks/heatmap_notebook.py.ipynb"
