rule read_lengths_plot:
    input:
        sample=f"{OUTDIR}/{{sample}}_non_host.fastq.gz",
    output:
        touch(f"{OUTDIR}/{{sample}}_read_lengths_hist.pdf"),
    params:
        num_of_bins=100,
        std_deviations=4,
    conda:
        "../envs/plots.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_rl_plot.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_rl_plot_snakemake.log",
    shell:
        "python workflow/scripts/plot_read_lengths.py "
        "-i {input} "
        "-s {params.std_deviations} "
        "-b {params.num_of_bins} "
        "-o {output} "
        "--title {input.sample}_deduplicated"
