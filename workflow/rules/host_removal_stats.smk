rule host_removal_stats:
    input:
        expand(f"{OUTDIR}/{{sample}}_remove_host_samtools.idxstats", sample=SAMPLES),
    output:
        temp(f"{OUTDIR}/host_removal_stats.tsv"),
    conda:
        "../envs/alignment.yaml"
    benchmark:
        f"{BENCHDIR}/host_removal_stats.benchmark"
    log:
        f"{LOGDIR}/host_removal_stats_snakemake.log",
    shell:
        "python workflow/scripts/samtools_idxstats.py -i {input} -o {output}"
