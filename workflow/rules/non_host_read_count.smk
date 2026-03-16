rule start_host_read_count:
    input:
        f"{OUTDIR}/{{sample}}_non_host.fastq.gz",
    output:
        f"{OUTDIR}/{{sample}}_non_host_read_count.txt",
    conda:
        "../envs/default.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_non_host_read_counts.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_non_host_read_count_snakemake.log",
    shell:
        "bash workflow/scripts/start_counts.sh {input} {output}"
