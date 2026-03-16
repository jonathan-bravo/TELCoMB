rule start_read_count:
    input:
        f"{OUTDIR}/{{sample}}_concat.fastq.gz",
    output:
        f"{OUTDIR}/{{sample}}_start_read_count.txt",
    conda:
        "../envs/default.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_read_counts.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_read_count_snakemake.log",
    shell:
        "bash workflow/scripts/start_counts.sh {input} {output}"
