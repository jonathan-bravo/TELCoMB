rule hard_trim_count:
    input:
        f"{OUTDIR}/{{sample}}_trimmed.log",
    output:
        f"{OUTDIR}/{{sample}}_hard_trim_count.txt",
    conda:
        "../envs/default.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_hard_trim_count.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_hard_trim_count_snakemake.log",
    shell:
        "cat {input} | "
        "awk '{{ if ($0 ~ /short/) print }}' | "
        "wc -l > {output}"


rule chimeric_count:
    input:
        f"{OUTDIR}/{{sample}}_trimmed.fastq.gz",
    output:
        f"{OUTDIR}/{{sample}}_chimeric_count.txt",
    conda:
        "../envs/default.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_chimer_count.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_chimer_count_snakemake.log",
    shell:
        "gunzip -c {input} | "
        "awk '{{ if ($0 ~ /_A/ || $0 ~ /_B/) print $0 }}' | "
        "wc -l > {output}"
