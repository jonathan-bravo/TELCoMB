rule concat_parts:
    input:
        f"{READS}/{{sample}}",
    output:
        temp(f"{OUTDIR}/{{sample}}_concat.fastq.gz"),
    conda:
        "../envs/default.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_concat_parts.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_concat_parts_snakemake.log",
    shell:
        "cat {input}/* > {output}"
