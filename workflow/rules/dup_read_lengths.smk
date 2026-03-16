rule pre_dedup_read_lengths:
    input:
        f"{OUTDIR}/{{sample}}_trimmed.fastq.gz",
    output:
        f"{OUTDIR}/{{sample}}_pre_dedup_rl.tsv",
    conda:
        "../envs/deduplication.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_pre_dedup_rl.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_pre_dedup_rl_snakemake.log",
    shell:
        "python workflow/scripts/read_lengths.py "
        "--infile {input} "
        "--outfile {output}"

rule post_dedup_read_lengths:
    input:
        f"{OUTDIR}/{{sample}}_dedup.fastq.gz",
    output:
        f"{OUTDIR}/{{sample}}_post_dedup_rl.tsv",
    conda:
        "../envs/deduplication.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_post_dedup_rl.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_post_dedup_rl_snakemake.log",
    shell:
        "python workflow/scripts/read_lengths.py "
        "--infile {input} "
        "--outfile {output}"
