rule raw_read_lengths_json:
    input:
        f"{OUTDIR}/{{sample}}_concat.fastq.gz",
    output:
         f"{OUTDIR}/{{sample}}_raw_reads_lengths.json",
    conda:
        "../envs/deduplication.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_raw_rlj.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_raw_rlj_snakemake.log", 
    shell:
        "echo {input}; "
        "python workflow/scripts/read_lengths_json.py "
        "{input} "
        "> {output}"

rule non_host_read_lengths_json:
    input:
        f"{OUTDIR}/{{sample}}_non_host.fastq.gz",
    output:
         f"{OUTDIR}/{{sample}}_non_host_reads_lengths.json",
    conda:
        "../envs/deduplication.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_non_host_rlj.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_non_host_rlj_snakemake.log", 
    shell:
        "echo {input}; "
        "python workflow/scripts/read_lengths_json.py "
        "{input} "
        "> {output}"

