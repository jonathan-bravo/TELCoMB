rule bin_reads_by_length:
    input:
        f"{OUTDIR}/{{sample}}_trimmed.fastq.gz",
    output:
        touch(f"{OUTDIR}/{{sample}}_bin_reads.done"),
    params:
        outdir=f"{OUTDIR}/{{sample}}_read_bins",
    conda:
        "../envs/deduplication.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_bin_reads.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_bin_reads_snakemake.log",
    shell:
        "mkdir -p {params.outdir}; "
        "python workflow/scripts/bin_reads.py "
        "--infile {input} "
        "--outdir {params.outdir}"


rule cluster_reads:
    input:
        f"{OUTDIR}/{{sample}}_bin_reads.done",
    output:
        touch(f"{OUTDIR}/{{sample}}_cluster_reads.done"),
    params:
        similarity_threshold=config["similarity_threshold"],
        indir=f"{OUTDIR}/{{sample}}_read_bins",
        outdir=f"{OUTDIR}/{{sample}}_read_clusters",
    conda:
        "../envs/deduplication.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_cluster_reads.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_cluster_reads_snakemake.log",
    shell:
        "mkdir -p {params.outdir}; "
        "python workflow/scripts/cluster_reads.py "
        "--threshold {params.similarity_threshold} "
        "--indir {params.indir} "
        "--outdir {params.outdir}; "
        "rm -rf {params.indir}; "
        "rm {input}"


rule blat_clustered_reads:
    input:
        f"{OUTDIR}/{{sample}}_cluster_reads.done",
    output:
        touch(f"{OUTDIR}/{{sample}}_blat.done"),
    params:
        rc=f"{OUTDIR}/{{sample}}_read_clusters/",
        o=f"{OUTDIR}/{{sample}}_psl_files/",
    threads: 10
    conda:
        "../envs/deduplication.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_blat.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_blat_reads_snakemake.log",
    shell:
        "mkdir -p {params.o}; "
        "python workflow/scripts/run_blat.py "
        "--outdir {params.o} "
        "--threads {threads} "
        "--read_clusters {params.rc}; "
        "rm -rf {params.rc}; "
        "rm {input}"


rule find_duplicates:
    input:
        f"{OUTDIR}/{{sample}}_blat.done",
    output:
        touch(f"{OUTDIR}/{{sample}}_find_duplcates.done"),
    params:
        similarity_threshold=SIMILARITY,
        pls_dir=f"{OUTDIR}/{{sample}}_psl_files/",
        outdir=f"{OUTDIR}/{{sample}}_duplicate_txts/",
    conda:
        "../envs/deduplication.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_find_dups.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_find_dups_snakemake.log",
    shell:
        "mkdir -p {params.outdir}; "
        "bash workflow/scripts/run_find_duplicates.sh "
        "{params.outdir} "
        "{params.pls_dir} "
        "{params.similarity_threshold}; "
        "rm -rf {params.pls_dir}; "
        "rm {input}"


rule merge_duplicates_lists:
    input:
        f"{OUTDIR}/{{sample}}_find_duplcates.done",
    output:
        f"{OUTDIR}/{{sample}}_duplicates.txt",
    params:
        indir=f"{OUTDIR}/{{sample}}_duplicate_txts",
    conda:
        "../envs/default.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_merge_dups_lists.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_merge_dup_lists_snakemake.log",
    shell:
        "cat {params.indir}/* > {output}; "
        "rm -rf {params.indir}"


rule deduplicate:
    input:
        reads=f"{OUTDIR}/{{sample}}_trimmed.fastq.gz",
        duplicates_list=f"{OUTDIR}/{{sample}}_duplicates.txt",
    output:
        reads=f"{OUTDIR}/{{sample}}_dedup.fastq.gz",
        dupes=f"{OUTDIR}/{{sample}}_dup_reads.fastq.gz",
    conda:
        "../envs/deduplication.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_dedup.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_dedup_snakemake.log",
    shell:
        "python workflow/scripts/deduplicate.py "
        "--reads {input.reads} "
        "--duplicates {input.duplicates_list} "
        "--out_reads {output.reads} "
        "--out_dupes {output.dupes}"
