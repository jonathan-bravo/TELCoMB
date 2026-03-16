rule align_reads_to_host:
    input:
        host=HOST_FILE,
        samples=f"{OUTDIR}/{{sample}}_dedup.fastq.gz",
    output:
        temp(f"{OUTDIR}/{{sample}}_host.sam"),
    conda:
        "../envs/alignment.yaml"
    threads: 32
    benchmark:
        f"{BENCHDIR}/{{sample}}_ato_host.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_ato_host_snakemake.log",
    shell:
        "minimap2 --secondary=no "
        "-t {threads} "
        "-o {output} "
        "-a {input.host} {input.samples}"

rule host_sam_to_bam:
    input:
        f"{OUTDIR}/{{sample}}_host.sam",
    output:
        temp(f"{OUTDIR}/{{sample}}_host_sorted.bam"),
    conda:
        "../envs/alignment.yaml"
    threads: 32
    benchmark:
        f"{BENCHDIR}/{{sample}}_host_sam2bam.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_host_sam2bam_snakemake.log",
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {output}"


rule remove_host_dna:
    input:
        f"{OUTDIR}/{{sample}}_host_sorted.bam",
    output:
        idx=temp(f"{OUTDIR}/{{sample}}_remove_host_samtools.idxstats"),
        bam=temp(f"{OUTDIR}/{{sample}}_host_removed_sorted.bam"),
    conda:
        "../envs/alignment.yaml"
    threads: 10
    benchmark:
        f"{BENCHDIR}/{{sample}}_host_removal.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_host_removal_snakemake.log",
    shell:
        "samtools index {input} && "
        "samtools view -h -F 2048 -b {input} | "
        "samtools idxstats - > {output.idx}; "
        "samtools view -h -f 4 -b {input} -o {output.bam}; "
        "rm {input}.bai"


rule non_host_reads:
    input:
        f"{OUTDIR}/{{sample}}_host_removed_sorted.bam",
    output:
        f"{OUTDIR}/{{sample}}_non_host.fastq.gz",
    conda:
        "../envs/alignment.yaml"
    threads: 4
    benchmark:
        f"{BENCHDIR}/{{sample}}_non_host_reads.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_non_host_reads_snakemake.log",
    shell:
        "samtools fastq {input} | bgzip -@ {threads} -c > {output}"
