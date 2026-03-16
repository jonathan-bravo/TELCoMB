rule viruses_alignment_stats:
    input:
        f"{OUTDIR}/{{sample}}_stats_viruses.sam",
    output:
        temp(f"{OUTDIR}/{{sample}}_all_viruses_samtools.idxstats"),
    params:
        cutoff=CUTOFF,
        bam=f"{OUTDIR}/{{sample}}_stats_viruses_sorted.bam",
        sftbam=f"{OUTDIR}/{{sample}}_stats_viruses_sorted_sftclp.bam",
    conda:
        "../envs/alignment.yaml"
    threads: 32
    benchmark:
        f"{BENCHDIR}/{{sample}}_stat_align.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_stat_align_snakemake.log",
    shell:
        "samtools view -h -F 2048 -Sb {input} | "
        "samtools sort -@ {threads} -o {params.bam}; "
        "samtools index -@ {threads} {params.bam}; "
        "python workflow/scripts/sftclipper.py "
        "--cutoff {params.cutoff} "
        "--bam {params.bam} "
        "--outfile {params.sftbam}; "
        "samtools index -@ {threads} {params.sftbam}; "
        "samtools idxstats {params.sftbam} > {output}; "
        "rm {params.bam} {params.bam}.bai; "
        "rm {params.sftbam} {params.sftbam}.bai"


rule merge_viral_alignment_stats:
    input:
        expand(f"{OUTDIR}/{{sample}}_all_viruses_samtools.idxstats", sample=SAMPLES),
    output:
        temp(f"{OUTDIR}/all_viruses_stats.tsv"),
    conda:
        "../envs/alignment.yaml"
    benchmark:
        f"{BENCHDIR}/merge_viral_stats.benchmark"
    log:
        f"{LOGDIR}/merge_viral_stats_snakemake.log",
    shell:
        "python workflow/scripts/samtools_idxstats.py -i {input} -o {output}"


rule on_target_stats:
    input:
        viral_stats=f"{OUTDIR}/all_viruses_stats.tsv",
        host_stats=f"{OUTDIR}/host_removal_stats.tsv",
    output:
        f"{OUTDIR}/on_target_stats.tsv",
    conda:
        "../envs/alignment.yaml"
    benchmark:
        f"{BENCHDIR}/on_target_stats.benchmark"
    log:
        f"{LOGDIR}/on_target_stats_snakemake.log",
    shell:
        "python workflow/scripts/merge_stats.py "
        "--host_stats {input.host_stats} "
        "--viral_stats {input.viral_stats} "
        "--outfile {output}"
