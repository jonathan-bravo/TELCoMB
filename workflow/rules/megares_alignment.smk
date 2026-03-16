rule align_to_megares:
    input:
        samples=f"{OUTDIR}/{{sample}}_non_host.fastq.gz",
        megares_seqs=ancient(f"{DATABASES}/megares_database_v3.00.fasta"),
    output:
        temp(f"{OUTDIR}/{{sample}}_ato_megares.sam"),
    conda:
        "../envs/alignment.yaml"
    threads: 32
    benchmark:
        f"{BENCHDIR}/{{sample}}_ato_megares.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_ato_megares.log",
    shell:
        "minimap2 -Y "
        "--secondary=no "
        "-t {threads} "
        "-ax map-ont "
        "{input.megares_seqs} "
        "{input.samples} "
        "-o {output}"

