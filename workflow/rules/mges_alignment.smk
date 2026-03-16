rule align_to_mges:
    input:
        samples=f"{OUTDIR}/{{sample}}_non_host.fastq.gz",
        mges_database=ancient(f"{DATABASES}/mges_combined.fasta"),
    output:
        temp(f"{OUTDIR}/{{sample}}_ato_mges.sam"),
    conda:
        "../envs/alignment.yaml"
    threads: 32
    benchmark:
        f"{BENCHDIR}/{{sample}}_ato_mges.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_ato_mges.log",
    shell:
        "minimap2 -Y "
        "--secondary=no "
        "-t {threads} "
        "-ax map-ont " 
        "{input.mges_database} "
        "{input.samples} "
        "-o {output}"
