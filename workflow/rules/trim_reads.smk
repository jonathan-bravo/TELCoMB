rule trim_reads:
    input:
        f"{OUTDIR}/{{sample}}_concat.fastq.gz",
    output:
        logfile=f"{OUTDIR}/{{sample}}_trimmed.log",
        trimmed_reads=f"{OUTDIR}/{{sample}}_trimmed.fastq.gz",
    params:
        crop=CROPLEN,
        barcodes=BARCODES,
    conda:
        "../envs/trimming.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_trim_reads.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_trim_reads_snakemake.log",
    shell:
        "python workflow/scripts/trim.py "
        "--infile {input} "
        "--logfile {output.logfile} "
        "--outfile {output.trimmed_reads} "
        "--barcodes {params.barcodes} "
        "--crop {params.crop}"
