rule resistome_and_mobilome:
    input:
        samples=f"{OUTDIR}/{{sample}}_non_host.fastq.gz", 
        megares_sam=f"{OUTDIR}/{{sample}}_ato_megares.sam",
        mges_sam=f"{OUTDIR}/{{sample}}_ato_mges.sam",
        reads_lengths=f"{OUTDIR}/{{sample}}_non_host_reads_lengths.json",
        overlap=f"{OUTDIR}/{{sample}}_deduplicated_overlapped_mges.csv",
        reads_lengths_json=f"{OUTDIR}/{{sample}}_non_host_reads_lengths.json",
        config_file=f"config/config.ini",
    output:
        resistome_richness=touch(f"{OUTDIR}/{{sample}}_deduplicated_amr_richness.csv"),
        resistome_diversity=touch(f"{OUTDIR}/{{sample}}_deduplicated_amr_features.csv"),
        mobilome=touch(f"{OUTDIR}/{{sample}}_deduplicated_mobilome.csv"),
    params:
        output_prefix=f"{OUTDIR}/{{sample}}_deduplicated",    
    conda:
        "../envs/pipeline.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_resistome_and_mobilome.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_resistome_and_mobilome_snakemake.log",
    shell:
        "python workflow/scripts/gen_resistome_and_mobilome.py "
        "-r {input.samples} "
        "-a {input.megares_sam} "
        "-m {input.mges_sam} "
        "-s {input.overlap} "
        "-c {input.config_file} "
        "-o {params.output_prefix} "
        "-rl {input.reads_lengths_json}"
