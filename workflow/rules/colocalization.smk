rule find_colocalizations:
    input:
        samples=f"{OUTDIR}/{{sample}}_non_host.fastq.gz", 
        megares_sam=f"{OUTDIR}/{{sample}}_ato_megares.sam",
        mges_sam=f"{OUTDIR}/{{sample}}_ato_mges.sam",
        overlap=f"{OUTDIR}/{{sample}}_deduplicated_overlapped_mges.csv",
        config_file=f"config/config.ini",
    output:
        colocalizations=f"{OUTDIR}/{{sample}}_deduplicated_colocalizations.csv",
        genes_list=temp(f"{OUTDIR}/{{sample}}_non_host.fastq.gz_genes_list.xlsx"),
    params:
        output_directory=f"{OUTDIR}", 
    conda:
        "../envs/pipeline.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_find_colocalizations.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_find_colocalizations_snakemake.log",
    shell:
        "python workflow/scripts/find_colocalizations.py "
        "-r {input.samples} "
        "--arg {input.megares_sam} "
        "--mge {input.mges_sam} "
        "-s {input.overlap} "
        "-c {input.config_file} "
        "-o {params.output_directory} "
        "> {output.colocalizations}"

rule colocalization_richness:
    input:
        colocalizations=f"{OUTDIR}/{{sample}}_deduplicated_colocalizations.csv",
        config_file=f"config/config.ini",
    output:
        f"{OUTDIR}/{{sample}}_deduplicated_colocalizations_richness.csv",
    conda:
        "../envs/pipeline.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_colocalizations_richness.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_colocalizations_richness_snakemake.log",
    shell:
        "python workflow/scripts/colocalization_richness.py "
        "-i {input.colocalizations} "
        "-c {input.config_file} "
        "> {output}"
