rule overlap:
    input:
        samples=f"{OUTDIR}/{{sample}}_non_host.fastq.gz", 
        megares_sam=f"{OUTDIR}/{{sample}}_ato_megares.sam",
        mges_sam=f"{OUTDIR}/{{sample}}_ato_mges.sam",
        reads_length=f"{OUTDIR}/{{sample}}_non_host_reads_lengths.json",
        config_file=f"config/config.ini",
    params:
        output_prefix=f"{OUTDIR}/{{sample}}_deduplicated",
    conda:
        "../envs/pipeline.yaml"
    output:
        overlap=f"{OUTDIR}/{{sample}}_deduplicated_overlapped_mges.csv",
    benchmark:
        f"{BENCHDIR}/{{sample}}_overlapped_mges.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_overlapped_mges.log",
    shell:
        "python workflow/scripts/find_overlap.py "
        "-r {input.samples} "
        "-a {input.megares_sam} "
        "-m {input.mges_sam} "
        "-c {input.config_file} "
        "-o {params.output_prefix}"

rule merge_overlap_info:
    input:
        #f"{OUTDIR}/{{sample}}_deduplicated_overlapped_mges.csv",
        data=expand("{output}" + "/" + "{sample_name}" + "_deduplicated_overlapped_mges.csv", sample_name=SAMPLES,output=OUTDIR),
    output: 
        merged_info=f"{OUTDIR}/merged_overlapped_mges_info.csv",
    benchmark:
        f"{BENCHDIR}/merge_overlap.benchmark"
    log:
        f"{LOGDIR}/merge_overlap.log",
    run:
        import csv
        merged_overlaped_mges_info = set()
        for overlap_info_filename in input:
            with open(overlap_info_filename) as overlap_info_file:
                overlap_reader = csv.reader(overlap_info_file, delimiter=',')
                for row in overlap_reader:
                    merged_overlaped_mges_info.add(row[0])
        with open(output[0], 'w') as merged:
            merged_writer = csv.writer(merged, delimiter=',')
            for mge in merged_overlaped_mges_info:
                merged_writer.writerow([mge])
