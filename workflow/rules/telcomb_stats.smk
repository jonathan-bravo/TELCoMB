rule amr_mges_on_target:
    input:
       sample_filename=f"{OUTDIR}/{{sample}}_start_read_count.txt",
       amr_file=f"{OUTDIR}/{{sample}}_deduplicated_amr_features.csv",
       mges_file=f"{OUTDIR}/{{sample}}_deduplicated_mobilome.csv",
       raw_read_file=f"{OUTDIR}/{{sample}}_start_read_count.txt",
       non_host_file=f"{OUTDIR}/{{sample}}_non_host_read_count.txt",
    output:
        f"{OUTDIR}/{{sample}}_telcomb_stats.csv",
    conda:
        "../envs/default.yaml",
    benchmark:
        f"{BENCHDIR}/{{sample}}_telcomb_stats.benchmark",
    log:
        f"{LOGDIR}/{{sample}}_telcomb_stats.log",
    shell:
        "bash workflow/scripts/telcomb_stats.sh "
        "{input.sample_filename} "
        "{input.amr_file} "
        "{input.mges_file} "
        "{input.raw_read_file} "
        "{input.non_host_file} "
        "{output}"

rule cat_telcomb_stats:
    input:
        data=expand("{output}" + "/" + "{sample_name}" + "_telcomb_stats.csv", sample_name=SAMPLES,output=OUTDIR),
    output:
        f"{OUTDIR}/no_header_telcomb_stats.csv",
    conda:
        "../envs/default.yaml",
    benchmark:
        f"{BENCHDIR}/concat_telcomb_stats.benchmark",
    log:
        f"{LOGDIR}/concat_telcomb_stats.log",
    shell:
        "cat {input.data} > {output}" 

rule add_header_telcomb_stats:
    input:
        f"{OUTDIR}/no_header_telcomb_stats.csv",
    output:
        f"{OUTDIR}/telcomb_stats.csv",
    conda:
        "../envs/default.yaml",
    benchmark:
        f"{BENCHDIR}/header_telcomb_stats.benchmark",
    log:
        f"{LOGDIR}/header_telcomb_stats.log",
    shell:
        "sed '1i barcode,ARGs,MGEs,raw_reads,non_host_reads,raw_read_on_target_amr,raw_read_on_target_mges,non_host_read_on_target_amr,non_host_read_on_target_mges' {input} > {output}"

