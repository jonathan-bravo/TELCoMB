rule amr_mges_on_target:
    input:
       sample_filename=f"{OUTDIR}/{{sample}}_start_read_count.txt",
       amr_file=f"{OUTDIR}/{{sample}}_deduplicated_amr_features.csv",
       mges_file=f"{OUTDIR}/{{sample}}_deduplicated_mobilome.csv",
       raw_read_file=f"{OUTDIR}/{{sample}}_start_read_count.txt",
       trim_file=f"{OUTDIR}/{{sample}}_pre_dedup_rl.tsv",
       hard_trim_file=f"{OUTDIR}/{{sample}}_hard_trim_count.txt",
       chimeric_file=f"{OUTDIR}/{{sample}}_chimeric_count.txt",
       dedup_file=f"{OUTDIR}/{{sample}}_post_dedup_rl.tsv",
       duplicate_file=f"{OUTDIR}/{{sample}}_duplicates.txt",
       non_host_file=f"{OUTDIR}/{{sample}}_non_host_read_count.txt",
       flowcell=FLOWCELL,
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
        "{input.trim_file} "
        "{input.hard_trim_file} "
        "{input.chimeric_file} "
        "{input.dedup_file} "
        "{input.duplicate_file} "
        "{input.non_host_file} "
        "{input.flowcell} "
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
        "sed '1i flowcell_id,barcode,ARGs,MGEs,raw_read_count,trim_read_count,hard_trim_read_count,chimeric_read_count,dedup_read_count,duplicate_read_count,non_host_read_rount,raw_read_on_target_amr,raw_read_on_target_mges,non_host_read_on_target_amr,non_host_read_on_target_mges' {input} > {output}"

