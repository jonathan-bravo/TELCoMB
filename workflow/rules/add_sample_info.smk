rule add_sample_info:
    input:
        f"{OUTDIR}/{{sample}}_pre_dedup_rl.tsv",
        f"{OUTDIR}/{{sample}}_post_dedup_rl.tsv",
        f"{OUTDIR}/{{sample}}_reads_per_strain.tsv",
        f"{OUTDIR}/{{sample}}_reads_per_strain_filtered.tsv",
        f"{OUTDIR}/{{sample}}_viral_targets.log",
        f"{OUTDIR}/{{sample}}_selected_viral_targets.log",
    output:
        touch(f"{OUTDIR}/{{sample}}_add_sample_info.done"),
    params:
        run_id=RUNID,
        sample_id="{sample}",
    conda:
        "../envs/default.yaml"
    benchmark:
        f"{BENCHDIR}/{{sample}}_add_sample_info.benchmark"
    log:
        f"{LOGDIR}/{{sample}}_sample_info_snakemake.log",
    shell:
        "python workflow/scripts/add_run_info.py "
        "--sample_id {params.sample_id} "
        "--run_id {params.run_id} "
        "--infiles {input}"
