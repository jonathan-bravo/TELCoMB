# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

################################################################################
#                                Imports                                       #
################################################################################

import os, ast, glob

################################################################################
#                             Long Read Steps                                  #
################################################################################

############################## Read Lengths ####################################

rule read_lengths:
    input:
        samples_dir + "/{sample_name}.fastq"
        # "{sample_name}"
    output:
         "{sample_name}" + config["EXTENSION"]["READS_LENGTH"]
    params:
        read_lengths_script = workflow.basedir + "/" + config["SCRIPTS"]["READS_LENGTH"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["DEDUP"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.rl.benchmark" 
    shell:
        "echo {input}; "
        "python {params.read_lengths_script} {input} > {output}"

########################### Deduplication ######################################
    
rule bin_reads_by_length:
    input:
        samples_dir + "/{sample_name}.fastq"
        # "{sample_name}"
    output:
        touch(tmp_dir + "/{sample_name}.bin.reads.done")
    params:
        outdir = tmp_dir + "/{sample_name}_read_bins",
        bin_script = workflow.basedir + "/" + config["SCRIPTS"]["BIN_READS"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["DEDUP"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.bin_reads.benchmark" 
    shell:
        "mkdir -p {params.outdir}; "
        "python {params.bin_script} "
        "--infile {input} "
        "--outdir {params.outdir}"

rule cluster_reads:
    input:
        tmp_dir + "/{sample_name}.bin.reads.done"
    output:
        touch(tmp_dir + "/{sample_name}.cluster.reads.done")
    params:
        indir = tmp_dir + "/{sample_name}_read_bins",
        outdir = tmp_dir + "/{sample_name}_read_clusters",
        cluster_script = workflow.basedir + "/" + config["SCRIPTS"]["CLUSTER_READS"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["DEDUP"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.cluster_reads.benchmark" 
    shell:
        "mkdir -p {params.outdir}; "
        "python {params.cluster_script} "
        "--indir {params.indir} "
        "--outdir {params.outdir}; "
        "rm -rf {params.indir}; "
        "rm {input}"

rule blat_clustered_reads:
    input:
        tmp_dir + "/{sample_name}.cluster.reads.done"
    output:
        touch(tmp_dir + "/{sample_name}.blat.done")
    params:
        rc = tmp_dir + "/{sample_name}_read_clusters/",
        o = tmp_dir + "/{sample_name}_psl_files/",
        blat_script = workflow.basedir + "/" + config["SCRIPTS"]["RUN_BLAT"]
    threads:
        config["BLAT"]["THREADS"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["DEDUP"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.blat.benchmark" 
    shell:
        "mkdir -p {params.o}; "
        "python {params.blat_script} "
        "--outdir {params.o} "
        "--threads {threads} "
        "--read_clusters {params.rc}; "
        "rm -rf {params.rc}; "
        "rm {input}"

rule find_duplicates:
    input:
        tmp_dir + "/{sample_name}.blat.done"
    output:
        touch(tmp_dir + "/{sample_name}.find.duplcates.done")
    params:
        similarity_threshold = config["MISC"]["DEDUPLICATION_SIMILARITY_THRESHOLD"],
        pls_dir = tmp_dir + "/{sample_name}_psl_files/",
        outdir = tmp_dir + "/{sample_name}_duplicate_txts/",
        run_find_dups_script = workflow.basedir + "/" + config["SCRIPTS"]["RUN_FIND_DUPLICATES"],
        find_dupes_script = workflow.basedir + "/" + config["SCRIPTS"]["FIND_DUPLICATES"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["DEDUP"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.find_dupes.benchmark" 
    shell:
        "mkdir -p {params.outdir}; "
        "{params.run_find_dups_script} "
        "{params.outdir} "
        "{params.pls_dir} "
        "{params.similarity_threshold} "
        "{params.find_dupes_script};"
        "rm -rf {params.pls_dir}; "
        "rm {input}"

rule merge_duplicates_lists:
    input:
        tmp_dir + "/{sample_name}.find.duplcates.done"
    output:
        tmp_dir + "/{sample_name}.duplicates.txt"
    params:
        indir = tmp_dir + "/{sample_name}_duplicate_txts/"
    conda:
        workflow.basedir + "/" + config["CONDA"]["DEDUP"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.merge_dupes.benchmark" 
    shell:
        "cat {params.indir}/* > {output}; "
        "rm -rf {params.indir}"

# change output so that it is NOT gzipped *********
rule deduplicate:
    input:
        reads = samples_dir + "/{sample_name}.fastq",
        duplicates_list = tmp_dir + "/{sample_name}.duplicates.txt"
    output:
        reads =  "{sample_name}" + config["EXTENSION"]["DEDUPLICATED"],
        dupes = "{sample_name}.dup.reads.fastq"
    params:
        dedup_script = workflow.basedir + "/" + config["SCRIPTS"]["DEDUPLICATE"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["DEDUP"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.dedup.benchmark" 
    shell:
        "python {params.dedup_script} "
        "--reads {input.reads} "
        "--duplicates {input.duplicates_list} "
        "--out_reads {output.reads} "
        "--out_dupes {output.dupes}"
        
rule not_deduplicated_reads:
    input:
        samples_dir + "/{sample_name}.fastq"
    output:
        "{sample_name}" + config["EXTENSION"]["NOT_DEDUPLICATED"]
    shell:
        "ln -s {input} {output}"