# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

################################################################################
#                                Imports                                       #
################################################################################

import os, ast, glob

################################################################################
#                           Short Reads Steps                                  #
################################################################################

############################## Assembly ########################################

# Should split the meta_spades step from the fasta2fastq step
rule meta_spades_assembly:
    input:
        f_reads = samples_dir + "/{sample_name}_R1.fastq",
        r_reads = samples_dir + "/{sample_name}_R2.fastq"
    output:
        touch("{sample_name}_assembly.done")
    params:
        kmers = config["SPADES"]["KMERS"],
        memory = config["SPADES"]["MEMORY"],
        phred = config["SPADES"]["PHRED"],
        outdir = tmp_dir + "/spades_files/{sample_name}/",
    conda:
        workflow.basedir + "/" + config["CONDA"]["ASSEMBLY"]
    threads:
        config["SPADES"]["THREADS"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.assembly.benchmark" 
    shell:
        "spades.py --meta "
        "--phred-offset {params.phred} "
        "-t {threads} "
        "-m {params.memory} "
        "-k {params.kmers} "
        "-1 {input.f_reads} "
        "-2 {input.r_reads} "
        "-o {params.outdir}"


rule fasta2fastq:
    input:
        "{sample_name}_assembly.done"
    output:
        "{sample_name}" + config["EXTENSION"]["ASSEMBLED"],
    params:
        reads = tmp_dir + "/spades_files/{sample_name}/scaffolds.fasta",
        fasta2fastq_script = workflow.basedir + "/" + config["SCRIPTS"]["FASTA2FASTQ"],
    conda:
        workflow.basedir + "/" + config["CONDA"]["ASSEMBLY"]
    shell: 
        "python {params.fasta2fastq_script} "
        "--fasta {params.reads} "
        "--fastq {output}"

############################## Read Lengths ####################################

rule read_lengths:
    input:
        "{sample_name}" + config["EXTENSION"]["ASSEMBLED"]
    output:
         "{sample_name}" + config["EXTENSION"]["READS_LENGTH"]
    params:
        read_lengths_script = workflow.basedir + "/" + config["SCRIPTS"]["READS_LENGTH"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["ASSEMBLY"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.rl.benchmark" 
    shell:
        "python {params.read_lengths_script} {input} > {output}"