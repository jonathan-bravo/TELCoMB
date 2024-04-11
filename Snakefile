# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

################################################################################
#                                Imports                                       #
################################################################################

import os, ast, glob

################################################################################
#                       Config File & Shorthands                               #
################################################################################

configfile: "config/config.json"
workdir: config["WORKFLOW"]["WORKDIR"]

databases_dir = workflow.basedir + "/" + config["WORKFLOW"]["DATABASES_DIR"]
samples_dir = workflow.basedir + "/" + config["WORKFLOW"]["SAMPLEDIR"] + "/samples"
tmp_dir = "tmp"
log_dir = "logs"

################################################################################
#                                 All Rule                                     #
################################################################################

if config["WORKFLOW"]["LONG_READS"].upper() == "TRUE":
    SAMPLES, = glob_wildcards(samples_dir + "/{sample_name}.fastq")
    if config["WORKFLOW"]["DEDUPLICATE"].upper() == "TRUE":
        DEDUP_STRING = config["EXTENSION"]["DEDUPLICATED"]
    else:
        DEDUP_STRING = config["EXTENSION"]["NOT_DEDUPLICATED"]
else:
    SAMPLES, = glob_wildcards(samples_dir + "/{sample_name}_R1.fastq") # CHANGE THIS TO FASTA, ALSO CONSIDER READ PAIRS
    DEDUP_STRING = config["EXTENSION"]["ASSEMBLED"]

EXTS = [
    DEDUP_STRING,
    config["EXTENSION"]["READS_LENGTH"],
    DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"],
    DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS_RICHNESS"],
    DEDUP_STRING + config['EXTENSION']['OVERLAP'],
    DEDUP_STRING + config["EXTENSION"]["RESISTOME_RICHNESS"],
    DEDUP_STRING + config["EXTENSION"]["RESISTOME_DIVERSITY"],
    DEDUP_STRING + config["EXTENSION"]["MOBILOME"],
    "_deduplicated_read_lengts_hist.pdf",
    "_colocalizations_plot.svg"]

rule all:
    input:
        expand("{sample_name}{ext}", sample_name=SAMPLES, ext=EXTS),
        "violin_plot_all_samples.svg",
        "heatmap_all_samples.svg"

if config["WORKFLOW"]["LONG_READS"].upper() == "TRUE":
    include: workflow.basedir + "/workflow/subworkflows/long_read.smk"
else:
    include: workflow.basedir + "/workflow/subworkflows/short_read.smk"

################################################################################
#                                Alignment                                     #
################################################################################

rule align_to_megares:
    input:
        reads = "{sample_name}" + DEDUP_STRING,
        megares_seqs = ancient(databases_dir + "/" + "megares_modified_database_v2.00.fasta")
    output:
        "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MEGARES"]
    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_OPTIONS"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["ALIGNMENT"]
    threads:
        config["MINIMAP2"]["THREADS"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.ato_megares.benchmark" 
    shell:
        "minimap2 -Y "
        "--secondary=no "
        "-t {threads} "
        "{params.minimap_flags} "
        "{input.megares_seqs} "
        "{input.reads} "
        "-o {output}"

rule align_to_mges:
    input:
        #reads = READS,
        reads = "{sample_name}" + DEDUP_STRING,
        mges_database = ancient(databases_dir + "/" + "mges_combined.fasta")
    output:
        "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MGES"]
    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_OPTIONS"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["ALIGNMENT"]
    threads:
        config["MINIMAP2"]["THREADS"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.ato_mges.benchmark" 
    shell:
        "minimap2 -Y "
        "--secondary=no "
        "-t {threads} "
        "{params.minimap_flags} "
        "{input.mges_database} "
        "{input.reads} "
        "-o {output}"

rule pass_config_file:
    output:
        out_config_file = "config.ini"
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".config.benchmark" 
    # conda:
    #     workflow.basedir + "/" + config["CONDA"]["PIPELINE"]
    run:
        import configparser
        with open(output.out_config_file,'w') as configfile_out:
            config_to_pass = dict(config)
            config_to_pass["DATABASE"] = dict()
            config_to_pass["DATABASE"]["MEGARES"] = databases_dir + "/" + "megares_modified_database_v2.00.fasta"
            config_to_pass["DATABASE"]["MEGARES_ONTOLOGY"] = databases_dir + "/" + "megares_modified_annotations_v2.00.csv"
            config_to_pass["DATABASE"]["MGES"] = databases_dir + "/" + "mges_combined.fasta"
            config_parser = configparser.ConfigParser()
            config_parser.read_dict(config_to_pass)
            config_parser.write(configfile_out)

rule overlap:
    input:
        megares_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MEGARES"],
        mges_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MGES"],
        reads_length = "{sample_name}" + config["EXTENSION"]["READS_LENGTH"],
        #dedup_reads_length = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"], # I have no idea where these come from
        config_file = "config.ini"
    params:
        overlap_script = workflow.basedir + "/" + config["SCRIPTS"]["FIND_OVERLAP"],
        output_prefix = "{sample_name}" + DEDUP_STRING
    conda:
        workflow.basedir + "/" + config["CONDA"]["PIPELINE"]
    output:
        overlap = "{sample_name}" + DEDUP_STRING + config['EXTENSION']['OVERLAP']
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.overlap.benchmark" 
    shell:
        "python {params.overlap_script} "
        "-r {wildcards.sample_name}.fastq "
        "-a {input.megares_sam} "
        "-m {input.mges_sam} "
        "-c {input.config_file} "
        "-o {params.output_prefix}"

rule merge_overlap_info:
    input:
        expand("{sample_name}" + DEDUP_STRING + config['EXTENSION']['OVERLAP'], sample_name = SAMPLES)
    output: 
        merged_info = 'merged_overlaped_mges_info.csv'
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".merge_overlap.benchmark" 
    # conda:
    #     workflow.basedir + "/" + config["CONDA"]["PIPELINE"]
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

rule resistome_and_mobilome:
    input:
        megares_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MEGARES"],
        mges_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MGES"],
        reads_length = "{sample_name}" + config["EXTENSION"]["READS_LENGTH"],
        #dedup_reads_length = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"],
        overlap = "{sample_name}" + DEDUP_STRING + config['EXTENSION']['OVERLAP'],
        config_file = "config.ini"
    output:
        resistome_richness = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["RESISTOME_RICHNESS"],
        resistome_diversity = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["RESISTOME_DIVERSITY"],
        mobilome = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["MOBILOME"]
    params:
        resistome_mobilome_script = workflow.basedir + "/" + config["SCRIPTS"]["GEN_RESISTOME_AND_MOBILOME"],
        output_prefix = "{sample_name}" + DEDUP_STRING
    conda:
        workflow.basedir + "/" + config["CONDA"]["PIPELINE"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.resistome_and_mobilome.benchmark" 
    shell:
        "python {params.resistome_mobilome_script} "
        "-r {wildcards.sample_name}.fastq "
        "-a {input.megares_sam} "
        "-m {input.mges_sam} "
        "-s {input.overlap} "
        "-c {input.config_file} "
        "-o {params.output_prefix}"

rule find_colocalizations:
    input:
        #reads = READS,
        reads = "{sample_name}" + DEDUP_STRING,
        megares_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MEGARES"],
        mges_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MGES"],
        reads_length = "{sample_name}" + config["EXTENSION"]["READS_LENGTH"],
        #dedup_reads_lenght = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"],
        overlap = "{sample_name}" + DEDUP_STRING + config['EXTENSION']['OVERLAP'],
        config_file = "config.ini"
    output:
        colocalizations = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"],
        genes_list = temp("{sample_name}" + DEDUP_STRING + config["EXTENSION"]["GENES_LIST"])
    params:
        find_colocalizations_script = workflow.basedir + "/" + config["SCRIPTS"]["FIND_COLOCALIZATIONS"],
        output_directory = os.getcwd()
    conda:
        workflow.basedir + "/" + config["CONDA"]["PIPELINE"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.find_colocalizations.benchmark" 
    shell:
        "python {params.find_colocalizations_script} "
        "-r {input.reads} "
        "--arg {input.megares_sam} "
        "--mge {input.mges_sam} "
        "-s {input.overlap} "
        "-c {input.config_file} "
        "-o {params.output_directory} "
        "> {output.colocalizations}"

rule colocalization_richness:
    input:
        colocalizations = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"],
        config_file = "config.ini"
    output:
        "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS_RICHNESS"]
    params:
        find_colocalizations_script = workflow.basedir + "/" + config["SCRIPTS"]["COLOCALIZATIONS_RICHNESS"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["PIPELINE"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.colocalization_richness.benchmark" 
    shell:
        "python {params.find_colocalizations_script} "
        "-i {input.colocalizations} "
        "-c {input.config_file} "
        "> {output}"


############################################################
## Plots
############################################################

rule read_lengths_plot:
    input:
        "{sample_name}" + DEDUP_STRING
    output:
        "{sample_name}_deduplicated_read_lengts_hist.pdf"
    params:
        num_of_bins = 100,
        std_deviations = 4,
        read_lengths_script = workflow.basedir + "/" + config["SCRIPTS"]["PLOT_RL"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["PLOTS"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.rl_plot.benchmark" 
    shell:
        "python {params.read_lengths_script} "
        "-i {input} "
        "-s {params.std_deviations} "
        "-b {params.num_of_bins} "
        "-o {output} "
        "--title {wildcards.sample_name}_deduplicated"

rule violin_plots_notebook:
    input:
        megares_db = databases_dir + "/megares_modified_database_v2.00.fasta",
        megares_annotation = databases_dir + "/megares_modified_annotations_v2.00.csv",
        config_file = "config.ini",
        data = expand("{sample_name}{ext}",sample_name=SAMPLES,ext=EXTS)
    output:
        "violin_plot_all_samples.svg"
    params:
        samples_list = SAMPLES,
        dedup_string = DEDUP_STRING,
        script = workflow.basedir + "/workflow/scripts/violin_notebook.py"
    conda:
        workflow.basedir + "/" + config["CONDA"]["PLOTS"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".violoin_plot.benchmark" 
    shell:
        "python {params.script} "
        "--dedup_string {params.dedup_string} "
        "--config_file {input.config_file} "
        "--output_plot {output} "
        "--samples_list {params.samples_list}"
    # notebook:
    #     "workflow/notebooks/violin_notebook.py.ipynb"


rule heatmap_notebook:
    input:
        megares_db = databases_dir + "/megares_modified_database_v2.00.fasta",
        megares_annotation = databases_dir + "/megares_modified_annotations_v2.00.csv",
        config_file = "config.ini",
        data = expand("{sample_name}{ext}",sample_name=SAMPLES,ext=EXTS)
    output:
        "heatmap_all_samples.svg"
    params:
        samples_list = SAMPLES,
        dedup_string = DEDUP_STRING,
        script = workflow.basedir +  "/workflow/scripts/heatmap_notebook.py"
    conda:
        workflow.basedir + "/" + config["CONDA"]["PLOTS"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".heatmap.benchmark" 
    shell:
        "python {params.script} "
        "--dedup_string {params.dedup_string} "
        "--config_file {input.config_file} "
        "--output_plot {output} "
        "--samples_list {params.samples_list}"
    # notebook:
    #     "workflow/notebooks/heatmap_notebook.py.ipynb"

# Going to have to change code here as well
rule colocalization_visualizations_notebook:
    input:
        megares_db = databases_dir + "/megares_modified_database_v2.00.fasta",
        megares_annotation = databases_dir + "/megares_modified_annotations_v2.00.csv",
        mges_db = databases_dir + "/mges_combined.fasta",
        read_lengths = "{sample_name}" + config["EXTENSION"]["READS_LENGTH"],
        colocalizations = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"],
        config_file = "config.ini"
    output:
        "{sample_name}_colocalizations_plot.svg"
    params:
        script = workflow.basedir + "/workflow/scripts/colocalizations_notebook.py"
    conda:
        workflow.basedir + "/" + config["CONDA"]["PLOTS"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".{sample_name}.colocalization_plot.benchmark" 
    shell:
        "python {params.script} "
        "--config_file {input.config_file} "
        "--read_lengths {input.read_lengths} "
        "--colocalizations {input.colocalizations} "
        "--output_plot {output}"
    # notebook:
    #     "workflow/notebooks/colocalizations_notebook.py.ipynb"

############################################################
## Databases
############################################################

rule get_megares:
    output:
        megares_seqs = os.path.join(databases_dir,"megares_modified_database_v2.00.fasta"),
        megares_ontology = os.path.join(databases_dir,"megares_modified_annotations_v2.00.csv")
    params:
        megares_seqs = config["DOWNLOADS"]["MEGARES_SEQS"],
        megares_ann = config["DOWNLOADS"]["MEGARES_ANN"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["DB"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".get_megares.benchmark" 
    shell:
        "mkdir -p {databases_dir}; "
        "wget {params.megares_seqs} -O {output.megares_seqs}; "
        "wget {params.megares_ann} -O {output.megares_ontology}"

# There has to be a better way to get the plasmid finder db...

rule get_plasmid_finder_db:
    params:
        git_repo = config["DOWNLOADS"]["PLASMID_FINDER_GIT"],
        commit = config["DOWNLOADS"]["PLASMID_FINDER_COMMIT"]
    output:
        temp(databases_dir + "/plasmid_finder_db.fasta")
    conda:
        workflow.basedir + "/" + config["CONDA"]["DB"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".get_pf.benchmark"
    shell:
        "mkdir -p {tmp_dir}; "
        "mkdir -p {databases_dir}; "
        "cd {tmp_dir}; "
        "git clone {params.git_repo}; "
        "cd plasmidfinder_db; "
        "git checkout {params.commit}; "
        "cd ../..; "
        "cat tmp/plasmidfinder_db/*.fsa > {output}; "
        "rm -rf tmp/plasmidfinder_db"

rule get_aclame_db:
    output:
        temp(databases_dir + "/aclame_db.fasta")
    params:
        aclame = config["DOWNLOADS"]["ACLAME"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["DB"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".get_aclame.benchmark"
    shell:
        "mkdir -p {databases_dir}; "
        "gdown {params.aclame} --fuzzy -O {output}"

rule get_iceberg_db:
    output:
        temp(databases_dir + "/iceberg_db.fasta")
    params:
        iceberg = config["DOWNLOADS"]["ICEBERG"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["DB"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".get_iceberg.benchmark"
    shell:
        "mkdir -p {databases_dir}; "
        "gdown {params.iceberg} --fuzzy -O {output}"

rule get_MGEs_DBs:
    input:
        databases_dir + "/plasmid_finder_db.fasta",
        databases_dir + "/aclame_db.fasta",
        databases_dir + "/iceberg_db.fasta"
    output:
        databases_dir + "/mges_combined.fasta"
    conda:
        workflow.basedir + "/" + config["CONDA"]["DB"]
    benchmark:
        workflow.basedir + "/benchmarks/" + config["WORKFLOW"]["WORKDIR"] + ".get_mges.benchmark"
    shell:
        "cat {input} > {output}"

############################################################
## Cleans
############################################################

# Need to run a test to see if this will work WILL UPDATE
# Also need to change some things with the config file before deleting all "tmp" files

# onsuccess:
#     shell("rm -f *.csv *.sam *.json *.pdf *_deduplicated.fastq")
#     shell("rm -rf {tmp_dir}")
#     shell("rm -f config.ini")
#     shell("rm -rf {databases_dir}")

# rule clean:
#     shell:
#         """
#         rm -f *.csv *.sam *.json *.pdf *_deduplicated.fastq
#         rm -rf {tmp_dir}
#         rm -f config.ini
#         """

# rule clean_databases:
#     shell:
#         """
#         rm -rf {databases_dir} {tmp_dir}
#         """