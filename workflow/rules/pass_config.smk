rule pass_config_file:
    input:
        f"config/config.json"
    output:
        out_config_file=f"{OUTDIR}/config.ini",
    benchmark:
        f"{BENCHDIR}/pass_config.benchmark"
    log:
        f"{LOGDIR}/pass_config.log",
    run:
        import configparser
        #config=input
        with open(output.out_config_file,'w') as configfile_out:
            config_to_pass = dict()
            config_to_pass["DATABASE"] = dict()
            config_to_pass["DATABASE"]["MEGARES"] = f"{DATABASES}/megares_database_v3.00.fasta"
            config_to_pass["DATABASE"]["MEGARES_ONTOLOGY"] = f"{DATABASES}/megares_annotations_v3.00.csv"
            config_to_pass["DATABASE"]["MGES"] = f"{DATABASES}/mges_combined.fasta"
            config_to_pass["MISC"] = dict()
            config_to_pass["MISC"]["GLOBAL_OVERLAP_THRESHOLD"] = 0.5 
            config_to_pass["MISC"]["GLOBAL_AMR_THRESHOLD"] = 0.8
            config_to_pass["MISC"]["GLOBAL_MGE_THRESHOLD"] = 0.5
            config_to_pass["MISC"]["DEDUPLICATION_SIMILARITY_THRESHOLD"] =0.9
            config_to_pass["MISC"]["OVERLAP_THRESHOLD_STRATEGY"] = "ONE"
            config_to_pass["MISC"]["MAX_BP_COLOCALIZATIONS_PLOT"] = 5000
            config_parser = configparser.ConfigParser()
            config_parser.read_dict(config_to_pass)
            config_parser.write(configfile_out)
