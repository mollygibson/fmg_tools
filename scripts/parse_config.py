#!/usr/bin/env python

import os

def main():
    config_open = open(os.environ['FMG_CONFIG_FP'], "U")

    config_hash = {}
    for item in config_open:
        pieces = item.split("\t")
        config_hash[pieces[0]] = pieces[1].rstrip()

    config_object = ParseConfig(config_hash['mgm_model_fp'], config_hash['resfam_database_fp'],config_hash['resfam_only_database_fp'], config_hash['pfam_database_fp'],config_hash['tigrfam_database_fp'])

    return(config_object)


class ParseConfig:
    def __init__(self,mgm,resfam,resfam_only,pfam,tigrfam):
        self.mgm_model_fp = mgm
        self.resfam_database_fp = resfam
        self.resfam_only_database_fp = resfam_only
        self.pfam_database_fp = pfam
        self.tigrfam_database_fp = tigrfam
