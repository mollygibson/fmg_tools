#!/usr/bin/env python

import os

def main(mapping_fp):
    mapping_open = open(mapping_fp, "U")

    abx_hash = {}
    for item in mapping_open:
        pieces = item.split("\t")
        abx_hash[pieces[0]] = pieces[1].rstrip()

    abx_object = ParseMapping(abx_hash)

    return(abx_object)


class ParseMapping:
    def __init__(self,id_to_antibiotic):
        self.abx = id_to_antibiotic
