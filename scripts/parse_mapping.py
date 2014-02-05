#!/usr/bin/env python

import os

def main(mapping_fp):
    mapping_open = open(mapping_fp, "U")

    abx_hash = {}
    id_hash = {}
    for item in mapping_open:
        pieces = item.split("\t")
        abx_hash[pieces[0]] = pieces[2].rstrip()
        id_hash[pieces[0]] = pieces[1].rstrip()

    map_object = ParseMapping(abx_hash, id_hash)

    return(map_object)


class ParseMapping:
    def __init__(self,id_to_antibiotic, id_to_id):
        self.abx = id_to_antibiotic
        self.id = id_to_id
