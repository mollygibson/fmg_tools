#!/usr/bin/env python

import os
import pandas

def main(mapping_fp):
    mapping = pandas.io.parsers.read_table(mapping_fp)

    abx_dict = mapping.set_index('id')['antibiotic'].to_dict()
    library_dict = mapping.set_index('id')['library'].to_dict()

    map_object = ParseMapping(abx_dict, library_dict)

    return(map_object)


class ParseMapping:
    def __init__(self,id_to_antibiotic, id_to_id):
        self.abx = id_to_antibiotic
        self.id = id_to_id
