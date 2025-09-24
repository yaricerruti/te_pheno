#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Jun 13 11:35:49 2024

@author: yari
"""

# from os import listdir

infile = snakemake.input[0]
outfile = snakemake.output[0]
# infile = "/mnt/filippo/prj/te_enhancer_conservation/age_dataset/v2/rmsk_age_raw.bed"
# outfile = "/mnt/filippo/prj/te_enhancer_conservation/age_dataset/v2/rmsk_age_processed.bed"

chrom_dict = {}

## reading file
with open(infile, "r") as f:
    for line in f:
        ## key: (chr, start, end, classname) ## TE copy
        key = (line.split("\t")[0].strip(), line.split("\t")[1].strip(), line.split("\t")[2].strip(), line.split("\t")[3].strip())
        ## value: (fragment_age, fragment_length) ## computed age for the fragment, length of the overlap
        fragment_age = int(line.split("\t")[4].strip())
        fragment_length = int(line.split("\t")[5].strip())
        if key not in chrom_dict.keys():
            chrom_dict[key] = [(fragment_age, fragment_length)]
        else:
            chrom_dict[key].append((fragment_age, fragment_length))

## selecting fragment_age with longest fragment_length, for each key (TE copy)
for key, value in chrom_dict.items():
    ## list of fragment_age values
    age_list = [element[0] for element in value]
    ## list of fragment length values
    fragment_list = [element[1] for element in value]
    ## assigning to each unique key (each unique TE copy) the age associated with the longest fragment
    chrom_dict[key] = int(age_list[fragment_list.index(max(fragment_list))])


## output
with open(outfile, "a") as f:
    for key, value in chrom_dict.items():
        f.write("{chrom:}\t{start:}\t{end:}\t{classname:}\t{age:}\n".format(chrom=key[0],
                                                                            start=key[1],
                                                                            end=key[2],
                                                                            classname=key[3],
                                                                            age=value))
