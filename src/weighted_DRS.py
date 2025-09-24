#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# from os import listdir

import numpy as np
# from sys import getsizeof

infile = snakemake.input[0]
outfile = snakemake.output[0]
# infile = "/mnt/filippo/prj/te_enhancer_conservation/DRS_dataset/v1/DRS/PLS_real_te_func_ovl_DRS.bed"
# outfile = "/mnt/filippo/prj/te_enhancer_conservation/DRS_dataset/v1/DRS/PLS_real_te_func_ovl_weighted_DRS.bed"

DRS_dict = {}

## reading file
with open(infile, "r") as f:
    for line in f:
        ## key: (chr, start, end, classname) ## TE copy
        key = (line.split("\t")[0].strip(), line.split("\t")[1].strip(), line.split("\t")[2].strip(), line.split("\t")[3].strip())
        ## value: (DRS, ovl_length) ## Depletion Rank Score for the TE overlapping a specific window, length of the overlap between TE and window
        DRS = float(line.split("\t")[4].strip())
        ovl_length = int(line.split("\t")[5].strip())
        if key not in DRS_dict.keys():
            DRS_dict[key] = [(DRS, ovl_length)]
        else:
            DRS_dict[key].append((DRS, ovl_length))


## computing DRS weighted average for each key (TE copy)
for key, value in DRS_dict.items():
    ## list of Depletion Rank Scores
    DRS_list = [element[0] for element in value]
    ## list of overlap length values
    ovl_list = [element[1] for element in value]
    ## assigning to each unique key (each unique TE copy) the average DRS weighted for overlap length
    DRS_dict[key] = float(np.average(a=DRS_list, weights=ovl_list))


## output
with open(outfile, "a") as f:
    for key, value in DRS_dict.items():
        f.write("{chrom:}\t{start:}\t{end:}\t{classname:}\t{DRS:}\n".format(chrom=key[0],
                                                                            start=key[1],
                                                                            end=key[2],
                                                                            classname=key[3],
                                                                            DRS=value))

##