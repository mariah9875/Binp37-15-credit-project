#!/usr/bin/env python3

import argparse


usage = "This program extracts gene type, gene name and gene id from GTF annotation file"
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('-i',metavar='InFile',help='A GTF input file',type=argparse.FileType('r'),required=1)
parser.add_argument('-t',metavar='gene_type',help='A txt file with gene_type',type=argparse.FileType('w'),required=1)
parser.add_argument('-o',metavar='gene_id',help='A txt file with gene_id',type=argparse.FileType('w'),required=1)
parser.add_argument('-n',metavar='gene_name',help='A txt file with gene_name',type=argparse.FileType('w'),required=1)
args=parser.parse_args()

gene_idlist = []
gene_typelist = []
gene_namelist = []


with args.i as fin:
    for line in fin:
        if not line.startswith("#"):
            line = line.split()
            gene_id1 = line[8] # select "gene_id"
            gene_id2 = line[9] # select gene id name
            transcrpipt1 = line[10] # select "transcript"
            transcrpipt2 = line[11] # select transcript name
            gene_type1 = line[12] # select "gene_type"
            gene_type2 = line[13] # select gene_type name
            gene_name1 = line[14] # select "gene_name"
            gene_name2 = line[15] # select gene name
            if gene_id1 == "gene_id":
                gene_idlist.append(gene_id2) # append gene_id to gene_id list
            if transcrpipt1 == "gene_type":
                gene_typelist.append(transcrpipt2) # append gene_type to gene_type list
            if gene_type1 == "gene_type":
                gene_typelist.append(gene_type2)
            if gene_type1 == "gene_name":
                gene_namelist.append(gene_type2) # append gene_name to gene_name list
            if gene_name1 == "gene_name":
                gene_namelist.append(gene_name2)


with args.o as fout:
    print("#gene_id", file=fout)
    for line in gene_idlist:
        print(line, file=fout) # write to gene_id file

with args.n as fout:
    print("#gene_name", file=fout)
    for line in gene_namelist:
        print(line, file=fout) # write to gene_name file

with args.t as fout:
    print("#gene_type", file=fout)
    for line in gene_typelist:
        print(line, file=fout) # write to gene_type file
