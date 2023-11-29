#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 10:01:42 2020

@author: u0125634
"""
import re
import os

#Define working directory
path = "/Users/u0125634/Documents/PhD/Projects/PTM_project_final/Data/Human/Alignment_glyco/Alignments/"
os.chdir(path)

f=open("/Users/u0125634/Documents/PhD/Projects/PTM_project_final/Data/Human/Alignment_glyco/knownCanonical.exonAA.fa", "r")

#Read proteins
with open("/Users/u0125634/Documents/PhD/Projects/PTM_project_final/Data/Human/Alignment_glyco/ensembl_prots.txt", 'r') as h:
    glyc_prot = h.read().splitlines()

#Read file to dictionary
dict_prots = {}
for line in f:
    line = line.rstrip()
    if line.startswith(">"):
        matchHeader = re.search(">(.+?)\..+?_(.+?)_.+", line)
        if matchHeader:
                ensembl = matchHeader.group(1)
                organism = matchHeader.group(2)
                if ensembl not in dict_prots:
                    dict_prots[ensembl] = {}
                if organism not in dict_prots[ensembl]:
                    dict_prots[ensembl][organism] = []
    elif line:
        sequence = line
        dict_prots[ensembl][organism].append(sequence)
        
    else:
        continue
    
#Format the dictionary and write to folder
for e in dict_prots.keys():
    if e in glyc_prot:
        if not os.path.exists(e):
            os.makedirs(e)
        a_file = path+e+'/alignment.txt'
        g = open(a_file,'w')
        for o in dict_prots[e].keys():
            separator = '' 
            seq = separator.join(dict_prots[e][o])
            g.write(o+'\t'+seq+'\n')
        g.close()
    else:
        continue
        