#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 10:51:17 2020

@author: u0125634
"""
import re

f=open("/Users/u0125634/Documents/PhD/Project/PTM_project/Data/Human/Alignment_glyco/knownCanonical.exonAA.fa", "r")
headers = open("/Users/u0125634/Documents/PhD/PTM_project/Data/Human/Alignment_glyco/headers.txt", "w")

all_titles = []
for line in f:
    line = line.strip('\n')
    if line.startswith(">"):
        title = line.rstrip()
        matchHeader = re.search(">(.+?)\..+", title)
        if matchHeader:
                title = matchHeader.group(1)
                all_titles.append(title)
                headers.write(title+"\n")