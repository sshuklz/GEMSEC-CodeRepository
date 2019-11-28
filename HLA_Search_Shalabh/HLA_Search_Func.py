#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Allele --> AAseq function

@author: Shalabh
"""

import csv

def Allele2AAseq(Allele_entry):
    
    Allele_read = Allele_entry[4:] # removes HLA- component which is redundant 
    Allele_read1 = Allele_read
    Allele_read2 = Allele_read
    
    ColNum = Allele_read.count(':') # distinguishes between short and long form by : seperators present
            
    if ColNum <= 2:
        if ColNum == 1:
            Allele_read1 = "".join([Allele_read, ":01"]) # assumes short version is short hand for longer i.e HLA-A*01:01 is HLA-A*01:01:01
            Allele_read2 = "".join([Allele_read, ":01:01"]) # assumes short version is short hand for longer i.e HLA-A*01:01 is HLA-A*01:01:01:01
        elif ColNum == 2:
            Allele_read2 = "".join([Allele_read, ":01"]) # assumes short version is short hand for longer i.e HLA-A*01:01:01 is HLA-A*01:01:01:01
    
    with open ('hla_ISOprot_DB.csv', 'r') as csv_file:
        csvread = csv.reader(csv_file, delimiter=' ')
        
        for row in csvread:
            
            if Allele_read == row[0] or Allele_read1 == row[0] or Allele_read2 == row[0]: #matched entered allele with stored allele, made arguments more lenient 
            
                print(Allele_entry) # print allele
                print(row[1]) # print AAseq
