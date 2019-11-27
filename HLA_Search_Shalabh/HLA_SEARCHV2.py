#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:06:49 2019

@author: Shalabh
"""

import numpy as np
import csv
import os

print ('Put raw .fasta files in current working directory')
print (os.getcwd())

Denovo = input("Enter (YES) for Denovo FASTA to CSV conversion or (NO) to search existing CSV file: ")
if Denovo == "YES":
    
    FASTAname = input("Enter .fasta filename from database (i.e hla_prot.fasta): ")
    
    ## DeNovo refrence table construction form FASTA file download
    
    Allele = np.array([]) # Allele and AAseq are array elements which contain those respective information
    AAseq = np.array([])
    CombStr = '' # fasta files split AAseq data by 60AA a line. Use .join to combine split strings
    n = 0 # to initiate first loop iteration, can use row number instead but this works
    
    HLA_ISOonly = input("CSV output for specific HLA isoforms only? (YES/NO): ") #refernce table of ISOrelevant alleles
    
    with open (FASTAname, 'r') as csv_file:             # Loop taking forever, current time taken is 267 seconds for full (V2 w/ ISOonly is 118 seconds)
        csvread = csv.reader(csv_file, delimiter=' ')
        
        for row in csvread:
            #row[0][0] == '>' denotes row as header row containing header information such as allele name in file
            if n >= 1 and row[0][0] == '>':
                
                if HLA_ISOonly == 'YES':
                    if row[1].count(':') <= 2 or row[1][-4:-1] == ':01' or row[1][-3:] == ':01': # conditions for ISOrelevent alleles
                        n = 1 # proceed as normal
                    else:
                        n = 2 # skip loop / AAseq collection for allele
                                        
                if n == 1:
                    print(CombStr)                      #CombStr needs at least one AA before it is entered in. On first run this condition is not satisfied
                    AAseq = np.append(AAseq,CombStr)    #store CombStr in AAseq corrosponding to header allele
                    CombStr = ''                        #clear CombStr for next header allele
                    print(CombStr)
                    Allele = np.append(Allele,row[1])   #store Header allele name in Allele
                    print(row[1])
                    
            elif n == 1:
                CombStr = "".join([CombStr, "".join(row)])          #Iterative joining of AA strings as long as it is not a header row in fasta file
           
            elif n == 0:                                            #first cylce allele name is entered, with no prior CombStr entry
                Allele = np.append(Allele,row[1])
                print(row[1])                                       #print function slows loop, remove all in final implemenation (10 second improvement)
                n = n + 1
                
    print(CombStr)
    AAseq = np.append(AAseq,CombStr) # last CombStr entry as loop has finished
    CombStr = ''
    print(CombStr)
    
    RefTable = np.c_[Allele,AAseq] #everything can be stored directly into RefTable in loop, so no need for Allele,AAseq in final implementaion
    
    CSVname = input("Save as .csv filename (i.e hla_prot_DB2.csv): ")
    np.savetxt(CSVname, RefTable, fmt='%s')
    
## Allele search against refrence table (made in session or previous entry)

## Allele search nomenclature (short form - HLA-id*0n:0n , id is string identifer for gene like 'A',...
## * is seperator, first 0n field is allele group, : is field seperator, second 0n field is specific HLA protein)
## Allele search nomenclature (long form - HLA-id*0n:0n:0n:0n , third 0n field shows it is silent varient, 
## fourth 0n field are changes in non coding region.

## In summary long nomenclature also contains genomic / expression information, where as short form is only interested in protein isoforms
## To convert between forms HLA-A*01:01:01:01 is equivalent to HLA-A*01:01

CSVname = input("Enter .csv filename from previous session (i.e hla_prot_DB2.csv): ")

Allele_entry = input("For allele search enter name of allele of interest i.e. HLA-id*0n:0n, or HLA-id*0n:0n:0n:0n's': ")
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

with open (CSVname, 'r') as csv_file:
    csvread = csv.reader(csv_file, delimiter=' ')
    
    for row in csvread:
        
        if Allele_read == row[0] or Allele_read1 == row[0] or Allele_read2 == row[0]: #matched entered allele with stored allele, made arguments more lenient 
        
            print('')
            print(Allele_entry) # print allele
            print(row[1]) # print AAseq
