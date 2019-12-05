#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:06:49 2019

@author: Shalabh
"""

import numpy as np
import csv
    
def FASTA2CSV_HLAiso(FASTAname): # input FASTA file name string
    
    ## DeNovo refrence table construction form FASTA file download
    
    Allele = np.array([]) # Allele and AAseq are array elements which contain those respective information
    AAseq = np.array([])
    CombStr = '' # fasta files split AAseq data by 60AA a line. Use .join to combine split strings
    n = 0 # to initiate first loop iteration, can use row number instead but this works
    
    with open (FASTAname, 'r') as csv_file:             # Loop taking forever, current time taken is 267 seconds for full (V2 w/ ISOonly is 118 seconds)
        csvread = csv.reader(csv_file, delimiter=' ')
        
        for row in csvread:
            #row[0][0] == '>' denotes row as header row containing header information such as allele name in file
            if n >= 1 and row[0][0] == '>':
                
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
                CombStr = "".join([CombStr, "".join(row)])  #Iterative joining of AA strings as long as it is not a header row in fasta file
           
            elif n == 0:                                    #first cylce allele name is entered, with no prior CombStr entry
                Allele = np.append(Allele,row[1])
                print(row[1])                               #print function slows loop, remove all in final implemenation (10 second improvement)
                n = n + 1
                
print(CombStr)
AAseq = np.append(AAseq,CombStr) # last CombStr entry as loop has finished
CombStr = ''
print(CombStr)

AAseq_ISO, AAseq_indices = np.unique(AAseq, return_index=True)
AAseq_indices = list(np.sort(AAseq_indices))

AAseq_ISO = AAseq[AAseq_indices]
Allele_ISO = Allele[AAseq_indices]

RefTable = np.c_[Allele,AAseq] # everything can be stored directly into RefTable in loop, so no need for Allele,AAseq in final implementaion
RefTable_ISO = np.c_[Allele_ISO,AAseq_ISO]

CSVname = input("Save as .csv filename (i.e hla_prot_DB2): ")
CSVname_ISO = "".join([CSVname, "_ISOonly.csv"])
CSVname = "".join([CSVname, ".csv"])
np.savetxt(CSVname, RefTable, fmt='%s')
np.savetxt(CSVname_ISO, RefTable_ISO, fmt='%s')
