# Shalabh Shukla
# The purpose is to construct different molecular configurations of the MHC binding slot for alleles of interest using RDkit, returns Mol data as alternate x value instead of AAsq
# The three working models that will be tested include: 'Full_Binding_Slot','AASeq_2Helices','AASeq_2Helices_Cyclic', can include residue models / other models if needed
# Match y values to data using Allele names (allele name nomencalature if finicky maybe use AAseq match instead)
# 50 is the alignment threshold for alpha helix 1, and 80 is the alignment threshold for alpha helix 2

THRESHOLD_HELIX1 = 50
THRESHOLD_HELIX2 = 80
MIN_LENGTH_SLOT = 180

def mhc_slot_dataset(directory = '../data/mhc_slot/',
                    data_files = ['MHC_Full_alleles.csv', 'MHC_Full_AAseq_DB.csv'],
                    mhc_slot_model = 'AASeq_2Helices_Cyclic', 
                    output_files = None):

    def MissingFile(file_name):
        
        try:
            file = open(file_name); file.close() 
        
        except FileNotFoundError:
            import os; os.getcwd()
            print(file_name + ' does not exist or not in wokring directory')
            raise
    
    def MissingModule(Module_name, Import_name):
        
        try:
            Import_name = __import__(Import_name) 
        
        except ImportError:
            
            try:
                print(Module_name + ' not installed, pip installing Modules')
                import sys, subprocess; subprocess.call([sys.executable, "-m", "pip", "install", Module_name]); print(Module_name + ' installed\n')
           
            except:
                print('failed module import')
                raise
                
    MissingFile(data_files[0]) # csv with allles of interest down a column
    MissingFile(data_files[1]) # csv containg MHC allele name in column 1 and MHC Sequences in column 2 w/ delim " "            
    MissingModule('biopython', 'Bio') # Sequence alligner to parse sequence data
    MissingModule('pandas', 'pandas') # Dataframe containing mol data
    
    from Bio import pairwise2; from Bio.Seq import Seq; from Bio.SubsMat.MatrixInfo import blosum62 # matrix optimizes AA reads in search algorithm (AA matrix scores matches and misses)
    RefAH1Seq = Seq("PWIEQEGPEYWDRNTQIFKTNTQTYRESLRNLRG") # Refrence Sequnece helix 1 for B:35:01 isoform AASeq[73:107] 
    RefAH2Seq = Seq("AQITQRKWEAARVAEQLRAYLEGLCVEWLRRYLENGKET") # Refrence Sequnece helix 2 for B:35:01 isoform AASeq[163:202]
    PrevAASeq = " "; Not_MatchedA1 = 0; Not_MatchedA2 = 0; Not_Long = 0; Not_Expressed = 0; Not_Unique = 0; Accepted  = 0; # intiliazing AA sequence read and counters
    MHC_Helm = [[],[],[]]; MHC_Mol = [[],[],[]] # helm anotation of peptides (more robust than smiles anotation for peptide configurations), and Rdkit Mol data initialization
    
    def FullSeq2HelixSeq(Seq_full,Seq_helix,AASeq): # helix match function
        global AASeq_Helix; global Alignment_Helix; global Helix_LastPos
        Alignment_Helix = list(sum((pairwise2.align.localds(Seq_full, Seq_helix, blosum62, -10, -0.5)), ())) # helix allignemnt, pairwise needs 5 args returns 5 elements in list  
        Helix_Dif = (len(Seq_helix)-(Alignment_Helix[4] - Alignment_Helix[3])) # alignment starts and ends on matches only
        
        if Helix_Dif > 0: # allignment not contigous to refrence helix if diffrence greater than 0
            
            for HelixPos in range(len(Seq_helix)):
                
                if AASeq[Alignment_Helix[3]] == Seq_helix[HelixPos]:
                    AASeq_Helix = AASeq[Alignment_Helix[3] - HelixPos : Alignment_Helix[4] + Helix_Dif - HelixPos] # allignment adjustments for mismatch
                    Helix_LastPos = Alignment_Helix[4] + Helix_Dif - HelixPos
                    break
                
        else:
            AASeq_Helix = AASeq[Alignment_Helix[3]:Alignment_Helix[4]]
            Helix_LastPos = Alignment_Helix[4] # perfect allignment case / insertion case
            
        return AASeq_Helix, Alignment_Helix, Helix_LastPos
    
    import pandas as pd 
    MHC_slot_df = pd.DataFrame() # final output file containg alleles as rows and MOL configurations as columns
    
    import csv
    from rdkit import Chem
                                
    with open (data_files[0], 'r') as Alleles_csv_file: # .csv or .txt doesn't really matter as names stored in single column
        Alleles = csv.reader(Alleles_csv_file) # no need to specify delimiter
        
        for Allele in Alleles:
                
            with open (data_files[1], 'r') as Sequnces_csv_file: # enter sequences csv file here (allele column 1 and sequences column 2)
                Sequnces = csv.reader(Sequnces_csv_file, delimiter=',') # check csv file for delimiter type
                
                for Sequnce in Sequnces :
                    
                    if Allele[0] == Sequnce[0][:len(Allele[0])]: # MHC sequence search leneince shown towards allele entry i.e A*01:01 accepted as A*01:01:01:01
                        print(Allele[0]); # print Allele
                        break
                    
    
            if  Allele[0][-1] not in ["Q","N"] and PrevAASeq != Sequnce[1] and len(Sequnce[1]) >= MIN_LENGTH_SLOT: # Q and N denote dubious expression status (substoichiometric detection)
                PrevAASeq = Sequnce[1];# Variants need have different AA seq (exonic variants only), Sequences that are small transcript variants are not Accepted (can't contain full binding slot)   
                AASeq_whole = Seq(Sequnce[1]) # MHC Sequence of interest
                
                Alignment_AH1 = FullSeq2HelixSeq(AASeq_whole,RefAH1Seq,Sequnce[1]) # search for alpha helix 1 alignment
                Alignment_AH2 = FullSeq2HelixSeq(AASeq_whole,RefAH2Seq,Sequnce[1]) # search for alpha helix 2 alignment
                
                if Alignment_AH1[1][2] >= THRESHOLD_HELIX1 and Alignment_AH2[1][2] >= THRESHOLD_HELIX2: # Arbitrary alignment thresholds for helix 1 and helix 2 given above 
                    print("Accepted \n"); Accepted  =  Accepted  + 1 # print acceptable AASeq
                    
                    AASeq_AH1 = Alignment_AH1[0] # helix 1 AAseq
                    AASeq_AH2 = Alignment_AH2[0] # helix 2 AAseq
                    AASeq_AH1_AH2_wG = AASeq_AH1 + ("G" * 4) + AASeq_AH2 # alpha helices linked w G linker, G linker length is 4
                    AASeq_AH1_AH2_wG_cyc = AASeq_AH1 + ("G" * 4) + AASeq_AH2 + ("G" * 4) # alpha helices linked w 2XG linker for cyclic peptide, G linker length is 4
                    AASeq_BS_AH1_AH2 = Sequnce[1][:Alignment_AH2[2]] # start of beta sheet to end of helix 2, full Sequence of binding slot
                    MHC_Seq = [AASeq_BS_AH1_AH2,AASeq_AH1_AH2_wG,AASeq_AH1_AH2_wG_cyc] # string Sequences of peptides
                    
                    for SeqTypes in range (len(MHC_Seq)-1):
                        MHC_Helm[SeqTypes] = "PEPTIDE1{" + ('.'.join(MHC_Seq[SeqTypes])) + "}$$$$V2.0" # assigning HELM anotation (more robust for describing peptide structures than smiles)
                    
                    MHC_Helm[-1] = "PEPTIDE1{" + '.'.join(AASeq_AH1_AH2_wG_cyc) + "}$PEPTIDE1,PEPTIDE1,1:R1-" + str(len(AASeq_AH1_AH2_wG_cyc)) + ":R2$$$V2.0" # cyclic HELM
                                                 
                    for MolTypes in range (len(MHC_Helm)):
                        MHC_Mol[MolTypes] = Chem.MolToSmiles(Chem.MolFromHELM(MHC_Helm[MolTypes])) # Mol data for helices and binding slot models
                    
                    #  Dataframe containing Allele names as index and MOL data (in smiles annotation) for different models of helices and binding slots
                    MHC_slot_df = MHC_slot_df.append(pd.DataFrame.from_dict(dict([(Allele[0],[MHC_Mol[0],MHC_Mol[1],MHC_Mol[2]])]), 
                                       orient= 'index',columns=['Full_Binding_Slot','AASeq_2Helices','AASeq_2Helices_Cyclic']))
                
                else: # rejected alleles that didn't allign with refernce helices for secondary quality control
                    Rejected = "Rejected | "
                    
                    if Alignment_AH1[1][2] < THRESHOLD_HELIX1:
                        Rejected = Rejected + "Failed helix 1 allignment match | "; Not_MatchedA1 = Not_MatchedA1 +1 
                        
                    if Alignment_AH2[1][2] < THRESHOLD_HELIX1:
                        Rejected = Rejected + "Failed helix 2 allignment match | "; Not_MatchedA2 = Not_MatchedA2 +1 
                        
                    print(Rejected + "\n")
                        
            else: # rejected alleles for various other parameters for preliminary quality control
                Rejected = "Rejected | "
                    
                if PrevAASeq == Sequnce[1]:
                    Rejected = Rejected + "Non unique exonic variant | "; Not_Unique = Not_Unique + 1
                
                if Allele[0][-1] in ["Q", "N"]:
                    Rejected = Rejected + "Dubious expression status | "; Not_Expressed = Not_Expressed + 1
                
                if len(Sequnce[1]) < MIN_LENGTH_SLOT:
                    Rejected = Rejected + "Sequence too short | "; Not_Long = Not_Long + 1
                    
                print(Rejected + "\n")
    
    print(MHC_slot_df)    
    print("\n" + "Alleles accepted = " + str(Accepted) + "\n") # quality control report
    print("Alleles not unique exonic variants = " + str(Not_Unique))
    print("Alleles with dubious expression status (N/Q status) = " + str(Not_Expressed))
    print("Alleles too short to contain MHC-1 binding slot = " + str(Not_Long) + "\n")
    print("Alleles don't allign to refrence MHC-1 helix 1 = " + str(Not_MatchedA1))
    print("Alleles don't allign to refrence MHC-1 helix 2 = " + str(Not_MatchedA2))
    
    if output_files:
        print('generating output files')
        pd.to_pickle(MHC_slot_df, "./MHC_MOLslot_file_FullDB.pkl")
        MHC_slot_df.to_csv(MHC_slot_df)
    
    return (MHC_slot_df[mhc_slot_model])