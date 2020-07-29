# mhc_slot_dataset

mhc_slot_dataset is a data parsing tool used to extract useful binding information from full protein amino acid sequences given as an input from a large database of protein alleles. The purpose is to construct representations of binding slots for alleles of binding proteins, returning x values as either strings (the amino acid sequence) or smiles (molecular configuration annototation). Currently there are four models to chose from: 'full_complex', 'binding_slot', 'linear', or 'cyclic' which captures varying degrees of accuracy as dummy representations of the protein binding slot. 

Below are some examples of mhc_slot_dataset usage:

```
x = mhc_slot_dataset(allele = ['HLA-A*02:01:01'], binding_model = 'cyclic', return_type = 'strings')

x = mhc_slot_dataset(allele = ['HLA-DRA*01:01:01:01','HLA-DRB1*01:01:01:01'], binding_model = 'full_complex', return_type = 'smiles') 

x = mhc_slot_dataset(allele = ['QJE37815_surface_glycoprotein'], binding_model = 'binding_slot', return_type = 'smiles') 
```
## Code flowchart

![Image of flowchart](https://github.com/sshuklz/GEMSEC-CodeRepository/blob/master/mhc_slot_dataset_flowchart.key.png)

## Hyper parameters

There are currently three hyper parameters which can be fine tuned to optimize search for specific proteins of interest:

1. BINDING_DOMAIN_REF - The given refrence sequences that are known to participate in binding (sourced from literature, crystal stuctures, known domains and more). Have to be changed for each new protein investigated. Furure iterations will remove the need for this hyperparameter by instead using a database as refrence. Below is an example of refrence sequences for MHC-I class proteins:

```
BINDING_DOMAIN_REF = [
    "PWIEQEGPEYWDRNTQIFKTNTQTYRESLRNLRG",      # Refrence Sequnece helix 1 for B:35:01 isoform AASeq[73:107] 
    "AQITQRKWEAARVAEQLRAYLEGLCVEWLRRYLENGKET"  # Refrence Sequnece helix 2 for B:35:01 isoform AASeq[163:202]
    ]
```

2. THRESHOLD_COVERAGE - A value between 0 and 1 which denotes percent sequence matched with refrence sequence using blosum scores obtained from the pairwise2 funtion of the BioPython package. This serves as a quality control measure to ensure sequences still contain binding features and are not too disimilar to refrence sequences. A low THRESHOLD_COVERAGE such as 0.3 should be reserved for highly variable binding domains such as those seen in MHC proteins, however a high THRESHOLD_COVERAGE such as 0.9 would be recomended for sequences with higly conserved domains such as ATPase. 
```
THRESHOLD_COVERAGE = 0.30 
```

Below is an example of a well conseverd sequence vs a varied sequence for refrence when deciding on a THRESHOLD_COVERAGE value: 

```
Highly conserved:                             Highly variable:
PWIEQEGPEYWDRNTQIFKTNTQTYRESLRNLRG            PWIEQEGPEYWDRNTQIFKTNTQTYRESLRNLRG
|||||||||||||.||||||||||||||||||||     VS     |||..|||.|.||.|.|||..|||..|...||||
PWIEQEGPEYWDRPTQIFKTNTQTYRESLRNLRG            PWILLEGPSYTDRPTSIFKSSTQTAAEGGSNLRG
```

3. G_LINKER_LENGTH - For certain binding models, linear and cyclic, a G-linker is used to attach disparate binding domains to each other. G-linkers are flexible in terms of length as all they consist of are series of G amino acids in succession. A single G amino acid is ~ 3.5 angstroms so use this to estimate what the appropriate distance between disprate binding domains should be for more effective results. An example is given below:

```
G_LINKER_LENGTH = 4 # Implies ... -G-G-G-G- ...
```
## User inputs in mhc_slot_dataset function fields

There are currently 5 user defined inputs that need to be specfied for mhc_slot_dataset to operate:

```
def mhc_slot_dataset(
    allele,
    binding_model,
    return_type,
    sequences_database,                              
    directory):
```

1. allele - Is the allele and corrosponding name of the allele that the user is interested in. Allele name must exist in allele dataset that is also specified by the user (see 4). If possible avoid using names with math operators although mhc_slot_dataset is able to accomadate such cases. Allele name must be given in square braces. For multicomplex cases and do the following: for homomultimeric simply repeat the allele name multiple times in accordance to the number of monomers that form the full binding pocket, for heteromultimeric complexes simply enter the serotype of the complex of interest.

![Image of complexes](https://raw.githubusercontent.com/sshuklz/GEMSEC-CodeRepository/master/complex_diagram.png)

```
allele = ['HLA-A*02:01:01'] # monomer
allele = ['P97VCP_ATPase','P97VCP_ATPase','P97VCP_ATPase'] # homomultimeric (trimer)
allele = ['HLA-DRA*01:01:01:01','HLA-DRB1*01:01:01:01'] # heteromultimeric (dimer)
```

2. binding_model - Various models that seek to represent the binding slot of binding proteins in varying degrees of resolution. The 'full_complex' is just that, all the compoenents of the binding complex. One order of complexity down is the 'binding_slot' which simply takes the tertiary compoents of the binding pocket, starting from the first alligend binding amino acid to the last containing everything inbetween. It may be the case that the quatenary / tertiary super structure is neccasay to orient the binding domain of the protein such that it complements the substrate / surface of interest. This would be the lock and key model where the slot is shaped spefically to accomadate its binding parter. However using the more lose and widely accepted model of the induced fit, all we need are the base amino acids of the binding domains (typically helices) and simply link disparate domains together in either a 'linear' or 'cyclic' fashion with the aid of G-linkers (specied in hyper parameters).  

![Image of models](https://github.com/sshuklz/GEMSEC-CodeRepository/blob/master/Binding_Model.png?raw=true)
```
binding_model = 'full_complex'
binding_model = 'binding_slot'
binding_model = 'linear'
binding_model = 'cyclic'
```

3. Return_type - Choice to either return the specified binding model as either a string or a smiles output. String ouputs are simply the amino acid sequences that form the respective binding models. The Smiles output however encodes the amino acid sequence into a molecular configuration that is widely understood and used. Below are some examples of what the output looks like. 
```
return_type = 'strings'
Loaded cyclic data for HLA-A*02:01:01 from HLA-A020101_strings.pickle 

PWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGGGGGAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETGGGG

return_type = 'smiles'
Loaded cyclic data for HLA-A*02:01:01 from HLA-A020101_smiles.pickle 

CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@@H]2CCCN2C(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCCN)NC(=O)CNC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CS)NC(=O)[C@H]([C@@H](C)O)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](CCCCN)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](CCCCN)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](C)NC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H]([C@@H](C)O)NC(=O)CNC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](CO)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CCC(=O)O)NC(=O)CNC(=O)[C@H](CC(=O)O)NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](CCC(=O)O)NC1=O
```

4. sequences_database - In order for mhc_slot_dataset to function you must provide it a dataset in the form of a csv file that contains the following two columns: allele name in the first column, and allele amino acid sequence in the second. These datasets can either be sourced from the literature, online, or can be your own personalized / customized mutant library.
```
sequences_database = 'MHC_Full_AAseq_DB.csv'
```

5. directory - is self explanatory, it is simply the location of sequences_database csv file.
```
directory = '../data/alleles/'
```
