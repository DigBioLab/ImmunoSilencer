import pandas as pd
from functions import*

# allotype_path = "A0101_probability.txt"
# binder_path = "test.fasta"
# with open(binder_path,'r') as file:
#     s=''
#     for line in file:
#         s+=line
    
    
# binder = Binder(s,epitopes=[[2,"QLVESGGGL",["HLA-A0101"],0.1]])

# logo=MHCI_Logo(allotype_path)

pdb_path = "data/pdbs/NY2-H04.pdb"

print(pdb_to_fasta(pdb_path))