# pdb is read inn
# Chain A is isolated
# Fasta is printed from chain A
pdb_path = input()

def split_binder_and_target(file_path, binder_chains=['A'], out_path=None, out_dir=''):
    """
    Takes a pdb and writes two new files one with binder chains and the other with non-binder chains.
    Returns file paths to new files (2), and also the target chain names as a list.
    """
    
    if type(binder_chains) == str:
        binder_chains = [binder_chains]
    
    if out_path is None: #create the output paths
        path_parts = file_path.split('.')
        binder_chain_path = f'{out_dir}{path_parts[0]}_binder.{path_parts[1]}'
        target_chain_path = f'{out_dir}{path_parts[0]}_target.{path_parts[1]}'

    binder_lines = []
    target_lines = {}
    
    with open(file_path,'r') as file:

        for line in file:
            if line[0:4] == 'ATOM':
                row = [x  for x in line.split(' ') if (x!='')]

                if row[4] in binder_chains:
                    binder_lines.append(line) 
                    
                elif row[4] in target_lines:
                    target_lines[row[4]].append(line)
                
                else:
                    target_lines[row[4]]=[]
                    
    
    with open(binder_chain_path,'w') as file:
        for line in binder_lines:
            file.write(line)
    
    with open(target_chain_path,'w') as file:
        for chain in target_lines:
            for line in target_lines[chain]:
                file.write(line)
 
    return binder_chain_path, target_chain_path
three_to_one = {
    "ALA": "A",    "ARG": "R",    "ASN": "N",    "ASP": "D",    "CYS": "C",    "GLN": "Q",    "GLU": "E",    "GLY": "G",
    "HIS": "H",    "ILE": "I",    "LEU": "L",    "LYS": "K",    "MET": "M",    "PHE": "F",    "PRO": "P",    "SER": "S",
    "THR": "T",    "TRP": "W",    "TYR": "Y",    "VAL": "V",}

def pdb_to_fasta(pdb_path:str) -> str:
    out_str=''
    binder_id = pdb_path.split('/')[-1].replace('.pdb','')
    residues = ''
    last_added_pos = 0
    current_chain = ''
    
    with open(pdb_path, 'r') as file:
        for line in file:
            if line[0:4]!='ATOM':
                continue
            
            line = (line.split(' '))
            line = [entry for entry in line if entry != '']
            
            if line[5] == str(last_added_pos):
                continue
            
            if line[4] != current_chain:
                
                if current_chain != '': #just to not add a line hsift as the first char
                    out_str+='\n'
                
                out_str+=f'>{binder_id}_chain_{line[4]}\n'
                current_chain = line[4]
                
            last_added_pos = line[5]
            out_str+=(three_to_one[line[3]])


    return(out_str)

binder_name = pdb_path.split('.')[0]
binder_path, target_path = split_binder_and_target(pdb_path, out_dir=binder_name+'/') #Writes new pdb file to disc and returns filename

with open (f'{binder_name}/{binder_name}_binder.fasta','w') as file:
    file.write(pdb_to_fasta(binder_path))

with open (f'{binder_name}/{binder_name}_target.fasta','w') as file:
    file.write(pdb_to_fasta(target_path))
