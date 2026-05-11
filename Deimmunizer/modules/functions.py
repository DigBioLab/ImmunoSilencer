import pandas as pd
import numpy as np
import random as rd
from Bio.PDB import PDBParser, NeighborSearch

#Load default settings
settings_dict={}
with open("Deimmunizer/input/default_settings.txt") as file:
    for line in file:
        if line[0]=='#' or line.strip()=='':
            continue
        row = line.split(',')
        settings_dict[row[0]] = float(row[1])

#overwrite default with user defined settings
with open("Deimmunizer/input/user_settings.txt") as file:
    for line in file:
        if line[0]=='#' or line.strip()=='':
            continue
        row = line.split(',')
        settings_dict[row[0]] = float(row[1])


three_to_one = {
    "ALA": "A",    "ARG": "R",    "ASN": "N",    "ASP": "D",    "CYS": "C",    "GLN": "Q",    "GLU": "E",    "GLY": "G",
    "HIS": "H",    "ILE": "I",    "LEU": "L",    "LYS": "K",    "MET": "M",    "PHE": "F",    "PRO": "P",    "SER": "S",
    "THR": "T",    "TRP": "W",    "TYR": "Y",    "VAL": "V",}

#One time run to create the database of theoretical human peptides
# def fastaread(filename:str):

#     headers = []
#     seqs = []
#     with open(filename, 'r') as file:
#         for line in file: # reading file line by line
            
#             if line.replace(' ','')[0] == '>': #Headers start with a '>'. I remove spaces from the check for robustness
#                 headers.append(line)
#                 seqs.append('') #A new header means a new sequence incomming. If not, an empty string for a sequence is also reasonable behaviour                              
                    
#             else:
#                 seqs[len(headers)-1] += line.replace(' ','').replace('\n','') #If the line is not a header it contains a sequence. I remove lineshifts and spaces.
#     return headers, seqs

# ref_path = "data/Human_Ref.fasta"
# headers, seqs = fastaread(ref_path)
# human_peptides = []
# for seq in seqs:

#     for i in range(len(seq)-8):
#         peptide = seq[i:i+9]
#         human_peptides.append(peptide)
    
# human_peptides=sorted(set(human_peptides))
        
# with open("data/Human_Ref.peptide",'w') as file:
#     for pep in human_peptides:        
#         file.write(pep+'\n')


def usage(msg:str):
    print(msg)
    return(1)



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

class MHCI_Logo:
    def __init__(self, matrix_path):
        self.name = matrix_path[-25:-16]
        self.logo = pd.read_csv(matrix_path,sep='\t')
        
    def get_mutations(self,aa,pos):
        vector=self.logo.iloc[int(pos)].drop([aa,"Position"]).copy()      
        return(vector)
    
    def get_mutation_weights(self,aa,pos,normalize=True,temp=settings_dict["residue_temperature_logo"]):
        vector = self.get_mutations(aa,pos)
        vector += 0.000001 #safety messure for the cases where x=0
        weights = ((vector)**(-1)) #This could also be some -log function or to another negative number
        weights = weights**(1/temp)
        if normalize:
            return(weights/sum(weights))
        return(weights)
        
    def get_allignment(self,core):
        allignment = pd.Series([row[core[i]] for i,row in self.logo.iterrows()])
        
        return(allignment)
    
    def get_mutation_spot(self,core,temp = settings_dict["position_temperature_logo"]): #temperature changes "how random" the selection should be. higher numbers mean more random. 0 is the non included minimum
        s=self.get_allignment(core)
        s=s**(1/temp)
        s=s/s.sum()
        pos=np.random.choice(len(s),p=s)
        return(pos,core[pos])
   
###-----define blossumSub metrics-------------####
with open("Deimmunizer/data/blosum62.txt",'r') as file:
    df=pd.read_csv(file,sep=r"\s+")
    
df.set_index("aa") 
aa = df.aa
mat= df.drop(columns="aa")

blossumsub_mutation_probs = {}
for i,X in enumerate(aa):

    series = mat.drop(columns=X).iloc[i]
    residue_list = series.index.tolist()
    scores = np.array(series.tolist())
    scores -= scores.min() #removes negative numbers. Important for the exponenent transformation to come
    scores_adj = (scores)**(1/settings_dict["residue_temperature_blossum"])
    p=(scores_adj-scores_adj.min())/(scores_adj-scores_adj.min()).sum()
    blossumsub_mutation_probs[X] = p,residue_list
###-----define blossumSub metrics-------------####

class Binder:
    def __init__(self, entry:str, n_mutations=0,epitopes=None):
        self.entry=entry.split('\n')
        self.name = self.entry[0][1:].strip()
        self.sequence = ''.join(self.entry[1:])
        self.n_mutations = n_mutations
        self.epitopes = [] if epitopes is None else epitopes
        self.pdb_path = self.name + '.pdb'
        self.bepi_result = None
        self.disco_result = None
        self.netmhci_result = None
        self.contact_positions = None
        self.in_silico_metrics = pd.DataFrame({"Pos":list(range(len(self.sequence))),
                                               "Residue": list(self.sequence),
                                               "Core": [(self.sequence+'--------')[i:i+9] for i in range(len(self.sequence))]})
        
    def __len__(self):
        return(len(self.sequence))
       
    def __str__(self):
        return(f'>{self.name} \n{self.sequence}')
    
    def __repr__(self):
        return(self.name)
    
    def n_epitopes(self):
        return(len(self.epitopes))

    
    def blossumsub(self, mut_pos, mutation_probs, report_to_path = None):
        
        res_from = self.sequence[mut_pos]
        res_to = np.random.choice(a=mutation_probs[res_from][1], p=mutation_probs[res_from][0])
        
        self.sequence = self.sequence[0:mut_pos] + res_to + self.sequence[mut_pos+1:]
        if report_to_path != None:
            
            print(f'BCR epitope detected: {res_from}{mut_pos} -> {res_to}{mut_pos}') 
                
        return(res_to)
            
    def logosub(self,allotype, mut_pos, report_to_path = None, logo_dir = 'Deimmunizer/data/MHCI_position_probabilities', bl_weight = settings_dict["blossum_weight_logo"]):
        logo = MHCI_Logo(f'{logo_dir}/{allotype[3:]}_probability.txt')
        global_pos = mut_pos
        core = self.in_silico_metrics.iloc[mut_pos].Core
        pos, aa_from = logo.get_mutation_spot(core)
        weights = logo.get_mutation_weights(aa_from, pos)
        #fetching and reformating the blossumsub metrics for this aminoacid
        weights_bl = pd.Series({aa: weight for aa,weight in zip(blossumsub_mutation_probs[aa_from][1],blossumsub_mutation_probs[aa_from][0])})
        weights_bl = weights_bl.reindex(weights.index)
        combined_weights = weights*(1-bl_weight) + weights_bl*bl_weight #applying weight
        combined_weights = combined_weights / combined_weights.sum()
            
        aa_to = np.random.choice(combined_weights.index, p=combined_weights.values)

        self.sequence = self.sequence[0:global_pos+pos] + aa_to + self.sequence[global_pos+pos+1:]
        
        if report_to_path != None:

            print(f'MHCI epitope detected: core at {global_pos} was {core}. {aa_from}{global_pos+pos}({pos})->{aa_to}') 
            
     
    def get_MHCI_epitopes_singles(self,MHCI_output_path):
        with open(MHCI_output_path ,'r') as file:
            epitopes=[]
            reading_data=0 #This needs to be equal to 2 because of mhcpans output
            for line in file:
                pickcc
                if line[0] == "#":
                    continue

                if line[0:5] == "-----":
                    reading_data+=1
                    
                    if reading_data == 4:
                        reading_data = 0 #This is also just to fit with the text file output

                
                if reading_data ==2: #Now we are in a data column
                    if "<= SB" in line:
                        line = line.split()
                        epitope = [line[0], line[2], [line[1]],line[-3]] #Position (1-indexed), sequence(Core), Allotype(s), score(s). 
                        self.epitopes.append(epitope)

                
        
    def get_MHCI_epitopes(self,df,allotpyes:list, peptide_lengths=[9,14], threshold=2.0):
        subset = df[df["ID"] == self.name]
        positive_subset = subset[subset["N_binders"] > 0]
        
        n = len(allotpyes)
        for i,row in positive_subset.iterrows():
            epitope=[row['Pos'],row["Peptide"],[],[]] #Position (1-indexed), sequence(Core), Alloty. 
            for j in range(n):
                
                epitope[2].append(allotpyes[j])
                epitope[3].append(row[f'Rank.{j}'])
            self.epitopes.append(epitope)  

    def report_epitopes(self,file_name,itteration):
        print(f'>{self.name}, {itteration}', file=file_name)
        for e in self.epitopes:
            print(f'{e[0]}\t{e[1]}\t{e[2]}\t{e[3]}', file=file_name)
        print(file=file_name)
        
    def load_BepiPred_results(self,result_path):
        self.bepi_result = pd.read_csv(result_path).drop(columns=["Residue","Accession","BepiPred-3.0 score"])
    
    def load_DiscoTope_results(self,result_path):
        self.disco_result = pd.read_csv(result_path).drop(columns=["pdb","chain","res_id","residue","rsa","length","alphafold_struc_flag"])
    
    def load_NetMHI_results_single(self,result_path):
        dfs=[]

        with open(result_path ,'r') as file:
            reading_data=0 #This needs to be equal to 2 because of mhcpans output
            
            for line in file:
                if line[0] == "#":
                    continue

                if line[0:5] == "-----":
                    reading_data+=1
                    
                    if reading_data == 3:
                        HLA= rows[-1][1] #This is a way to fetch the current HLA
                        df = pd.DataFrame(rows, columns=header)
                        df = df[["pos","Core","1-log50k(aff)","Affinity(nM)","%Rank"]]
                        df = df.rename(columns={name: f'{name}_{HLA}' for name in ["Core","1-log50k(aff)","Affinity(nM)","%Rank"]})
                        dfs.append(df)
                        print(df)
                    if reading_data == 4:
                        reading_data = 0 #This is also just to fit with the text file output
                    continue
                
                if reading_data == 1:
                    header=line.split()
                    header.pop()
                    rows=[]
                
                if reading_data == 2: #Now we are in a data column
                    line = line.split()
                    if "<=" in line:
                        line = line[:-2]
                    rows.append(line)
       
        out = dfs[0]
        for df in dfs[1:]:
            out = pd.merge(out, df, on="pos")
           
        out=out.filter(regex=r'^%')
        self.netmhci_result = out

    def merge_software_results(self):
        df = self.in_silico_metrics
        if self.bepi_result is not None:
            df = pd.concat([df,self.bepi_result.add_prefix("BP_")],axis=1)
        if self.netmhci_result is not None:
            df = pd.concat([df,self.netmhci_result.add_prefix("NM_")],axis=1)
        if self.disco_result is not None:
            df = pd.concat([df,self.disco_result.add_prefix("DT_")],axis=1)    
            
        return(df)
    
    def scan_humanness(self): #Returns a list of tuples containing epitopes [Position (1-indexed), sequence(Core), Allotype(s), score(s)] and humanness [min_distance, [peptides]]
        return [(epitope, get_humanness(epitope[1])) for epitope in epitopes]

    def get_epitopes_from_metrics(self,softwares:list, NM_threshold=2, BP_threshold=0.1512, DT_threshold=1.5, MHCI_allotypes=["HLA-A0201","HLA-A0301"],filter_epitopes=False): #
        #Allowed softwares [NM,BP,DT,IG] short for [NetMHCI-pan, BepiPred, DiscoTope, Immunogen]
        df = self.in_silico_metrics
        df["Epitopes"]=''
        if "NM" in softwares:
            for allotype in MHCI_allotypes:
                df[f'NM_%Rank_{allotype}'] = pd.to_numeric(df[f'NM_%Rank_{allotype}'])
                df.loc[df[f'NM_%Rank_{allotype}'] < NM_threshold, "Epitopes"] += f'NM_{allotype},'

                
        if "BP" in softwares:
            df["BP_BepiPred-3.0 linear epitope score"] = pd.to_numeric(df["BP_BepiPred-3.0 linear epitope score"])
            df.loc[df["BP_BepiPred-3.0 linear epitope score"] > BP_threshold, "Epitopes"] += "BP,"
            
        if "DT" in softwares:
            df["DT_calibrated_score"] = pd.to_numeric(df["DT_calibrated_score"])
            df.loc[df["DT_calibrated_score"] > DT_threshold, "Epitopes"] += "DT,"
            
            
        #2nd part is mesuring "humanness"
        df["Epitope_Humanness"]=np.nan
        df.loc[df["Epitopes"] != '', "Epitope_Humanness"] = df.loc[df["Epitopes"] != '', "Core"].apply(get_humanness)        
        
        # 3rd part load epitopes into the binder instance
        for row in df.iterrows():
           if row[1].Epitopes != '':
                epitope=[row[1].Pos, row[1].Core, row[1].Epitopes, get_humanness(row[1].Core)]
                self.epitopes.append(epitope)
                
        # 4th part: Mark contact residues
        if self.contact_positions != []:
            df["Contact Residue"] = df["Pos"].isin(self.contact_positions)
        
        self.in_silico_metrics=df
        
        if filter_epitopes:
            df = self.in_silico_metrics
            return df[df["Epitopes"] != '']
        
    def pick_mutation_spots_from_metrics(self,n_mutations = settings_dict["mutations_per_cycle"],BP_weight = settings_dict["BP_weight"],NM_weight = settings_dict["NM_weight"], DT_weight = settings_dict["DT_weight"], active_site_punish = -10):
        
        scoring_weights={"BP" : BP_weight,
                         "NM" : NM_weight,
                         "DT" : DT_weight}

        #Scan wich softwares are pressent
        allowed_softwares = ["BP","NM","DT"]
        softwares = [] 
        for colname in self.in_silico_metrics.columns:
            pre = colname[0:2]
            if pre in allowed_softwares:
                if pre not in softwares:
                    softwares.append(pre)

        epitope_positions= [] 
        mutation_weights = []

        for i,row in self.in_silico_metrics.iterrows():

            if type(row["Epitopes"]) is not str: 
                continue #skip non flagged residues
            
            epitope_positions.append(i)
            weight = 0 #weight determines likelihood to be picked for mutation
            
            for kind in row["Epitopes"].split(',')[:-1]:
                weight += scoring_weights[kind[0:2]] #
                
            weight += row["Contact Residue"]*active_site_punish
        
            mutation_weights.append(weight)


        
        if mutation_weights == []:
            return([])

        #softmax
        mutation_weights = np.array(mutation_weights)
        mutation_weights-=mutation_weights.max()
        mutation_weights = np.exp(mutation_weights) / np.exp(mutation_weights).sum()

        size = int(n_mutations)

        if len(mutation_weights) < size:
            size=len(mutation_weights)

        
        samples = np.random.choice(epitope_positions, size=size, replace=False, p=mutation_weights)
        return(samples)

    def mutate(self,positions:list, report_to_path=None,suffix=None): #sort priority
        rows = self.in_silico_metrics.iloc[positions]
        
        for i,row in rows.iterrows():

            kind = row.Epitopes.split(',')[0]
            if kind[0:3] == 'NM_':
                self.logosub(allotype = kind, mut_pos=row.Pos, report_to_path = report_to_path)     
                   
            elif kind == 'BP':
                self.blossumsub(mut_pos = row.Pos, mutation_probs=blossumsub_mutation_probs , report_to_path = report_to_path)
            
            elif kind == 'DT':
                self.blossumsub(mut_pos = row.Pos, mutation_probs=blossumsub_mutation_probs , report_to_path = report_to_path)
        
        if suffix != None:
            self.name += suffix
    
def load_NetMHCI_results_poly(file_name):
    df=pd.read_csv(file_name,sep ='\t',skiprows=1)
    df=df.dropna()
    df["ID"]=df["ID"].str[0:10].str.replace('_','|')
    
    with open(file_name, 'r') as file:
        allotpyes=[a for a in file.readline().split('\t') if 'HLA' in a]
    
    return(df.rename(columns={"Rank":"Rank.0"}),allotpyes)

def load_binders(fasta_path):
    with open(fasta_path, 'r') as file:
        fastas=[]
        i=-1
        for row in file:
            if row[0]==">":
                fastas.append(row)
                i+=1
            elif row=='\n':
                continue
            else:
                fastas[i]+=row
                
    binders=[Binder(fasta) for fasta in fastas]
    return(binders)
     
def load_epitope_report(file_name):
    epitopes={}
    
    with open(file_name, 'r') as file:
        for line in file:
            if line[0]==">":
                binder_name=line[1:11]
                epitopes[binder_name]=[]

            elif line == '\n':
                continue 
            else:
                line = line.replace('\n','')
                parts = line.split('\t') #[Pos,Seq,[Allotypes],[scores]]
                pos = int(parts[0])
                seq = parts[1]
                allotpyes = parts[2].strip("[]").split(',')
                scores = [float(x) for x in parts[3].strip("[]").split(',')]
                
                epitopes[binder_name].append([pos,seq,allotpyes,scores])
                
    return(epitopes)            
            
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

def Load_disco(file):
    df = pd.read_csv(file)
    return df

def get_peptide_dist(seq1,seq2):
    dist=0
    for a1,a2 in zip(seq1,seq2):
        if a1!=a2:
            dist+=1
    return(int(dist))

def get_humanness(query,human_peptides = '../../../data/Human_Ref.peptide',with_kandidates=False): # Returns [min_distance, [human peptides]]
    kandidates = []
    min_dist = 9
    
    with open(human_peptides, 'r') as file:
        for line in file:
            line=line[0:9]
            dist=get_peptide_dist(query,line)
            if dist < min_dist:
                kandidates=[line]
                min_dist=dist
            elif dist == min_dist:
                kandidates.append(line)
    if with_kandidates:
        return(min_dist,kandidates)
    return min_dist

def get_contact_residues(pdb_path, distance_threshold=4.0):
    contact_pos=[]
    
    parser=PDBParser(QUIET=True)
    structure = parser.get_structure("complex",pdb_path)
    
    binder_atoms = [atom for atom in structure[0]["A"].get_atoms()]
    target_atoms = [atom for atom in structure[0]["B"].get_atoms()]
    
    ns = NeighborSearch(target_atoms)
    
    for atom in binder_atoms:
        neighbors = ns.search(atom.coord, distance_threshold)
        if neighbors:
            residue = atom.get_parent()
            contact_pos.append(residue.get_id()[1])
    return(list(set(contact_pos)))

