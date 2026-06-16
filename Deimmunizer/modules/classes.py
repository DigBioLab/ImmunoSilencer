import pandas as pd
import numpy as np
import random as rd
import os
from Bio.PDB import PDBParser, NeighborSearch

def get_settings(run_name):
    settings_dict={}
    with open(f"{run_name}/input/default_settings.txt") as file:
        for line in file:
            if line[0]=='#' or line.strip()=='':
                continue
            row = line.split(',')
            settings_dict[row[0]] = float(row[1])

    #overwrite default with user defined settings
    with open(f"{run_name}/input/user_settings.txt") as file:
        for line in file:
            if line[0]=='#' or line.strip()=='':
                continue
            row = line.split(',')
            settings_dict[row[0]] = float(row[1])
    return settings_dict

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



class MHCI_Logo:
    def __init__(self, matrix_path):
        self.name = matrix_path[-25:-16]
        self.logo = pd.read_csv(matrix_path,sep='\t')
        
    def get_mutations(self,aa,pos):
        vector=self.logo.iloc[int(pos)].drop([aa,"Position"]).copy()      
        return(vector)
    
    def get_mutation_weights(self,aa,pos,normalize=True,temp=None):
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
    
    def get_mutation_spot(self,core,temp): #temperature changes "how random" the selection should be. higher numbers mean more random. 0 is the non included minimum
        s=self.get_allignment(core)
        s=s**(1/temp)
        s=s/s.sum()
        pos=np.random.choice(len(s),p=s)
        return(pos,core[pos])
   


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
        self.ig_results = None
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

    
    def blossumsub(self, mut_pos, mutation_probs, report_to_path = None, epitope_type="BCR"):
        
        res_from = self.sequence[mut_pos]
        res_to = np.random.choice(a=mutation_probs[res_from][1], p=mutation_probs[res_from][0])
        
        self.sequence = self.sequence[0:mut_pos] + res_to + self.sequence[mut_pos+1:]
        if report_to_path != None:
            
            print(f'{epitope_type} epitope detected. Mutation with blossum66: {res_from}{mut_pos} -> {res_to}{mut_pos}') 
                
        return(res_to)
            
    def logosub(self,allotype,
                mut_pos,
                report_to_path = None,
                logo_dir = 'data/MHCI_position_probabilities',
                position_temperature_logo = None,
                residue_temperature_logo = None,
                blossum_weight_logo = None,
                blossumsub_mutation_probs=None):
        
        logo = MHCI_Logo(f'{logo_dir}/{allotype[3:]}_probability.txt')
        global_pos = mut_pos
        core = self.in_silico_metrics.iloc[mut_pos].Core
        pos, aa_from = logo.get_mutation_spot(core, temp=position_temperature_logo)
        weights = logo.get_mutation_weights(aa_from, pos, temp=residue_temperature_logo)
        #fetching and reformating the blossumsub metrics for this aminoacid
        weights_bl = pd.Series({aa: weight for aa,weight in zip(blossumsub_mutation_probs[aa_from][1],blossumsub_mutation_probs[aa_from][0])})
        weights_bl = weights_bl.reindex(weights.index)
        combined_weights = weights*(1-blossum_weight_logo) + weights_bl*blossum_weight_logo #applying weight
        combined_weights = combined_weights / combined_weights.sum()
            
        aa_to = np.random.choice(combined_weights.index, p=combined_weights.values)

        self.sequence = self.sequence[0:global_pos+pos] + aa_to + self.sequence[global_pos+pos+1:]
        
        if report_to_path != None:

            print(f'MHCI epitope detected: core at {global_pos} was {core}. {aa_from}{global_pos+pos}({pos})->{aa_to}') 
            
        
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
    
    def load_ImmunoGen_results(self, result_path):
        df = pd.read_csv(result_path)
        df["core"] = df["core_seq"]
        self.ig_results = df[["core","core_pos","pIRS_rank"]]
    
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
                        HLA= rows[-1][1].replace("*","") #This is a way to fetch the current HLA
                        df = pd.DataFrame(rows, columns=header)
                        df = df[["Pos","Core","%Rank_EL"]]
                        df = df.rename(columns={name: f'{name}_{HLA}' for name in ["Core","%Rank_EL"]})
                        dfs.append(df)

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
            out = pd.merge(out, df, on="Pos")
           
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
        if self.ig_results is not None:
            df = pd.concat([df,self.ig_results.add_prefix("IG_")],axis=1)    
            
        return(df)
    
    def get_contact_residues2(self, pdb_path, distance_threshold=4.0):

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("complex", pdb_path)
        model = structure[0]


        target_atoms = [atom for chain in model if chain.id != "A" for atom in chain.get_atoms()]
        ns = NeighborSearch(target_atoms)

        # Build a mapping: residue_number -> min distance to any target atom
        min_dist: dict[int, float] = {}

        for residue in model["A"].get_residues():
            res_num = residue.get_id()[1] - 1
            res_min = float("inf")

            for atom in residue.get_atoms():
                # Search with a generous radius, then take the true minimum
                neighbors = ns.search(atom.coord, distance_threshold * 10, level="A")
                for target_atom in neighbors:
                    d = atom - target_atom          # Bio.PDB overloads __sub__ → Å distance
                    if d < res_min:
                        res_min = d

            min_dist[res_num] = res_min if res_min < float("inf") else float("nan")

        # Write results back into the dataframe
        
        df = self.in_silico_metrics
        df["Target distance"] = df["Pos"].map(min_dist)
        df["Contact residue"] = df["Target distance"] <= distance_threshold
        self.in_silico_metrics = df

    
    # def scan_humanness(self): #Returns a list of tuples containing epitopes [Position (1-indexed), sequence(Core), Allotype(s), score(s)] and humanness [min_distance, [peptides]]
    #     return [(epitope, get_humanness(epitope[1])) for epitope in epitopes]

    def get_epitopes_from_metrics(self,softwares:list, NM_threshold=2, BP_threshold=0.1512, maximum_bepi_length=None, DT_threshold=1.5, MHCI_allotypes=["HLA-A0201","HLA-A0301"],filter_epitopes=False,IG_threshold=83): #
        #Allowed softwares [NM,BP,DT] short for [NetMHCI-pan, BepiPred, DiscoTope, Immunogen]
        df = self.in_silico_metrics
        df["Epitopes"]=''
        if "NM" in softwares:
            for allotype in MHCI_allotypes:
                print(df)
                df[f'NM_%Rank_EL_{allotype}'] = pd.to_numeric(df[f'NM_%Rank_EL_{allotype}'])
                df.loc[df[f'NM_%Rank_EL_{allotype}'] < NM_threshold, "Epitopes"] += f'NM_{allotype},'
                
                
        if "BP" in softwares:
            df["BP_BepiPred-3.0 linear epitope score"] = pd.to_numeric(df["BP_BepiPred-3.0 linear epitope score"])
            current_window = []
            
            for i,row in df.iterrows():            
                if len(current_window)==maximum_bepi_length:
                    df.loc[int(current_window[0]), "Epitopes"] = f"BP_{';'.join([current_window[0], current_window[-1]])},"
                    current_window = []
                
                if row["BP_BepiPred-3.0 linear epitope score"]>BP_threshold or current_window != []:
                    current_window.append(str(row["Pos"]))

            # For safty. wrap any unclosed epitopes
            if current_window != []:
                df.loc[int(current_window[0]), "Epitopes"] = f"BP_{';'.join(current_window)},"
                current_window = []
            
        if "DT" in softwares:
            df["DT_calibrated_score"] = pd.to_numeric(df["DT_calibrated_score"])
            df.loc[df["DT_calibrated_score"] > DT_threshold, "Epitopes"] += "DT,"
            
        if "IG" in softwares:
            df["IG_pIRS_rank"] = pd.to_numeric(df["IG_pIRS_rank"])
            ig_epitopes={} #dict of {indicies of epitope cores: scores}
            
            for i,row in df.iterrows():
                ig_score = row["IG_pIRS_rank"]
                if ig_score > IG_threshold:
                    print(f'{ig_score} over threshold ({i})')
                    ig_epitope_position = i + int(row["IG_core_pos"])
                    
                    if ig_epitope_position in ig_epitopes.keys():
                        ig_epitopes[ig_epitope_position] = max(ig_score, ig_epitopes[ig_epitope_position])
                    else:
                        ig_epitopes[ig_epitope_position] = ig_score
            print(ig_epitopes)     
            for i,row in df.iterrows():
                if i in ig_epitopes.keys():
                    df["Epitopes"].iloc[i] += 'IG,'
            
        #2nd part is mesuring "humanness"
        df["Epitope_Humanness"]=np.nan
        df.loc[df["Epitopes"] != '', "Epitope_Humanness"] = df.loc[df["Epitopes"] != '', "Core"].apply(get_humanness)        
        
        # 3rd part load epitopes into the binder instance
        for row in df.iterrows():
           if row[1].Epitopes != '':
                epitope=[row[1].Pos, row[1].Core, row[1].Epitopes, get_humanness(row[1].Core)]
                self.epitopes.append(epitope)
       
        
        self.in_silico_metrics=df
        
        if filter_epitopes:
            df = self.in_silico_metrics
            return df[df["Epitopes"] != '']
        
    def pick_mutation_spots_from_metrics(self,n_mutations, BP_weight,NM_weight,DT_weight,IG_weight,active_site_punish = -10):
        
        scoring_weights={"BP" : BP_weight,
                         "NM" : NM_weight,
                         "DT" : DT_weight,
                         "IG" : IG_weight}

        #Scan wich softwares are pressent
        allowed_softwares = ["BP","NM","DT","IG"]
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
                
            weight += row["Contact residue"]*active_site_punish
        
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

    def mutate(self,positions:list,
               report_to_path=None,suffix=None,
               blossumsub_mutation_probs=None,
               bepipred_position_temperature = None,
               BP_threshold=0.1512,
               BP_position_temperature=1,
               position_temperature_logo=None,
               residue_temperature_logo=None,
               blossum_weight_logo=None):

        rows = self.in_silico_metrics.iloc[positions]
        
        for i,row in rows.iterrows():

            kind = row.Epitopes.split(',')[0]
            if kind[0:3] == 'NM_':
                if os.path.exists(f'data/MHCI_position_probabilities/{kind[3:]}_probability.txt'):
                    self.logosub(allotype = kind,
                             mut_pos=row.Pos,
                             report_to_path = report_to_path,
                             position_temperature_logo = position_temperature_logo,
                             residue_temperature_logo =  residue_temperature_logo,
                             blossum_weight_logo = blossum_weight_logo,
                             blossumsub_mutation_probs=blossumsub_mutation_probs)
                else:
                    self.blossumsub(mut_pos = row.Pos, mutation_probs=blossumsub_mutation_probs , report_to_path = report_to_path, epitope_type="MHCI")
                   
            elif kind[0:3] == 'BP_':
                core_pos = kind[3:].split(';')
                core = self.in_silico_metrics.iloc[int(core_pos[0]): int(core_pos[1])+1]
                scores = []
                positions = []
                for i,aa in core.iterrows():
                    positions.append(int(aa["Pos"]))
                    scores.append(float(aa["BP_BepiPred-3.0 linear epitope score"])-BP_threshold - 100*bool(aa["Contact residue"]))
                
                scores = np.array(scores)
                scores-=scores.min()
                scores += 0.0001
                scores = scores**BP_position_temperature
                scores = scores/scores.sum()
                
                mut_pos = np.random.choice(positions, p=scores)
                self.blossumsub(mut_pos = mut_pos, mutation_probs=blossumsub_mutation_probs , report_to_path = report_to_path)
            
            elif kind == 'DT':
                self.blossumsub(mut_pos = row.Pos, mutation_probs=blossumsub_mutation_probs , report_to_path = report_to_path)
            elif kind == 'IG':
                self.blossumsub(mut_pos = np.random.randint(row.Pos, row.Pos+9), mutation_probs=blossumsub_mutation_probs , report_to_path = report_to_path, epitope_type = "MHC-II")
        
        if suffix != None:
            self.name += suffix