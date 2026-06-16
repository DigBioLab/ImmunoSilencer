from classes import*
from sys import argv
import pandas as pd

run_name = input()
settings_dict=get_settings(f'../../../{run_name}')


argv.pop(0)

while argv != []:
    arg = argv.pop(0)
    if arg == '-f':
        fasta_path = argv.pop(0)
        binder = load_binders(fasta_path)[0]
        binder.name = binder.name.split('_')[0]
        
    if arg == '-m':
        metrics = pd.read_csv(argv.pop(0),sep='\t')
        
    if arg == '-t':
        target_path = argv.pop(0)
        target = load_binders(target_path)[0]
        
        
###-----define blossumSub metrics-------------####
with open("data/blosum62.txt",'r') as file:
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

binder.in_silico_metrics = metrics

spots = binder.pick_mutation_spots_from_metrics(n_mutations=settings_dict["mutations_per_cycle"],
                                                BP_weight=settings_dict["BP_weight"],
                                                NM_weight=settings_dict["NM_weight"],
                                                DT_weight=settings_dict["DT_weight"],
                                                IG_weight=settings_dict["IG_weight"])

rows = binder.in_silico_metrics.iloc[spots]

binder.mutate(spots,
              report_to_path='info.txt',
              suffix='I',
              blossumsub_mutation_probs=blossumsub_mutation_probs,
              BP_threshold=settings_dict["BP_threshold"],
              BP_position_temperature=settings_dict["BP_position_temperature"],
              position_temperature_logo=settings_dict["position_temperature_logo"],
              residue_temperature_logo=settings_dict["residue_temperature_logo"],
              blossum_weight_logo=settings_dict["blossum_weight_logo"])

with open(f"{binder.name}/{binder.name}_chain_A",'w') as file:
    print(binder, file=file)

