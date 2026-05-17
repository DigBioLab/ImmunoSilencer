from ImmunoSilencer.Deimmunizer.modules.classes import*
from sys import argv
run_name = input()
settings_dict=get_settings(f'../../../{run_name}')

argv.pop(0) # I don't need this scripts file name

while argv != []:
    arg = argv.pop(0)
    if arg == '-t':
        target_path = argv.pop(0)
        target = load_binders(target_path) #
        target_chains=[chain.name[-1] for chain in target] #contains eg. ['B','C','D'] or ['B']
        
    if arg == '-a':
        chain_A = argv.pop(0)
        binder = load_binders(chain_A)[0]
        binder.name = binder.name.split('_')[0] #Always just one

if len(target) == 1:
    target = target[0] #unpack the single target chain
    
    df = pd.DataFrame({"binder_id": [binder.name],
                "binder_chain": ['A'], # hardcoded
                "target_id": [target.name+'_target'], #change later
                "A_seq": [binder.sequence],
                "target_chains": ['["B"]'], #expand this to more chains
                "target_chain_range": [f'["1:{len(target)}"]'],
                "msa_info": ['["A:run_msa", "B:run_msa"]'],
                "target_subchain_B_seq": [target.sequence]})

else:
    df = pd.DataFrame({"binder_id": [binder.name],
                       "binder_chain":['A'],
                       "target_id":[binder.name+'_target'],
                       "A_seq":[binder.sequence],
                       "target_chains":['[' + ', '.join(f'"{chain}"' for chain in target_chains) + ']'],
                       "target_chain_range":['[' + ', '.join(f'"1:{len(chain.sequence)}"' for chain in target) + ']'],
                       "msa_info":['[' + '"A:run_msa", ' +', '.join(f'"{chain}:run_msa"' for chain in target_chains) + ']']      
    })
    
    for chain_name, chain in zip(target_chains,target):
        df[f'target_subchain_{chain_name}_seq'] = chain.sequence

print(f'scoring_input.csv')
df.to_csv(f'scoring_input.csv',index=False)
