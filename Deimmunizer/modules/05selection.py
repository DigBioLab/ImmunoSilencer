from sys import argv
from classes import*
settings_dict = get_settings('.')

ipae_limit = 8
ipsae_limit = 0.8
n_binders=int(input())

argv.pop(0)
candidates = argv[:]


def score_epitopes(candidate_dir,
                   NM_weight_selection=settings_dict["NM_weight_selection"],
                   BP_weight_selection=settings_dict["BP_weight_selection"],
                   DT_weight_selection=settings_dict["DT_weight_selection"],
                   IG_weight_selection=settings_dict["IG_weight_selection"]):
    
    epitope_df = pd.read_csv(f'{candidate_dir}/metrics_epitopes.csv', sep = '\t')
    weights = {'NM': NM_weight_selection,
               'BP': BP_weight_selection,
               'DT': DT_weight_selection,
               'IG': DT_weight_selection}
    score=0
    for flag,humanness in zip(list(epitope_df["Epitopes"]), list(epitope_df["Epitope_Humanness"])):
        if float(humanness) == 0:
            continue
        s=set([word[0:2] for word in flag.split(',') if word != ''])

        for prefix in s:
            score+=weights[prefix]
    return(score)

def score_binding(binbder_dir):
    df = pd.read_csv(f"{binbder_dir}/ipsae_and_ipae.csv")
    return df["af3_ipSAE_min"].iloc[0]
    

def candidate_evaluation(candidate_dir, original_dir="archive/round_0/CAN-1", ipsae_hard_limit = None):
    '''
    Evaluates the performance of a binder with mutaions to the original binder.
    The score becomes 1, if epitope score increased thus always loosing to the default.
    The score is then a function of %epitope cahnge and %ipsea change
    
    '''
    
    score_can = score_epitopes(candidate_dir)
    score_original = score_epitopes(original_dir)
    ipsae_min_can = score_binding(candidate_dir)
    ipsae_min_original = score_binding(original_dir)

    
    if score_can - score_original > 0:
        return 1
    
    delta_epitope_ratio = (score_can - score_original)/(score_original + 0.00001)
    delta_binding_ratio = (ipsae_min_can - ipsae_min_original)/(ipsae_min_original+0.0001)

    
    
    if ipsae_hard_limit != None:
        if ipsae_min_can < ipsae_hard_limit:
            return 1
        
    return (delta_epitope_ratio*100)**3 - (delta_binding_ratio*100)**5 

candidates.sort(key = candidate_evaluation)


next_generation = []

for can in candidates:
    for _ in range(int(settings_dict["top_coppies"])):
        next_generation.append(can)


for binder in next_generation[0:n_binders]:
    print(binder)


