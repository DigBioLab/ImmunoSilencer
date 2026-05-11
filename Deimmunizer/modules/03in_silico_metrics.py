from classes import*
from sys import argv
run_name=input()
#Load default settings
settings_dict = get_settings(f'../../../{run_name}')
softwares_all = []
out_file = argv.pop() #This is where I put my final csv
argv.pop(0) # I don't need this scripts file name

while argv != []:
    arg = argv.pop(0)
    if arg == '-b':
        bepi_path = argv.pop(0)
        softwares_all.append("BP")
    if arg == '-d':
        disco_path = argv.pop(0)
        softwares_all.append("DT")    
    if arg == '-n':
        netmhci_path = argv.pop(0)
        softwares_all.append("NM")
    if arg == '-f':
        fasta_path = argv.pop(0)
    if arg =='-p':
        full_pdb = argv.pop(0)

softwares=[]
for prefix in softwares_all:
    if float(settings_dict[f'{prefix}_weight']) != 0:
        softwares.append(prefix)

binder = load_binders(fasta_path)[0]
with open("other.dat", 'r') as file:
    file.readline()#scratch the first line
    MHCI_allotypes = file.readline().replace('\n','').split(',')
    print(MHCI_allotypes)

binder.get_contact_residues2(full_pdb, distance_threshold=settings_dict["contact_residue_threshold"])


if "BP" in softwares:
    binder.load_BepiPred_results(bepi_path)
     
if "NM" in softwares:
    binder.load_NetMHI_results_single(netmhci_path)

if "DT" in softwares:
    binder.load_DiscoTope_results(disco_path)

df = binder.merge_software_results()
binder.in_silico_metrics = df
metrics_in_epitopes = binder.get_epitopes_from_metrics(softwares=softwares,
                                                       MHCI_allotypes = MHCI_allotypes,
                                                       filter_epitopes=True,
                                                       maximum_bepi_length=settings_dict["maximum_bepi_length"],
                                                       BP_threshold = settings_dict["BP_threshold"],
                                                       DT_threshold = settings_dict["DT_threshold"],
                                                       NM_threshold = settings_dict["NM_threshold"])
binder.in_silico_metrics.to_csv(out_file,sep='\t')
metrics_in_epitopes.to_csv(f'{out_file[0:-4]}_epitopes.csv' ,sep='\t')
