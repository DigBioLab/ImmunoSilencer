#Load default settings
from classes import*
import os
run_name = input()

if os.path.exists(f'{run_name}/input/HLA.txt'):
    if os.path.getsize(f'{run_name}/input/HLA.txt') > 0:
        with open(f'{run_name}/input/HLA.txt', 'r') as file:
            for line in file:
                print(line)
                exit()

settings_dict=get_settings(run_name)

allele_pupulation = "data/HLA/MHCI_distribution.csv"
available_HLA_path =  "data/HLA/avilable_alleles"
percentile = float(settings_dict["HLA_percentiles"])

avilable_alleles = []

with open(available_HLA_path) as file:
    for line in file:
        row=line.split('\t')
        for allele in row:
            avilable_alleles.append(allele)

avilable_alleles=set(avilable_alleles)

df = pd.read_csv(allele_pupulation, sep= ";").to_numpy()[1:]

def get_top_percentile(df_locus, percentile):
    df_locus.sort(key= lambda c: -float(c[1].replace(',','.')))
    df_locus = np.array(df_locus)
    fractile=0
    actual_fractile=0
    alleles=[]
    for allele in df_locus:
        fractile+=float(allele[1].replace(',','.'))
        if allele[0].replace(":",'') in avilable_alleles:
            alleles.append(allele)
            actual_fractile += float(allele[1].replace(',','.'))
        else:
            pass
#            print(f"warning {allele[0]} not modeled in NetMHCIpan-4.2. Frequency is {allele[1]}")
            
        if fractile>=percentile:
#            print(f"Total fractile: {actual_fractile}")
            return alleles

df_A=[]
df_B=[]
df_C=[]

for row in df:
    if row[0][4]=='A':
        df_A.append(row)
    elif row[0][4]=='B':
        df_B.append(row)
    elif row[0][4]=='C':
        df_C.append(row)


HLA_As = list(get_top_percentile(df_A,percentile))
HLA_Bs = list(get_top_percentile(df_B,percentile))
HLA_Cs = list(get_top_percentile(df_C,percentile))

all_HLAs = [ele[0].replace(":",'') for ele in (HLA_As + HLA_Bs + HLA_Cs)]

print(",".join(list(all_HLAs)))
        