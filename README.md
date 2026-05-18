# ImmunoSilencer module
This repository contains the code for a module that takes a pdb of a binder and target and looks for mutations that remove predicted epitopes
for t-cells and/or b-cells by trying and testing different spot mutations through multiple iterations, keeping the best ones (based on epitope removal and af3 ipSAE_min).

# Input structure
Input folder may be located anywhere, but for safty should be kept outside of Deimmunizer/ in case of naming overlaps.
The Input folder must conatin the following (optinal subfolders and files are in paranthesis):

input_dir/
|
-- binder_and_target.pdb
-- default_settings.txt
-- (user_settings.txt)
-- (HLA.txt)

## binder_and_target.pdb
This file must contain a pdb with no more than one model. The A chain must be the binder, and subsequent chains are combined as the target.
This is the only file where the naming is irelevant (.pdb is mandetory), but there must be exactly one file ending in .pdb in the input folder.

## default_settings.txt
This file contains all settings for the run. All paramters must be defined. See exhaustive list in Parameters.

## (user_settings.txt)
This file contains any parameter you wish to override from the default_settings.txt file. This is useful to store relevant paramters in for easier overeview, but is equivlant to having the same definitions in default_settings.txt

## (HLA.txt)
Force a certain set of HLA-alleles to be evaluated instead of the usual population based approach. This text file must contain a single line containg all allels seperated with a comma and no spaces eg. HLA-A0101,HLA-B0702,HLA-C0602.


# Pipeline descirption

## Initial_scoring.nf
- Bepipred-3.0, DiscoTope-3.0 and netMHCpan-4.2 is run on the A-chain of the input binder
- de novo binder scoring is run (using only alphafold) on the sequences of the input pdb.
- the putput is copied into CAN-1, CAN-2,...,CAN-n, with n being the $treewidth$

## main_loop.nf
- Mutations are picked at random (guided by parameters)
- de novo binder scoring is run (using only alphafold) on new sequneces, yielding both foldings and binding prediction of new candidates.
new candidates are given the suffix 'I'. eg CAN-1I
- Selection of best performing mutations. Selected binders go to next loop until $treedepth$ is met. new seed candidates are given the suffix '-x' with x being 1-$treewidth$ eg. CAN-1I-1, in the case where CAN-1I won the selection.

## Mutation cases
Mutations are defined by 2 variables: position and residue. These are picked depending on what software has flagged the selected epitope.

- DiscoTope-3.0: Position is always the position of the epitope flag. The residue is picked from distribution weighted by the BLOSUM62 substition scores of the current residue (Blosmsub).

- BepiPred-3.0: Position is picked from a distribution weigthed by BepiPred scores upstream of the flagged position. Residue is picked with blosumsub

- netMHCpan-4.2: Position is picked based on the highest log-probabilty of pressent residues in the 9-mer core. Residue is picked based on the lowest log-probabilty of the position. These are specific to the HLA-allele that flagged the postion. (Logosub)

# Requirements 
This module only works when running from dtu HPC and by having:
- Execute permissions in: /dtu/projects/RFdiffusion/closed-loop/
- python on PATH, with Bio.PDB and pandas.
- Access to GPU

# How to run
Ensure you have a folder with the structure explain in *{Input structure}. Define the 4 lines in the top of main.sh. All dir-paths must be from root and end with '/'
 - working_dir: path of the Deimmunizer dir
 - input_dir:  path to input dir
 - output_dir: path to output dir. This dir must not contain a subfolder with the same run name as the newly initiated run
 - run_name: Simply for structure in your output. Must be unqie, must not be left empty. "/" not alloved.
Then execute main.sh. 

To run with, you need GPU, so either be on interactive node, or queue a job using bsub < main.sh. 


# Parameters
## Complexity
- Treedepth: Number of iterations in the main loop
- Treewidth: Number of candidates in each main loop
- mutations_per_cycle: Number of mutations done in in each main loop per candidate
- top_coppies: Number of coppies made of each top placement in selection. Copying halts when $treewidth$ is met.
## thresholds
- BP_threshold: Threshold for BepiPred-3.0 epitope flag
- NM_threshold: Threshold for netMHCpan-4.2 epitope flag
- DT_threshold: Threshold for DiscoTope-3.0 epitope flag
## Temperatures
All temperatures control how biased some selections are made. All selections are based on probabilty distributions generated from different scores. Temperatures of 0 force the maximum score to be picked, as temperature rises, the distributions aproach uniform.

- residue_temperature_blossum: Control the varivance of residue selection in $Blosumsub$
- BP_position_temperature: control the variance in selecting position in Bepipred-3.0 predicted epitopes
- position_temperature_logo: control variance in selecting position in logosub
- residue_temperature_logo: control variance in selecting residue in logosub

## weights
Each software has two weights assosiated with them. A mutation weight and a selectino weight.
- NM_weight, BP_weight, DT_weight: Weight of selecting the an epitope flag with the given software. If set to 0 The software results are completly ignored in all steps in the entire pipeline.

- NM_weight_selection, BP_weight_selection, DT_weight_selection: How much the selection step punishes the pressens of epitope flags of different softwares.

- blossum_weight_logo: must be between (or equal to) 0 and 1. How much residue selection in Logosub is based  on blosum62. 1 is equivilant to blosumsub and 0 is purely based on log-probability of residues.


## other
- Maximum_bepi_length: length of BepiPred epitope cores.
- contact_residue_threshold: The distance thrshold distance between a binder residue and the target before you flag the residue as a contact residue. Contact residues are hevily penalized in mutation selection.


# Output structure
The initial module produces a round_0 folder containing the metrics of the input binder.

Each iteration produces new mutatinos and the metrics of these are stored in round_x. The final generation is stored in final_gen, with the winner being the CAN with the final suffix '-1'. A copy of the input is also saved to the output folder. 