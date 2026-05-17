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

## user_settings.txt
This file contains any parameter you wish to override from the default_settings.txt file. This is useful to store relevant paramters in for easier overeview, but is equivlant to having the same definitions in default_settings.txt

## HLA.txt
Force a certain set of HLA-alleles to be evaluated instead of the usual population based approach. This text file must contain a single line containg all allels seperated with a comma and no spaces eg. HLA-A0101,HLA-B0702,HLA-C0602.


# Pipeline descirption

## Initial_scoring.nf
- Bepipred-3.0, DiscoTope-3.0 and netMHCpan-4.2 is run on the A-chain of the input binder
- de novo binder scoring is run (using only alphafold) on the sequences of the input pdb.
- the putput is copied into CAN-1, CAN-2,...,CAN-n, with n being the $treewidth$

## main_loop.nf
- Mutations are picked at random (guided by parameters)
- de novo binder scoring is run (using only alphafold) on new sequneces
- Selection of best performing mutations. Selected binders go to next loop until $treedepth$ is met



