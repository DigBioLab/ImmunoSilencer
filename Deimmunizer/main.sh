#!/bin/sh
### General options
### –- specify queue --
#BSUB -q gpuv100
### -- set the job Name --
#BSUB -J binder_scoring_Deimunizer
### -- ask for number of cores (default: 1) --
#BSUB -n 8
### -- specify that the cores must be on the same host --
#BSUB -R "span[hosts=1]"
### -- Select the resources: 1 gpu in exclusive process mode --
#BSUB -gpu "num=1:mode=exclusive_process"
### -- set walltime limit: hh:mm --  maximum 24 hours for GPU-queues right now
#BSUB -W 24:00
# request 5GB of system-memory
#BSUB -R "rusage[mem=5GB]"
### -- Specify the output and error file. %J is the job-id --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o /zhome/c8/e/204914/DBL054_RasmusB/scoring/%J.out
#BSUB -e /zhome/c8/e/204914/DBL054_RasmusB/scoring/%J.err
# -- end of LSF options --

### --- specify directories --- ### Must be from root
working_dir="/zhome/c8/e/204914/DBL054_RasmusB/Deimmunizer/" #location of the Deimmunizer dir
input_dir="/zhome/c8/e/204914/DBL054_RasmusB/input2/" #Dir containing pdb and Default_settings.txt
output_dir="/zhome/c8/e/204914/DBL054_RasmusB/output/" #Dir must be empty
run_name="DT_full_retest" #must be unqie if you want to run multiple runs simutaniosly. "/" not alloved. May not be left empty
### --- ------------------- --- ###

cd $working_dir

## fetching leaf names ##
input_leaf="${input_dir%/}"
input_leaf="${input_leaf##*/}/"
mkdir -p $output_dir
output_leaf="${output_dir%/}"
output_leaf="${output_leaf##*/}/"

## copying and naming input
mkdir -p "${input_dir}data"
mkdir -p "${working_dir}${run_name}"
cp -r "${input_dir}" "${working_dir}${run_name}"
mv "${working_dir}${run_name}/${input_leaf}" "${working_dir}${run_name}/input"
input_run="${working_dir}${run_name}/input"
## -- ##



treewidth=$(grep "^treewidth" "${input_run}/default_settings.txt" | cut -d',' -f2)
treedepth=$(grep "^treedepth" "${input_run}/default_settings.txt" | cut -d',' -f2)


### Overrides from user_settings ###
val=$(grep '^treewidth,' "${input_run}/user_settings.txt" | cut -d',' -f2)
[ -n "$val" ] && treewidth="$val"

val=$(grep '^treedepth,' "${input_run}/user_settings.txt" | cut -d',' -f2)
[ -n "$val" ] && treedepth="$val"

### ---------------------------- ###

input_pdb=("$input_run"/*.pdb)

# Check if exactly one file exists
if [ ${#input_pdb[@]} -ne 1 ] || [ ! -e "${input_pdb[0]}" ]; then
    echo "Error: Expected exactly one .pdb file in $input_dir"
    exit 1
fi

input_pdb="${input_pdb##*/}"
input_pdb="${input_pdb%.*}" #isolate the name of the binder

echo "${run_name}" | python modules/00load_settings.py > "${input_run}/data/HLA_alleles.txt" #Load the HLA alleles

mkdir -p "${run_name}/tmp"
mkdir -p "${run_name}/alive_binders"

nextflow run Initial_scoring.nf --run_name "${run_name}"


## ----------------- Setting up width ----------------## 
for i in $(seq 1 $treewidth)  
do
    mkdir "${run_name}/alive_binders/CAN-${i}"
    cp -r "${run_name}"/input/data/${input_pdb}/* "${run_name}/alive_binders/CAN-${i}"
    mv "${run_name}/alive_binders/CAN-${i}/${input_pdb}".pdb "${run_name}/alive_binders/CAN-${i}/CAN-${i}.pdb"
done

cd "${run_name}"

mkdir -p archive
mkdir archive/round_0
cp -r alive_binders/* archive/round_0

## iteraring through the pipe line
for i in $(seq 1 $treedepth)  
do
    ls alive_binders
    cd ../
    nextflow run main_loop.nf --run_name "${run_name}"
    cd "${run_name}"

    #select what binders get to go to next loop
    rm -rf tmp/*
    mkdir archive/round_$i
    cp -r new_gen/* archive/round_$i
    j=1    
    echo $treewidth | python ../modules/05selection.py alive_binders/* new_gen/* | while read f; do
        base=$(basename "$f")
        cp -r "$f" "tmp/${base}-$j"
        mv tmp/${base}-$j/*.pdb tmp/${base}-$j/${base}-$j.pdb
        ((j++))
    done

    ls tmp
    rm -rf alive_binders/*
    mv tmp/* alive_binders/
    rm -rf tmp/* new_gen/*
done

cd "${working_dir}${run_name}"
mkdir archive/final_gen
cp -r alive_binders/* archive/final_gen/
mkdir archive/input
cp -r input/* archive/input/


##-- send to outputdir --##
mkdir "${output_dir}${run_name}"
mv archive "${output_dir}${run_name}"

## -- Clean up -- ##
cd ..
rm -rf "${run_name}"