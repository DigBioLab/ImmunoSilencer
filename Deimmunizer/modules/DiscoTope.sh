#!/bin/bash
input_pdb=$1
source /dtu/projects/RFdiffusion/scripts_collection/binder_desing_metrics/miniconda3/bin/activate
conda activate /dtu/projects/RFdiffusion/closed-loop/discotope3_web/Discotope_env
 
python /dtu/projects/RFdiffusion/closed-loop/discotope3_web/src/predict_webserver.py \
--pdb_or_zip_file "$input_pdb" \
--struc_type solved \
--out_dir disco_out \
--cpu_only \
--models_dir /dtu/projects/RFdiffusion/closed-loop/discotope3_web/models

conda deactivate