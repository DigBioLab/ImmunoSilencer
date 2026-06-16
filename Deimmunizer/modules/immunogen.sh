#!/bin/sh
(
    source /dtu/projects/RFdiffusion/closed-loop/ImmunoGeNN/env_ImmunoGeNN/bin/activate
    #python /dtu/projects/RFdiffusion/closed-loop/ImmunoGeNN/run.py --fasta_file $1 --outdir IG_out --skip_plots true --model_dir ../../../data/IG_model --human_references_pkl ../../../data/IG_data_record/human_references_9mers.pkl.lz4
    python /dtu/projects/RFdiffusion/closed-loop/ImmunoGeNN/run.py --fasta_file $1 --outdir IG_out --skip_plots true --human_references_pkl ../../../data/IG_data_record/human_references_9mers.pkl.lz4
)