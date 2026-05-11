#!/bin/sh
### General options
### –- specify queue --
#BSUB -q gpua100
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
#BSUB -R "rusage[mem=10GB]"
### -- set the email address --
### -- Specify the output and error file. %J is the job-id --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o /zhome/c8/e/204914/DBL054_RasmusB/scoring/%J.out
#BSUB -e /zhome/c8/e/204914/DBL054_RasmusB/scoring/%J.err
# -- end of LSF options --


metrics_csv="/zhome/c8/e/204914/DBL054_RasmusB/nextflow_stuff/work/e0/fb8b1769c5b3d3b298ec7587c6ab04/CAN-1/scoring_input.csv" /zhome/c8/e/204914/DBL054_RasmusB/nextflow_stuff/modules/binder_scoring.sh
