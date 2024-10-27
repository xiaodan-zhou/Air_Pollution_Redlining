#!/bin/tcsh
    #BSUB -J simulation[1-20]
    #BSUB -q stat
    #BSUB -W 60:00
    #BSUB -n 1
    ##BSUB -R "rusage[mem=5GB]"
    #BSUB -o out.%J
    #BSUB -e err.%J
    /usr/local/usrapps/$GROUP/$USER/env_redlining/bin/Rscript ./src/redlining_model_latent_u1w3.R $LSB_JOBINDEX
