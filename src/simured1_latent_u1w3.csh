#!/bin/tcsh
    #BSUB -J simulation[1-2000]
    #BSUB -q stat
    #BSUB -W 2:00
    #BSUB -n 1
    ##BSUB -R "rusage[mem=5GB]"
    #BSUB -o out.%J
    #BSUB -e err.%J
    /usr/local/usrapps/$GROUP/$USER/env_redlining/bin/Rscript ./src/simured1_latent_u1w3.R $LSB_JOBINDEX