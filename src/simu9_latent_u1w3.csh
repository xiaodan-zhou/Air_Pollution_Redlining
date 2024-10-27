#!/bin/tcsh
    #BSUB -J simulation[1-1000]
    #BSUB -q stat
    #BSUB -W 00:30
    #BSUB -n 1
    ##BSUB -R "rusage[mem=5GB]"
    #BSUB -o out.%J
    #BSUB -e err.%J
    /usr/local/usrapps/$GROUP/$USER/env_redlining/bin/Rscript ./src/simu9_latent_u1w3.R $LSB_JOBINDEX
