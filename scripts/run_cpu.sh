#!/bin/bash
#BSUB -W 120
#BSUB -n 16
#BSUB -R span[hosts=1] 
#BSUB -o tmp/out.%J
#BSUB -e tmp/err.%J
source /usr/share/Modules/init/bash
export PATH=/usr/local/usrapps/tmschaef/jkott/julia-1.10.3/bin:$PATH
export LD_LIBRARY_PATH=;
