#!/bin/bash
#BSUB -W 2880
#BSUB -n 1
#BSUB -q gpu
#BSUB -R "select[a10 || a30 || a100 || h100 || l40 || l40s]"
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -o tmp/out.%J
#BSUB -e tmp/err.%J
source /usr/share/Modules/init/bash
export PATH=/usr/local/usrapps/tmschaef/jkott/julia-1.10.3/bin:$PATH
export LD_LIBRARY_PATH=;
module load cuda/12.3
