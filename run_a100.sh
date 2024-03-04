#!/bin/bash
#BSUB -W 480
#BSUB -n 1
#BSUB -q new_gpu
#BSUB -R "select[a100]"
#BSUB -gpu "num=1:mode=shared:mps=yes"
#BSUB -o tmp/out.%J
#BSUB -e tmp/err.%J
source /usr/share/Modules/init/bash
export JULIA_DEPOT_PATH=~/perm/julia
module load julia/1.8.0
module load cuda/12.0
