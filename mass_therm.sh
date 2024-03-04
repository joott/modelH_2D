#!/usr/bin/env bash
L=64

for mass in `seq -2.0 -0.1 -4.0`; do

    id=${RANDOM}
    TMPFILE=`mktemp --tmpdir=/home/jkott/perm/tmp`
    cp run_cpu.sh $TMPFILE

    echo "julia -t 16 mass_thermalize.jl --cpu --fp64 --mass=$mass --rng=$id $L 1.0" >> $TMPFILE
    echo "rm $TMPFILE" >> $TMPFILE

    echo $TMPFILE

    chmod +x $TMPFILE
    bsub < $TMPFILE

done
