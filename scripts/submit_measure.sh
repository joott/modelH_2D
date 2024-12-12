#!/usr/bin/env bash
L=64

for mass in `seq -3.86 -0.01 -3.89`; do

    id=${RANDOM}
    TMPFILE=`mktemp --tmpdir=/home/jkott/perm/tmp`
    cp run_arb.sh $TMPFILE

    init=`ls /home/jkott/perm/modelH_2D/thermalized/thermalized_L_64_mass_-3.85_id_* | head -1`
    echo "julia mass_measure.jl --fp64 --init=$init --mass=$mass --rng=$id $L 0.1" >> $TMPFILE
    echo "rm $TMPFILE" >> $TMPFILE

    echo $TMPFILE

    chmod +x $TMPFILE
    bsub < $TMPFILE

done
