#!/usr/bin/env bash

for L in 72; do
for mass in -0.045 -0.04 -0.035; do
for i in {1..8}; do

    id=${RANDOM}
    TMPFILE=`mktemp --tmpdir=/home/jkott/perm/tmp`
    cp run_arb.sh $TMPFILE

    echo "julia thermalize.jl --mass=$mass --fp64 --rng=$id $L 0.01" >> $TMPFILE
    echo "rm $TMPFILE" >> $TMPFILE

    echo $TMPFILE

    chmod +x $TMPFILE
    bsub < $TMPFILE

done
done
done
