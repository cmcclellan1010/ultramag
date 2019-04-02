#!/bin/bash

pushd ./
cd /media/connor/Lore/ultramag/PRESTO/
for f in *.dat
do
    echo "FFTing file $f..."
    realfft "$f"
done
popd

