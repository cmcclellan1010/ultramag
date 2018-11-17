#!/bin/bash

pushd ./
cd /media/connor/Lore/ultramag/PRESTO/TESTO/
find . -name '*.dat' -exec realfft {} +
popd

