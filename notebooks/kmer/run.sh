#!/bin/bash

python -m jcvi.assembly.kmer histogram reads.histo Nepal 21 \
    --pdf \
    --vmin=2 \
    --vmax=200
