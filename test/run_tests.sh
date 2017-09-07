#! /bin/bash

root -l -b -q generate_files.C &> /dev/null

[[ -d tmp ]] || mkdir tmp
export TMPDIR=$(pwd)/tmp

python tests.py
