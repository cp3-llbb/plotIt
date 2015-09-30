#!/bin/bash
# Sort all plots ending with different suffixes in corresponding directories

# Suffixes for lepton categories
leptons="mumu elel elmu muel"
# Other suffixes
categs="jj bb bb_MET"
# Plot file extension (can be wildcard)
ext=png

# Default it current directory
basedir=${1=.}

pushd $basedir

for lep in $leptons; do
  for categ in $categs; do

    dir=categ_${lep}_${categ}
    
    if [ ! -d $dir ]; then mkdir $dir ; fi
    
    for file in *_${lep}_${categ}.${ext}; do mv $file $dir ; done

  done
done

popd

