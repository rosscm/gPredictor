#!/bin/bash
SCRIPT_DIR=$(pwd)
FASTA_DIR=$SCRIPT_DIR/fasta

cd $FASTA_DIR
for i in {1..22} {X..Y}
do
  wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr'$i'.fa.gz'
done

gunzip *.fa.gz
