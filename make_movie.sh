#!/bin/bash
# Template script for running visualizer.py
# Written Aug 5 2015 by Alex DeGrave

################################### VARIABLES ##################################
TOPOLOGY=../../template/movie_reference.pdb
ROOT_ITER=520
ROOT_SEG=99
WESTH5=../../west.h5 # use either this or the fort23
NUMATOMS=113
TPLEN=3
RAW_PDB_DIR=raw_pdbs
BBQ_OUT_DIR=BBQ_pdbs
OUTPUT_MODE=movie
TP_ID=2
FORT23=../../traj_segs/000522/000000/fort.23
################################################################################


rm -r $RAW_PDB_DIR
rm -r $BBQ_OUT_DIR

mkdir $BBQ_OUT_DIR #|| exit 1
mkdir $RAW_PDB_DIR #|| exit 1

python visualizer.py --topo $TOPOLOGY\
    --iter $ROOT_ITER\
    --seg-id $ROOT_SEG\
    --raw-outdir $RAW_PDB_DIR\
    --bbq-outdir $BBQ_OUT_DIR\
    --output-mode $OUTPUT_MODE\
    --westh5 $WESTH5
#    --timepoint-id $TP_ID
#    --fort23 $FORT23\
