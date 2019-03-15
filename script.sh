#!/bin/bash
#PBS -N calcStats_container1.1
#PBS -W group_list=jdpellet
#PBS -q standard
#PBS -l place=pack:shared
#PBS -l select=1:ncpus=2:mem=11gb:pcmem=6gb
#PBS -l walltime=24:00:00
#PBS -l cput=48:00:00
#PBS -M alexprescott@email.arizona.edu
#PBS -m bea

cd ~/ClimateReanalysis/

cp /extra/alexprescott/prec20c/*.tif /tmp/

echo "tif files copied to /tmp/"

module load singularity
date
singularity pull docker://alexprescott/climatereanalysis/1.1
date
singularity run climatereanalysis_1.1.sif
date


