#!/bin/bash
#PBS -N wget_CHELSA
#PBS -W group_list=jdpellet
#PBS -q standard
#PBS -l place=pack:shared
#PBS -l select=1:ncpus=1:mem=5gb:pcmem=6gb
#PBS -l walltime=8:00:00
#PBS -l cput=8:00:00
#PBS -M alexprescott@email.arizona.edu
#PBS -m bea

cd /extra/alexprescott/prec20c/
date
wget -nv -nd -r -N https://www.wsl.ch/lud/chelsa/data/timeseries20c/prec/ -a log2.txt
date

