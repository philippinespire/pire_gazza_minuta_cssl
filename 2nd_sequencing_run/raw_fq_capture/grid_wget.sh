#!/bin/bash

#SBATCH --job-name=Wget
#SBATCH -o Wget-%j.out
#SBATCH --time=00:00:00
#SBATCH --exclusive
#SBATCH --ntasks=40
#SBATCH -p main

# this script download files (meant for fq.gz files) from TAMUCC "whebshare/grid"

# EXECUTE
# place this script where you want to download the files 
# get the grid parent directory (which is just a website) where the files will be downloaded from. Usually Sharon will post this in the species slack channel.
# visually check in your browser that the is a file name tamucc_files.txt, which list the rest of the files in that directory
# execute as "sbatch grid_wget.sh <grid parent directory>

module load parallel

# save outdir and  parent grid directory into variables
OUTDIR=$1
GRIDDIR=$2

# downlownd the list of lifes
wget $GRIDDIR/tamucc_files.txt

# use that list to download all files in parallel
cat tamucc_files.txt | parallel --no-notice -kj40 wget -P $OURDIR $GRIDDIR/{} 
