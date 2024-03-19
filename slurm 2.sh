#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --partition=batch
#SBATCH -J gmap
#SBATCH -o gmap.%J.out
#SBATCH -e gmap.%J.err
#SBATCH --time=04:00:00
#SBATCH --mem=100G

$ module av nextflow
nextflow/20.04.1 nextflow/21.04.2 