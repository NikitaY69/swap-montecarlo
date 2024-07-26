#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=test
#SBATCH --mem=32GB
#SBATCH --output=test.out
#SBATCH --error=test.err

./test.exe --input /home/allaglo/production/configs/T0.025_N10000/run_81.xy --outdir /home/allaglo/test/ --tau 100000 --lin 50 --log 50 --T 0.025 --N 10000