#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=test
#SBATCH --mem=32GB
#SBATCH --output=test.out
#SBATCH --error=test.err


./test.exe --input /home/allaglo/production/configs/T0.04/run_601.xy --outdir /home/allaglo/test/