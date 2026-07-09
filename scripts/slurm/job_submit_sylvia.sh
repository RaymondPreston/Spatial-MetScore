#!/bin/bash
#SBATCH --job-name=Tumor_Purity
#SBATCH --output=/scratch/hhuan40/Spatial-MetScore/logs/Tumor_Purity_%j.out
#SBATCH --error=/scratch/hhuan40/Spatial-MetScore/logs/Tumor_Purity_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --partition=c64-m512

# Send email notifications
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hhuan40@emory.edu 

echo "Starting Tumor Purity Calculation..."

# Initialize Conda for Bash shell and source it
conda init bash > /dev/null 2>&1
source ~/.bashrc

# Activate your specific R environment
conda activate spatial_R

# Run the R script using Rscript command
Rscript /scratch/hhuan40/Spatial-MetScore/notebooks/tumor_purity.R

echo "Tumor Purity Calculation complete."