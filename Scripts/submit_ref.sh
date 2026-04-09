#!/bin/bash
#SBATCH --job-name=ref_prep          # Name of the job [cite: 560]
#SBATCH --account=general            # Use 'general' for standard CPU partitions [cite: 561]
#SBATCH --partition=c64-m512         # Default general-purpose CPU partition [cite: 1359]
#SBATCH --nodes=1                    # Number of nodes requested [cite: 562]
#SBATCH --ntasks=1                   # Number of tasks [cite: 563]
#SBATCH --cpus-per-task=4            # Number of CPU cores for this task (adjust as needed) [cite: 595]
#SBATCH --mem=32G                    # Amount of RAM required (adjust based on your dataset size) [cite: 576]
#SBATCH --time=04:00:00              # Max time limit (Hours:Minutes:Seconds) [cite: 572]
#SBATCH --output=%x_%j.out           # Standard output log (creates ref_prep_JOBID.out) [cite: 573]
#SBATCH --error=%x_%j.err            # Standard error log (creates ref_prep_JOBID.err) [cite: 574]

# Initialize Conda for the Bash shell [cite: 588]
conda init bash > /dev/null 2>&1     # [cite: 586]
source ~/.bashrc                     # [cite: 587]

# Activate your specific environment 
conda activate Spatial_bio        # [cite: 589]

# Execute the python script
# Assuming $USER resolves to your netID, otherwise replace it explicitly
python /scratch/rprest2/Spatial-MetScore/Scripts/reference_preparation.py