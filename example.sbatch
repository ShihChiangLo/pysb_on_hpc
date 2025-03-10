#!/bin/bash
#SBATCH --job-name= NAME_of_the_job        # Job name
#SBATCH --output=NAME_of_the_output_file.log     # Output file
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=1              # Tasks per node
#SBATCH --cpus-per-task=1                # Number of CPU cores per task
#SBATCH --mem=8G                        # Memory allocation
#SBATCH --time=2:00:00                   # Time limit hrs:min:sec

module purge                             # Purge loaded modules
# Run the Python script inside the Singularity container
singularity exec \
    /path_to_the_repository_of_sif_file/pysb_oct6_2024_latest.sif \
    python /path_to_the_repository_of_py_file/example.py