#!/usr/local_rwth/bin/zsh

#SBATCH --job-name=envi_15
#SBATCH --output=logs/df_15_genes_all-%j.out
#SBATCH --error=logs/df_15_genes_all-%j.err
#SBATCH --time=08:00:00
#SBATCH --gres=gpu:1
#SBATCH --mem=180G
#SBATCH --partition=c18g
#SBATCH --cpus-per-task=30
#SBATCH --account=rwth1209 
#SBATCH --signal=2
#SBATCH --nodes=1
#SBATCH --export=ALL
#SBATCH --no-requeue

# Make conda work with jupyter from slurm
# Insert this after any #SLURM commands
export CONDA_ROOT=$HOME/anaconda3
. $CONDA_ROOT/etc/profile.d/conda.sh
export PATH="$CONDA_ROOT/bin:$PATH"
# but naturally before using any python scripts
​
echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "Running on nodes: $SLURM_NODELIST"
echo "------------------------------------------------------------"


# Activate the appropriate conda environment
conda activate /work/rwth1209/enviroments/envi

# Run the command
python /work/rwth1209/dana_projects/test_new_imputation_tools/enVI/save_as_df/df_test.py