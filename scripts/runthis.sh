#SBATCH --job-name=tough_run
#SBATCH --verbose
#SBATCH --time=50:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70GB
#SBATCH --array=1-3
#SBATCH --mail-type=END
#SBATCH --mail-user=xm595@nyu.edu
#SBATCH --output=sample%A_%a.out
#SBATCH --error=sample%A_%a.err


echo "$(pwd)"
bash ./pipeline.sh sample_$SLURM_ARRAY_TASK_ID.txt
