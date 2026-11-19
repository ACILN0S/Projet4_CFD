#!/bin/bash
# Partie B sur Ruche (mésocentre Paris-Saclay) : 48 tâches indépendantes, un cas par tâche
# (indice = SLURM_ARRAY_TASK_ID, cf. CASES dans scripts/partie_B.py ; « python scripts/partie_B.py liste »).
# Doc : https://mesocentre.pages.centralesupelec.fr/user_doc/ruche/08_slurm_examples/
# Soumission, depuis $WORKDIR/Projet4_CFD : sbatch scripts/ruche_partie_B.slurm
# Durée : quelques minutes aux bas Re ; à Re = 5000, ~1 h (N = 128) et plus d'une journée (M2, N = 256) ; d'où 72 h.
#SBATCH --job-name=B_sweep
#SBATCH --output=results/logs/B_%a.o%A
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=72:00:00
#SBATCH --partition=cpu_long
#SBATCH --array=0-47

set -x
cd "${SLURM_SUBMIT_DIR}"
export MKL_THREADING_LAYER=SEQUENTIAL OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
"${SLURM_SUBMIT_DIR}/env/bin/python" scripts/partie_B.py "${SLURM_ARRAY_TASK_ID}"
