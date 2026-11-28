#!/bin/bash
# Partie C sur le Mac (Apple M1, 4 cœurs performance) : les cas de scripts/partie_C.py, 4 à la fois,
# un processus par cas et un seul thread par processus (comme en Partie B, CLAUDE.md §3).
# Le 5ᵉ cas démarre dès qu'un cœur se libère (xargs -P 4). caffeinate -i empêche la mise en veille
# tant que les calculs tournent (laisser le Mac branché).
# Usage, depuis la racine du projet : bash scripts/lancer_partie_C.sh   (ou avec des indices : ... 5 6)
# Journaux : results/logs/C_<indice>.log ; suivi : tail -f results/logs/C_*.log
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results/logs
export MKL_THREADING_LAYER=SEQUENTIAL OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONUNBUFFERED=1

CASES="${*:-0 1 2 3 4}"     # par défaut tous les cas, les plus longs d'abord (ordre de CASES)
echo "$(date '+%F %T') lancement des cas : ${CASES}"
echo ${CASES} | tr ' ' '\n' | caffeinate -i xargs -P 4 -I{} \
    sh -c '.venv/bin/python scripts/partie_C.py {} > results/logs/C_{}.log 2>&1 \
           && echo "$(date "+%F %T") cas {} terminé" || echo "$(date "+%F %T") cas {} ÉCHEC (voir results/logs/C_{}.log)"'
echo "$(date '+%F %T') tous les cas sont finis"
