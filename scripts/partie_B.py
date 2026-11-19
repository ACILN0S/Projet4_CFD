"""Calculs de la Partie B (balayage en Re, schéma centré) : un cas par appel, pour le calculateur Ruche.

48 cas : 11 Re de 100 à 5000 (ceux de Ghia 1982 et d'Erturk 2005, plus des valeurs intermédiaires pour P(Re)),
M1 et M2, N = 64 et 128 ; plus N = 256 à Re = 5000 (grille de Ghia) et à Re = 3200 (convergence en maillage).
Critère d'arrêt : ||R|| < ε ||R(t_ref)||, ε = 1e-3, t_ref = 30 (choisi à Re = 5000, docs/adr/0003 §2) ;
tolérance du gradient conjugué 1e-8 (Partie A).

Usage : python scripts/partie_B.py <indice du cas>     (python scripts/partie_B.py liste : les cas)
Sur Ruche : sbatch scripts/ruche_partie_B.slurm (un cas par tâche du tableau SLURM).
Sortie : results/B_<méthode>_Re<Re>_N<N>_tol1e-08_res1e-03_tref30.npz, relu par le notebook B (run_path).
"""
import os
os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")   # cf. CLAUDE.md, avant numpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import numpy as np

from cavity.projection import run_projection
from cavity.psi_omega import run_psi_omega

RE_LIST = (100, 200, 400, 700, 1000, 1500, 2000, 2500, 3200, 4000, 5000)
CG_TOL, EPS, T_REF, T_MAX = 1e-8, 1e-3, 30.0, 1000.0
RUNNERS = {"M1": run_projection, "M2": run_psi_omega}
# Les plus longs d'abord (Re élevé, N grand, M2) : ils démarrent en premier si les cœurs manquent
CASES = ([(m, re, n) for re in reversed(RE_LIST) for n in (128, 64) for m in ("M2", "M1")]
         + [(m, re, 256) for re in (5000, 3200) for m in ("M2", "M1")])


def run_path(method, re, n):
    """Même nom que run_path du notebook B."""
    return ROOT / "results" / f"B_{method}_Re{re}_N{n}_tol{CG_TOL:.0e}_res{EPS:.0e}_tref{T_REF:g}.npz"


if __name__ == "__main__":
    if sys.argv[1] == "liste":
        for k, case in enumerate(CASES):
            print(k, *case)
        sys.exit(0)
    method, re, n = CASES[int(sys.argv[1])]
    run = RUNNERS[method](n, re, "centered", cg_tol=CG_TOL, residual_tol=EPS, t_ref=T_REF, t_max=T_MAX)
    np.savez(run_path(method, re, n), **{key: np.asarray(val) for key, val in run.items()})
    print(f"{method} Re = {re} N = {n} : t = {run['t_final']:.1f}, {run['n_steps']} pas, "
          f"{run['wall_time'] / 60:.1f} min, arrêt par le critère : {run['converged']}", flush=True)
