"""Calculs de la Partie C (Re = 10 000 > Re_c ≈ 8000, perte de stationnarité) : un cas par appel.

Choix de l'équipe (PLAN.md, étape 5) : Re = 10 000, N = 128 puis 256, M1 et M2, Euler explicite, départ du repos
(même condition initiale qu'en Parties A et B). 11 cas :
    - M1 et M2, advection centrée, Δt de stabilité (common.time_step) : calculs principaux ;
    - M2 centré avec Δt/2 : erreur d'Euler sur la période et l'amplitude (l'état de M2 ne dépend pas de Δt,
      contrairement à M1 dont la projection non incrémentale mêle une erreur de projection, docs/adr/0002 §7) ;
    - M1 et M2 upwind : choix « malheureux » (diffusion numérique D = ½ U h (1 - C) ≫ ν, cours §5.4.2) ;
    - M1 et M2 centrés à N = 256 (grille de Ghia à Re = 10 000) : sensibilité au maillage ; Erturk et al. (2005)
      trouvent un calcul périodique sur 129² mais stationnaire sur 257² à Re = 7500 (ref_data/README.md).
    - M1 et M2 centrés à N = 128 jusqu'à t = 1000 : à t = 300, l'écoulement moyen n'est pas établi (option A du 24/11) ;
    - M2 puis M1 centrés à N = 128 jusqu'à t = 2000 (demandes du 25/11).
Pas d'arrêt sur le résidu (residual_tol = 0) : on calcule jusqu'à t_max, le régime pouvant être périodique.
Tolérance du gradient conjugué 1e-8, comme en Parties A et B (docs/adr/0003 §3).

Usage (depuis la racine du projet) :
    python scripts/partie_C.py liste            les cas et leur indice
    python scripts/partie_C.py bench <indice>   3000 pas, temps par pas (estimation du coût, rien n'est sauvé)
    python scripts/partie_C.py <indice>         calcul complet
Lancement de tous les cas en parallèle sur le Mac : scripts/lancer_partie_C.sh.
Sortie : results/C_<méthode>_<schéma>_Re10000_N<N>_dt<Δt>.npz, relu par notebooks/C_unsteady.ipynb.
"""
import os
os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")   # cf. CLAUDE.md, avant numpy
import platform
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import numba
import numpy as np

from cavity.common import time_step
from cavity.projection import run_projection
from cavity.psi_omega import run_psi_omega

RE = 10000
CG_TOL, T_MAX = 1e-8, 300.0
RUNNERS = {"M1": run_projection, "M2": run_psi_omega}
# (méthode, schéma, facteur sur le Δt de stabilité, N, t_max) ; indices 0-4 : N = 128, les plus longs d'abord ;
# indices 5-6 : N = 256, ajoutés ensuite ; indices 7-8 : M1 et M2 centrés jusqu'à t = 1000 (option A du 24/11 :
# à t = 300, l'écoulement moyen évolue encore et l'oscillation est amortie ; s'éteint-elle ou s'entretient-elle ?).
# Indices 9 et 10 : M2 et M1 centrés jusqu'à t = 2000 (demandes du 25/11).
# On ajoute toujours à la fin : un indice désigne toujours le même cas.
CASES = [("M2", "centered", 0.5, 128, T_MAX), ("M1", "centered", 1.0, 128, T_MAX), ("M2", "centered", 1.0, 128, T_MAX),
         ("M1", "upwind", 1.0, 128, T_MAX), ("M2", "upwind", 1.0, 128, T_MAX),
         ("M2", "centered", 1.0, 256, T_MAX), ("M1", "centered", 1.0, 256, T_MAX),
         ("M1", "centered", 1.0, 128, 1000.0), ("M2", "centered", 1.0, 128, 1000.0),
         ("M2", "centered", 1.0, 128, 2000.0), ("M1", "centered", 1.0, 128, 2000.0)]

# Sondes (u, v enregistrés à chaque pas), en coordonnées (x, y) ; toutes tombent sur des nœuds des grilles 128 et 256
# (mêmes points physiques pour comparer les deux grilles).
#   - à 2 mailles (N = 128) des coins du couvercle : là où le Reynolds de maille est maximal, des zigzags nœud à nœud y
#     signaleraient une instationnarité numérique et non la bifurcation (PLAN.md §2, « Singularités aux coins ») ;
#   - centre (tourbillon principal) et près des trois autres coins (tourbillons secondaires BL, BR, TL).
PROBES = {"coin_couvercle_gauche": (2 / 128, 126 / 128), "coin_couvercle_droit": (126 / 128, 126 / 128),
          "centre": (0.5, 0.5), "bas_gauche": (0.125, 0.125), "bas_droit": (0.875, 0.125),
          "haut_gauche": (0.125, 0.875)}

def snapshot_times(t_max):
    """Instantanés : tous les 5 (évolution sur tout le calcul), puis tous les 0.2 sur les 20 derniers
    (plusieurs cycles, ~20 images par cycle pour les périodes mesurées, 3.6 à 5)."""
    return np.concatenate([np.arange(5.0, t_max - 20.0 + 1e-9, 5.0),
                           np.arange(t_max - 20.0 + 0.2, t_max + 1e-9, 0.2)])


def case_dt(scheme, factor, n):
    """Δt du cas : limite de stabilité du cours (common.time_step) multipliée par factor (≤ 1)."""
    return factor * time_step(1.0 / n, 1.0 / RE, scheme)


def run_path(method, scheme, factor, n, t_max=T_MAX):
    """Nom du fichier de résultats, relu par le notebook C (suffixe _t<t_max> seulement si t_max ≠ 300)."""
    suffix = "" if t_max == T_MAX else f"_t{t_max:g}"
    return ROOT / "results" / f"C_{method}_{scheme}_Re{RE}_N{n}_dt{case_dt(scheme, factor, n):.2e}{suffix}.npz"


def run_case(method, scheme, factor, n, t_max, with_records):
    """Lance un cas ; with_records : sondes, instantanés et journal (calcul complet seulement)."""
    probe_ij = [(round(x * n), round(y * n)) for x, y in PROBES.values()]
    extra = {"probes": probe_ij, "snapshot_times": snapshot_times(t_max), "log_every": 5.0} if with_records else {}
    return RUNNERS[method](n, RE, scheme, cg_tol=CG_TOL, residual_tol=0.0, t_max=t_max,
                           dt=case_dt(scheme, factor, n), **extra)


if __name__ == "__main__":
    if sys.argv[1] == "liste":
        for k, (method, scheme, factor, n, t_max) in enumerate(CASES):
            print(k, method, scheme, f"N = {n}", f"Δt = {case_dt(scheme, factor, n):.3e}", f"t_max = {t_max:g}",
                  f"pas ≈ {t_max / case_dt(scheme, factor, n):.2e}")
        sys.exit(0)

    if sys.argv[1] == "bench":
        method, scheme, factor, n, t_max = CASES[int(sys.argv[2])]
        RUNNERS[method](16, 100, scheme, t_max=0.05)      # compilation Numba hors mesure
        n_bench = 3000
        run = run_case(method, scheme, factor, n, n_bench * case_dt(scheme, factor, n) * (1 - 1e-9), with_records=False)
        ms = 1e3 * run["wall_time"] / run["n_steps"]
        print(f"{method} {scheme} Δt×{factor:g} N = {n} : {ms:.2f} ms/pas, {run['cg_iterations'].mean():.0f} it. GC en moyenne ; "
              f"calcul complet ≈ {ms * t_max / case_dt(scheme, factor, n) / 3.6e6:.1f} h", flush=True)
        sys.exit(0)

    method, scheme, factor, n, t_max = CASES[int(sys.argv[1])]
    run = run_case(method, scheme, factor, n, t_max, with_records=True)
    run["probe_names"] = np.array(list(PROBES))
    run["probe_xy"] = np.array(list(PROBES.values()))
    run["dt_factor"] = factor
    run["machine"] = f"{platform.platform()} | {platform.processor()} | Python {platform.python_version()} | " \
                     f"numpy {np.__version__} | numba {numba.__version__}"
    run["date"] = datetime.now().isoformat(timespec="seconds")
    path = run_path(method, scheme, factor, n, t_max)
    np.savez(path, **{key: np.asarray(val) for key, val in run.items()})
    print(f"{method} {scheme} Δt×{factor:g} N = {n} : t = {run['t_final']:.1f}, {run['n_steps']} pas, "
          f"{run['wall_time'] / 3600:.2f} h -> {path.name}", flush=True)
