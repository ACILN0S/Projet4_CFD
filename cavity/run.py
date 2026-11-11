"""Boucle en temps de la méthode 1 (lot L5) : `run(config) -> dict de résultats`.

Structure reprise du Projet 3 (`run_simulation` → dict). On part du repos (ou d'un état fourni),
on avance avec `mac_solver.step`, et tous les `n_diag` pas on enregistre :
    t, K, P, ε, le résidu stationnaire ‖u^{n+1} − u^n‖∞ / Δt et max |∇·u|.
Arrêt : t ≥ t_final, ou résidu < tol_steady (si tol_steady > 0).
(Sauvegarde/reprise .npz : à ajouter pour les longs runs de la Partie C.)
"""
import time

import numpy as np

from cavity.diagnostics import dissipation, kinetic_energy, lid_power
from cavity.mac_solver import apply_bc, divergence, new_fields, step
from cavity.params import lid_velocity


def run(cfg, u0=None, v0=None):
    """Lance une simulation M1 décrite par `cfg` (cavity.params.Config).

    u0, v0 : état initial optionnel (continuation depuis un autre run) ; sinon repos.
    Retourne un dict : séries "t", "K", "P", "eps", "residual", "div_max" (tableaux numpy),
    champs finaux "u", "v", "n_steps", "steady" (bool), "wall_time" (s), et la config "cfg".
    """
    N, h, nu, dt = cfg.N, cfg.h, cfg.nu, cfg.dt
    upwind = {"centered": False, "upwind1": True}[cfg.scheme]
    if u0 is None:
        u, v = new_fields(N)
    else:
        u, v = u0.copy(), v0.copy()
    apply_bc(u, v, lid_velocity(0.0))

    series = {name: [] for name in ("t", "K", "P", "eps", "residual", "div_max")}
    div = np.empty((N, N))

    def record(t, residual):
        """Enregistre les diagnostics de l'état courant (u, v) au temps t."""
        divergence(u, v, h, div)
        series["t"].append(t)
        series["K"].append(kinetic_energy(u, v, h))
        series["P"].append(lid_power(u, nu, h, lid_velocity(t)))
        series["eps"].append(dissipation(u, v, nu, h))
        series["residual"].append(residual)
        series["div_max"].append(np.max(np.abs(div)))

    t, n, steady = 0.0, 0, False
    record(t, np.nan)
    start = time.perf_counter()
    n_steps_total = int(round(cfg.t_final / dt))
    while n < n_steps_total:
        u_old = u
        u, v = step(u, v, t, dt, nu, h, upwind)
        n += 1
        t = n * dt                 # pas d'accumulation d'erreur d'arrondi sur t
        if n % cfg.n_diag == 0 or n == n_steps_total:
            # résidu stationnaire sur les faces u intérieures (‖∂u/∂t‖∞ ≈ 0 en régime établi)
            residual = np.max(np.abs(u[1:-1, 1:-1] - u_old[1:-1, 1:-1])) / dt
            record(t, residual)
            if cfg.tol_steady > 0 and residual < cfg.tol_steady:
                steady = True
                break

    result = {name: np.array(values) for name, values in series.items()}
    result.update(u=u, v=v, n_steps=n, steady=steady,
                  wall_time=time.perf_counter() - start, cfg=cfg)
    return result
