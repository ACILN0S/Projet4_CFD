"""Méthode M1 : projection / pas fractionnaire (cours §7.3.2, workshop 15) — Projet 4, Parties A, B, C.

Inconnues : u, v, p aux nœuds de la grille colocalisée (formes (N+1, N+1), [i, j] ↔ (x, y)).
Un pas de temps (Euler explicite, cours §3.1) :
    1) prédiction sans pression : u* = u + Δt [ -(u·∇)u + ν Δu ]          (intérieur)
    2) pression                 : Δp = (1/Δt) ∇·u*,  ∂p/∂n = 0, moyenne nulle
    3) correction               : u^{n+1} = u* - Δt ∇p                    (intérieur)
Conditions aux limites : u = v = 0 sur les parois fixes, u = U_t(t), v = 0 au couvercle,
u = 0 aux deux nœuds de coin du couvercle (ils appartiennent aux parois fixes).

Toutes les briques (opérateurs, Poisson, pas de temps, K, P) viennent de cavity/common.py.
"""
import time

import numpy as np

from cavity.common import (advection, kinetic_energy, laplacian, lid_power, lid_velocity,
                           make_grid, solve_poisson, time_step)


def apply_velocity_bc(u, v, u_lid):
    """Impose (en place) les vitesses de paroi.

    Parois fixes : u = v = 0. Couvercle (j = N) : u = U_t sur les nœuds intérieurs i = 1..N-1,
    v = 0. Les deux coins du couvercle gardent u = 0 (choix du PLAN, §2 « Nœud de coin »).
    """
    u[0, :] = 0.0
    u[-1, :] = 0.0
    u[:, 0] = 0.0
    u[1:-1, -1] = u_lid
    v[0, :] = 0.0
    v[-1, :] = 0.0
    v[:, 0] = 0.0
    v[:, -1] = 0.0


def projection_step(u, v, p, t, dt, h, nu, scheme, cg_tol):
    """Avance (u, v, p) de t à t + Δt par la méthode de projection.

    Entrées : u, v, p de forme (N+1, N+1) au temps t (non modifiés) ; Δt, h, ν ;
              scheme = "centered" ou "upwind" (advection) ; cg_tol = tolérance du gradient conjugué.
    Sorties : u, v, p au temps t + Δt (nouveaux tableaux) et le nombre d'itérations du gradient conjugué.
    """
    # 1) Prédiction (Euler explicite, cours §3.1) : advection + diffusion, sans pression.
    #    Advection centrée ou upwind (cours §2.2, §5.6.2), diffusion centrée ordre 2 (§4.2.1).
    u_star = u.copy()
    v_star = v.copy()
    u_star[1:-1, 1:-1] += dt * (-advection(u, v, u, h, scheme) + nu * laplacian(u, h))
    v_star[1:-1, 1:-1] += dt * (-advection(u, v, v, h, scheme) + nu * laplacian(v, h))
    # Les vitesses de paroi sont imposées au nouveau temps t + Δt (u* = U_t(t+Δt) au couvercle).
    apply_velocity_bc(u_star, v_star, lid_velocity(t + dt))

    # 2) Pression : Δp = (1/Δt) ∇·u* (cours §7.3.2, « discrete pressure equation »).
    #    Divergence centrée d'ordre 2 aux nœuds intérieurs (cours §2.2, éq. 2.9) ; elle utilise
    #    les valeurs de paroi de u* (en particulier U_t au couvercle).
    div_star = ((u_star[2:, 1:-1] - u_star[:-2, 1:-1])
                + (v_star[1:-1, 2:] - v_star[1:-1, :-2])) / (2.0 * h)
    # Gradient conjugué Neumann, départ à chaud depuis la pression du pas précédent (cours §4.7.3).
    p_new, n_iter = solve_poisson(div_star / dt, h, "neumann", phi0=p, tol=cg_tol)

    # 3) Correction : u^{n+1} = u* - Δt ∇p, gradient centré aux nœuds intérieurs seulement
    #    (les valeurs de paroi sont imposées, on ne les corrige pas).
    u_new = u_star
    v_new = v_star
    u_new[1:-1, 1:-1] -= dt * (p_new[2:, 1:-1] - p_new[:-2, 1:-1]) / (2.0 * h)
    v_new[1:-1, 1:-1] -= dt * (p_new[1:-1, 2:] - p_new[1:-1, :-2]) / (2.0 * h)
    return u_new, v_new, p_new, n_iter


def run_projection(n_cells, reynolds, scheme="centered", cg_tol=1e-8,
                   steady_tol=1e-6, t_min=10.0, t_max=200.0):
    """Lance M1 depuis le repos jusqu'au stationnaire (ou jusqu'à t_max).

    Entrées :
        n_cells    : N (grille (N+1) × (N+1)) ; reynolds : Re, donc ν = 1/Re ;
        scheme     : advection "centered" (défaut) ou "upwind" ;
        cg_tol     : tolérance (résidu relatif) du gradient conjugué, cf. docs/adr/0003 ;
        steady_tol : arrêt quand max|u^{n+1} - u^n| et max|v^{n+1} - v^n| < steady_tol ...
        t_min      : ... et t > t_min (la rampe du couvercle 1 - e^{-t} doit être terminée) ;
        t_max      : garde-fou ; si atteint, le résultat est marqué converged = False.
    Sortie : dictionnaire avec la grille, les champs finaux, l'historique (t, K, P, itérations du
             gradient conjugué) et les informations d'exécution (Δt, nombre de pas, temps de calcul,
             max|∂u/∂t| à l'arrêt, converged).
    """
    nu = 1.0 / reynolds
    x, y, h = make_grid(n_cells)
    dt = time_step(h, nu, scheme)                 # limites de stabilité du cours §5.3

    # Départ du repos : u = v = p = 0 (et U_t(0) = 0, cohérent avec la rampe).
    u = np.zeros((n_cells + 1, n_cells + 1))
    v = np.zeros_like(u)
    p = np.zeros_like(u)

    hist_t, hist_k, hist_p, hist_iter, hist_change = [], [], [], [], []
    t = 0.0
    n_steps = 0
    change = np.inf
    wall_start = time.perf_counter()

    while t < t_max:
        u_new, v_new, p, n_iter = projection_step(u, v, p, t, dt, h, nu, scheme, cg_tol)
        t += dt
        n_steps += 1
        # Variation sur un pas (même idée que le résidu des méthodes itératives, cours §4.3.2)
        change = max(np.abs(u_new - u).max(), np.abs(v_new - v).max())
        u, v = u_new, v_new

        # Diagnostics imposés par l'énoncé, à chaque pas (coût négligeable devant le Poisson)
        hist_t.append(t)
        hist_k.append(kinetic_energy(u, v, h))
        hist_p.append(lid_power(u, h, nu, lid_velocity(t)))
        hist_iter.append(n_iter)
        hist_change.append(change)            # pour suivre l'approche du stationnaire

        if not np.isfinite(change):
            raise FloatingPointError(f"M1 diverge à t = {t:.4f} (N = {n_cells}, Re = {reynolds}, Δt = {dt:.3e})")
        if change < steady_tol and t > t_min:
            break

    converged = change < steady_tol and t > t_min
    return {
        "method": "M1", "scheme": scheme, "N": n_cells, "Re": reynolds, "nu": nu,
        "x": x, "y": y, "h": h, "dt": dt, "u": u, "v": v, "p": p,
        "t": np.array(hist_t), "K": np.array(hist_k), "P": np.array(hist_p),
        "cg_iterations": np.array(hist_iter), "change": np.array(hist_change),
        "n_steps": n_steps, "t_final": t, "wall_time": time.perf_counter() - wall_start,
        "max_change": change, "max_dudt": change / dt, "converged": converged, "cg_tol": cg_tol,
    }
