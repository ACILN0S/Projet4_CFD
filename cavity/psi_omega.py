"""Méthode M2 : fonction de courant / vorticité (cours §7.3.4) — Projet 4, Parties A, B, C.

Inconnues : ω (vorticité) et ψ (fonction de courant) aux nœuds de la grille colocalisée
(formes (N+1, N+1), [i, j] ↔ (x, y)) ; u, v en sont déduites.
Convention (docs/pseudocode.md) : u = ∂ψ/∂y, v = -∂ψ/∂x, ω = ∂v/∂x - ∂u/∂y, donc Δψ = -ω.

Un pas de temps (Euler explicite, cours §3.1) :
    1) transport de ω     : ω^{n+1} = ω^n + Δt [ -(u ∂ω/∂x + v ∂ω/∂y) + ν Δω ]   (intérieur)
    2) fonction de courant : Δψ = -ω^{n+1},  ψ = 0 sur les parois (Dirichlet)
    3) vitesses            : u = ∂ψ/∂y, v = -∂ψ/∂x (centré), vitesses de paroi imposées
    4) ω de paroi          : formule de Thom (développement de Taylor de ψ, outil du cours §2.1)
Pas de pression, donc pas de condition de Neumann : c'est ce qui rend la comparaison M1 / M2
utile pour la question de l'ordre de convergence (docs/adr/0002).

Toutes les briques (opérateurs, Poisson, pas de temps, boucle en temps, critère d'arrêt, K, P)
viennent de cavity/common.py, comme pour M1 ; ce fichier ne contient que le pas de temps propre à M2.
"""
import numpy as np

from cavity.common import (advection, apply_velocity_bc, laplacian, lid_velocity, make_grid,
                           march_to_steady, solve_poisson, time_step)


def velocity_from_streamfunction(psi, h, u_lid):
    """Vitesses à partir de ψ : u = ∂ψ/∂y, v = -∂ψ/∂x (différences centrées d'ordre 2, cours §2.2).

    Entrées : ψ de forme (N+1, N+1), h, U_t. Sorties : u, v de forme (N+1, N+1).
    Aux nœuds intérieurs : u[i,j] = (ψ[i,j+1] - ψ[i,j-1]) / 2h,  v[i,j] = -(ψ[i+1,j] - ψ[i-1,j]) / 2h.
    Aux parois : vitesses imposées (u = U_t au couvercle, 0 ailleurs), comme pour M1.
    """
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)
    u[1:-1, 1:-1] = (psi[1:-1, 2:] - psi[1:-1, :-2]) / (2.0 * h)
    v[1:-1, 1:-1] = -(psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2.0 * h)
    apply_velocity_bc(u, v, u_lid)
    return u, v


def apply_thom_wall_vorticity(omega, psi, h, u_lid):
    """Vorticité de paroi par la formule de Thom (en place). Dérivation : docs/pseudocode.md.

    Taylor (cours §2.1) vers l'intérieur, avec ψ = 0 et ∂ψ/∂n = vitesse tangentielle à la paroi,
    et ω = -∂²ψ/∂n² sur la paroi (∂²ψ/∂t² = 0 le long d'une paroi où ψ = 0) :
        couvercle (j = N) : ω[i,N] = -2 ψ[i,N-1] / h² - 2 U_t / h
        bas       (j = 0) : ω[i,0] = -2 ψ[i,1]   / h²
        gauche    (i = 0) : ω[0,j] = -2 ψ[1,j]   / h²
        droite    (i = N) : ω[N,j] = -2 ψ[N-1,j] / h²
    Formules appliquées aux nœuds de paroi hors coins (i ou j = 1..N-1). Les quatre coins ne servent
    à aucun stencil 5 points des nœuds intérieurs : on les laisse à 0 (aux coins du couvercle, ω est
    de toute façon singulière, cf. PLAN.md §2).
    Thom est d'ordre 1 localement mais préserve l'ordre 2 global (Erturk et al. 2005, §3 ; ref_data/README.md).
    """
    h2 = h * h
    omega[1:-1, -1] = -2.0 * psi[1:-1, -2] / h2 - 2.0 * u_lid / h     # couvercle
    omega[1:-1, 0] = -2.0 * psi[1:-1, 1] / h2                         # bas
    omega[0, 1:-1] = -2.0 * psi[1, 1:-1] / h2                         # gauche
    omega[-1, 1:-1] = -2.0 * psi[-2, 1:-1] / h2                       # droite
    for i, j in ((0, 0), (0, -1), (-1, 0), (-1, -1)):
        omega[i, j] = 0.0


def psi_omega_step(omega, psi, u, v, t, dt, h, nu, scheme, cg_tol):
    """Avance (ω, ψ, u, v) de t à t + Δt par la méthode ψ–ω.

    Entrées : ω, ψ, u, v de forme (N+1, N+1) au temps t (non modifiés) ; Δt, h, ν ;
              scheme = "centered" ou "upwind" (advection de ω) ; cg_tol = tolérance du gradient conjugué.
    Sorties : ω, ψ, u, v au temps t + Δt (nouveaux tableaux) et le nombre d'itérations du gradient conjugué.
    """
    # 1) Transport de ω (Euler explicite, cours §3.1 ; équation de §7.3.4) aux nœuds intérieurs.
    #    Advection centrée ou upwind (cours §2.2, §5.6.2), diffusion centrée ordre 2 (§4.2.1),
    #    avec les vitesses et la vorticité de paroi du temps t.
    omega_new = omega.copy()
    omega_new[1:-1, 1:-1] += dt * (-advection(u, v, omega, h, scheme) + nu * laplacian(omega, h))

    # 2) Δψ = -ω : gradient conjugué Dirichlet (ψ = 0 sur les parois), départ à chaud (cours §4.7.3).
    psi_new, n_iter = solve_poisson(-omega_new[1:-1, 1:-1], h, "dirichlet", phi0=psi, tol=cg_tol)

    # 3) Vitesses au temps t + Δt (couvercle à U_t(t+Δt)).
    u_lid = lid_velocity(t + dt)
    u_new, v_new = velocity_from_streamfunction(psi_new, h, u_lid)

    # 4) Vorticité de paroi (Thom) à partir du nouveau ψ et de U_t(t+Δt) : elle servira au pas suivant.
    apply_thom_wall_vorticity(omega_new, psi_new, h, u_lid)
    return omega_new, psi_new, u_new, v_new, n_iter


def run_psi_omega(n_cells, reynolds, scheme="centered", cg_tol=1e-8,
                  residual_tol=1e-2, t_ref=10.0, t_max=200.0, dt=None):
    """Lance M2 depuis le repos jusqu'au stationnaire (ou jusqu'à t_max).

    Mêmes entrées, même pas de temps (common.time_step : ω obéit à une équation d'advection-diffusion
    de même ν), même critère d'arrêt (résidu stationnaire sur u et v, docs/adr/0003 §2) et mêmes
    diagnostics que M1 (run_projection) : seule la méthode change. dt : None ⇒ Δt de stabilité ;
    sinon Δt imposé, plus petit que cette limite (témoin de l'étude « convergence en Δt »).
    Sortie : dictionnaire avec la grille, les champs finaux (ω, ψ, u, v), l'historique et les
             informations d'exécution (mêmes clés que M1).
    """
    nu = 1.0 / reynolds
    x, y, h = make_grid(n_cells)
    dt_stable = time_step(h, nu, scheme)          # limites de stabilité du cours §5.3
    if dt is None:
        dt = dt_stable
    elif dt > dt_stable:
        raise ValueError(f"Δt = {dt:.3e} dépasse la limite de stabilité {dt_stable:.3e} (N = {n_cells}, Re = {reynolds})")

    # Départ du repos : ω = ψ = u = v = 0.
    zeros = np.zeros((n_cells + 1, n_cells + 1))
    state = {"omega": zeros.copy(), "psi": zeros.copy(), "u": zeros.copy(), "v": zeros.copy()}

    def step(s, t):
        omega, psi, u, v, n_iter = psi_omega_step(s["omega"], s["psi"], s["u"], s["v"],
                                                  t, dt, h, nu, scheme, cg_tol)
        return {"omega": omega, "psi": psi, "u": u, "v": v}, n_iter

    state, info = march_to_steady(step, state, h, dt, nu, residual_tol, t_ref, t_max,
                                  label=f"M2 (N = {n_cells}, Re = {reynolds})")
    return {"method": "M2", "scheme": scheme, "N": n_cells, "Re": reynolds, "nu": nu,
            "x": x, "y": y, "h": h, "dt": dt, "cg_tol": cg_tol, **state, **info}
