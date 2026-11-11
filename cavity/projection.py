"""Méthode M1 : projection / pas fractionnaire (cours §7.3.2, workshop 15) — Projet 4, Parties A, B, C.

Inconnues : u, v, p aux nœuds de la grille colocalisée (formes (N+1, N+1), [i, j] ↔ (x, y)).
Un pas de temps (Euler explicite, cours §3.1) :
    1) prédiction sans pression : u* = u + Δt [ -(u·∇)u + ν Δu ]          (intérieur)
    2) pression                 : Δp = (1/Δt) ∇·u*,  ∂p/∂n = 0, moyenne nulle
                                  (∂p/∂n = 0 d'ordre 2 par nœuds fantômes : docs/adr/0002, §6)
    3) correction               : u^{n+1} = u* - Δt ∇p                    (intérieur)
Conditions aux limites : u = v = 0 sur les parois fixes, u = U_t(t), v = 0 au couvercle,
u = 0 aux deux nœuds de coin du couvercle (ils appartiennent aux parois fixes).

Toutes les briques (opérateurs, Poisson, pas de temps, boucle en temps, critère d'arrêt, K, P)
viennent de cavity/common.py ; ce fichier ne contient que le pas de temps propre à M1.
"""
import numpy as np

from cavity.common import (advection, apply_velocity_bc, laplacian, lid_velocity, make_grid,
                           march_to_steady, solve_poisson, time_step)


def divergence_all_nodes(u, v, h):
    """Divergence ∂u/∂x + ∂v/∂y à TOUS les nœuds, parois comprises (forme (N+1, N+1)).

    Nécessaire avec la Neumann à nœuds fantômes : les nœuds de paroi sont des inconnues de la pression,
    il y faut donc le second membre (docs/adr/0002, §6). Chaque dérivée est :
      - centrée d'ordre 2 aux nœuds où les deux voisins existent (cours §2.2, éq. 2.9) ;
      - décentrée d'ordre 2 aux deux extrémités, (-3 f_0 + 4 f_1 - f_2)/2h et (3 f_N - 4 f_{N-1} + f_{N-2})/2h
        (workshop 11.2, éq. 11.10 ; même formule que pour ∂u/∂y dans P).
    C'est exactement np.gradient(..., edge_order=2). Le long d'une paroi, la dérivée tangentielle est donc
    centrée (valeurs de paroi connues) et la dérivée normale décentrée ; aux coins, les deux sont décentrées.
    Aux nœuds intérieurs, on retrouve la divergence centrée de la version d'ordre 1.
    Aux deux coins du couvercle, u saute de U_t à 0 : la divergence y vaut ~U_t/h (singularité de coin).
    """
    return np.gradient(u, h, axis=0, edge_order=2) + np.gradient(v, h, axis=1, edge_order=2)


def projection_step(u, v, p, t, dt, h, nu, scheme, cg_tol, pressure_bc="neumann_ghost"):
    """Avance (u, v, p) de t à t + Δt par la méthode de projection.

    Entrées : u, v, p de forme (N+1, N+1) au temps t (non modifiés) ; Δt, h, ν ;
              scheme = "centered" ou "upwind" (advection) ; cg_tol = tolérance du gradient conjugué ;
              pressure_bc = "neumann_ghost" (ordre 2, nœuds fantômes, retenu : docs/adr/0002) ou
              "neumann" (ordre 1, φ_0 = φ_1, version d'origine gardée pour la comparaison).
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
    #    Ordre 2 (fantômes) : second membre à tous les nœuds, divergence d'ordre 2 aux parois.
    #    Ordre 1 : second membre aux nœuds intérieurs seulement (divergence centrée).
    #    Dans les deux cas, la divergence utilise les valeurs de paroi de u* (U_t au couvercle).
    div_star = divergence_all_nodes(u_star, v_star, h)
    if pressure_bc == "neumann":
        div_star = div_star[1:-1, 1:-1]
    # Gradient conjugué Neumann, départ à chaud depuis la pression du pas précédent (cours §4.7.3).
    p_new, n_iter = solve_poisson(div_star / dt, h, pressure_bc, phi0=p, tol=cg_tol)

    # 3) Correction : u^{n+1} = u* - Δt ∇p, gradient centré aux nœuds intérieurs seulement
    #    (les valeurs de paroi sont imposées, on ne les corrige pas). Avec les fantômes, p aux nœuds de
    #    paroi est une vraie inconnue d'ordre 2, utilisée ici par le gradient des nœuds voisins des parois.
    u_new = u_star
    v_new = v_star
    u_new[1:-1, 1:-1] -= dt * (p_new[2:, 1:-1] - p_new[:-2, 1:-1]) / (2.0 * h)
    v_new[1:-1, 1:-1] -= dt * (p_new[1:-1, 2:] - p_new[1:-1, :-2]) / (2.0 * h)
    return u_new, v_new, p_new, n_iter


def run_projection(n_cells, reynolds, scheme="centered", cg_tol=1e-8,
                   residual_tol=1e-2, t_ref=10.0, t_max=200.0, pressure_bc="neumann_ghost", dt=None):
    """Lance M1 depuis le repos jusqu'au stationnaire (ou jusqu'à t_max).

    Entrées :
        n_cells    : N (grille (N+1) × (N+1)) ; reynolds : Re, donc ν = 1/Re ;
        scheme     : advection "centered" (défaut) ou "upwind" ;
        pressure_bc : "neumann_ghost" (défaut, ordre 2) ou "neumann" (ordre 1, pour comparaison) ;
        cg_tol     : tolérance (résidu relatif) du gradient conjugué, cf. docs/adr/0003 §3 ;
        residual_tol : ε du critère d'arrêt sur le résidu stationnaire (comme au Projet 1, question b) :
                       arrêt dès que ||R(t)|| < ε ||R(t_ref)||, avec R = (u^{n+1} - u^n)/Δt
                       (common.steady_residual) ; ε = 1e-2 fixé par le test de l'ADR 0003 §2 ;
        t_ref        : instant de référence et début de la surveillance du résidu (t_ref = 10 :
                       rampe du couvercle 1 - e^{-t} terminée) ; on ne prend pas t = 0 car le fluide
                       part du repos (résidu initial nul) ;
        t_max        : garde-fou ; si atteint, le résultat est marqué converged = False ;
        dt           : None (défaut) ⇒ Δt des limites de stabilité (common.time_step) ; sinon Δt imposé,
                       PLUS PETIT que cette limite (études « Δt fixe » et « convergence en Δt » du
                       notebook A : l'état stationnaire de M1 dépend de Δt, docs/adr/0002 §7.3).
    Sortie : dictionnaire avec la grille, les champs finaux (u, v, p), l'historique (t, K, P, itérations
             du gradient conjugué, ||R||) et les informations d'exécution (Δt, nombre de pas, temps de
             calcul, ||R(t_ref)||, converged = arrêt par le critère et non par t_max).
    """
    nu = 1.0 / reynolds
    x, y, h = make_grid(n_cells)
    dt_stable = time_step(h, nu, scheme)          # limites de stabilité du cours §5.3
    if dt is None:
        dt = dt_stable
    elif dt > dt_stable:
        raise ValueError(f"Δt = {dt:.3e} dépasse la limite de stabilité {dt_stable:.3e} (N = {n_cells}, Re = {reynolds})")

    # Départ du repos : u = v = p = 0 (et U_t(0) = 0, cohérent avec la rampe).
    zeros = np.zeros((n_cells + 1, n_cells + 1))
    state = {"u": zeros.copy(), "v": zeros.copy(), "p": zeros.copy()}

    def step(s, t):
        u, v, p, n_iter = projection_step(s["u"], s["v"], s["p"], t, dt, h, nu, scheme, cg_tol,
                                          pressure_bc)
        return {"u": u, "v": v, "p": p}, n_iter

    state, info = march_to_steady(step, state, h, dt, nu, residual_tol, t_ref, t_max,
                                  label=f"M1 (N = {n_cells}, Re = {reynolds})")
    return {"method": "M1", "scheme": scheme, "pressure_bc": pressure_bc, "N": n_cells, "Re": reynolds,
            "nu": nu, "x": x, "y": y, "h": h, "dt": dt, "cg_tol": cg_tol, **state, **info}
