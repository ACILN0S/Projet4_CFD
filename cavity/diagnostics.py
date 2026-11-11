"""Diagnostics globaux de la méthode 1 (grille MAC) : K, P, ε (lot L4).

Fonctions pures : elles lisent (u, v) — CL déjà appliquées, fantômes à jour — et renvoient un nombre.
Positions MAC et décalage +1 des fantômes : voir l'en-tête de cavity/mac_solver.py.

Choix des formules (docs/pseudocode.md §1.7) : ce sont celles pour lesquelles le bilan DISCRET
    dK/dt = P − ε
est exact pour le schéma centré. On les obtient en multipliant l'équation discrète de u par u et en
sommant par parties (« produit scalaire » discret Σ u² h²). Ce sont aussi des approximations d'ordre 2
des intégrales continues de l'énoncé.
"""
import numpy as np


def kinetic_energy(u, v, h):
    """K = ½ ∫ (u² + v²) dS, chaque composante à sa place sur la grille MAC.

        K = ½ h² [ Σ u² sur les faces u intérieures + Σ v² sur les faces v intérieures ]
    (les faces de paroi ont une vitesse normale nulle et ne contribuent pas).
    """
    return 0.5 * h**2 * (np.sum(u[1:-1, 1:-1] ** 2) + np.sum(v[1:-1, 1:-1] ** 2))


def lid_power(u, nu, h, U_lid):
    """Puissance injectée par le couvercle P = ∫₀¹ ν ∂u/∂y|_{y=1} U_t dx.

    ∂u/∂y au couvercle = (u_fantôme − u_intérieur) / h : différence centrée sur la paroi (ordre 2),
    puisque la paroi est à mi-chemin entre les deux. Intégrale : rectangles de largeur h centrés
    sur les N−1 faces u intérieures (i = 1..N−1). Les deux faces de coin (i = 0, N) sont sur les
    parois latérales (u = 0) et ne contribuent pas au bilan discret.
    Attention (CLAUDE.md §3.1) : à cause de la singularité des coins, P croît comme log N.
    """
    du_dy = (u[1:-1, -1] - u[1:-1, -2]) / h          # fantôme du haut − dernière face intérieure
    return nu * U_lid * np.sum(du_dy) * h


def dissipation(u, v, nu, h):
    """Dissipation visqueuse ε = ν ∫ |∇u|² dS  (= ν ∫ ω² dS en continu pour ces conditions aux limites).

    Somme des carrés de TOUTES les différences entre voisins (Δ/h)² × h² = Δ², avec un poids ½ pour
    les différences qui traversent une paroi (fantôme ↔ première face intérieure) : c'est
    exactement le terme de bord qui sort de la sommation par parties du laplacien.
    """
    # u : différences en x entre faces voisines (paroi u = 0 incluse), i = 0..N−1, rangées intérieures
    du_x = u[1:, 1:-1] - u[:-1, 1:-1]
    # u : différences en y, faces intérieures en x, de la rangée fantôme du bas à celle du haut
    du_y = u[1:-1, 1:] - u[1:-1, :-1]
    # v : différences en y (paroi v = 0 incluse), colonnes intérieures
    dv_y = v[1:-1, 1:] - v[1:-1, :-1]
    # v : différences en x, de la colonne fantôme de gauche à celle de droite
    dv_x = v[1:, 1:-1] - v[:-1, 1:-1]

    def half_ends(d2, axis):
        """Somme de d² avec un poids ½ sur la première et la dernière différence (celles de paroi)."""
        first = np.take(d2, 0, axis=axis)
        last = np.take(d2, -1, axis=axis)
        return np.sum(d2) - 0.5 * (np.sum(first) + np.sum(last))

    total = (np.sum(du_x**2) + half_ends(du_y**2, axis=1)
             + np.sum(dv_y**2) + half_ends(dv_x**2, axis=0))
    return nu * total
