"""Tests du module cavity/poisson.py (lot L1, partagé par M1 et M2).

On vérifie, sur des solutions analytiques :
  1. que la DCT (Neumann) et la DST (Dirichlet) convergent à l'ordre 2 vers la solution continue
     (l'erreur vient seulement du laplacien 5 points, pas du solveur) ;
  2. que le résidu discret est au niveau de l'erreur machine (le solveur est « exact ») ;
  3. que le SOR du Projet 1 redonne la même solution que la DCT.
"""
import numpy as np

from cavity.poisson import solve_poisson_dirichlet, solve_poisson_neumann, sor_neumann

GRIDS = [16, 32, 64]


def observed_order(hs, errors):
    """Pente de log(erreur) en fonction de log(h) : l'ordre de convergence observé."""
    return np.polyfit(np.log(hs), np.log(errors), 1)[0]


def neumann_case(N):
    """φ = cos(πx) cos(πy) aux centres des cellules : ∂φ/∂n = 0 sur les 4 parois, Δφ = −2π² φ."""
    h = 1.0 / N
    xc = (np.arange(N) + 0.5) * h
    X, Y = np.meshgrid(xc, xc, indexing="ij")          # [i, j] ↔ (x, y)
    phi = np.cos(np.pi * X) * np.cos(np.pi * Y)
    return h, phi, -2.0 * np.pi**2 * phi


def neumann_laplacian(phi, h):
    """Laplacien 5 points aux centres avec fantômes miroirs (φ_fantôme = φ_intérieur ⇔ ∂φ/∂n = 0)."""
    p = np.pad(phi, 1, mode="edge")
    return (p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2] - 4.0 * phi) / h**2


def test_neumann_dct_order2_and_exact_residual():
    errors, hs = [], []
    for N in GRIDS:
        h, phi_exact, b = neumann_case(N)
        phi = solve_poisson_neumann(b, h)
        # φ n'est défini qu'à une constante près : on compare des champs de moyenne nulle
        errors.append(np.max(np.abs(phi - (phi_exact - phi_exact.mean()))))
        hs.append(h)
        # le résidu discret doit être nul à la précision machine
        assert np.max(np.abs(neumann_laplacian(phi, h) - b)) < 1e-9 * np.max(np.abs(b))
    assert abs(observed_order(hs, errors) - 2.0) < 0.1


def test_dirichlet_dst_order2_and_exact_residual():
    errors, hs = [], []
    for N in GRIDS:
        h = 1.0 / N
        x = np.arange(1, N) * h                             # nœuds intérieurs seulement
        X, Y = np.meshgrid(x, x, indexing="ij")
        psi_exact = np.sin(np.pi * X) * np.sin(2.0 * np.pi * Y)   # nul sur les 4 parois
        b = -5.0 * np.pi**2 * psi_exact                     # Δψ = −(1 + 4) π² ψ
        psi = solve_poisson_dirichlet(b, h)
        errors.append(np.max(np.abs(psi - psi_exact)))
        hs.append(h)
        p = np.pad(psi, 1)                                  # ψ = 0 sur les parois
        lap = (p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2] - 4.0 * psi) / h**2
        assert np.max(np.abs(lap - b)) < 1e-9 * np.max(np.abs(b))
    assert abs(observed_order(hs, errors) - 2.0) < 0.1


def test_sor_matches_dct():
    h, _, b = neumann_case(32)
    phi_sor, n_iter = sor_neumann(np.zeros_like(b), b, h, omega=1.9, tol=1e-10, max_iter=20000)
    phi_dct = solve_poisson_neumann(b, h)
    assert n_iter < 20000
    assert np.max(np.abs(phi_sor - phi_dct)) < 1e-6
