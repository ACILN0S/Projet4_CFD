"""Solveurs de l'équation de Poisson Δφ = b sur le carré unité (lot L1, partagé par M1 et M2).

Trois solveurs, tous pour le laplacien 5 points d'ordre 2 :

- `solve_poisson_neumann` (DCT)   : Neumann homogène, inconnues aux CENTRES des cellules.
                                   Utilisé par M1 (projection de la pression, grille MAC).
- `solve_poisson_dirichlet` (DST) : Dirichlet homogène, inconnues aux NŒUDS intérieurs.
                                   Utilisé par M2 (Δψ = −ω avec ψ = 0 sur les parois).
- `sor_neumann` (SOR)             : même problème que la DCT, par la méthode itérative du cours
                                   (Projet 1, cellules 13 et 17). Sert à comparer coût et précision.

Pourquoi une transformée en cosinus/sinus résout le problème exactement :
sur une grille uniforme, les vecteurs cos(π k (i + ½)/N) (resp. sin(π k i/N)) sont des vecteurs
propres du laplacien discret 1D avec Neumann (resp. Dirichlet), de valeur propre
    λ_k = (2 cos(π k / N) − 2) / h².
En 2D, λ_{kl} = λ_k + λ_l. On passe donc b dans la base propre (DCT/DST), on divise par λ_{kl},
et on revient : c'est une résolution DIRECTE (pas d'itération) en O(N² log N).

Convention d'indices : [i, j] ↔ (x, y).
"""
import numba
import numpy as np
from scipy.fft import dctn, dstn, idctn, idstn


def solve_poisson_neumann(b, h):
    """Résout Δφ = b, ∂φ/∂n = 0, par DCT-II (inconnues aux centres des cellules).

    Discrétisation : (φ[i+1,j] + φ[i-1,j] + φ[i,j+1] + φ[i,j-1] − 4 φ[i,j]) / h² = b[i,j],
    avec fantômes miroirs φ[−1, j] = φ[0, j], etc. (dérivée normale nulle sur la paroi).

    Entrée  : b, forme (N, N), aux centres ((i + ½) h, (j + ½) h). Doit être de somme nulle
              (condition de compatibilité de Neumann) ; c'est le cas pour b = ∇·u* en M1.
    Sortie  : φ, forme (N, N), de moyenne nulle (φ n'est défini qu'à une constante près).
    """
    N = b.shape[0]
    k = np.arange(N)
    lam_1d = (2.0 * np.cos(np.pi * k / N) - 2.0) / h**2       # valeurs propres 1D (Neumann)
    lam = lam_1d[:, None] + lam_1d[None, :]                     # λ_{kl} = λ_k + λ_l
    b_hat = dctn(b, type=2, norm="ortho")                       # b dans la base des cosinus
    lam[0, 0] = 1.0                 # mode constant (λ = 0) : on évite la division par zéro…
    phi_hat = b_hat / lam
    phi_hat[0, 0] = 0.0             # … et on le met à 0 ⇔ φ de moyenne nulle
    return idctn(phi_hat, type=2, norm="ortho")


def solve_poisson_dirichlet(b, h):
    """Résout Δψ = b, ψ = 0 sur les parois, par DST-I (inconnues aux nœuds intérieurs).

    Discrétisation : laplacien 5 points aux nœuds (i h, j h), i, j = 1..N−1, avec ψ = 0 aux nœuds
    de paroi (i ou j = 0 ou N).

    Entrée  : b, forme (N−1, N−1), aux nœuds intérieurs.
    Sortie  : ψ, forme (N−1, N−1), aux mêmes nœuds (on ajoutera les zéros de paroi dans M2).
    """
    N = b.shape[0] + 1
    k = np.arange(1, N)
    lam_1d = (2.0 * np.cos(np.pi * k / N) - 2.0) / h**2       # valeurs propres 1D (Dirichlet)
    lam = lam_1d[:, None] + lam_1d[None, :]                     # toutes < 0 : pas de mode nul
    return idstn(dstn(b, type=1, norm="ortho") / lam, type=1, norm="ortho")


@numba.njit(cache=True)
def sor_neumann(phi, b, h, omega, tol, max_iter):
    """Résout Δφ = b, ∂φ/∂n = 0 aux centres par SOR (repris du Projet 1, cellules 13 et 17).

    Même équation discrète que `solve_poisson_neumann`. Au lieu de fantômes, on écrit directement
    la ligne de la matrice : avec φ_fantôme = φ_intérieur, le voisin manquant disparaît de
    l'équation, et le coefficient diagonal vaut le NOMBRE de voisins réels (4 à l'intérieur,
    3 le long d'une paroi, 2 dans un coin) :
        φ_GS[i,j] = (Σ_voisins φ − h² b[i,j]) / n_voisins
        φ[i,j]    = (1 − ω) φ[i,j] + ω φ_GS[i,j]          (ω = 1 : Gauss–Seidel)

    Critère d'arrêt (comme au Projet 1) : ‖r‖₂ / ‖b‖₂ < tol, avec r = b − Δ_h φ.
    Entrées : phi (N, N) estimation initiale, modifiée EN PLACE ; b (N, N) de somme nulle.
    Sorties : (φ de moyenne nulle, nombre d'itérations).
    """
    N = phi.shape[0]
    b_norm = np.sqrt(np.sum(b * b))
    for it in range(max_iter):
        # --- un balayage SOR (lexicographique) ---
        for i in range(N):
            for j in range(N):
                s = 0.0
                n = 0
                if i > 0:
                    s += phi[i - 1, j]
                    n += 1
                if i < N - 1:
                    s += phi[i + 1, j]
                    n += 1
                if j > 0:
                    s += phi[i, j - 1]
                    n += 1
                if j < N - 1:
                    s += phi[i, j + 1]
                    n += 1
                phi_gs = (s - h * h * b[i, j]) / n
                phi[i, j] = (1.0 - omega) * phi[i, j] + omega * phi_gs
        # --- résidu r = b − Δ_h φ (même stencil), calculé sans tableau intermédiaire ---
        r2 = 0.0
        for i in range(N):
            for j in range(N):
                s = 0.0
                n = 0
                if i > 0:
                    s += phi[i - 1, j]
                    n += 1
                if i < N - 1:
                    s += phi[i + 1, j]
                    n += 1
                if j > 0:
                    s += phi[i, j - 1]
                    n += 1
                if j < N - 1:
                    s += phi[i, j + 1]
                    n += 1
                r = b[i, j] - (s - n * phi[i, j]) / (h * h)
                r2 += r * r
        if np.sqrt(r2) < tol * b_norm:
            break
    # la constante additive est arbitraire : on la fixe comme la DCT (moyenne nulle)
    return phi - np.mean(phi), it + 1
