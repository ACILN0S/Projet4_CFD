"""Briques communes aux deux méthodes (M1 projection, M2 ψ–ω) — Projet 4, cavité entraînée.

Rôle : tout ce que M1 et M2 partagent, écrit une seule fois.
    - grille colocalisée aux nœuds et vitesse du couvercle U_t(t) ;
    - laplacien centré d'ordre 2 et advection (centrée ordre 2 ou upwind ordre 1) ;
    - solveur de Poisson par gradient conjugué SANS MATRICE (Dirichlet pour ψ, Neumann pour p) ;
    - pas de temps tiré des limites de stabilité du cours ;
    - diagnostics imposés par l'énoncé : K(t) et P(t).

Conventions (cf. docs/pseudocode.md) :
    - adimensionnement L = U = ρ = 1, donc ν = 1/Re ;
    - nœuds (x_i, y_j) = (i h, j h), i, j = 0..N, h = 1/N ; les nœuds de bord sont SUR les parois ;
    - tableaux numpy de forme (N+1, N+1), indexés [i, j] ↔ (x, y) ;
    - « intérieur » = nœuds i, j = 1..N-1, tableau de forme (N-1, N-1) ; l'élément [i-1, j-1] d'un
      tableau intérieur correspond au nœud (i, j).

Numba (recommandé par le prof) : les noyaux coûteux sont décorés par @njit, qui les compile en code
machine à leur premier appel (cache=True : la compilation est gardée sur disque d'une session à l'autre).
On peut alors écrire des boucles explicites sur i et j, qui reprennent mot pour mot les formules du
cours, sans le surcoût de numpy à chaque opération. Pas de fastmath : on garde l'arithmétique exacte
IEEE, donc les mêmes résultats (aux arrondis près) qu'une version numpy. Les fonctions publiques
(sans tiret bas) vérifient leurs arguments en Python pur, puis appellent le noyau compilé.
"""
import numpy as np
from numba import njit


# ---------------------------------------------------------------------------------------------
# Grille et condition au couvercle
# ---------------------------------------------------------------------------------------------

def make_grid(n_cells):
    """Grille colocalisée du carré unité.

    Entrée : n_cells = N, nombre d'intervalles par direction.
    Sorties : x, y de forme (N+1,) et h = 1/N.
    Comme aux workshops 8 et 11 : N+1 points, extrémités comprises (les parois).
    """
    h = 1.0 / n_cells
    x = np.linspace(0.0, 1.0, n_cells + 1)    # x_i = i h, i = 0..N
    y = np.linspace(0.0, 1.0, n_cells + 1)    # y_j = j h, j = 0..N
    return x, y, h


def lid_velocity(t):
    """Vitesse du couvercle U_t(t) = 1 - exp(-t) (énoncé, éq. 2, avec τ_lid = L/U = 1)."""
    return 1.0 - np.exp(-t)


# ---------------------------------------------------------------------------------------------
# Opérateurs aux différences finies (nœuds intérieurs uniquement)
# ---------------------------------------------------------------------------------------------

@njit(cache=True)
def laplacian(f, h):
    """Laplacien centré d'ordre 2, 5 points (cours §2.3 et §4.2.1).

        lap(f)[i,j] = (f[i+1,j] + f[i-1,j] + f[i,j+1] + f[i,j-1] - 4 f[i,j]) / h²

    Entrée : f de forme (N+1, N+1). Sortie : forme (N-1, N-1), nœuds intérieurs
    (l'élément [i-1, j-1] de la sortie est le nœud (i, j)).
    """
    n = f.shape[0] - 1
    out = np.empty((n - 1, n - 1))
    inv_h2 = 1.0 / (h * h)
    for i in range(1, n):
        for j in range(1, n):
            out[i - 1, j - 1] = (f[i + 1, j] + f[i - 1, j] + f[i, j + 1] + f[i, j - 1]
                                 - 4.0 * f[i, j]) * inv_h2
    return out


@njit(cache=True)
def _advection_kernel(a, b, f, h, upwind):
    """Noyau compilé de advection() ; upwind = True (upwind ordre 1) ou False (centré ordre 2)."""
    n = f.shape[0] - 1
    out = np.empty((n - 1, n - 1))
    for i in range(1, n):
        for j in range(1, n):
            if upwind:
                # décentré amont : on regarde d'où vient le fluide (signe de la vitesse)
                if a[i, j] > 0.0:
                    dfdx = (f[i, j] - f[i - 1, j]) / h          # différence arrière (éq. 2.6)
                else:
                    dfdx = (f[i + 1, j] - f[i, j]) / h          # différence avant (éq. 2.4)
                if b[i, j] > 0.0:
                    dfdy = (f[i, j] - f[i, j - 1]) / h
                else:
                    dfdy = (f[i, j + 1] - f[i, j]) / h
            else:
                # centré d'ordre 2 (éq. 2.9)
                dfdx = (f[i + 1, j] - f[i - 1, j]) / (2.0 * h)
                dfdy = (f[i, j + 1] - f[i, j - 1]) / (2.0 * h)
            out[i - 1, j - 1] = a[i, j] * dfdx + b[i, j] * dfdy
    return out


def advection(a, b, f, h, scheme):
    """Terme d'advection a ∂f/∂x + b ∂f/∂y aux nœuds intérieurs.

    Entrées : a, b (champs advectants) et f (champ advecté), formes (N+1, N+1) ;
              scheme = "centered" ou "upwind".
    Sortie : forme (N-1, N-1).

    - "centered" : différence centrée d'ordre 2 (cours §2.2, éq. 2.9)
          ∂f/∂x ≈ (f[i+1,j] - f[i-1,j]) / 2h
    - "upwind"   : décentré amont d'ordre 1 (cours §2.2, éq. 2.4/2.6 ; §5.6.2)
          ∂f/∂x ≈ (f[i,j] - f[i-1,j]) / h  si a > 0   (on regarde d'où vient le fluide)
          ∂f/∂x ≈ (f[i+1,j] - f[i,j]) / h  si a < 0
      même chose en y avec le signe de b.
    C'est LE choix « sage / malheureux » de l'énoncé (équation modifiée, cours §5.4.2).
    """
    if scheme not in ("centered", "upwind"):
        raise ValueError(f"schéma d'advection inconnu : {scheme!r} (attendu 'centered' ou 'upwind')")
    return _advection_kernel(a, b, f, h, scheme == "upwind")


# ---------------------------------------------------------------------------------------------
# Poisson Δφ = rhs par gradient conjugué sans matrice (cours §4.7.3, workshop 11.3)
# ---------------------------------------------------------------------------------------------

@njit(cache=True)
def _fill_boundary_kernel(phi, neumann):
    """Noyau compilé de fill_boundary() ; neumann = True (p) ou False (Dirichlet, ψ)."""
    n = phi.shape[0] - 1
    if neumann:
        # bord = voisin intérieur ; lignes i = 0, N puis colonnes j = 0, N (les coins reçoivent
        # ainsi la valeur du nœud intérieur en diagonale ; ils ne servent à aucun stencil)
        for j in range(n + 1):
            phi[0, j] = phi[1, j]          # paroi gauche  x = 0
            phi[n, j] = phi[n - 1, j]      # paroi droite  x = 1
        for i in range(n + 1):
            phi[i, 0] = phi[i, 1]          # paroi du bas  y = 0
            phi[i, n] = phi[i, n - 1]      # couvercle     y = 1
    else:
        for k in range(n + 1):
            phi[0, k] = 0.0
            phi[n, k] = 0.0
            phi[k, 0] = 0.0
            phi[k, n] = 0.0


def fill_boundary(phi, bc):
    """Remplit (en place) les nœuds de bord de phi selon la condition aux limites.

    - "dirichlet" (ψ) : φ = 0 sur les parois (cours §4.2.2.1).
    - "neumann"   (p) : ∂φ/∂n = 0 par la formule décentrée d'ordre 1 (cours §4.2.2.2, éq. 4.27-4.28,
                        avec g = 0) : φ_0 = φ_1, c'est-à-dire « bord = voisin intérieur ».
                        Ordre 1 : voir docs/adr/0002-neumann-pression.md.
    """
    if bc not in ("dirichlet", "neumann"):
        raise ValueError(f"condition aux limites inconnue : {bc!r} (attendu 'dirichlet' ou 'neumann')")
    _fill_boundary_kernel(phi, bc == "neumann")
    return phi


@njit(cache=True)
def _apply_minus_laplacian(d, h, neumann, out):
    """Produit A·d sans stocker A, avec A = -Δ_h (laplacien 5 points changé de signe).

    Entrée : d de forme (N+1, N+1) ; seul son intérieur compte, son bord est réécrit ici.
    Sortie : écrite dans out, forme (N-1, N-1).

    Pourquoi -Δ : c'est -A qui est symétrique définie positive dans le problème de Poisson
    (cours §4.7.3, note 2), et le gradient conjugué exige une matrice SDP.
    En Neumann, la règle « bord = voisin » garde A symétrique (le coefficient diagonal passe de 4 à 3
    le long d'une paroi, à 2 dans un coin) mais A·constante = 0 : A est seulement semi-définie positive
    et φ n'est défini qu'à une constante près.
    """
    _fill_boundary_kernel(d, neumann)
    n = d.shape[0] - 1
    inv_h2 = 1.0 / (h * h)
    for i in range(1, n):
        for j in range(1, n):
            out[i - 1, j - 1] = -(d[i + 1, j] + d[i - 1, j] + d[i, j + 1] + d[i, j - 1]
                                  - 4.0 * d[i, j]) * inv_h2


@njit(cache=True)
def _remove_interior_mean(phi):
    """φ -= moyenne(φ) sur les nœuds intérieurs (normalisation Neumann, slides du prof)."""
    n = phi.shape[0] - 1
    total = 0.0
    for i in range(1, n):
        for j in range(1, n):
            total += phi[i, j]
    mean = total / ((n - 1) * (n - 1))
    for i in range(1, n):
        for j in range(1, n):
            phi[i, j] -= mean


@njit(cache=True)
def _conjugate_gradient(b, phi, h, neumann, tol, max_iter):
    """Noyau compilé du gradient conjugué (cours §4.7.3), sur A φ = b avec A = -Δ_h.

    Entrées : b de forme (N-1, N-1) (moyenne déjà retirée en Neumann) ; phi de forme (N+1, N+1),
              estimation initiale, MODIFIÉE en place (elle devient la solution).
    Sorties : (nombre d'itérations, résidu relatif final). Itérations = -1 si max_iter est atteint.
    """
    n = phi.shape[0] - 1
    if neumann:
        _remove_interior_mean(phi)                 # même normalisation que la solution cherchée

    # Résidu initial r0 = b - A φ0 et direction d0 = r0 (d stockée en (N+1, N+1) pour son bord)
    r = np.empty_like(b)
    ad = np.empty_like(b)
    d = np.zeros_like(phi)
    _apply_minus_laplacian(phi, h, neumann, ad)
    rr = 0.0
    b_norm2 = 0.0
    for i in range(n - 1):
        for j in range(n - 1):
            r[i, j] = b[i, j] - ad[i, j]
            d[i + 1, j + 1] = r[i, j]
            rr += r[i, j] * r[i, j]
            b_norm2 += b[i, j] * b[i, j]

    b_norm = np.sqrt(b_norm2)
    if b_norm == 0.0:                              # second membre nul : la solution est φ = 0
        phi[:, :] = 0.0
        return 0, 0.0
    target = tol * b_norm

    iterations = 0
    while np.sqrt(rr) > target:
        if iterations >= max_iter:
            return -1, np.sqrt(rr) / b_norm
        _apply_minus_laplacian(d, h, neumann, ad)  # produit matrice-vecteur sans matrice
        d_ad = 0.0
        for i in range(n - 1):
            for j in range(n - 1):
                d_ad += d[i + 1, j + 1] * ad[i, j]
        alpha = rr / d_ad                          # pas optimal le long de d
        rr_new = 0.0
        for i in range(n - 1):
            for j in range(n - 1):
                phi[i + 1, j + 1] += alpha * d[i + 1, j + 1]   # mise à jour de la solution
                r[i, j] -= alpha * ad[i, j]                    # mise à jour du résidu
                rr_new += r[i, j] * r[i, j]
        beta = rr_new / rr                         # nouvelle direction A-conjuguée aux précédentes
        for i in range(n - 1):
            for j in range(n - 1):
                d[i + 1, j + 1] = r[i, j] + beta * d[i + 1, j + 1]
        rr = rr_new
        if neumann:
            _remove_interior_mean(phi)             # moyenne nulle à chaque itération (slides)
        iterations += 1

    _fill_boundary_kernel(phi, neumann)            # bord de la solution : même règle que dans A
    return iterations, np.sqrt(rr) / b_norm


def solve_poisson(rhs, h, bc, phi0=None, tol=1e-10, max_iter=None):
    """Résout Δφ = rhs par gradient conjugué sans matrice (cours §4.7.3, workshop 11.3).

    Entrées :
        rhs  : second membre aux nœuds intérieurs, forme (N-1, N-1) ;
        h    : pas d'espace ;
        bc   : "dirichlet" (φ = 0 aux parois, pour ψ) ou "neumann" (∂φ/∂n = 0, pour p) ;
        phi0 : estimation initiale de forme (N+1, N+1), typiquement la solution du pas de temps
               précédent (départ « à chaud ») ; None ⇒ départ de zéro ; phi0 n'est pas modifié ;
        tol  : arrêt quand ||r|| ≤ tol·||b|| (résidu relatif, cours §4.3.2, cf. docs/adr/0003) ;
        max_iter : garde-fou ; None ⇒ 10·(N-1)² (bien au-delà de la théorie, cf. ci-dessous).
    Sorties : phi de forme (N+1, N+1) (bord rempli selon bc) et le nombre d'itérations.

    On résout le système équivalent A φ = b avec A = -Δ_h et b = -rhs (A SDP, cf. _apply_minus_laplacian).
    Algorithme du cours §4.7.3 (sans préconditionnement), noté ici avec r (résidu) et d (direction) :
        r0 = b - A x0, d0 = r0
        α  = (r·r)/(d·A d) ; x += α d ; r_new = r - α A d
        β  = (r_new·r_new)/(r·r) ; d = r_new + β d
    En exactes, le gradient conjugué converge en au plus (N-1)² itérations ; en pratique en O(N)
    itérations pour le laplacien (complexité O(n^{3/2}), cours §4.8, Tab. 4.1).

    Spécificités Neumann (slides du prof : existence et unicité) :
        - existence : la somme du second membre doit être nulle (A symétrique et A·1 = 0 ⇒ b ⊥ 1),
          on retire donc sa moyenne : b -= mean(b) ;
        - unicité : φ est défini à une constante près, on impose une moyenne nulle à chaque itération,
          φ -= mean(φ) ; cela ne change pas le résidu puisque A·constante = 0.
    """
    if bc not in ("dirichlet", "neumann"):
        raise ValueError(f"condition aux limites inconnue : {bc!r} (attendu 'dirichlet' ou 'neumann')")
    n_int = rhs.shape[0]                     # N - 1 nœuds intérieurs par direction
    if rhs.shape != (n_int, n_int):
        raise ValueError(f"rhs doit être carré (N-1, N-1), reçu {rhs.shape}")
    if max_iter is None:
        max_iter = 10 * n_int**2

    # Estimation initiale (copie : phi0 n'est pas modifié) et second membre b = -rhs
    phi = np.zeros((n_int + 2, n_int + 2)) if phi0 is None else np.array(phi0, dtype=float)
    b = -np.asarray(rhs, dtype=float)
    if bc == "neumann":
        b = b - b.mean()                     # condition de compatibilité (existence)

    iterations, rel_res = _conjugate_gradient(b, phi, h, bc == "neumann", tol, max_iter)
    if iterations < 0:
        raise RuntimeError(f"gradient conjugué non convergé en {max_iter} itérations "
                           f"(||r||/||b|| = {rel_res:.3e}, tol = {tol:.1e}, bc = {bc})")
    return phi, iterations


# ---------------------------------------------------------------------------------------------
# Pas de temps (cours §5.3)
# ---------------------------------------------------------------------------------------------

def time_step(h, nu, scheme, safety=0.8):
    """Pas de temps d'Euler explicite donné par les limites de stabilité du cours, marge 0.8.

    - Diffusion, Euler + centré 2D : Fo = ν Δt / h² ≤ 1/4            (cours §5.3.3, éq. 5.127)
    - Advection upwind 2D          : C_x + C_y ≤ 1, avec |u|, |v| ≤ U = 1, donc Δt ≤ h/2
                                                                    (cours §5.3.3, éq. 5.131)
    - Advection centrée + Euler    : instable seule ; stable grâce à la diffusion si C² ≤ 2 Fo ≤ 1
                                     (cours §5.3.1.4, advection-diffusion). Avec C = U Δt / h et
                                     Fo = ν Δt / h², cela donne Δt ≤ 2ν / U² = 2ν.
    Remarque : la limite C² ≤ 2 Fo est établie en 1D dans le cours ; on l'applique telle quelle.
    Laquelle des limites est active selon N et Re : docs/adr/0003-pas-de-temps-et-arrets.md.
    """
    dt = safety * min(h**2 / (4.0 * nu), h / 2.0)
    if scheme == "centered":
        dt = min(dt, safety * 2.0 * nu)
    elif scheme != "upwind":
        raise ValueError(f"schéma d'advection inconnu : {scheme!r}")
    return dt


# ---------------------------------------------------------------------------------------------
# Diagnostics imposés : K(t) et P(t)  (appelés une fois par pas : numpy suffit)
# ---------------------------------------------------------------------------------------------

def kinetic_energy(u, v, h):
    """K = ½ ∫∫ (u² + v²) dS par la méthode des trapèzes en 2D (workshop 8.2).

    Trapèzes appliqués en x puis en y : poids h aux nœuds intérieurs, h/2 sur les parois
    (donc h²/4 aux coins). Seul le couvercle a une vitesse de paroi non nulle.
    Entrées : u, v de forme (N+1, N+1). Sortie : un réel.
    """
    e = 0.5 * (u**2 + v**2)
    return np.trapezoid(np.trapezoid(e, dx=h, axis=1), dx=h)


def lid_power(u, h, nu, u_lid):
    """P = ∫₀¹ ν ∂u/∂y|_{y=1} U_t dx (puissance injectée par le couvercle, énoncé).

    - ∂u/∂y au couvercle : formule décentrée amont d'ordre 2 (workshop 11.2, éq. 11.10 ;
      cours §2.2, fin : « forward or backward formula with a larger stencil »)
          dudy[i] = (3 u[i,N] - 4 u[i,N-1] + u[i,N-2]) / 2h
    - intégrale en x par les trapèzes (workshop 8.2), poids ½ aux deux coins.
    Aux coins, u = 0 (nœud des parois fixes) : dudy y vaut 0, ils ne contribuent pas.
    Entrées : u de forme (N+1, N+1), h, ν, U_t. Sortie : un réel.
    """
    dudy = (3.0 * u[:, -1] - 4.0 * u[:, -2] + u[:, -3]) / (2.0 * h)
    return np.trapezoid(nu * dudy * u_lid, dx=h)
