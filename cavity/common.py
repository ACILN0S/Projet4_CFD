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
import time

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


def apply_velocity_bc(u, v, u_lid):
    """Impose (en place) les vitesses de paroi (M1 et M2).

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
# Trois conditions aux limites (codes entiers pour les noyaux Numba) :
#   DIRICHLET     φ = 0 aux parois (ψ, M2) ; inconnues = nœuds intérieurs.
#   NEUMANN       ∂φ/∂n = 0 d'ordre 1, φ_0 = φ_1 (ancienne pression de M1, gardée pour la comparaison
#                 du notebook A §1.c) ; inconnues = nœuds intérieurs.
#   NEUMANN_GHOST ∂φ/∂n = 0 d'ordre 2 par nœuds fantômes, miroir pair φ_{-1} = φ_1 (pression de M1,
#                 docs/adr/0002 option 2) ; inconnues = TOUS les nœuds, parois comprises.
DIRICHLET, NEUMANN, NEUMANN_GHOST = 0, 1, 2
BC_CODES = {"dirichlet": DIRICHLET, "neumann": NEUMANN, "neumann_ghost": NEUMANN_GHOST}


@njit(cache=True)
def _fill_boundary_kernel(phi, mode):
    """Noyau compilé de fill_boundary() (modes DIRICHLET et NEUMANN ; rien à faire en NEUMANN_GHOST,
    où les nœuds de paroi sont des inconnues)."""
    n = phi.shape[0] - 1
    if mode == NEUMANN:
        # bord = voisin intérieur ; lignes i = 0, N puis colonnes j = 0, N (les coins reçoivent
        # ainsi la valeur du nœud intérieur en diagonale ; ils ne servent à aucun stencil)
        for j in range(n + 1):
            phi[0, j] = phi[1, j]          # paroi gauche  x = 0
            phi[n, j] = phi[n - 1, j]      # paroi droite  x = 1
        for i in range(n + 1):
            phi[i, 0] = phi[i, 1]          # paroi du bas  y = 0
            phi[i, n] = phi[i, n - 1]      # couvercle     y = 1
    elif mode == DIRICHLET:
        for k in range(n + 1):
            phi[0, k] = 0.0
            phi[n, k] = 0.0
            phi[k, 0] = 0.0
            phi[k, n] = 0.0


def fill_boundary(phi, bc):
    """Remplit (en place) les nœuds de bord de phi selon la condition aux limites.

    - "dirichlet" (ψ) : φ = 0 sur les parois (cours §4.2.2.1).
    - "neumann"       : ∂φ/∂n = 0 par la formule décentrée d'ordre 1 (cours §4.2.2.2, éq. 4.27-4.28,
                        avec g = 0) : φ_0 = φ_1, « bord = voisin intérieur ». Ordre 1 : docs/adr/0002.
    - "neumann_ghost" : rien à remplir, les nœuds de paroi sont des inconnues du système.
    """
    if bc not in BC_CODES:
        raise ValueError(f"condition aux limites inconnue : {bc!r} (attendu {list(BC_CODES)})")
    _fill_boundary_kernel(phi, BC_CODES[bc])
    return phi


@njit(cache=True)
def _trapezoid_weight(k, n):
    """Poids des trapèzes en 1D : ½ aux deux extrémités (parois), 1 ailleurs (workshop 8.2)."""
    return 0.5 if (k == 0 or k == n) else 1.0


@njit(cache=True)
def _apply_matrix(d, h, mode, out):
    """Produit A·d sans stocker A (tableaux (N+1, N+1) ; out n'est écrit que sur les nœuds actifs).

    - DIRICHLET, NEUMANN : A = -Δ_h aux nœuds intérieurs, après avoir rempli le bord de d.
      Pourquoi -Δ : c'est -A qui est symétrique définie positive (cours §4.7.3, note 2), et le gradient
      conjugué exige une matrice SDP. En NEUMANN, la règle « bord = voisin » garde A symétrique
      (coefficient diagonal 3 le long d'une paroi, 2 dans un coin), mais A·constante = 0.
    - NEUMANN_GHOST : à TOUS les nœuds, laplacien 5 points avec fantômes miroirs φ_{-1} = φ_1,
      φ_{N+1} = φ_{N-1} (dérivée normale centrée d'ordre 2 nulle, cours §2.2 éq. 2.9). Ligne de paroi
      (gauche) : (2 φ_1 + φ_{0,j+1} + φ_{0,j-1} - 4 φ_0) / h² ; le facteur 2 casse la symétrie. On la rétablit
      en multipliant chaque ligne par le poids des trapèzes w_i w_j (½ sur une paroi, ¼ dans un coin) :
      A = -W Δ_h est symétrique (ex. : coefficient paroi→intérieur ½·2 = 1 = coefficient intérieur→paroi).
      A·constante = 0 aussi : φ défini à une constante près.
    """
    n = d.shape[0] - 1
    inv_h2 = 1.0 / (h * h)
    if mode == NEUMANN_GHOST:
        for i in range(n + 1):
            im = 1 if i == 0 else i - 1            # fantôme miroir à gauche : φ_{-1} = φ_1
            ip = n - 1 if i == n else i + 1        # fantôme miroir à droite : φ_{N+1} = φ_{N-1}
            for j in range(n + 1):
                jm = 1 if j == 0 else j - 1
                jp = n - 1 if j == n else j + 1
                lap = (d[ip, j] + d[im, j] + d[i, jp] + d[i, jm] - 4.0 * d[i, j]) * inv_h2
                out[i, j] = -_trapezoid_weight(i, n) * _trapezoid_weight(j, n) * lap
    else:
        _fill_boundary_kernel(d, mode)
        for i in range(1, n):
            for j in range(1, n):
                out[i, j] = -(d[i + 1, j] + d[i - 1, j] + d[i, j + 1] + d[i, j - 1]
                              - 4.0 * d[i, j]) * inv_h2


@njit(cache=True)
def _remove_mean(phi, mode):
    """φ -= moyenne(φ) (normalisation Neumann, slides du prof).

    NEUMANN : moyenne simple sur les nœuds intérieurs (les inconnues).
    NEUMANN_GHOST : moyenne pondérée par les poids des trapèzes sur tous les nœuds, c'est-à-dire
                    ∫∫ φ dS = 0 (workshop 8.2).
    """
    n = phi.shape[0] - 1
    lo, hi = (0, n + 1) if mode == NEUMANN_GHOST else (1, n)
    total = 0.0
    weight = 0.0
    for i in range(lo, hi):
        for j in range(lo, hi):
            w = _trapezoid_weight(i, n) * _trapezoid_weight(j, n) if mode == NEUMANN_GHOST else 1.0
            total += w * phi[i, j]
            weight += w
    mean = total / weight
    for i in range(lo, hi):
        for j in range(lo, hi):
            phi[i, j] -= mean


@njit(cache=True)
def _conjugate_gradient(b, phi, h, mode, tol, max_iter):
    """Noyau compilé du gradient conjugué (cours §4.7.3), sur A φ = b (A : cf. _apply_matrix).

    Entrées : b de forme (N+1, N+1), nul hors des nœuds actifs (compatibilité déjà assurée en Neumann) ;
              phi de forme (N+1, N+1), estimation initiale, MODIFIÉE en place (elle devient la solution).
    Nœuds actifs (les inconnues) : intérieurs en DIRICHLET et NEUMANN, tous en NEUMANN_GHOST.
    Sorties : (nombre d'itérations, résidu relatif final). Itérations = -1 si max_iter est atteint.
    """
    n = phi.shape[0] - 1
    lo, hi = (0, n + 1) if mode == NEUMANN_GHOST else (1, n)
    neumann = mode != DIRICHLET
    if neumann:
        _remove_mean(phi, mode)                    # même normalisation que la solution cherchée

    # Résidu initial r0 = b - A φ0 et direction d0 = r0
    r = np.zeros_like(phi)
    ad = np.zeros_like(phi)
    d = np.zeros_like(phi)
    _apply_matrix(phi, h, mode, ad)
    rr = 0.0
    b_norm2 = 0.0
    for i in range(lo, hi):
        for j in range(lo, hi):
            r[i, j] = b[i, j] - ad[i, j]
            d[i, j] = r[i, j]
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
        _apply_matrix(d, h, mode, ad)              # produit matrice-vecteur sans matrice
        d_ad = 0.0
        for i in range(lo, hi):
            for j in range(lo, hi):
                d_ad += d[i, j] * ad[i, j]
        alpha = rr / d_ad                          # pas optimal le long de d
        rr_new = 0.0
        for i in range(lo, hi):
            for j in range(lo, hi):
                phi[i, j] += alpha * d[i, j]       # mise à jour de la solution
                r[i, j] -= alpha * ad[i, j]        # mise à jour du résidu
                rr_new += r[i, j] * r[i, j]
        beta = rr_new / rr                         # nouvelle direction A-conjuguée aux précédentes
        for i in range(lo, hi):
            for j in range(lo, hi):
                d[i, j] = r[i, j] + beta * d[i, j]
        rr = rr_new
        if neumann:
            _remove_mean(phi, mode)                # moyenne nulle à chaque itération (slides)
        iterations += 1

    _fill_boundary_kernel(phi, mode)               # bord de la solution : même règle que dans A
    return iterations, np.sqrt(rr) / b_norm


def solve_poisson(rhs, h, bc, phi0=None, tol=1e-10, max_iter=None):
    """Résout Δφ = rhs par gradient conjugué sans matrice (cours §4.7.3, workshop 11.3).

    Entrées :
        rhs  : second membre ; forme (N-1, N-1) (nœuds intérieurs) pour "dirichlet" et "neumann",
               forme (N+1, N+1) (tous les nœuds, parois comprises) pour "neumann_ghost" ;
        h    : pas d'espace ;
        bc   : "dirichlet" (φ = 0 aux parois, pour ψ), "neumann" (∂φ/∂n = 0 d'ordre 1) ou
               "neumann_ghost" (∂φ/∂n = 0 d'ordre 2 par nœuds fantômes, pour p ; docs/adr/0002) ;
        phi0 : estimation initiale de forme (N+1, N+1), typiquement la solution du pas de temps
               précédent (départ « à chaud ») ; None ⇒ départ de zéro ; phi0 n'est pas modifié ;
        tol  : arrêt quand ||r|| ≤ tol·||b|| (résidu relatif, cours §4.3.2, cf. docs/adr/0003) ;
        max_iter : garde-fou ; None ⇒ 10 × (nombre d'inconnues).
    Sorties : phi de forme (N+1, N+1) et le nombre d'itérations.

    On résout le système équivalent A φ = b, avec A = -Δ_h et b = -rhs (Dirichlet, Neumann), ou
    A = -W Δ_h et b = -W rhs (fantômes, W = poids des trapèzes), A symétrique (cf. _apply_matrix).
    Algorithme du cours §4.7.3 (sans préconditionnement), noté ici avec r (résidu) et d (direction) :
        r0 = b - A x0, d0 = r0
        α  = (r·r)/(d·A d) ; x += α d ; r_new = r - α A d
        β  = (r_new·r_new)/(r·r) ; d = r_new + β d
    En pratique O(N) itérations pour le laplacien (complexité O(n^{3/2}), cours §4.8, Tab. 4.1).

    Spécificités Neumann (slides du prof : existence et unicité) :
        - existence : A symétrique et A·1 = 0 ⇒ b doit être orthogonal aux constantes. En "neumann",
          on retire la moyenne de rhs ; en "neumann_ghost", sa moyenne PONDÉRÉE par les trapèzes, soit
          ∫∫ rhs dS = 0, la condition de compatibilité du problème continu (∫ Δφ = ∮ ∂φ/∂n = 0) ;
        - unicité : φ est défini à une constante près, on impose une moyenne nulle (pondérée pour les
          fantômes) à chaque itération ; cela ne change pas le résidu puisque A·constante = 0.
    """
    if bc not in BC_CODES:
        raise ValueError(f"condition aux limites inconnue : {bc!r} (attendu {list(BC_CODES)})")
    mode = BC_CODES[bc]
    rhs = np.asarray(rhs, dtype=float)
    n_nodes = rhs.shape[0] if mode == NEUMANN_GHOST else rhs.shape[0] + 2     # N + 1
    expected = (n_nodes, n_nodes) if mode == NEUMANN_GHOST else (n_nodes - 2, n_nodes - 2)
    if rhs.shape != expected:
        raise ValueError(f"rhs de forme {rhs.shape} incompatible avec bc = {bc!r} (attendu {expected})")
    n_unknowns = rhs.size
    if max_iter is None:
        max_iter = 10 * n_unknowns

    # Estimation initiale (copie : phi0 n'est pas modifié)
    phi = np.zeros((n_nodes, n_nodes)) if phi0 is None else np.array(phi0, dtype=float)

    # Second membre b du système A φ = b, nul hors des nœuds actifs
    b = np.zeros((n_nodes, n_nodes))
    if mode == NEUMANN_GHOST:
        w1d = np.ones(n_nodes)
        w1d[0] = w1d[-1] = 0.5                        # poids des trapèzes (workshop 8.2)
        weights = np.outer(w1d, w1d)
        rhs = rhs - np.sum(weights * rhs) / np.sum(weights)   # ∫∫ rhs dS = 0 (compatibilité)
        b[:, :] = -weights * rhs
    else:
        if mode == NEUMANN:
            rhs = rhs - rhs.mean()                    # somme nulle (compatibilité)
        b[1:-1, 1:-1] = -rhs

    iterations, rel_res = _conjugate_gradient(b, phi, h, mode, tol, max_iter)
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
# Critère d'arrêt : résidu stationnaire (docs/adr/0003-pas-de-temps-et-arrets.md §2)
# ---------------------------------------------------------------------------------------------

def steady_residual(u_new, u_old, v_new, v_old, dt):
    """Norme L2 du résidu stationnaire R = (u^{n+1} - u^n) / Δt, composantes u et v, nœuds intérieurs.

    Pourquoi c'est un résidu : avec Euler explicite, (u^{n+1} - u^n)/Δt vaut exactement le membre de
    droite discret (-advection + ν Δu - ∇p) évalué sur le champ courant, c'est-à-dire de combien ce
    champ viole les équations STATIONNAIRES discrètes. R → 0 au stationnaire (même idée que le résidu
    r = b - A x des méthodes itératives, cours §4.3.2, et que la question b du Projet 1).
    Contrairement à max|u^{n+1} - u^n|, R ne dépend pas de Δt.
    Norme L2 comme au Projet 1 : ||R|| = sqrt( Σ (R_u² + R_v²) ) sur les nœuds intérieurs (on exclut les
    parois, où les vitesses sont imposées et suivent la rampe du couvercle).
    """
    r_u = (u_new[1:-1, 1:-1] - u_old[1:-1, 1:-1]) / dt
    r_v = (v_new[1:-1, 1:-1] - v_old[1:-1, 1:-1]) / dt
    return np.sqrt(np.sum(r_u**2) + np.sum(r_v**2))


# ---------------------------------------------------------------------------------------------
# Boucle en temps commune à M1 et M2 (Euler explicite, cours §3.1)
# ---------------------------------------------------------------------------------------------

def march_to_steady(step, state, h, dt, nu, residual_tol, t_ref, t_max, label):
    """Avance un état jusqu'au stationnaire (ou jusqu'à t_max), avec la même règle pour M1 et M2.

    Entrées :
        step  : fonction step(state, t) -> (nouvel état, itérations du gradient conjugué), qui avance
                d'un pas Δt ; l'état est un dictionnaire contenant au moins "u" et "v" (N+1, N+1) ;
        state : état initial (fluide au repos) ; h, Δt, ν ;
        residual_tol, t_ref : critère d'arrêt ||R|| < ε ||R(t_ref)|| sur le résidu stationnaire
                (steady_residual ; docs/adr/0003 §2), surveillé à partir de t_ref ;
        t_max : garde-fou ; si atteint, converged = False ; label : nom pour les messages d'erreur.
    Sortie : (état final, dictionnaire d'historique et d'informations d'exécution).
    Diagnostics K et P (énoncé) enregistrés à chaque pas : leur coût est négligeable devant le Poisson.
    """
    hist_t, hist_k, hist_p, hist_iter, hist_res = [], [], [], [], []
    t = 0.0
    n_steps = 0
    residual_ref = None                           # ||R(t_ref)||, mesuré au premier pas où t ≥ t_ref
    converged = False
    wall_start = time.perf_counter()

    while t < t_max:
        new_state, n_iter = step(state, t)
        t += dt
        n_steps += 1
        # Résidu stationnaire ||R|| = ||(u^{n+1} - u^n)/Δt|| (cours §4.3.2, Projet 1 question b)
        residual = steady_residual(new_state["u"], state["u"], new_state["v"], state["v"], dt)
        state = new_state

        hist_t.append(t)
        hist_k.append(kinetic_energy(state["u"], state["v"], h))
        hist_p.append(lid_power(state["u"], h, nu, lid_velocity(t)))
        hist_iter.append(n_iter)
        hist_res.append(residual)

        if not np.isfinite(residual):
            raise FloatingPointError(f"{label} diverge à t = {t:.4f} (Δt = {dt:.3e})")
        # Critère d'arrêt : résidu réduit d'un facteur ε par rapport à sa valeur en t_ref
        if t >= t_ref:
            if residual_ref is None:
                residual_ref = residual
            elif residual < residual_tol * residual_ref:
                converged = True
                break

    info = {
        "t": np.array(hist_t), "K": np.array(hist_k), "P": np.array(hist_p),
        "cg_iterations": np.array(hist_iter), "residual": np.array(hist_res),
        "n_steps": n_steps, "t_final": t, "wall_time": time.perf_counter() - wall_start,
        "residual_ref": residual_ref, "converged": converged,
        "residual_tol": residual_tol, "t_ref": t_ref,
    }
    return state, info


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
