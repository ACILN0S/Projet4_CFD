"""Méthode 1 (« cours ») : variables primitives (u, v), grille décalée MAC, projection de Chorin.

Pseudo-code de référence : docs/pseudocode.md §1. Un pas de temps =
    RK3-SSP sur  ∂u/∂t = −A(u) + ν Δu,  avec une PROJECTION après chaque étage
    (projection = on retire le gradient ∇φ qui rend le champ à divergence nulle).

Grille MAC, N × N cellules de côté h = 1/N, convention [i, j] ↔ (x, y)
--------------------------------------------------------------------------
    p[i, j]   au centre   ((i + ½) h, (j + ½) h)        i, j  = 0..N−1        forme (N, N)
    u[i, jj]  face vert.  (i h, (jj − ½) h)             i = 0..N, jj = 0..N+1  forme (N+1, N+2)
    v[ii, j]  face horiz. ((ii − ½) h, j h)             ii = 0..N+1, j = 0..N  forme (N+2, N+1)

Les indices `jj` (pour u) et `ii` (pour v) sont DÉCALÉS DE +1 à cause des fantômes :
    u[:, 0] et u[:, N+1]  = rangées fantômes sous la paroi du bas / au-dessus du couvercle,
    v[0, :] et v[N+1, :]  = colonnes fantômes à gauche / à droite.
Valeurs « inconnues » (intérieures) : u[1..N−1, 1..N] et v[1..N, 1..N−1].
Les faces de paroi u[0, :], u[N, :], v[:, 0], v[:, N] valent toujours 0 (vitesse normale nulle).

Advection, forme conservative (flux) : A_u = ∂(uu)/∂x + ∂(vu)/∂y, A_v = ∂(uv)/∂x + ∂(vv)/∂y.
Avec le schéma centré, cette forme conserve exactement l'énergie cinétique discrète : c'est ce qui
rend le bilan dK/dt = P − ε exact (aux erreurs en temps près). L'upwind 1 ajoute une dissipation
numérique ν_num ≈ |U| h / 2 : le bilan ne ferme plus (contre-exemple volontaire, Partie C).
"""
import numba
import numpy as np

from cavity.params import lid_velocity
from cavity.poisson import solve_poisson_neumann


def new_fields(N):
    """Champs de vitesse au repos : u de forme (N+1, N+2), v de forme (N+2, N+1) (fantômes compris)."""
    return np.zeros((N + 1, N + 2)), np.zeros((N + 2, N + 1))


@numba.njit(cache=True)
def apply_bc(u, v, U_lid):
    """Conditions aux limites de la cavité, écrites EN PLACE dans u et v.

    - Vitesse NORMALE à la paroi : posée directement sur la face de paroi (= 0, parois imperméables).
    - Vitesse TANGENTIELLE : la paroi passe au milieu entre la face fantôme et la première face
      intérieure, donc on impose la MOYENNE : (u_fant + u_int) / 2 = u_paroi
      ⇒ u_fant = 2 u_paroi − u_int  (ordre 2 au point de paroi).
    Paroi du haut (couvercle) : u_paroi = U_lid = U_t(t) ; les trois autres : 0.
    """
    nu_x, nu_y = u.shape          # (N+1, N+2)
    nv_x, nv_y = v.shape          # (N+2, N+1)
    # u normal aux parois gauche (i = 0) et droite (i = N)
    for jj in range(nu_y):
        u[0, jj] = 0.0
        u[nu_x - 1, jj] = 0.0
    # u tangentiel en bas (u_paroi = 0) et en haut (u_paroi = U_lid) : fantômes
    for i in range(nu_x):
        u[i, 0] = -u[i, 1]
        u[i, nu_y - 1] = 2.0 * U_lid - u[i, nu_y - 2]
    # v normal aux parois du bas (j = 0) et du haut (j = N)
    for ii in range(nv_x):
        v[ii, 0] = 0.0
        v[ii, nv_y - 1] = 0.0
    # v tangentiel à gauche et à droite (v_paroi = 0) : fantômes
    for j in range(nv_y):
        v[0, j] = -v[1, j]
        v[nv_x - 1, j] = -v[nv_x - 2, j]


@numba.njit(cache=True)
def laplacian(u, v, h, Lu, Lv):
    """Laplacien 5 points (ordre 2) aux points intérieurs de u et de v, écrit dans Lu, Lv.

        Δu[i,jj] ≈ (u[i+1,jj] + u[i−1,jj] + u[i,jj+1] + u[i,jj−1] − 4 u[i,jj]) / h²   (idem pour v)

    Près des parois, le stencil utilise les fantômes : il faut donc avoir appelé `apply_bc` avant.
    Lu, Lv : mêmes formes que u, v ; seules les cases intérieures sont écrites.
    """
    N = u.shape[0] - 1
    for i in range(1, N):                 # u intérieurs : i = 1..N−1, jj = 1..N
        for jj in range(1, N + 1):
            Lu[i, jj] = (u[i + 1, jj] + u[i - 1, jj] + u[i, jj + 1] + u[i, jj - 1]
                         - 4.0 * u[i, jj]) / h**2
    for ii in range(1, N + 1):            # v intérieurs : ii = 1..N, j = 1..N−1
        for j in range(1, N):
            Lv[ii, j] = (v[ii + 1, j] + v[ii - 1, j] + v[ii, j + 1] + v[ii, j - 1]
                         - 4.0 * v[ii, j]) / h**2


@numba.njit(cache=True)
def face_value(q_minus, q_plus, a, upwind):
    """Valeur transportée q sur une face, entre les deux voisins q_minus (côté −) et q_plus (côté +).

    - centré  : moyenne (q_minus + q_plus) / 2                     → ordre 2, non dissipatif
    - upwind1 : valeur amont selon le signe de la vitesse a        → ordre 1, dissipatif
    """
    if upwind:
        return q_minus if a > 0.0 else q_plus
    return 0.5 * (q_minus + q_plus)


@numba.njit(cache=True)
def advection(u, v, h, upwind, Au, Av):
    """Terme d'advection conservatif aux points intérieurs de u et v, écrit dans Au, Av.

    Pour u au point (i h, (j + ½) h) :
        A_u = [F(i+½) − F(i−½)] / h + [G(j+1) − G(j)] / h
        F = flux de u en x, aux centres des cellules : F = ū · q,  ū = moyenne de 2 u voisins
        G = flux de u en y, aux nœuds (coins)        : G = v̄ · q,  v̄ = moyenne de 2 v voisins
    où q est la valeur de u transportée sur la face (`face_value` : centré ou upwind 1).
    Symétriquement pour v. Sur les parois, la vitesse transportante (normale) est nulle,
    donc les flux de paroi sont nuls : aucune quantité de mouvement ne traverse les parois.
    Il faut avoir appelé `apply_bc` avant (les flux près des parois lisent les fantômes).
    """
    N = u.shape[0] - 1
    # ---- équation de u : points i = 1..N−1, jj = 1..N ----
    for i in range(1, N):
        for jj in range(1, N + 1):
            # flux en x aux centres (i+½) et (i−½) : vitesse transportante = moyenne des u voisins
            a_e = 0.5 * (u[i, jj] + u[i + 1, jj])
            a_w = 0.5 * (u[i - 1, jj] + u[i, jj])
            F_e = a_e * face_value(u[i, jj], u[i + 1, jj], a_e, upwind)
            F_w = a_w * face_value(u[i - 1, jj], u[i, jj], a_w, upwind)
            # flux en y aux nœuds du haut et du bas : v aux nœuds = moyenne des 2 v voisins en x
            # (les v voisins du nœud (i h, y) sont v[i, ·] et v[i+1, ·] à cause du décalage +1)
            a_n = 0.5 * (v[i, jj] + v[i + 1, jj])
            a_s = 0.5 * (v[i, jj - 1] + v[i + 1, jj - 1])
            G_n = a_n * face_value(u[i, jj], u[i, jj + 1], a_n, upwind)
            G_s = a_s * face_value(u[i, jj - 1], u[i, jj], a_s, upwind)
            Au[i, jj] = (F_e - F_w) / h + (G_n - G_s) / h
    # ---- équation de v : points ii = 1..N, j = 1..N−1 ----
    for ii in range(1, N + 1):
        for j in range(1, N):
            # flux en y aux centres (j+½) et (j−½) : vitesse transportante = moyenne des v voisins
            a_n = 0.5 * (v[ii, j] + v[ii, j + 1])
            a_s = 0.5 * (v[ii, j - 1] + v[ii, j])
            G_n = a_n * face_value(v[ii, j], v[ii, j + 1], a_n, upwind)
            G_s = a_s * face_value(v[ii, j - 1], v[ii, j], a_s, upwind)
            # flux en x aux nœuds de droite et de gauche : u aux nœuds = moyenne des 2 u voisins en y
            a_e = 0.5 * (u[ii, j] + u[ii, j + 1])
            a_w = 0.5 * (u[ii - 1, j] + u[ii - 1, j + 1])
            F_e = a_e * face_value(v[ii, j], v[ii + 1, j], a_e, upwind)
            F_w = a_w * face_value(v[ii - 1, j], v[ii, j], a_w, upwind)
            Av[ii, j] = (F_e - F_w) / h + (G_n - G_s) / h


@numba.njit(cache=True)
def divergence(u, v, h, div):
    """Divergence discrète au centre de chaque cellule, écrite dans div (forme (N, N)).

        (∇·u)[i,j] = (u[i+1, j+1] − u[i, j+1]) / h + (v[i+1, j+1] − v[i+1, j]) / h
    (les « +1 » viennent du décalage des fantômes : u[·, jj = j+1] est dans la rangée de la cellule j).
    """
    N = div.shape[0]
    for i in range(N):
        for j in range(N):
            div[i, j] = (u[i + 1, j + 1] - u[i, j + 1]) / h + (v[i + 1, j + 1] - v[i + 1, j]) / h


@numba.njit(cache=True)
def subtract_gradient(u, v, phi, h):
    """u ← u − ∇φ sur les faces INTÉRIEURES (en place). Les faces de paroi ne sont jamais corrigées.

        u[i, jj] −= (φ[i, jj−1] − φ[i−1, jj−1]) / h     (différence entre les 2 cellules voisines en x)
        v[ii, j] −= (φ[ii−1, j] − φ[ii−1, j−1]) / h     (différence entre les 2 cellules voisines en y)
    """
    N = phi.shape[0]
    for i in range(1, N):
        for jj in range(1, N + 1):
            u[i, jj] -= (phi[i, jj - 1] - phi[i - 1, jj - 1]) / h
    for ii in range(1, N + 1):
        for j in range(1, N):
            v[ii, j] -= (phi[ii - 1, j] - phi[ii - 1, j - 1]) / h


def project(u, v, h):
    """Projection de Chorin : rend (u, v) à divergence nulle, en place.

    On résout Δφ = ∇·u* (Neumann : ∂φ/∂n = 0 car on ne touche pas aux faces de paroi), puis
    u = u* − ∇φ. Comme le laplacien discret est exactement « divergence du gradient » sur la
    grille MAC, on obtient ∇·u = 0 à l'erreur d'arrondi près (≈ 1e−13).
    On ne divise pas par Δt : φ vaut alors Δt × pression, mais la vitesse projetée est la même
    (cf. docs/pseudocode.md §1.4) ; la pression n'est pas utilisée par nos diagnostics.
    """
    N = u.shape[0] - 1
    div = np.empty((N, N))
    divergence(u, v, h, div)
    phi = solve_poisson_neumann(div, h)
    subtract_gradient(u, v, phi, h)


def rhs(u, v, t, nu, h, upwind):
    """Second membre SANS pression : R = −A(u) + ν Δu, aux points intérieurs (0 ailleurs).

    Commence par appliquer les CL avec la vitesse du couvercle au temps t de l'étage.
    Retourne (Ru, Rv), de mêmes formes que u, v.
    """
    apply_bc(u, v, lid_velocity(t))
    Au, Av = np.zeros_like(u), np.zeros_like(v)
    Lu, Lv = np.zeros_like(u), np.zeros_like(v)
    advection(u, v, h, upwind, Au, Av)
    laplacian(u, v, h, Lu, Lv)
    return -Au + nu * Lu, -Av + nu * Lv


def step(u, v, t, dt, nu, h, upwind):
    """Un pas de temps de t à t + Δt : RK3-SSP (Shu–Osher) avec projection à chaque étage.

    Mêmes coefficients qu'au Projet 3 (cellule 19, `update_T_method_3`) :
        u1      = P[ u^n + Δt R(u^n, t) ]
        u2      = P[ ¾ u^n + ¼ (u1 + Δt R(u1, t + Δt)) ]
        u^{n+1} = P[ ⅓ u^n + ⅔ (u2 + Δt R(u2, t + Δt/2)) ]
    où P = projection. Projeter chaque étage garde chaque étage à divergence nulle ; c'est
    équivalent à appliquer RK3 à l'équation projetée ∂u/∂t = P R(u).
    Retourne le nouveau couple (u, v), CL appliquées au temps t + Δt.
    """
    Ru, Rv = rhs(u, v, t, nu, h, upwind)
    u1, v1 = u + dt * Ru, v + dt * Rv
    project(u1, v1, h)

    Ru, Rv = rhs(u1, v1, t + dt, nu, h, upwind)
    u2, v2 = 0.75 * u + 0.25 * (u1 + dt * Ru), 0.75 * v + 0.25 * (v1 + dt * Rv)
    project(u2, v2, h)

    Ru, Rv = rhs(u2, v2, t + 0.5 * dt, nu, h, upwind)
    u_new = u / 3.0 + 2.0 / 3.0 * (u2 + dt * Ru)
    v_new = v / 3.0 + 2.0 / 3.0 * (v2 + dt * Rv)
    project(u_new, v_new, h)
    apply_bc(u_new, v_new, lid_velocity(t + dt))
    return u_new, v_new
