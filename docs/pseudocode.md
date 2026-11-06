# Pseudo-code des deux méthodes (Phase 0)

> **Brouillon à relire et valider (ou réécrire à la main)** avant la Phase 1, comme recommandé aux
> Projets 1 et 3. Il fixe les conventions d'indices, les conditions aux limites et l'ordre des
> opérations. Les choix encore ouverts (D1–D8 de [PLAN.md](../PLAN.md) §8) sont signalés par **[Dx]**.

Notations communes : $L = U = 1$, $\nu = 1/Re$, $h = 1/N$, $U_t(t) = 1 - e^{-t}$.
Indices : `[i, j]` ↔ $(x, y)$ (`i` le long de $x$, `j` le long de $y$), tableaux numpy de forme `(nx, ny)`.
Convention de signe (identique à Erturk 2005, cf. [ref_data/README.md](../ref_data/README.md)) :
$u = \partial_y\psi$, $v = -\partial_x\psi$, $\omega = \partial_x v - \partial_y u$ ; le tourbillon principal
a $\psi < 0$, $\omega < 0$.

---

## Méthode 1 — variables primitives, grille MAC, projection de Chorin

### 1.1 Grille (N × N cellules)

| Champ | Position | Indices utiles | Avec fantômes |
|---|---|---|---|
| `p[i, j]` | centre $((i+\tfrac12)h, (j+\tfrac12)h)$ | `i, j = 0..N-1` | — (Neumann géré dans le solveur) |
| `u[i, j]` | face verticale $(i h, (j+\tfrac12)h)$ | `i = 0..N`, `j = 0..N-1` | rangées fantômes `j = -1` et `j = N` → forme `(N+1, N+2)` |
| `v[i, j]` | face horizontale $((i+\tfrac12)h, j h)$ | `i = 0..N-1`, `j = 0..N` | colonnes fantômes `i = -1` et `i = N` → forme `(N+2, N+1)` |

(En Python on décalera les indices fantômes de +1 ; la docstring de chaque fonction le rappellera.)

### 1.2 Conditions aux limites — `apply_bc(u, v, t)`

```
# vitesses NORMALES aux parois : posées directement sur la face de paroi
u[0, :]  = 0 ; u[N, :]  = 0          # parois gauche/droite
v[:, 0]  = 0 ; v[:, N]  = 0          # parois bas/haut
# vitesses TANGENTIELLES : par fantômes (la moyenne fantôme/intérieur = valeur de paroi)
u[:, -1] = -u[:, 0]                  # bas  : u_paroi = 0
u[:,  N] = 2*U_t(t) - u[:, N-1]      # haut : u_paroi = U_t(t)   (couvercle)
v[-1, :] = -v[0, :]                  # gauche : v_paroi = 0
v[ N, :] = -v[N-1, :]                # droite : v_paroi = 0
```
Coins du couvercle : la valeur de paroi saute de $U_t$ à 0 — c'est la singularité (PLAN §3.3.1) ;
avec ce traitement, les faces de coin `u[0, N-1]` et `u[N, N-1]` restent à 0.

### 1.3 Opérateurs discrets (ordre 2, aux points de `u` et de `v`)

```
# Laplacien (5 points) au point u[i,j] (même chose pour v)
L(u)[i,j] = (u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] - 4 u[i,j]) / h²

# Advection, forme conservative : A_u = ∂(uu)/∂x + ∂(vu)/∂y au point u[i,j]
#   flux en x évalués aux centres (i±½, j), flux en y aux nœuds (i, j±½)
#   vitesse transportante = moyenne de 2 voisins ; quantité transportée = selon le schéma [D2]
#     centré  : moyenne des 2 voisins           (ordre 2, non dissipatif)
#     upwind1 : valeur amont (signe de la vitesse transportante) (ordre 1, ν_num ≈ |U| h / 2)
#     QUICK   : interpolation quadratique amont (ordre 3, repris du Projet 3)
A_u[i,j] = (F(i+½,j) - F(i-½,j)) / h + (G(i,j+½) - G(i,j-½)) / h
# A_v symétrique.

# Divergence au centre p[i,j] — exactement nulle après projection
div[i,j] = (u[i+1,j] - u[i,j]) / h + (v[i,j+1] - v[i,j]) / h

# Gradient de φ aux faces intérieures
(∂φ/∂x)[i,j] = (φ[i,j] - φ[i-1,j]) / h       pour i = 1..N-1  (faces de paroi : jamais corrigées)
(∂φ/∂y)[i,j] = (φ[i,j] - φ[i,j-1]) / h       pour j = 1..N-1
```

### 1.4 Un pas de temps — `step_m1(state, dt, t)`

Intégration en temps **[D2]** : RK3-SSP (Projet 3, cellule 19) **avec projection à chaque étage**
(chaque étage reste à divergence nulle ; c'est la variante la plus simple à justifier).

```
fonction rhs(u, v, t):                      # second membre SANS pression
    apply_bc(u, v, t)
    return -A_u(u,v) + ν L(u),  -A_v(u,v) + ν L(v)

fonction project(u*, v*, dt_eff, t):
    apply_bc(u*, v*, t)
    b = div(u*, v*) / dt_eff                # N×N, somme(b) = 0 exactement (flux de paroi nuls)
    φ = solve_poisson_neumann(b)            # ∇²φ = b, ∂φ/∂n = 0, moyenne de φ fixée à 0  [D3]
    u* -= dt_eff ∂φ/∂x ; v* -= dt_eff ∂φ/∂y  (faces intérieures)
    apply_bc(u*, v*, t)
    return u*, v*                            # max|div| ≈ 1e-13

# RK3-SSP (Shu–Osher)
u1 = project(u  + dt rhs(u, t),                          dt,      t + dt)
u2 = project(¾u + ¼(u1 + dt rhs(u1, t+dt)),              dt/4,    t + dt/2)
u  = project(⅓u + ⅔(u2 + dt rhs(u2, t+dt/2)),            2dt/3,   t + dt)
```
Remarque : le `dt_eff` de chaque projection n'influence que l'échelle de φ (la vitesse projetée
est la même) ; la pression physique, si on en a besoin, vaut $p \approx \varphi$ à $O(\Delta t)$ près.

### 1.5 Solveur de Poisson Neumann — `solve_poisson_neumann(b)` **[D3]**

```
# Variante rapide (production) : DCT-II 2D, exacte pour le laplacien 5 points MAC avec Neumann homogène
b̂ = dctn(b, type=2)
λ[k,l] = (2 cos(π k / N) - 2)/h² + (2 cos(π l / N) - 2)/h²
φ̂[k,l] = b̂[k,l] / λ[k,l]   pour (k,l) ≠ (0,0) ;  φ̂[0,0] = 0   (constante arbitraire)
φ = idctn(φ̂, type=2)

# Variante « cours » : SOR du Projet 1 (cellules 13, 17), fantômes φ[-1,j] = φ[0,j] etc.
#   critère ‖r‖/‖r0‖ < tol ; à comparer en coût et en effet de tol sur la divergence résiduelle.
```

### 1.6 Pas de temps et boucle — `run(config)`

```
dt = min(CFL · h / U_max,  Fo · h² / ν)       # CFL ≈ 0.5 (RK3), Fo ≈ 0.25 (2D explicite)
t = 0 ; u = v = 0
tant que t < t_final et non stationnaire :
    u, v = step_m1(u, v, dt, t) ; t += dt
    tous les n_diag pas : enregistrer K, P, ε, résidu R = ‖u^{n+1} - u^n‖∞ / dt
    tous les n_save pas : sauvegarde .npz (reprise possible)
stationnaire  ⇔  R < tol_R  et  |dK/dt| / K < tol_K   (et t ≫ 1 pour que la rampe soit finie)
```

### 1.7 Diagnostics (M1)

```
# K : vitesses ramenées aux centres
K = ½ Σ_{i,j} [ ((u[i,j]+u[i+1,j])/2)² + ((v[i,j]+v[i,j+1])/2)² ] h²

# P : ∂u/∂y au couvercle par différence fantôme–intérieur (ordre 2 au point de paroi)
P = Σ_{i=0..N} w_i · ν (u[i,N] - u[i,N-1]) / h · U_t · h     # w_i : poids des trapèzes (½ aux coins)
#   = Σ ν · 2(U_t - u[i,N-1])/h · U_t · w_i h    ; coins i=0, N : contribution singulière (log N)

# ε : vorticité aux nœuds (coins de cellules), intérieur + parois
ω[i,j] = (v[i,j] - v[i-1,j]) / h - (u[i,j] - u[i,j-1]) / h
ε = ν Σ ω² w_ij h²

# Contrôle : dK/dt ≈ P - ε (dK/dt par différences centrées sur la série temporelle)
```

---

## Méthode 2 — fonction de courant–vorticité (ψ, ω), grille colocalisée

### 2.1 Grille

Nœuds $x_i = i h$, $y_j = j h$, `i, j = 0..N` → tableaux `(N+1, N+1)`. Les nœuds de bord sont
sur les parois (pas de fantômes).

### 2.2 Équations

$$\partial_t\omega + u\,\partial_x\omega + v\,\partial_y\omega = \nu\,\Delta\omega,\qquad
\Delta\psi = -\omega,\qquad u = \partial_y\psi,\ v = -\partial_x\psi$$

### 2.3 Conditions aux limites

```
ψ = 0 sur les 4 parois (parois imperméables, débit total nul)

# Vorticité de paroi (cœur du problème M2) — à comparer [Projet 1 : ordre des CL]
# Couvercle (j = N), avec ψ[i,N] = 0 et ∂ψ/∂y = U_t :
Thom  (ordre 1) : ω[i,N] = -2 ψ[i,N-1] / h² - 2 U_t / h
Woods (ordre 2) : ω[i,N] = -3 ψ[i,N-1] / h² - 3 U_t / h - ω[i,N-1] / 2
# Paroi du bas (j = 0), ∂ψ/∂y = 0 :
Thom  : ω[i,0] = -2 ψ[i,1] / h²
Woods : ω[i,0] = -3 ψ[i,1] / h² - ω[i,1] / 2
# Parois gauche/droite : idem avec ψ[1,j] et ψ[N-1,j] (vitesse de paroi nulle).
# Coins : ω non défini (singularité) ; valeur conventionnelle (ex. 0 ou moyenne) — sans effet sur
# l'intérieur avec un stencil 5 points, mais à préciser pour ε et P.
```
Dérivation de Thom (couvercle) : Taylor $\psi_{N-1} = \psi_N - h\,\partial_y\psi + \tfrac{h^2}{2}\partial_{yy}\psi$,
avec $\psi_N = 0$, $\partial_y\psi = U_t$ et, à la paroi, $\omega = -\partial_{yy}\psi$ (car $\partial_{xx}\psi = 0$ le long
de la paroi). Woods ajoute le terme en $h^3$ en estimant $\partial_{yyy}\psi$ par $(\omega_{N-1} - \omega_N)/h$.

### 2.4 Un pas de temps — `step_m2(state, dt, t)`

```
fonction rhs(ω, t):
    ψ = solve_poisson_dirichlet(-ω intérieur)        # Δψ = -ω, ψ = 0 au bord  [D3 : DST-I ou SOR]
    u[i,j] =  (ψ[i,j+1] - ψ[i,j-1]) / 2h            # intérieur
    v[i,j] = -(ψ[i+1,j] - ψ[i-1,j]) / 2h
    u, v de paroi = valeurs imposées (U_t au couvercle, 0 ailleurs)
    ω de paroi = Thom ou Woods (à partir de ψ et de U_t(t))
    return -(u ∂xω + v ∂yω) + ν Δω   aux nœuds intérieurs   # advection [D2] : centré / upwind1 / QUICK

# RK3-SSP sur ω (intérieur uniquement), CL de paroi recalculées à chaque étage
ω1 = ω + dt rhs(ω, t)
ω2 = ¾ω + ¼(ω1 + dt rhs(ω1, t+dt))
ω  = ⅓ω + ⅔(ω2 + dt rhs(ω2, t+dt/2))
```

### 2.5 Solveur de Poisson Dirichlet — `solve_poisson_dirichlet(b)` **[D3]**

```
# DST-I 2D sur les (N-1)² nœuds intérieurs, exacte pour le laplacien 5 points
b̂ = dstn(b, type=1)
λ[k,l] = (2 cos(π k / N) - 2)/h² + (2 cos(π l / N) - 2)/h²     k, l = 1..N-1
ψ = idstn(b̂ / λ, type=1)
```
(Interface commune avec M1 : `solve_poisson(rhs, bc)` avec `bc ∈ {"neumann", "dirichlet"}`.)

### 2.6 Diagnostics (M2)

```
K = ½ Σ (u² + v²) w_ij h²                         # trapèzes 2D sur les nœuds
P = Σ_i ν (∂u/∂y)[i,N] U_t w_i h  avec  ∂u/∂y = -ω  au couvercle  (car ∂v/∂x = 0 le long du couvercle)
  = -ν U_t Σ_i ω[i,N] w_i h                       # ω[i,N] < 0  ⇒  P > 0
ε = ν Σ ω² w_ij h²
```

---

## Points communs à vérifier en Phase 1 (tests à écrire AVANT les noyaux)

1. Poisson Neumann (DCT) et Dirichlet (DST) sur solution analytique → erreur ≈ ordre 2 en $h$ par rapport
   à la solution continue, et résidu discret ≈ 1e-12.
2. Opérateurs (laplacien, advection centrée) sur champs analytiques → pente 2 sur ≥ 3 grilles ;
   upwind1 → pente 1.
3. M1 : `max |div u| < 1e-12` après `project`.
4. M2 : ordre global Thom vs Woods sur une solution manufacturée (attendu : Thom dégrade, cf. Projet 1 q. g).
5. Bilan `dK/dt = P - ε` sur un run court à Re = 100.
