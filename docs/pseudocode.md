# Pseudo-code des deux méthodes (version minimale, tout vient du cours)

Notations : $L = U = 1$, $\nu = 1/Re$, $h = 1/N$, $U_t(t) = 1 - e^{-t}$.
Grille **colocalisée** : nœuds $(x_i, y_j) = (ih, jh)$, $i, j = 0..N$, tableaux numpy de forme `(N+1, N+1)`,
indexés `[i, j]` ↔ $(x, y)$. Les nœuds de bord sont sur les parois.
Convention (cours §7.3.4, signes corrigés) : $u = \partial_y\psi$, $v = -\partial_x\psi$,
$\omega = \partial_x v - \partial_y u$, d'où $\Delta\psi = -\omega$.

---

## Briques communes (`common.py`)

```
# Laplacien centré ordre 2, nœuds intérieurs (cours §4.2.1)
lap(f)[i,j] = (f[i+1,j] + f[i-1,j] + f[i,j+1] + f[i,j-1] - 4 f[i,j]) / h²

# Advection a·∂f/∂x + b·∂f/∂y, nœuds intérieurs (cours §2.2, §5.6.2)
#   centré  : ∂f/∂x ≈ (f[i+1,j] - f[i-1,j]) / 2h
#   upwind  : ∂f/∂x ≈ (f[i,j] - f[i-1,j]) / h  si a > 0,   (f[i+1,j] - f[i,j]) / h  si a < 0
#   (même chose en y avec b)

# Poisson Δφ = rhs : gradient conjugué SANS MATRICE (cours §4.7.3, workshop 11.3)
#   Inconnues : nœuds intérieurs. On résout -Δφ = -rhs (car -A est symétrique définie positive, §4.7.3).
#   apply_A(φ) : 1) remplir le bord  2) renvoyer -lap(φ) aux nœuds intérieurs
#     Dirichlet (ψ) : bord = 0
#     Neumann  (p)  : bord = voisin intérieur (p_0 = p_1, ordre 1, cours §4.2.2.2)
#                     → l'opérateur reste symétrique ; il annule les constantes (p défini à C près)
#   Neumann seulement :
#     rhs -= mean(rhs)   # condition de compatibilité : sinon pas de solution (slides : « existence »)
#
#   r = b - apply_A(φ0) ; d = r          # φ0 = solution du pas de temps précédent (départ à chaud)
#   répéter jusqu'à ||r|| ≤ tol·||b|| :
#       Ad = apply_A(d)
#       α  = (r·r) / (d·Ad)
#       φ += α d ;  r_new = r - α Ad
#       β  = (r_new·r_new) / (r·r) ;  d = r_new + β d ;  r = r_new
#       Neumann : φ -= mean(φ)            # moyenne nulle, demandé par les slides ;
#                                         # sans effet sur r car apply_A(constante) = 0
#   puis remplir le bord de φ (même règle que dans apply_A)

# Pas de temps (cours §5.3) — marge 0.8 sur chaque limite
dt = 0.8 * min( h²/(4ν) ,  h/2 )            # Fourier 2D, CFL upwind 2D (C_x + C_y ≤ 1)
si advection centrée : dt = min(dt, 0.8 * 2ν)   # C² ≤ 2 Fo (advection-diffusion, §5.3.1.4)

# Diagnostics
K = ½ Σ (u² + v²) h²                         # somme des trapèzes 2D (workshop 8)
dudy[i] = (3 u[i,N] - 4 u[i,N-1] + u[i,N-2]) / 2h      # décentré ordre 2 (aide-mémoire §2.6)
P = Σ_i ν dudy[i] U_t h                      # trapèzes en x (poids ½ aux deux coins)
```

Coins du couvercle : on prend $u = 0$ aux deux nœuds de coin (ils appartiennent aux parois fixes).

---

## M1 — Projection (cours §7.3.2, workshop 15)

```
u = v = 0 ; t = 0
tant que non stationnaire (ou t < t_final) :
    # 1) prédiction sans pression (Euler explicite)
    u* = u + dt * ( -(u ∂x u + v ∂y u) + ν lap(u) )      # intérieur
    v* = v + dt * ( -(u ∂x v + v ∂y v) + ν lap(v) )
    CL : u*, v* = valeurs de paroi (0 partout, u* = U_t(t+dt) au couvercle)

    # 2) pression :  Δp = (1/dt) ∇·u*      (∇· centré ordre 2 aux nœuds intérieurs)
    p = solve_neumann( (∂x u* + ∂y v*) / dt , p0 = p )   # gradient conjugué, moyenne de p nulle

    # 3) correction :  u = u* - dt ∇p      (∇p centré, nœuds intérieurs seulement)
    u = u* - dt ∂x p ;  v = v* - dt ∂y p
    t += dt
    stationnaire  ⇔  max|u_new - u_old| < 1e-6 et t > 10   (rampe terminée)
    tous les n pas : enregistrer t, K, P
```

---

## M2 — Fonction de courant / vorticité (cours §7.3.4)

```
ω = ψ = 0 ; t = 0
tant que ... :
    # 1) transport de ω (Euler explicite), intérieur
    ω_int = ω + dt * ( -(u ∂x ω + v ∂y ω) + ν lap(ω) )

    # 2) Δψ = -ω, ψ = 0 aux parois
    ψ = solve_dirichlet(-ω)

    # 3) vitesses (centré) à l'intérieur ; parois : valeurs imposées
    u = (ψ[i,j+1] - ψ[i,j-1]) / 2h ;  v = -(ψ[i+1,j] - ψ[i-1,j]) / 2h
    u au couvercle = U_t(t+dt), 0 ailleurs ; v = 0 sur les parois

    # 4) ω de paroi (formule de Thom, voir dérivation ci-dessous)
    couvercle : ω[i,N] = -2 ψ[i,N-1] / h² - 2 U_t / h
    bas       : ω[i,0] = -2 ψ[i,1]   / h²
    gauche    : ω[0,j] = -2 ψ[1,j]   / h²
    droite    : ω[N,j] = -2 ψ[N-1,j] / h²
    t += dt ;  mêmes diagnostics et critère d'arrêt que M1
```

**Dérivation de Thom (couvercle)**, avec le seul outil du cours §2.1 (Taylor) :
$\psi_{N-1} = \psi_N - h\,\partial_y\psi + \tfrac{h^2}{2}\partial_{yy}\psi + O(h^3)$, avec $\psi_N = 0$,
$\partial_y\psi = u = U_t$, et à la paroi $\omega = -\partial_{yy}\psi$ (car $\partial_{xx}\psi = 0$ le long de la
paroi). D'où $\omega_N = -2\psi_{N-1}/h^2 - 2U_t/h$. Les autres parois : même calcul avec une vitesse de paroi nulle.
Remarque (oral) : c'est la méthode des nœuds fantômes du Projet 3 écrite autrement. Avec un fantôme
$\psi_{N+1} = \psi_{N-1} + 2hU_t$ (dérivée $\partial_y\psi = U_t$ centrée) et le laplacien centré au couvercle,
on retrouve exactement la même formule.

Pour les lignes de courant de M1 (Partie B) : $\omega = \partial_x v - \partial_y u$ (centré) puis
`ψ = solve_dirichlet(-ω)` — même fonction que M2.
