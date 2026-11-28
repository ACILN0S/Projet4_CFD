# Structure du code — rôle de chaque fichier et de chaque fonction

Ce document décrit le code du Projet 4 tel qu'il est sur la branche `Partie-C` le 24 novembre 2026 : le code de Bastian (commit `ca71e32`), plus les ajouts de la Partie C. Il sert à retrouver vite où se trouve chaque étape des méthodes, et à pouvoir réexpliquer chaque fonction à l'oral.

## 1. Vue d'ensemble

**Il n'y a aucune classe dans le projet** : uniquement des fonctions, et des dictionnaires pour transporter l'état et les résultats. C'est un choix de simplicité : chaque fonction reprend une formule ou une étape du cours, et se lit de haut en bas.

Le code est organisé en trois couches :

| Couche | Fichiers | Rôle |
| --- | --- | --- |
| **Bibliothèque de calcul** | `cavity/common.py`, `cavity/projection.py`, `cavity/psi_omega.py` | les deux méthodes (M1, M2) et toutes leurs briques |
| **Scripts de lancement** | `scripts/partie_B.py`, `scripts/ruche_partie_B.slurm`, `scripts/partie_C.py`, `scripts/lancer_partie_C.sh` | calculs longs, un processus par cas, sauvegarde en `.npz` |
| **Analyse** | `notebooks/A_Re100.ipynb`, `B_sweep_Re.ipynb`, `C_unsteady.ipynb` | vérifications, figures, comparaisons, interprétation |

Volume : `common.py` 557 lignes, `projection.py` et `psi_omega.py` 122 lignes chacun, scripts 18 à 113 lignes. Plus de la moitié des lignes sont des commentaires et des docstrings (règle de `CLAUDE.md` : chaque choix cite sa section du cours).

**Conventions communes** (en-tête de `common.py`, `docs/pseudocode.md`) :

- adimensionnement $L = U = \rho = 1$, donc $\nu = 1/Re$ ;
- grille **colocalisée** aux nœuds $(x_i, y_j) = (ih, jh)$, $i, j = 0..N$, $h = 1/N$ ; les nœuds de bord sont **sur** les parois ;
- tous les champs sont des tableaux numpy de forme $(N+1, N+1)$, indexés `[i, j]` ↔ $(x, y)$ ;
- « intérieur » = nœuds $1..N-1$, tableau de forme $(N-1, N-1)$ : l'élément `[i-1, j-1]` est le nœud $(i, j)$.

## 2. Qui appelle qui

Le chemin d'un calcul, du lancement à la sauvegarde :

    scripts/partie_C.py <indice>                      (ou un notebook, ou scripts/partie_B.py)
     └─ run_projection(...)  [M1]   ou   run_psi_omega(...)  [M2]
         ├─ make_grid, time_step                      grille et Δt de stabilité
         ├─ état initial : fluide au repos (dictionnaire de tableaux nuls)
         └─ march_to_steady(step, state, ...)         boucle en temps commune
             ├─ step(state, t)  →  projection_step   ou   psi_omega_step
             │    ├─ advection, laplacian            opérateurs (Numba)
             │    ├─ solve_poisson                   gradient conjugué (Numba)
             │    └─ apply_velocity_bc, lid_velocity (+ Thom pour M2)
             ├─ steady_residual                       critère d'arrêt
             ├─ kinetic_energy, lid_power             diagnostics K et P
             └─ sondes, instantanés, journal          (Partie C, désactivés par défaut)
     └─ np.savez(results/...npz)                      relu par le notebook

**Point clé de l'architecture** : M1 et M2 ne diffèrent que par leur fonction « un pas de temps » (`projection_step`, `psi_omega_step`). Tout le reste est écrit une seule fois dans `common.py` : la boucle en temps, le critère d'arrêt, K, P, le pas de temps, le solveur de Poisson. La comparaison M1/M2 est donc honnête : seule la méthode change.

`march_to_steady` ne connaît pas la méthode : elle reçoit une fonction `step(state, t)` qui avance l'état d'un pas, et un dictionnaire `state` qui contient au moins `"u"` et `"v"`. C'est ce qui permet de la partager.

## 3. cavity/common.py : les briques communes

### 3.1 Grille et couvercle

| Fonction | Rôle | Cours |
| --- | --- | --- |
| `make_grid(n_cells)` | renvoie `x`, `y` (forme $(N+1,)$) et $h = 1/N$ | workshops 8, 11 |
| `lid_velocity(t)` | $U_t(t) = 1 - e^{-t}$ | énoncé, éq. 2 |
| `apply_velocity_bc(u, v, u_lid)` | impose en place les vitesses de paroi : 0 sur les parois fixes, $u = U_t$, $v = 0$ au couvercle ; les deux coins du couvercle gardent $u = 0$ | PLAN.md §2 |

### 3.2 Opérateurs aux différences finies (Numba)

| Fonction | Rôle | Cours |
| --- | --- | --- |
| `laplacian(f, h)` | laplacien centré 5 points, d'ordre 2, aux nœuds intérieurs ; sortie $(N-1, N-1)$ | §2.3, §4.2.1 |
| `advection(a, b, f, h, scheme)` | terme $a\,\partial_x f + b\,\partial_y f$ aux nœuds intérieurs ; `"centered"` (ordre 2) ou `"upwind"` (ordre 1, décentré selon le signe de la vitesse) ; vérifie `scheme` puis appelle le noyau | §2.2, §5.6.2, §5.4.2 |
| `_advection_kernel(...)` | noyau compilé de `advection` : la double boucle `for i, j` | — |

Le même opérateur `advection` sert à M1 ($f = u$ puis $f = v$) et à M2 ($f = \omega$). Le choix centré/upwind est le choix « sage / malheureux » de l'énoncé.

### 3.3 Solveur de Poisson : gradient conjugué sans matrice

C'est la partie la plus longue (environ 230 lignes), et la plus coûteuse à l'exécution. Elle résout $\Delta\phi = \text{rhs}$ pour trois conditions aux limites, repérées par un code entier (les noyaux Numba n'acceptent pas les chaînes de caractères) :

| Code | Condition | Utilisée pour | Inconnues |
| --- | --- | --- | --- |
| `DIRICHLET` | $\phi = 0$ aux parois | $\psi$ (M2) | nœuds intérieurs |
| `NEUMANN` | $\partial_n\phi = 0$ d'ordre 1, $\phi_0 = \phi_1$ | ancienne pression de M1 (comparaison, notebook A) | nœuds intérieurs |
| `NEUMANN_GHOST` | $\partial_n\phi = 0$ d'ordre 2 par nœuds fantômes miroirs | **pression de M1** (ADR 0002) | tous les nœuds |

| Fonction | Rôle |
| --- | --- |
| `solve_poisson(rhs, h, bc, phi0, tol, max_iter)` | **fonction publique** : vérifie les formes, construit le second membre $b$ du système symétrique $A\phi = b$ (avec la condition de compatibilité en Neumann : moyenne de rhs retirée), appelle le gradient conjugué, lève une erreur s'il ne converge pas en `max_iter` itérations ; renvoie $\phi$ et le nombre d'itérations |
| `_conjugate_gradient(b, phi, h, mode, tol, max_iter)` | l'algorithme du cours §4.7.3, sans préconditionnement ; arrêt sur le résidu relatif $\lVert r\rVert \le \text{tol}\,\lVert b\rVert$ ; en Neumann, retire la moyenne de $\phi$ à chaque itération (exigé par les slides du prof) |
| `_apply_matrix(d, h, mode, out)` | produit $A\,d$ **sans stocker $A$** : $A = -\Delta_h$ (SDP, ce qu'exige le gradient conjugué) ; pour les fantômes, chaque ligne est multipliée par les poids des trapèzes pour rester symétrique |
| `_remove_mean(phi, mode)` | $\phi \mathrel{-}= \text{moyenne}(\phi)$ : moyenne simple, ou pondérée par les trapèzes (fantômes) |
| `fill_boundary(phi, bc)` / `_fill_boundary_kernel` | remplit les nœuds de bord selon la condition (Dirichlet : 0 ; Neumann : bord = voisin) |
| `_trapezoid_weight(k, n)` | poids ½ aux extrémités, 1 ailleurs (workshop 8.2) |
| `BC_CODES` | dictionnaire nom → code entier |

**Départ à chaud** : `phi0` est la solution du pas de temps précédent. C'est ce qui rend le gradient conjugué rapide une fois le transitoire passé (en Partie C, M1 : 21 itérations par pas à la fin, contre plusieurs centaines au démarrage).

### 3.4 Pas de temps et critère d'arrêt

| Fonction | Rôle | Cours |
| --- | --- | --- |
| `time_step(h, nu, scheme, safety=0.8)` | $\Delta t$ = 0.8 × la plus contraignante des limites : Fourier $h^2/4\nu$, CFL $h/2$, et pour le centré $2\nu$ | §5.3 |
| `steady_residual(new, old, dt)` | norme L2 de $R = (f^{n+1} - f^n)/\Delta t$ aux nœuds intérieurs ; sur $u, v$ pour M1, sur $\omega$ pour M2 | §4.3.2, ADR 0003 §2 |

### 3.5 Boucle en temps commune

`march_to_steady(step, state, h, dt, nu, residual_tol, t_ref, t_max, label, residual_fields, probes, snapshot_times, log_every)`

À chaque pas :

1. `step(state, t)` avance l'état de $\Delta t$ et renvoie le nombre d'itérations du gradient conjugué ;
2. calcul du résidu stationnaire ;
3. enregistrement de $t$, $K$, $P$, itérations, résidu ;
4. (Partie C) $u$, $v$ aux sondes, copie des champs aux instants demandés, ligne de journal ;
5. arrêt avec erreur si le résidu n'est plus fini (divergence) ;
6. critère d'arrêt : $\lVert R\rVert < \varepsilon\,\lVert R(t_\text{ref})\rVert$ (avec $\varepsilon$ = `residual_tol`) ; `residual_tol = 0` désactive l'arrêt (Partie C : calcul jusqu'à $t_\text{max}$).

Elle renvoie l'état final et un dictionnaire `info` : historiques (`t`, `K`, `P`, `cg_iterations`, `residual`), `n_steps`, `t_final`, `wall_time`, `converged` (arrêt par le critère et non par $t_\text{max}$), et en Partie C `probe_u`, `probe_v`, `snap_t`, `snap_<champ>`.

Les trois options de la Partie C (`probes`, `snapshot_times`, `log_every`) valent `None` par défaut : le calcul est alors identique bit à bit à celui des Parties A et B (vérifié). L'historique des sondes est alloué une seule fois (`t_max/dt + 1` lignes) pour ne pas ralentir la boucle.

### 3.6 Diagnostics imposés

| Fonction | Formule | Cours |
| --- | --- | --- |
| `kinetic_energy(u, v, h)` | $K = \tfrac12 \iint (u^2 + v^2)\,dS$, trapèzes en x puis en y | workshop 8.2 |
| `lid_power(u, h, nu, u_lid)` | $P = \int_0^1 \nu\,\partial_y u|_{y=1}\,U_t\,dx$ ; $\partial_y u$ décentrée d'ordre 2 $(3u_N - 4u_{N-1} + u_{N-2})/2h$, puis trapèzes | workshop 11.2 |

Elles sont appelées une fois par pas, en numpy pur : leur coût est négligeable devant le Poisson.

## 4. cavity/projection.py : méthode M1

Méthode de projection (pas fractionnaire, cours §7.3.2, workshop 15). Inconnues : $u$, $v$, $p$.

| Fonction | Rôle |
| --- | --- |
| `projection_step(u, v, p, t, dt, h, nu, scheme, cg_tol, pressure_bc)` | **un pas de temps de M1**, en trois étapes (voir ci-dessous) ; renvoie les nouveaux $u, v, p$ et le nombre d'itérations |
| `divergence_all_nodes(u, v, h)` | $\nabla\cdot u$ à tous les nœuds, parois comprises (centrée à l'intérieur, décentrée d'ordre 2 aux bords : `np.gradient(..., edge_order=2)`) ; nécessaire pour la Neumann à nœuds fantômes |
| `run_projection(n_cells, reynolds, scheme, cg_tol, residual_tol, t_ref, t_max, pressure_bc, dt, probes, snapshot_times, log_every)` | **point d'entrée** : grille, $\Delta t$ (refuse un $\Delta t$ au-dessus de la limite de stabilité), état au repos, définit `step` et appelle `march_to_steady` ; renvoie un dictionnaire avec paramètres, champs finaux et historiques |

**Un pas de M1** (`projection_step`) :

1. **prédiction** sans pression : $u^* = u + \Delta t\,[-(u\cdot\nabla)u + \nu\Delta u]$ aux nœuds intérieurs, puis vitesses de paroi imposées à $t + \Delta t$ ;
2. **pression** : $\Delta p = \nabla\cdot u^*/\Delta t$, Neumann d'ordre 2 (fantômes), départ à chaud depuis le $p$ précédent ;
3. **correction** : $u^{n+1} = u^* - \Delta t\,\nabla p$, gradient centré aux nœuds intérieurs.

## 5. cavity/psi_omega.py : méthode M2

Fonction de courant – vorticité (cours §7.3.4). Inconnues : $\omega$ et $\psi$ ; $u$, $v$ en sont déduites. Convention : $u = \partial_y\psi$, $v = -\partial_x\psi$, $\Delta\psi = -\omega$.

| Fonction | Rôle |
| --- | --- |
| `psi_omega_step(omega, psi, u, v, t, dt, h, nu, scheme, cg_tol)` | **un pas de temps de M2**, en quatre étapes (voir ci-dessous) |
| `velocity_from_streamfunction(psi, h, u_lid)` | $u = \partial_y\psi$, $v = -\partial_x\psi$ par différences centrées, puis vitesses de paroi imposées |
| `apply_thom_wall_vorticity(omega, psi, h, u_lid)` | **formule de Thom** : $\omega$ de paroi par Taylor sur $\psi$ ; au couvercle $\omega = -2\psi_{N-1}/h^2 - 2U_t/h$ ; coins à 0 |
| `run_psi_omega(...)` | **point d'entrée**, mêmes paramètres et même sortie que `run_projection` (sans `pressure_bc`) |

**Un pas de M2** (`psi_omega_step`) :

1. **transport de $\omega$** : $\omega^{n+1} = \omega + \Delta t\,[-(u\cdot\nabla)\omega + \nu\Delta\omega]$ aux nœuds intérieurs ;
2. **Poisson** $\Delta\psi = -\omega^{n+1}$, Dirichlet $\psi = 0$, départ à chaud ;
3. **vitesses** à partir du nouveau $\psi$ ;
4. **vorticité de paroi** par Thom, qui servira au pas suivant.

Pas de pression, donc pas de condition de Neumann. Mais la formule de Thom, d'ordre 1 localement, sous-estime $|\omega|$ à la paroi d'autant plus que $Re$ est grand (Partie B, `Doc/partie_B_limites.md` §3) : c'est la principale faiblesse de M2 à haut Reynolds.

## 6. Numba : pourquoi certaines fonctions commencent par un tiret bas

Les noyaux coûteux (laplacien, advection, gradient conjugué et ses auxiliaires) sont décorés par `@njit(cache=True)` : Numba les compile en code machine au premier appel, et garde la compilation sur disque. On peut ainsi écrire des boucles `for i, j` explicites, qui reprennent mot pour mot les formules du cours, sans perdre la vitesse.

Règle de nommage :

- **sans tiret bas** (`advection`, `solve_poisson`, `fill_boundary`) : fonction publique, en Python pur ; elle vérifie ses arguments (schéma, condition aux limites, formes) et donne des messages d'erreur clairs ;
- **avec tiret bas** (`_advection_kernel`, `_conjugate_gradient`, `_apply_matrix`…) : noyau compilé, appelé seulement par la fonction publique.

Exception : `laplacian` est compilée directement, sans fonction publique, car elle n'a aucun argument à vérifier.

Pas de `fastmath` : arithmétique IEEE exacte, mêmes résultats qu'une version numpy aux arrondis près. Chaque calcul utilise un seul cœur ; on parallélise en lançant plusieurs cas à la fois.

## 7. Scripts de lancement

| Fichier | Rôle |
| --- | --- |
| `scripts/partie_B.py` | 48 cas de la Partie B (11 $Re$ de 100 à 5000, M1 et M2, $N = 64$, 128, plus 256 à $Re$ = 3200 et 5000) ; `python scripts/partie_B.py <indice>` lance un cas et sauve `results/B_...npz` ; `liste` affiche les cas |
| `scripts/ruche_partie_B.slurm` | soumission des 48 cas au calculateur Ruche (tableau SLURM, un cas par tâche) |
| `scripts/partie_C.py` | 9 cas de la Partie C à $Re = 10\,000$ ; modes `liste`, `bench <i>` (coût par pas sur 3000 pas), `<i>` (calcul complet avec sondes, instantanés, journal) |
| `scripts/lancer_partie_C.sh` | lance des cas de la Partie C en parallèle sur le Mac (4 à la fois, un thread par processus, `caffeinate` contre la veille), un journal par cas dans `results/logs/` |

Fonctions de `scripts/partie_C.py` :

| Fonction / constante | Rôle |
| --- | --- |
| `CASES` | liste des cas `(méthode, schéma, facteur sur Δt, N, t_max)` ; un indice désigne toujours le même cas (on ajoute à la fin) |
| `PROBES` | les 6 sondes en coordonnées $(x, y)$, communes aux grilles 128 et 256 |
| `snapshot_times(t_max)` | instants des instantanés : tous les 5, puis tous les 0.2 sur les 20 dernières unités |
| `case_dt(scheme, factor, n)` | $\Delta t$ du cas = `factor` × `time_step` |
| `run_path(method, scheme, factor, n, t_max)` | nom du fichier de résultats (une seule source de vérité, importée par le notebook C) |
| `run_case(...)` | appelle `run_projection` ou `run_psi_omega` avec `residual_tol = 0` (pas d'arrêt avant $t_\text{max}$) |

Chaque script fixe `MKL_THREADING_LAYER=SEQUENTIAL` **avant** d'importer numpy (sinon l'algèbre linéaire de numpy, donc matplotlib, fait planter le processus sur une des machines de l'équipe).

## 8. Notebooks : analyse et vérifications

Les notebooks ne lancent pas les calculs longs : ils relisent les `.npz` (et, pour la Partie A, lancent les petits calculs eux-mêmes). Les fonctions qu'ils définissent sont des outils d'analyse, pas de calcul.

| Notebook | Sections | Fonctions propres |
| --- | --- | --- |
| `A_Re100.ipynb` | 1. vérification du Poisson sur solution connue ; 2. écoulement à $Re = 100$ ; 3. validation M1, M2 contre Ghia ; 4. convergence en maillage | `run_path`, `get_run` (calcule ou relit), `centerlines` (profils sur les lignes médianes), `order` (ordre de convergence), `streamfunction`, `ghia_error` |
| `B_sweep_Re.ipynb` | 1. calculs ; 2. évolution de l'écoulement avec $Re$ ; 3. comparaison Ghia et Erturk ; 4. $P(Re)$ et loi d'échelle | `run_path`, `load_run`, `get_run`, `centerlines`, `slope` (pente log-log), `streamfunction`, `vortex_centers`, `reference`, `max_error` |
| `C_unsteady.ipynb` | 1. budget ; 2. K, P, résidu ; 3. période et amplitude ; 4. physique ou numérique ; 5. cycle de vorticité ; 6. sensibilité ; 7. littérature ; 8. bilan | `case_name`, `oscillation` (période, amplitude, taux $\sigma$), `regime`, `vorticity`, `snap_vorticity`, `reference`, `mean_profiles` |

Le notebook C importe `CASES`, `PROBES`, `run_path`… depuis `scripts/partie_C.py` : les noms de fichiers ne sont écrits qu'à un seul endroit.

## 9. Données et documentation

| Dossier / fichier | Contenu |
| --- | --- |
| `ref_data/` | profils de Ghia (1982) et Erturk (2005) en CSV, tourbillon principal ; `README.md` : ce que chaque article a fait, coquilles connues, $Re_c$ |
| `results/` | sorties `.npz` des calculs (non versionnées, trop lourdes), petits résumés CSV de la Partie A, journaux `logs/` |
| `figures/` | figures produites par les notebooks |
| `docs/adr/` | décisions argumentées : 0001 méthodes et outils, 0002 Neumann de la pression, 0003 pas de temps, critère d'arrêt et tolérance du gradient conjugué, 0004 optimisation du code |
| `docs/pseudocode.md` | pseudo-code des deux méthodes et conventions |
| `PLAN.md`, `CLAUDE.md` | plan du projet par étapes ; règles de travail (strict minimum, idées du cours seulement) |
| `Doc/` | énoncé, poly du cours, documents d'explication (limites des Parties A et B, choix de la Partie C, ce document) |

## 10. Pour retrouver une étape du cours dans le code

| Étape du cours | Où |
| --- | --- |
| Euler explicite (§3.1) | `projection_step`, `psi_omega_step` (`+= dt * ...`) |
| Différences centrées / upwind (§2.2) | `_advection_kernel` |
| Laplacien 5 points (§4.2.1) | `laplacian`, `_apply_matrix` |
| Gradient conjugué (§4.7.3) | `_conjugate_gradient` |
| Neumann et compatibilité (slides) | `solve_poisson`, `_remove_mean` |
| Limites de stabilité (§5.3) | `time_step` |
| Équation modifiée, diffusion numérique (§5.4.2) | notebook C §6 (calcul de $D$ et $Re_\text{eff}$) |
| Projection (§7.3.2) | `projection_step` |
| $\psi$–$\omega$ et Thom (§7.3.4, §2.1) | `psi_omega_step`, `apply_thom_wall_vorticity` |
| Trapèzes (workshop 8.2) | `kinetic_energy`, `lid_power`, `_trapezoid_weight` |
