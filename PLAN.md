# Projet 4 — Écoulement dans une cavité entraînée (lid-driven cavity) à haut Reynolds

> Document de cadrage : ce qu'on cherche à faire, comment on va y arriver, comment on découpe le travail,
> et les décisions de design qu'il reste à trancher (section 8).

---

## 1. Le problème en résumé

### 1.1 La configuration physique

Une boîte carrée 2D de côté $L$, remplie d'un fluide newtonien incompressible. Trois parois sont fixes
(adhérence : $u=v=0$). Le **couvercle** (paroi du haut, $y=L$) glisse horizontalement à la vitesse $U_t(t)$
et entraîne le fluide par frottement visqueux. Le fluide se met à tourner : c'est la **cavité entraînée**,
le cas test le plus célèbre de la CFD incompressible.

```
          U_t(t)  ──────►
        ┌──────────────────┐ y = L   (u = U_t, v = 0)
        │                  │
 u=v=0  │   ↻ tourbillon   │  u=v=0
        │     principal    │
        │                  │
        └──────────────────┘ y = 0   (u = v = 0)
       x = 0              x = L
```

On résout **uniquement Navier–Stokes incompressible** (pas de scalaire transporté) :

$$\nabla\cdot\mathbf u = 0,\qquad \partial_t \mathbf u + (\mathbf u\cdot\nabla)\mathbf u = -\tfrac{1}{\rho}\nabla P + \nu\,\Delta\mathbf u$$

Le couvercle démarre en douceur : $U_t(t) = U\,(1-e^{-t/\tau})$ avec $\tau = L/U$.
Un seul paramètre gouverne tout : $Re = UL/\nu$.

### 1.2 Adimensionnement (ce qu'on codera réellement)

On prend $L = U = \rho = 1$, donc $\nu = 1/Re$, temps en unités de $L/U$, et $U_t(t) = 1 - e^{-t}$.
Toutes les grandeurs sont alors directement comparables à la littérature (Ghia et al. et suivants
travaillent ainsi).

### 1.3 Ce qu'on doit mesurer en permanence

| Diagnostic | Formule (adimensionnée) | Sens physique |
|---|---|---|
| Énergie cinétique | $K(t) = \tfrac12\int_\Omega (u^2+v^2)\,dS$ | « quantité de mouvement de rotation » stockée |
| Puissance injectée | $P(t) = \int_0^1 \nu\,\partial_y u\big|_{y=1}\,U_t(t)\,dx$ | travail du couvercle sur le fluide |
| *(ajout proposé)* Dissipation | $\varepsilon(t) = \nu\int_\Omega \omega^2\,dS$ | énergie détruite par viscosité |

**Bilan d'énergie** : $\dfrac{dK}{dt} = P - \varepsilon$. En régime stationnaire, $P = \varepsilon$.
C'est un **contrôle de cohérence gratuit et très convaincant** du solveur (même esprit que le bilan
d'énergie thermique du Projet 3).

### 1.4 Ce que la physique va nous montrer (à quoi s'attendre)

| Régime | Ce qui se passe | Référence |
|---|---|---|
| $Re \sim 100$ | Un tourbillon principal décalé vers le haut-droite, deux minuscules tourbillons secondaires dans les coins du bas. Écoulement **stationnaire**. | Ghia, Ghia & Shin (1982) |
| $Re \sim 1000$ | Le tourbillon principal se recentre, les tourbillons de coin bas grossissent. | Ghia (1982), Botella & Peyret (1998), Erturk (2005) |
| $Re \sim 3000$–$7500$ | Apparition d'un tourbillon dans le coin haut-gauche, puis de tourbillons **tertiaires** dans les coins. Couches limites de plus en plus fines ($\delta \sim Re^{-1/2}$). Toujours stationnaire. | Ghia (1982), Erturk (2005) |
| $Re > Re_c \approx 8000$ | **Bifurcation de Hopf** : la solution stationnaire devient instable, l'écoulement devient **périodique** en temps. | Auteri et al. (2002), Bruneau & Saad (2006), Poliashenko & Aidun (1995)… |

> ⚠️ Les valeurs de $Re_c$ rapportées varient selon les auteurs (de ~7400 à ~8050 selon méthodes et
> maillages). **Cette dispersion est elle-même un résultat** : elle montre la sensibilité de la
> bifurcation aux choix numériques — c'est exactement ce que la Partie C nous demande d'étudier.
> (Références et valeurs exactes à revérifier lors de la phase bibliographie.)

### 1.5 Ce qu'on nous demande concrètement

- **Au moins deux méthodes numériques**, dont au moins une issue du cours. L'objectif explicite est de
  **comparer** et de montrer l'impact d'un choix « sage » ou « malheureux ».
- **Partie A** : écoulement stationnaire à $Re=100$, **validé** (littérature + autre solveur) avec une
  **étude de maillage** qui justifie la grille.
- **Partie B** : plusieurs $Re < 7500$, évolution de la **topologie** (tourbillons), comparaison
  littérature, **$P(Re)$ et loi d'échelle**.
- **Partie C** : $Re > Re_c$, le plus haut $Re$ qu'on peut se permettre (justifié), **sensibilité aux
  choix numériques**, comparaison littérature.
- Esprit de l'énoncé : étude « réaliste » et ouverte. **Un échec est acceptable s'il est expliqué**
  (diagnostic de la cause). Et, comme dans les projets précédents : on doit **maîtriser chaque ligne du
  code** et pouvoir l'expliquer/l'adapter à l'oral.

---

## 2. Ce qu'on récupère des projets précédents

| Acquis | Où | Utilisation dans le Projet 4 |
|---|---|---|
| Solveurs itératifs Jacobi / GS / **SOR**, critère sur résidu $\|r\|/\|r_0\|<\varepsilon$, version Numba optimisée | Projet 1, cellules 13, 17 | Résolution de l'équation de Poisson (pression ou fonction de courant). Sert de **méthode « cours »** et de **contre-exemple de coût** face à un solveur direct. |
| **Leçon : une CL d'ordre 1 détruit l'ordre global** (Neumann naïf → pente 1 au lieu de 2) | Projet 1, q. g), cellule 26 | Directement transposable : conditions de paroi par **cellules fantômes**, vorticité de paroi (formule de Thom vs ordre 2). On vérifiera l'ordre global, pas seulement l'ordre intérieur. |
| Extrapolation de **Richardson / GCI**, valeur ± incertitude | Projet 1 q. h), Projet 3 §5 | Justification de la grille (Partie A), barres d'erreur sur $P$, positions des tourbillons, $Re_c$. |
| Schémas d'**advection upwind 1 / upwind 2 / upwind-biaisé 3**, RK2, **RK3-SSP** | Projet 3, cellules 15, 17, 19 | Termes non linéaires $(\mathbf u\cdot\nabla)\mathbf u$ ou $(\mathbf u\cdot\nabla)\omega$. **Comparaison des schémas = cœur de la Partie C** (la diffusion numérique de l'upwind 1 agit comme un $Re$ effectif plus faible → peut *tuer* la bifurcation). |
| Analyse de **von Neumann** (dissipation/dispersion), limites CFL + Fourier | Projet 3 §2.6 | Choix du $\Delta t$, explication quantitative de la viscosité numérique : $\nu_{num}\approx U\Delta x/2$ pour l'upwind 1, soit $Re_{eff}$ bien inférieur au $Re$ visé. |
| Cellules fantômes, `apply_boundary_conditions`, structure `run_simulation` → dict de résultats, style de figures | Projet 3 cellules 1, 12, 28 | Squelette de code réutilisé presque tel quel. |
| Environnement conda **`CFD`** (numpy 2.5, scipy 1.18, numba 0.68, matplotlib 3.11) | `~/miniconda3/envs/CFD` | Déjà fonctionnel, on le garde. Machine : 32 cœurs logiques → on peut lancer plusieurs $Re$ en parallèle. |

---

## 3. Stratégie numérique proposée

### 3.1 Méthode 1 (« cours ») — Variables primitives $(u,v,p)$, grille décalée MAC, projection de Chorin

**Grille décalée (MAC)** : $p$ au centre des cellules, $u$ sur les faces verticales, $v$ sur les faces
horizontales. Avantages : pas de découplage pair/impair de la pression (« damier »), divergence discrète
nulle à la précision machine, conservation discrète de l'énergie avec l'advection centrée.

```
   +----v----+----v----+
   |         |         |
   u    p    u    p    u        N×N pressions, (N+1)×N valeurs de u, N×(N+1) valeurs de v
   |         |         |        + 1 rangée de cellules fantômes pour les CL de paroi
   +----v----+----v----+
```

**Un pas de temps (projection / pas fractionnaire)** :

```
1. u* = u^n + Δt · [ -A(u^n) + ν L(u^n) ]          (prédiction, sans pression ; RK ou AB2)
2. Résoudre  ∇² φ = (∇·u*) / Δt                     (Poisson, Neumann homogène sur les parois)
3. u^{n+1} = u* - Δt ∇φ                              (projection → ∇·u^{n+1} = 0)
4. Appliquer les CL (fantômes) avec U_t(t^{n+1})
```

- **Conditions de paroi** par fantômes : par ex. en haut $u_{fant} = 2U_t - u_{int}$ (la moyenne vaut
  $U_t$ sur la paroi), $v = 0$ directement sur la face de paroi.
- **Advection** : plusieurs variantes comparables dans le même code (choix de design, §8) :
  centrée ordre 2 (conservative, non dissipative), upwind 1 (le « mauvais choix » volontaire),
  QUICK / upwind-biaisé 3 (repris du Projet 3).
- **Diffusion** : explicite centrée (simple). Contrainte $\Delta t \le Re\,\Delta x^2/4$ — limitante à bas
  $Re$ seulement.
- **Pas de temps** : $\Delta t = \min(\text{CFL}\cdot\Delta x/U,\ \text{Fo}\cdot Re\,\Delta x^2)$.
- **Poisson** : c'est l'étape la plus coûteuse, faite à **chaque pas de temps**. Trois options (§8) :
  SOR du Projet 1 (lent, mais « cours »), **DCT rapide** (exact et en $O(N^2\log N)$ : sur un carré avec
  Neumann, le laplacien discret MAC est diagonalisé par la transformée en cosinus), ou LU creux
  factorisé une fois (`scipy.sparse.linalg.splu`).

### 3.2 Méthode 2 (proposée) — Formulation fonction de courant–vorticité $(\psi,\omega)$

En 2D, on peut éliminer la pression :

$$\partial_t\omega + u\,\partial_x\omega + v\,\partial_y\omega = \nu\Delta\omega,\qquad
\Delta\psi = -\omega,\qquad u = \partial_y\psi,\ v = -\partial_x\psi$$

- Grille **colocalisée** aux nœuds (comme les Projets 1 et 3), $\psi = 0$ sur toutes les parois.
- **Vorticité de paroi** : formule de **Thom** (ordre 1 formellement) ou de **Woods / Briley** (ordre 2) —
  à comparer, c'est exactement la leçon du Projet 1 sur l'ordre des CL.
- Poisson pour $\psi$ avec Dirichlet → SOR (Projet 1) ou **DST rapide**.
- Bonus : $P$ s'obtient directement ($\partial_y u|_{lid} = -\omega|_{lid}$ car $v=0$ le long du couvercle).

**Pourquoi cette paire de méthodes ?** Deux formulations réellement indépendantes (variables, grille,
traitement des parois, solveur de Poisson différents) : si elles s'accordent, c'est une validation
croisée très forte (« autre solveur validé » demandé par l'énoncé). Et chacune a ses faiblesses connues
(CL de vorticité pour $\psi$–$\omega$, erreur de splitting pour la projection), ce qui nourrit la
discussion.

Alternatives possibles pour la méthode 2 (cf. §8, décision D1) : solveur **stationnaire de Newton**
(permet même de suivre la branche stationnaire instable au-delà de $Re_c$), **Lattice Boltzmann**
(paradigme totalement différent), **spectral Chebyshev** (référence très précise, mais plus difficile).

### 3.3 Points durs identifiés à l'avance

1. **Singularité des coins supérieurs.** La vitesse est discontinue aux coins du couvercle
   ($U$ sur le couvercle, $0$ sur la paroi). La vorticité y diverge en $1/r$, donc
   $\int \nu\,\partial_y u\,dx$ diverge **logarithmiquement** quand on raffine : la puissance $P$
   « exacte » est infinie et la valeur discrète croît lentement avec $N$ (~$\log N$).
   → Il faudra le montrer, puis choisir comment rapporter $P(Re)$ (décision D5). Ce point est un
   bon exemple de « limitation à expliquer ».
2. **Coût à haut $Re$.** Couches limites en $Re^{-1/2}$ → à $Re = 10^4$, il faut typiquement $N\ge 256$.
   Avec $\Delta t \approx 0.5/N$, et un transitoire de plusieurs centaines à milliers d'unités $L/U$
   (près de $Re_c$ la perturbation croît très lentement), on parle de $10^5$–$10^6$ pas de temps.
   → Numba + Poisson rapide **indispensables** ; budget de calcul à estimer en Phase 5.
3. **Ralentissement critique près de $Re_c$.** Le taux de croissance de l'instabilité tend vers 0 :
   à $Re = 8100$, l'oscillation peut mettre très longtemps à apparaître. → Simuler plutôt à
   $Re = 8500$–$10000$ et extrapoler $Re_c$ (voir Phase 5).
4. **Diffusion numérique.** Avec l'upwind 1, $\nu_{num}\approx U\Delta x/2$. À $N=256$ :
   $\nu_{num}\approx 2\cdot10^{-3}$, alors que $\nu = 10^{-4}$ à $Re=10^4$ → $Re_{eff}\sim 500$ !
   Le schéma « voit » un écoulement à bas Reynolds et ne bifurquera jamais. **Excellent contre-exemple
   pour la Partie C.**

---

## 4. Plan de travail par phases

Chaque phase a un **critère de sortie** : on ne passe à la suivante que s'il est rempli.

### Phase 0 — Mise en place (½ journée)
- [ ] Structure du dépôt (§6), `environment.yml`, `AGENTS.md`/`CLAUDE.md` (style de code, conventions).
- [ ] Bibliographie : récupérer et **saisir les tableaux de Ghia (1982)** (profils $u(x=0.5,y)$ et
      $v(x,y=0.5)$, positions et intensités des tourbillons pour $Re$ = 100, 400, 1000, 3200, 5000, 7500,
      10000), Botella & Peyret (1998) pour $Re=1000$, Erturk (2005), et 2–3 articles sur la bifurcation
      ($Re_c$, période/fréquence, amplitude).
- [ ] Rédiger à la main le **pseudo-code** de la méthode 1 et de la méthode 2 (recommandation des projets
      précédents : comprendre avant de coder).

**Sortie** : fichiers de référence `ref_data/*.csv` + pseudo-code validé.

### Phase 1 — Méthode 1 : briques de base vérifiées (1–2 jours)
- [ ] Grille MAC, indexation, cellules fantômes, CL (avec $U_t(t)$).
- [ ] Opérateurs discrets : divergence, gradient, laplacien, advection (au moins centrée + upwind 1).
- [ ] Solveur de Poisson (au moins deux variantes, cf. D3) + test sur solution analytique
      (ordre 2 vérifié, comme au Projet 1 q. g).
- [ ] Tests unitaires : $\nabla\cdot\mathbf u^{n+1} \approx 10^{-12}$ après projection ; opérateurs
      appliqués à des champs analytiques → ordre 2 mesuré.
- [ ] *(Optionnel mais très convaincant)* **Solutions manufacturées (MMS)** pour vérifier l'ordre du
      solveur complet en espace et en temps.
- [ ] Diagnostics : $K$, $P$, $\varepsilon$, résidu stationnaire $\|\mathbf u^{n+1}-\mathbf u^n\|/\Delta t$.

**Sortie** : solveur M1 qui tourne, divergence nulle, ordres vérifiés.

### Phase 2 — Partie A avec la méthode 1 : $Re = 100$ (1 jour)
- [ ] Simulation depuis le repos avec la rampe jusqu'au stationnaire (critère sur le résidu et $dK/dt$).
- [ ] Courbes $K(t)$, $P(t)$, $\varepsilon(t)$, bilan $dK/dt = P - \varepsilon$.
- [ ] Comparaison Ghia : profils de vitesse sur les lignes médianes, position et $\psi_{min}$ du
      tourbillon principal, tourbillons secondaires.
- [ ] **Étude de maillage** $N = 32, 64, 128, 256$ : ordre observé, Richardson/GCI sur des grandeurs
      scalaires ($\psi_{min}$, $u_{min}$ sur l'axe, $K$, $P$) → **justifier la grille retenue**.
- [ ] Mise en évidence de la croissance en $\log N$ de $P$ (singularité de coin).

**Sortie** : figure « nos profils vs Ghia » + tableau de convergence + grille justifiée.

### Phase 3 — Méthode 2 et validation croisée (1–2 jours)
- [ ] Implémentation $\psi$–$\omega$ (réutilise les schémas d'advection et le Poisson).
- [ ] Comparaison CL de vorticité Thom vs ordre 2 (ordre global mesuré).
- [ ] Partie A refaite avec M2 → comparaison **M1 vs M2 vs Ghia** (écarts quantifiés).
- [ ] Comparaison des coûts (temps CPU à précision égale, comme au Projet 3 §6).

**Sortie** : deux solveurs indépendants qui s'accordent à $Re=100$ (et idéalement $Re=1000$, cas
Botella & Peyret, beaucoup plus discriminant).

### Phase 4 — Partie B : $Re < 7500$ (2 jours)
- [ ] Balayage $Re \in \{100, 400, 1000, 2000, 3200, 5000, 7500^-\}$ avec grille adaptée (étude de maillage
      allégée à chaque $Re$ ou aux $Re$ extrêmes).
- [ ] Accélérer : initialiser $Re_k$ avec la solution de $Re_{k-1}$ (continuation) — à condition d'en
      discuter (D6) — et/ou lancer les $Re$ en parallèle sur les 32 cœurs.
- [ ] **Topologie** : lignes de courant, détection automatique des tourbillons (extrema de $\psi$,
      changement de signe), positions/intensités vs Ghia, $Re$ d'apparition des tourbillons secondaires
      et tertiaires.
- [ ] **$P(Re)$** en log-log : on attend une transition entre $P^* \sim Re^{-1}$ (régime visqueux, couche
      limite d'épaisseur $\sim L$) et $P^* \sim Re^{-1/2}$ (couche limite de Blasius $\delta\sim L\,Re^{-1/2}$).
      Ajustement de pente par morceaux, avec barres d'erreur (maillage) et discussion de la partie
      « coin » en $Re^{-1}\log N$.
- [ ] Comparer M1 et M2 sur au moins 2 valeurs de $Re$ (ex. 1000 et 5000).

**Sortie** : planche de topologies + graphe $P(Re)$ + loi d'échelle justifiée physiquement.

### Phase 5 — Partie C : perte de stationnarité (3–4 jours, la plus coûteuse)
- [ ] **Budget de calcul** : mesurer le coût par pas à $N=256$ et $384/512$, en déduire le $Re$ et la durée
      simulée accessibles → c'est la **justification** demandée du « plus haut $Re$ ».
- [ ] Simulations longues à $Re$ = 8500, 9000, 10000 (+ plus haut si budget) : séries temporelles
      $K(t)$, $P(t)$, sondes ponctuelles (ex. près du coin bas-droit, positions utilisées dans la
      littérature).
- [ ] Analyse : FFT → **fréquence** du cycle limite, amplitude, portrait de phase ($K$ vs $P$).
- [ ] **Estimation de $Re_c$** : près d'une bifurcation de Hopf, $A^2 \propto (Re - Re_c)$ → régression
      linéaire de $A^2(Re)$, extrapolation à $A=0$. Alternative : taux de croissance/décroissance d'une
      perturbation de part et d'autre de $Re_c$.
- [ ] **Sensibilité aux choix numériques** (au $Re$ le plus élevé) : maillage, schéma d'advection
      (centrée / QUICK / **upwind 1 qui doit échouer**), $\Delta t$/CFL, méthode M1 vs M2, éventuellement
      précision du Poisson. Pour chacun : stationnaire ou périodique ? fréquence ? amplitude ?
- [ ] Comparaison littérature : $Re_c$, fréquence, forme des structures oscillantes.

**Sortie** : diagramme de bifurcation $A(Re)$, $Re_c$ estimé ± incertitude, tableau de sensibilité.

### Phase 6 — Synthèse et rendu (1–2 jours)
- [ ] Figures finales homogènes (style du Projet 3 : couleur fixe par méthode/schéma).
- [ ] Rédaction : choix, validations, résultats, **explication des échecs éventuels**.
- [ ] Relecture « oral » : chaque fonction expliquable ligne à ligne.

---

## 5. Découpage du travail (lots indépendants)

Si on travaille à plusieurs (ou pour paralléliser les sessions avec l'IA), voici des lots aux interfaces
claires :

| Lot | Contenu | Dépend de | Interface |
|---|---|---|---|
| **L0** Biblio & données de référence | Ghia, Botella, Erturk, articles Hopf | — | `ref_data/*.csv` + fiche résumé |
| **L1** Poisson | SOR (Projet 1 adapté), DCT/DST, LU creux ; tests | — | `solve_poisson(rhs, bc) -> champ` |
| **L2** Noyau M1 (MAC) | grille, CL, opérateurs, pas de temps | L1 | `step_m1(state, dt, t)` |
| **L3** Noyau M2 ($\psi$–$\omega$) | grille, CL de vorticité, pas de temps | L1 | `step_m2(state, dt, t)` |
| **L4** Diagnostics & post-traitement | $K$, $P$, $\varepsilon$, $\psi$ depuis $(u,v)$, détection de tourbillons, FFT, Richardson | — (testable sur champs analytiques) | fonctions pures sur champs |
| **L5** Pilotage | `run_simulation`, critère stationnaire, sauvegarde/reprise `.npz`, runs parallèles | L2/L3, L4 | `run(config) -> results` |
| **L6** Études | notebooks Parties A, B, C | tout | figures + tableaux |

Le noyau commun **L1 + L4** est partagé par les deux méthodes : c'est ce qui garantit que les
comparaisons M1/M2 portent bien sur la méthode et pas sur le post-traitement.

---

## 6. Architecture de code proposée

```
Projet4_CFD/
├── PLAN.md                  ← ce document
├── AGENTS.md                ← conventions pour l'IA (style, numba, docstrings…)
├── environment.yml
├── cavity/                  ← code source (modules .py, importés par les notebooks)
│   ├── params.py            ← paramètres, adimensionnement, U_t(t)
│   ├── poisson.py           ← SOR, DCT/DST, LU creux
│   ├── mac_solver.py        ← méthode 1
│   ├── psi_omega_solver.py  ← méthode 2
│   ├── advection.py         ← schémas (centré, upwind1, QUICK/upwind3)
│   ├── diagnostics.py       ← K, P, ε, ψ, vortex, FFT
│   └── run.py               ← boucle en temps, critères d'arrêt, sauvegarde/reprise
├── tests/                   ← tests des briques (ordre, divergence, Poisson…)
├── ref_data/                ← données littérature (Ghia…)
├── results/                 ← sorties .npz (non versionnées si trop lourdes)
├── notebooks/
│   ├── A_Re100.ipynb
│   ├── B_sweep_Re.ipynb
│   └── C_unsteady.ipynb
└── figures/
```

Principes (dans la continuité des Projets 1 et 3) :
- Noyaux de calcul en **boucles `i, j` explicites + `@numba.njit`**, formules reconnaissables.
- Docstrings, commentaires en français.
- **Sauvegarde/reprise** des simulations longues (Partie C) : indispensable pour ne pas perdre des
  heures de calcul.
- Chaque résultat chiffré du rapport doit être reproductible par une cellule de notebook.

---

## 7. Calendrier indicatif

| Phase | Durée estimée | Cumul |
|---|---|---|
| 0 — Mise en place, biblio | 0.5 j | 0.5 j |
| 1 — Briques M1 | 1.5 j | 2 j |
| 2 — Partie A (M1) | 1 j | 3 j |
| 3 — M2 + validation croisée | 1.5 j | 4.5 j |
| 4 — Partie B | 2 j | 6.5 j |
| 5 — Partie C | 3.5 j | 10 j |
| 6 — Synthèse | 1.5 j | 11.5 j |

À recaler sur la date de rendu (question Q1). Les simulations longues de la Partie C peuvent tourner
en tâche de fond pendant les Phases 4 et 6.

---

## 8. Décisions de design à prendre

Pour chaque décision : options, et ma **recommandation** (★).

### D1 — Choix de la méthode 2
| Option | + | − |
|---|---|---|
| ★ **$\psi$–$\omega$ différences finies** | Très proche du cours, indépendante de M1, pas de pression, $P$ direct | CL de vorticité délicate (mais c'est un sujet de discussion intéressant) ; propre au 2D |
| Newton stationnaire (sur M1 discrétisé) | Convergence quadratique, suit la branche instable au-delà de $Re_c$ | Pas d'instationnaire ; jacobienne à assembler ; plus avancé |
| Lattice Boltzmann (D2Q9) | Paradigme très différent, simple à coder | Faiblement compressible, hors cours, $Re$ élevé délicat (stabilité) |
| Spectral Chebyshev | Référence ultra-précise | Difficile, singularité de coin pénalisante |

### D2 — Variantes à comparer *à l'intérieur* d'une méthode (pour le « choix sage vs malheureux »)
★ Advection **centrée ordre 2** (référence), **upwind 1** (contre-exemple de diffusion numérique),
**QUICK / upwind-biaisé 3** (repris du Projet 3). Intégration en temps : **RK3-SSP** (déjà codé) ou AB2.

### D3 — Solveur de Poisson
★ **DCT/DST rapide** comme solveur de production (exact, rapide) **+ SOR** du Projet 1 gardé pour la
comparaison de coût et parce qu'il est « du cours ». LU creux comme alternative si la DCT pose
problème (ex. si l'on change de géométrie).

### D4 — Format du code
★ **Modules `.py` + notebooks d'analyse** (les simulations longues de la Partie C se pilotent mal
depuis un notebook unique). Alternative : un notebook par partie comme aux Projets 1 et 3.

### D5 — Traitement de la puissance $P$ face à la singularité de coin
Options : (a) rapporter $P$ à grille fixée en montrant explicitement la dépendance en $\log N$ ;
(b) vérifier via $\varepsilon$ (même divergence, mais calcul indépendant) ; (c) en plus, isoler la
contribution « hors coins » (intégrale sur $x\in[\delta, 1-\delta]$) qui, elle, converge.
★ **(a) + (b) + (c)** : la loi d'échelle est tirée de la partie convergente, la partie singulière est
expliquée. Le couvercle reste **uniforme** comme l'impose l'énoncé (pas de couvercle régularisé, sauf
éventuellement en annexe).

### D6 — Démarrage des simulations
L'énoncé impose la rampe $U_t(t)$. Pour les états **stationnaires** (Parties A, B), le résultat final ne
devrait pas dépendre de l'histoire. ★ Rampe depuis le repos pour A et pour les runs « officiels »,
**continuation** depuis un $Re$ voisin autorisée en B pour gagner du temps (en le vérifiant sur un cas).
Pour C, toujours la rampe depuis le repos (le transitoire fait partie du problème), plus éventuellement
des redémarrages depuis un état convergé pour étudier la sensibilité.

### D7 — Grilles visées
★ A : $N = 32 \to 256$. B : $N = 128$ ($Re \le 1000$), $256$ ($Re \ge 3200$). C : $256$ et $384$ ou $512$
selon le budget mesuré en Phase 5.

### D8 — Environnement
L'environnement conda Windows `CFD` existant fonctionne (numba OK). ★ On le garde plutôt que de
passer sous WSL, sauf contrainte de l'enseignant.

---

## 9. Questions ouvertes (à me renseigner)

- **Q1.** Date de rendu, et format attendu (rapport PDF ? notebook ? soutenance orale ?) ;
  y a-t-il une limite de pages ?
- **Q2.** Travailles-tu **seul ou en groupe** ? (impacte le découpage §5)
- **Q3.** D1 : es-tu d'accord avec $\psi$–$\omega$ comme méthode 2 ? Y a-t-il une méthode que l'enseignant
  a particulièrement mise en avant en cours (projection, SIMPLE, $\psi$–$\omega$…) ?
- **Q4.** D4 : modules `.py` + notebooks, ou tout en notebooks comme avant ?
- **Q5.** Langue du code (commentaires) et du rapport : français comme aux Projets 1 et 3 ?
- **Q6.** Budget de calcul : peut-on laisser tourner des simulations plusieurs heures (voire une nuit)
  sur ta machine ?
- **Q7.** As-tu accès aux articles de référence (Ghia 1982, Bruneau & Saad 2006, Auteri 2002…) via la
  bibliothèque de l'école, ou faut-il partir de ce qui est disponible librement ?
- **Q8.** Niveau d'implication souhaité : veux-tu écrire toi-même certaines briques (comme le
  recommandait le Projet 1 pour SOR) et que je relise, ou que je propose le code et que tu le revoies ?
