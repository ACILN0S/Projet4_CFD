# CLAUDE.md — Projet 4 CFD : cavité entraînée à haut Reynolds

Projet de CFD (CentraleSupélec, 2026-27). Énoncé : [Doc/FinalProject_LidDrivenCavity_HighRe.pdf](Doc/FinalProject_LidDrivenCavity_HighRe.pdf).
Cadrage complet, phases et décisions : [PLAN.md](PLAN.md) — **le lire avant toute tâche non triviale**.
Ce fichier résume ce qu'il faut savoir en permanence ; en cas de conflit, l'énoncé prime, puis PLAN.md.

## 1. Le problème

- Cavité carrée 2D, NS incompressible seul (pas de scalaire). 3 parois fixes (no-slip), couvercle
  `y = L` à vitesse `U_t(t) = U (1 - exp(-t/τ))`, `τ = L/U`. Paramètre unique : `Re = UL/ν`.
- **Adimensionnement codé** : `L = U = ρ = 1` → `ν = 1/Re`, `U_t(t) = 1 - exp(-t)`.
- **Diagnostics obligatoires** (à suivre en continu) :
  - `K(t) = ½ ∫ (u² + v²) dS`
  - `P(t) = ∫₀¹ ν ∂u/∂y|_{y=1} U_t(t) dx`
  - ajouté : `ε(t) = ν ∫ ω² dS` ; **bilan `dK/dt = P − ε`** = contrôle de cohérence systématique.
- **Exigences de l'énoncé** :
  - ≥ 2 méthodes numériques, dont ≥ 1 issue du cours ; objectif = **comparer** et montrer l'effet d'un
    choix « sage » vs « malheureux ».
  - **Partie A** : stationnaire à `Re = 100`, validé (littérature + autre solveur), **grille justifiée**.
  - **Partie B** : plusieurs `Re < 7500`, évolution de la topologie, comparaison littérature, **`P(Re)` + loi d'échelle**.
  - **Partie C** : `Re > Re_c ≈ 8000`, plus haut `Re` abordable (**justifié** par un budget mesuré),
    **sensibilité aux choix numériques**, comparaison littérature.
  - Un échec est acceptable **s'il est diagnostiqué et expliqué**.

## 2. Méthodes retenues (propositions de PLAN.md §3, §8 — à confirmer par l'utilisateur)

| | Méthode 1 (« cours ») | Méthode 2 |
|---|---|---|
| Variables | `(u, v, p)` | `(ψ, ω)` |
| Grille | décalée **MAC** + cellules fantômes | colocalisée aux nœuds |
| Pas de temps | projection de Chorin (prédiction → Poisson → correction) | transport de ω puis `Δψ = −ω` |
| Parois | fantômes, ex. `u_fant = 2U_t − u_int` | ω de paroi : **Thom** (ordre 1) vs **Woods/Briley** (ordre 2) |
| Poisson | Neumann → **DCT** (production), SOR (comparaison/cours), `splu` en secours | Dirichlet → **DST**, SOR |

- **Advection** (comparées dans le même code) : centrée ordre 2 (référence), **upwind 1** (contre-exemple
  volontaire : `ν_num ≈ UΔx/2` ⇒ `Re_eff` ≪ `Re`, tue la bifurcation), QUICK / upwind-biaisé 3.
- **Temps** : RK3-SSP (déjà codé au Projet 3) ou AB2. `Δt = min(CFL·Δx/U, Fo·Re·Δx²)`.
- **Grilles visées** : A `N = 32 → 256` ; B `N = 128` (Re ≤ 1000), `256` (Re ≥ 3200) ; C `256`, `384`/`512` selon budget.
- **Démarrage** : rampe depuis le repos pour A, C et les runs « officiels » ; continuation entre Re voisins
  autorisée en B seulement (vérifiée sur un cas).

Ne pas trancher seul une décision D1–D8 ouverte de PLAN.md §8 : proposer, puis demander.

## 3. Pièges physiques/numériques connus (à garder en tête à chaque résultat)

1. **Singularité des coins du couvercle** : ω ~ 1/r ⇒ `P` diverge en `log N`. Ne jamais présenter `P` comme
   « convergé en maillage » sans le dire. Rapporter : `P` à N fixé + dépendance en log N, vérif via ε, et
   contribution hors coins `x ∈ [δ, 1−δ]` (qui converge) pour la loi d'échelle.
2. **Ordre des CL** : une CL d'ordre 1 détruit l'ordre global (leçon du Projet 1). Toujours mesurer l'ordre
   **global**, pas seulement intérieur.
3. **Coût à haut Re** : couches limites en `Re^{-1/2}`, 10⁵–10⁶ pas en Partie C ⇒ Numba + Poisson rapide
   indispensables, **sauvegarde/reprise `.npz`** obligatoire pour les runs longs.
4. **Ralentissement critique près de Re_c** : simuler plutôt 8500–10000 et extrapoler `Re_c` via `A² ∝ (Re − Re_c)`.
5. Les valeurs de `Re_c` de la littérature divergent (~7400–8050) : c'est un résultat à discuter, pas un bug.

## 4. Références à utiliser pour la validation

Ghia, Ghia & Shin (1982) — profils médians + tourbillons, Re = 100…10000 ; Botella & Peyret (1998) — Re = 1000
de référence ; Erturk (2005) ; bifurcation : Auteri et al. (2002), Bruneau & Saad (2006), Poliashenko & Aidun (1995).
Données saisies dans `ref_data/*.csv` (source + tableau cités en en-tête). **Ne jamais inventer une valeur de
référence** : si elle n'est pas vérifiée, le signaler explicitement.

## 5. Architecture du code (PLAN.md §6)

```
cavity/          params.py, poisson.py, mac_solver.py, psi_omega_solver.py,
                 advection.py, diagnostics.py, run.py
tests/           tests des briques (ordre, divergence, Poisson, bilan d'énergie…)
ref_data/        données littérature (CSV)
results/         sorties .npz (lourdes → non versionnées)
notebooks/       A_Re100.ipynb, B_sweep_Re.ipynb, C_unsteady.ipynb
figures/
Projet_1_Pressure/, Projet_3_advection/   anciens TP à réutiliser (lecture seule)
```

Interfaces entre lots (PLAN.md §5) : `solve_poisson(rhs, bc)`, `step_m1(state, dt, t)`, `step_m2(state, dt, t)`,
diagnostics = fonctions pures sur champs, `run(config) -> results`. **L1 (Poisson) et L4 (diagnostics) sont
partagés** par M1 et M2 : la comparaison doit porter sur la méthode, pas sur le post-traitement.

### Réutilisation des projets précédents
- Projet 1 (`Projet_1_Pressure/TP.ipynb`) : SOR/Jacobi/GS + version Numba (cellules 13, 17), ordre des CL (cellule 26), Richardson/GCI.
- Projet 3 (`Projet_3_advection/TP.ipynb`) : upwind 1/2/3, RK2, RK3-SSP (cellules 15, 17, 19), von Neumann,
  `apply_boundary_conditions`, `run_simulation → dict`, style de figures (cellules 1, 12, 28).
- Partir de ce code plutôt que réécrire ; citer la cellule d'origine en commentaire.

## 6. Conventions de code

- **Langue** : docstrings, commentaires et rapport en **français** ; identifiants en anglais.
- **Noyaux de calcul** : boucles `i, j` explicites sous `@numba.njit`, formules reconnaissables (lien
  direct avec les équations du rapport). Pas de vectorisation « astucieuse » illisible dans les noyaux.
- **Exception assumée à la règle d'immutabilité** des règles globales : les noyaux numériques écrivent
  **en place** dans des tableaux préalloués (performance, Numba). L'immutabilité s'applique au reste
  (configs, dicts de résultats, post-traitement).
- Indexation : documenter dans la docstring la convention `[i, j]` ↔ `(x, y)` et la position MAC de chaque
  champ (centres / faces / nœuds, avec fantômes).
- Paramètres physiques et numériques regroupés dans une config (`params.py`), pas de nombres magiques.
- Chaque résultat chiffré du rapport doit être **reproductible par une cellule de notebook**.
- L'utilisateur doit **maîtriser chaque ligne** (oral) : privilégier la clarté, expliquer les choix non
  évidents, éviter les dépendances exotiques. Bibliothèques autorisées : numpy, scipy, numba, matplotlib.
- Figures : couleur fixe par méthode/schéma (style du Projet 3), axes et unités adimensionnées étiquetés.
- Fichiers < 800 lignes, fonctions courtes ; pas de `print` de debug laissé dans le code.

## 7. Vérification et tests

Ici, la « couverture » qui compte est **physique et numérique**, pas un pourcentage de lignes :
- Opérateurs sur champs analytiques → **ordre 2 mesuré** (pente sur 3+ grilles).
- Poisson sur solution analytique (Neumann et Dirichlet) → ordre 2 + erreur machine pour la DCT/DST.
- Après projection : `max |∇·u| ≲ 1e-12`.
- Bilan `dK/dt = P − ε` vérifié sur chaque run (résidu rapporté).
- MMS (solutions manufacturées) pour l'ordre espace/temps du solveur complet si le temps le permet.
- Validation Ghia/Botella avec écarts **quantifiés** ; étude de maillage avec Richardson/GCI.
- TDD pour les briques (`tests/`, pytest) : écrire le test de convergence avant le noyau.

**Exécution** (environnement conda Windows `CFD` : numpy 2.5, scipy 1.18, numba 0.68, matplotlib 3.11) :

```bash
conda run -n CFD python -m pytest tests -q
```

pytest 9.2 est installé dans `CFD` (conda-forge).
Machine 32 cœurs logiques → balayages en Re lancés en parallèle (processus séparés).
Runs longs (Partie C) : lancer en arrière-plan, avec sauvegarde périodique et reprise.

## 8. Workflow avec Claude

- Respecter les **phases et critères de sortie** de PLAN.md §4 ; cocher les cases au fur et à mesure.
- Proposer le pseudo-code avant d'implémenter une méthode ; l'utilisateur peut vouloir écrire certaines
  briques lui-même (Q8) — demander avant de tout coder.
- Pour chaque résultat : comparer à l'attendu physique (PLAN.md §1.4) et signaler tout écart plutôt que
  l'ajuster silencieusement.
- Consigner les décisions tranchées (D1–D8) en ADR dans `docs/adr/`.
- Questions ouvertes (PLAN.md §9) : date/format de rendu, seul ou en groupe, choix M2, budget de calcul,
  accès aux articles, niveau d'implication.

### Skills ECC disponibles localement (`.claude/skills/`)

| Skill | Usage dans ce projet |
|---|---|
| `python-patterns`, `coding-standards` | style Python (hors noyaux Numba : voir §6) |
| `python-testing`, `tdd-workflow` | tests pytest des briques, tests de convergence paramétrés |
| `verification-loop` | vérification avant de clore une phase (tests + diff + bilan) |
| `benchmark-optimization-loop` | optimiser Poisson / pas de temps avec garde-fou de correction (Phase 5 : budget) |
| `search-first`, `documentation-lookup` | vérifier l'API scipy (`dctn`, `splu`…), numba avant de coder |
| `scientific-thinking-literature-review`, `deep-research` | bibliographie Ghia, Botella, Erturk, Re_c |
| `scientific-thinking-scholar-evaluation` | relecture critique du rapport et des sources |
| `architecture-decision-records` | tracer les décisions D1–D8 |
| `strategic-compact` | compacter le contexte entre phases longues |
