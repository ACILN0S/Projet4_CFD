# CLAUDE.md — Projet 4 CFD : cavité entraînée à haut Reynolds

Énoncé : [Doc/FinalProject_LidDrivenCavity_HighRe.pdf](Doc/FinalProject_LidDrivenCavity_HighRe.pdf).
Plan : [PLAN.md](PLAN.md). Pseudo-code : [docs/pseudocode.md](docs/pseudocode.md).
Ordre de priorité : énoncé > PLAN.md > ce fichier.

## 0. Deux règles qui priment sur tout

1. **Strict minimum.** Faire seulement ce que l'énoncé exige : 2 méthodes, Parties A/B/C, K et P.
   Pas de bonus (tests unitaires, bilan d'énergie, Richardson, FFT, sauvegarde/reprise…)
   sauf demande explicite. Avant d'ajouter quoi que ce soit : « l'énoncé l'exige-t-il ? »
2. **Uniquement des idées du cours** (poly *Numerical Methods*, Vicquelin & Vié 2025-26).
   Pas de grille décalée MAC, pas de cellules fantômes, pas de DCT/DST, pas de RK3-SSP/QUICK.
   Chaque choix numérique cite sa section du cours dans un commentaire. En cas de doute : demander.

L'utilisateur (et son équipe) doit pouvoir réexpliquer chaque ligne à l'oral.

## 1. Le problème

- Cavité carrée 2D, NS incompressible, couvercle `U_t(t) = 1 - exp(-t)` (adimensionné : `L = U = ρ = 1`, `ν = 1/Re`).
- Diagnostics : `K(t) = ½ ∫ (u² + v²) dS` et `P(t) = ∫₀¹ ν ∂u/∂y|_{y=1} U_t dx`. Rien d'autre.

## 2. Méthodes retenues (cf. PLAN.md §2)

- Grille **colocalisée** aux nœuds `(N+1, N+1)`, `[i, j]` ↔ `(x, y)`, bords sur les parois.
- **M1** : projection / pas fractionnaire (cours §7.3.2, workshop 15). **M2** : ψ–ω (cours §7.3.4),
  ω de paroi par Taylor (Thom).
- Euler explicite (§3.1) ; diffusion centrée ordre 2 ; advection **centrée** ou **upwind 1**
  (choix sage / malheureux, expliqué par l'équation modifiée §5.4.2) ; repli RK3 (§3.9.4) en Partie C si besoin.
- Poisson : gradient conjugué sans matrice (§4.7.3, workshop 11.3), départ = solution du pas précédent.
  Pression : Neumann partout, moyenne du second membre retirée (compatibilité) et `p -= p.mean()` à chaque
  itération (exigé par les slides du prof).
- `Δt` = limites de stabilité du cours (§5.3) : Fourier 2D, CFL upwind, `C² ≤ 2Fo` pour le centré.

## 3. Code

- Fichiers : `cavity/common.py`, `cavity/projection.py`, `cavity/psi_omega.py` + 3 notebooks (A, B, C).
- **Numba** (`@njit(cache=True)`, recommandé par le prof ; décision du 2026-11-11) pour les noyaux coûteux
  (opérateurs, gradient conjugué), en boucles `for i, j` explicites qui reprennent les formules du cours ;
  pas de `fastmath` (mêmes résultats qu'en numpy aux arrondis près). numpy pour le reste.
  Bibliothèques : numpy, numba, scipy, matplotlib. Les calculs indépendants (plusieurs N, plusieurs Re)
  tournent en parallèle, un processus par cas.
- Environnement CFD : `MKL_THREADING_LAYER=SEQUENTIAL` à fixer avant d'importer numpy (sinon toute
  algèbre linéaire numpy, donc matplotlib, fait planter le processus sur cette machine).
- Français pour docstrings, commentaires et rapport ; identifiants en anglais.
- **Commentaires abondants** : en-tête de fichier (rôle, méthode, partie) ; docstring par fonction
  (équation discrétisée, formes des tableaux, section du cours) ; un commentaire *quoi/pourquoi* par bloc.
  Dans les notebooks : cellule markdown avant chaque calcul (ce qu'on fait, ce qu'on attend) et
  interprétation courte après chaque figure.
- Vérifications : directement dans les notebooks (Poisson sur solution connue, comparaison Ghia),
  pas de dossier de tests.
- Lancer avec l'environnement conda `CFD` (`~/miniconda3/envs/CFD/python.exe` sous Windows).

## 4. Références (strict minimum)

- Ghia, Ghia & Shin (1982) : `ref_data/ghia1982_*.csv`, lus avec `np.loadtxt(..., delimiter=",")`.
  Grilles 129² (Re ≤ 3200) / 257² (Re ≥ 5000), ordre 2, sans étude de convergence chiffrée ; ni K ni P.
- Erturk et al. (2005) : étude de maillage (401/513/601 + Richardson), Thom ⇒ ordre 2 global.
  Articles dans `doc_ref/`, résumé dans `ref_data/README.md`.
- `Re_c` ≈ 8000 (fourchette publiée 7400–8050) : `ref_data/README.md`.
- Ne jamais inventer une valeur de référence.

## 5. À garder en tête

- Coins du couvercle : saut de vitesse ⇒ `P` diverge en continu et croît en `ln N` en discret ; comparer
  les `P` à `N` fixé ; pente de `P(Re)` biaisée vers plus fort que −1/2 (terme de coin ∝ 1/Re). Détail : PLAN.md §2.
- Upwind à haut `Re` : `D = ½ U h (1 − C) ≫ ν` ⇒ `Re` effectif faible ⇒ pas d'oscillation attendue.
- Pour chaque résultat, comparer à l'attendu physique et signaler un écart plutôt que l'ajuster.
