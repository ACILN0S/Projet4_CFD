# ADR 0004 — Optimisation du code de calcul (version de travail, non livrée)

- **Statut** : accepté (2026-11-13). Version **de travail** pour accélérer les calculs des Parties B et C ;
  **la version livrée est celle d'avant l'optimisation**.
- **Sauvegarde de la version d'avant** : branche git `sauvegarde/avant-optim-2026-11-13` (et tag du même nom),
  commit `6bc3206` « partie B » ; archive `Projet4_CFD_sauvegarde_2026-11-13.zip` (dossier parent du projet,
  avec les `.npz` de `results/`). Revenir à cette version : `git checkout sauvegarde/avant-optim-2026-11-13 -- cavity/`.
- **Code concerné** : `cavity/common.py`, `cavity/projection.py`, `cavity/psi_omega.py`.

## 1. Pourquoi

À Re = 5000, N = 128, le schéma centré impose Δt = 0.8 × 2ν = 3.2·10⁻⁴ (ADR 0003 §1) : il faut ~10⁶ pas pour
atteindre le stationnaire, et le gradient conjugué fait de 50 à 270 itérations par pas. Dans l'ancienne version,
M2 à Re = 5000, N = 128 a demandé ~18 h (17 calculs en parallèle sur le PC). La Partie B étoffée (11 Re, 46 cas)
et la Partie C ne sont pas faisables à ce coût.

## 2. Règle suivie : aucune méthode ni formule ne change

Même grille, mêmes schémas (Euler explicite, diffusion centrée, advection centrée / upwind, Thom, Neumann à nœuds
fantômes), même gradient conjugué (même algorithme, même tolérance, même départ à chaud, moyenne nulle de p à chaque
itération), même Δt, même critère d'arrêt, mêmes K et P. Seule change **la manière d'exécuter** ces calculs.

## 3. Ce qui a été fait

1. **Un pas de temps = un seul noyau compilé** (`_projection_step_kernel`, `_psi_omega_step_kernel`).
   Avant : chaque pas enchaînait en Python ~15 opérations numpy (copies, tranches, `np.gradient`, `np.trapezoid`…),
   chacune allouant un tableau temporaire. Maintenant : boucles `for i, j` compilées par Numba, qui reprennent les
   mêmes formules (fonctions `laplacian_at`, `advection_at`, écrites une fois et partagées par M1, M2 et les
   fonctions publiques `laplacian`, `advection`).
2. **Double tampon, aucune allocation par pas.** Les champs « temps t » et « temps t + Δt » sont deux jeux de tableaux
   alloués une fois, échangés après chaque pas (`u, u_new = u_new, u`). Les tableaux de travail du gradient conjugué
   (second membre, résidu, direction, A·direction) sont alloués une fois par calcul (`PoissonWorkspace`).
3. **Gradient conjugué en 3 balayages de la grille par itération au lieu de 6** :
   - le produit A·d et le produit scalaire d·(A d) dans la même boucle (`_apply_matrix_dot`) ;
   - mise à jour de φ et de r, r·r et somme pondérée de φ (pour sa moyenne) dans la même boucle ;
   - nouvelle direction d.
   La moyenne de φ calculée à une itération est retirée au début de la mise à jour suivante (et après la dernière) :
   c'est exactement le même calcul, (φ − m) + α d, puisque rien ne lit φ entre les deux.
4. **Sommes vectorisées** : `fastmath={"reassoc"}` sur le gradient conjugué **seulement**, et ce seul drapeau
   (pas de `fastmath=True`). Il autorise le compilateur à faire les sommes (r·r, d·Ad, moyenne) dans un autre ordre
   que i, j croissants, donc à additionner plusieurs termes à la fois (instructions vectorielles) au lieu d'attendre
   chaque addition. Chaque opération reste exacte au sens IEEE ; seul l'ordre des additions change.
5. **Boucles de mise à jour sur les tableaux aplatis** (vues 1D, sans copie) : le compilateur ne vectorise ces boucles
   qu'en 1D. Elles parcourent tous les nœuds ; aux nœuds inactifs (parois en Dirichlet), le poids de la moyenne vaut 0
   et r, A·d y restent nuls, donc rien ne change.
6. **Diagnostics compilés** : résidu stationnaire, K et P en boucles Numba (plus de temporaires numpy).
7. **Historique** (t, K, P, ||R||, itérations) dans des tableaux numpy agrandis par doublement, au lieu de listes
   Python (×4 en mémoire pour 10⁶ pas). Options ajoutées à `run_projection` / `run_psi_omega` (désactivées par
   défaut) : `record_every` (profils médians enregistrés au cours du temps, pour l'étude du critère d'arrêt) et
   `log_every` (ligne d'avancement).

## 4. Ce qui a été essayé et écarté

- **Multithreading Numba (`prange`)** : à N ≤ 128, synchroniser les threads coûte plus que le calcul (45 µs par
  boucle contre 5 µs en séquentiel avec la couche par défaut ; gain ×1.75 seulement avec OpenMP à 4 threads).
  On garde un processus par cas, en parallèle.
- **Préconditionnement du gradient conjugué, multigrille** : changeraient la méthode (hors cours) ; non faits.
- **LU creuse (SuperLU) à la place du gradient conjugué** : utilisée **temporairement** pour le seul calcul long de
  calibration du critère d'arrêt (Re = 5000, `scripts/calibration_Re5000.py`, ~5 fois plus rapide ; écart au gradient
  conjugué ~10⁻⁹ sur u, v, vérifié), puis **retirée du code**.

## 5. Piège rencontré : le cache de Numba

`cache=True` garde le code compilé sur disque, mais Numba ne surveille que le fichier où la fonction est écrite.
Un noyau de `projection.py` qui appelle `common.py` réutilisait donc l'ancien code compilé après une modification de
`common.py`, **sans avertissement** (constaté pendant les mesures). Les noyaux de `projection.py` et `psi_omega.py`
qui appellent `common.py` sont donc en `cache=False` : 1 à 3 s de compilation au lancement de chaque processus.

## 6. Vérification : mêmes résultats que la version d'avant

Comparaison directe, version d'avant contre version optimisée, mêmes paramètres :

| Cas (3000 pas, Re = 5000 sauf mention) | max \|Δu\|, \|Δv\| | max \|Δp\| ou \|Δω\| | itérations du gradient conjugué |
| --- | --- | --- | --- |
| M1, N = 128, centré | 2·10⁻¹⁶ | p : 6·10⁻¹⁵ | identiques à chaque pas |
| M2, N = 128, centré | 2·10⁻¹⁵ | ω : 2·10⁻¹³ | identiques |
| M1, N = 64, centré / upwind | 3·10⁻¹⁶ | p : 2·10⁻¹⁴ | identiques |
| M1 Neumann d'ordre 1 (Re = 100, N = 32) ; M2 upwind (Re = 1000, N = 32) | 1·10⁻¹⁵ | ω : 6·10⁻¹⁴ | identiques |
| `solve_poisson` (Dirichlet, Neumann d'ordre 1, fantômes) | φ : 2·10⁻¹⁵ | — | identiques |
| Calcul complet M1, Re = 100, N = 32 (critère ε = 10⁻²) | 2·10⁻¹⁰ | — | arrêt au **même pas** (t = 18.45) |

Les écarts sont des écarts d'arrondi (10⁻¹⁶ à 10⁻¹³), qui grossissent lentement au fil des pas jusqu'au niveau de la
tolérance du gradient conjugué (10⁻¹⁰ après 1500 pas), très loin des écarts à Ghia (~10⁻²). Les résultats calculés
avec la version optimisée valent donc pour la version livrée. Même vérification refaite sur un nœud de Ruche
(Intel Xeon Gold 6230) : mêmes écarts, mêmes itérations.

## 7. Gain mesuré

PC portable (AMD Ryzen 9 8940HX, sur batterie), Re = 5000, 3000 premiers pas (transitoire : beaucoup d'itérations) :

| Cas | avant (ms / pas) | après (ms / pas) | accélération |
| --- | --- | --- | --- |
| M1, N = 128 | 11.5 | 3.4 – 3.8 | ×3.1 – 3.4 |
| M2, N = 128 | 7.1 | 3.3 | ×2.2 |
| M1, N = 64 | 1.7 | 0.58 | ×3.0 |

Les mesures sur batterie varient de ±10 %. Sur Ruche, un cœur est ~3.8 fois plus lent que ce PC (12.7 ms / pas pour
M2, N = 128) : Ruche sert à lancer beaucoup de cas en même temps, pas à accélérer un cas.

## 8. Ce qu'il faut savoir pour la suite

- Calculs faits avec la version optimisée : calibration du critère d'arrêt (Re = 5000, `results/C5000_*.npz`, avec
  la LU temporaire), Partie B étoffée sur Ruche (`results/B_*_res1e-03_tref30.npz`, 11 Re, N = 64, 128 et 256 à
  Re = 5000).
- Pour livrer : remettre `cavity/` dans l'état de la sauvegarde (commande en tête de ce fichier). Les notebooks
  n'utilisent que `run_projection`, `run_psi_omega`, `solve_poisson` et `time_step`, dont l'interface n'a pas changé
  (seules options ajoutées : `record_every`, `log_every`, à ne pas utiliser dans les notebooks livrés).
- `CLAUDE.md` §3 décrit les règles Numba de la version optimisée (`reassoc`, `cache=False`) : à remettre comme avant
  avec le code livré.

## 9. Retour à la version simple (2026-11-19)

- `cavity/` et `CLAUDE.md` remis dans l'état de la sauvegarde (`git checkout HEAD -- cavity/ CLAUDE.md`), sur décision de
  l'équipe : on livre et on présente le code simple, y compris pour la Partie C (2 à 3 fois plus lent sur Ruche).
- Seul changement reporté depuis la version optimisée : le résidu de M2 mesuré sur ω (`common.steady_residual` et
  `march_to_steady(..., residual_fields)`, ADR 0003 §2), et les valeurs par défaut du critère (ε = 1e-3, t_ref = 30, t_max = 1000).
- Vérifié : le code simple redonne les calculs de la Partie B (faits avec la version optimisée) à 2·10⁻⁹ près sur u, v
  (M1 et M2, Re = 400, N = 64 ; arrêt à 1 ou 2 pas près).
- La version optimisée et les scripts de diagnostic (calibration, Woods, Δt) sont archivés dans
  `Projet4_CFD_version_optimisee_2026-11-19.zip` (dossier parent du projet), avec l'ancienne version longue du notebook B.
