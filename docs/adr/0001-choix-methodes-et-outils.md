# ADR 0001 — Choix des méthodes et des outils (décisions D1–D8)

- **Statut** : accepté (2026-11-06, fin de la Phase 0)
- **Contexte** : PLAN.md §8 proposait une recommandation (★) pour chaque décision D1–D8.
  L'utilisateur a demandé de clore la Phase 0 et de lancer la Phase 1 : les recommandations sont
  retenues telles quelles, dans l'esprit « simplicité » de CLAUDE.md §0.

| # | Décision retenue | Pourquoi (en une ligne) |
|---|---|---|
| D1 | Méthode 2 = **ψ–ω différences finies**, grille colocalisée | Proche du cours, indépendante de M1 (variables, grille, CL, Poisson) → validation croisée crédible |
| D2 | Advection **centrée ordre 2** (référence) vs **upwind 1** (contre-exemple) ; temps **RK3-SSP** | Le contraste « sage / malheureux » demandé par l'énoncé ; RK3 déjà codé au Projet 3. QUICK seulement si besoin |
| D3 | Poisson : **DCT (Neumann) / DST (Dirichlet)** en production + **SOR** du Projet 1 pour comparaison | Exact et rapide ; SOR = méthode du cours et contre-exemple de coût |
| D4 | **Modules `.py`** (`cavity/`) + **notebooks** d'analyse | Les longs runs de la Partie C se pilotent mal depuis un notebook |
| D5 | P rapporté à N fixé + dépendance en log N, vérifié par ε, plus la contribution hors coins | La singularité de coin rend P non convergé en maillage ; on l'explique au lieu de le cacher |
| D6 | Rampe depuis le repos (A, C, runs officiels) ; continuation entre Re voisins autorisée en B | Gain de temps en B, vérifié sur un cas |
| D7 | Grilles : A `N = 32 → 256` ; B `128` / `256` ; C `256` + `384`/`512` selon budget | Couches limites en Re^(-1/2) |
| D8 | Environnement conda Windows (pas de WSL) | Fonctionne avec Numba |

**Note D8 (constat du 2026-11-06)** : l'environnement `CFD` décrit dans `environment.yml` n'existe pas
sur cette machine ; on utilise l'environnement `base` d'Anaconda (numpy 1.24, scipy 1.10, numba 0.57,
pytest 7.4). Le code n'utilise que des API présentes dans les deux jeux de versions.

**Conséquences** : la Phase 1 code M1 (MAC + projection + RK3) avec les deux schémas d'advection,
le module Poisson commun (DCT, DST, SOR) et les diagnostics K, P, ε.
