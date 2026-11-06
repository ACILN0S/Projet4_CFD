# Données de référence (littérature)

Toutes les grandeurs sont adimensionnées avec $L = U = 1$ (comme notre code).
Chargement : [`cavity/ref_data.py`](../cavity/ref_data.py) (`load_profile`, `load_vortices`).
Contrôles de saisie : [`tests/test_ref_data.py`](../tests/test_ref_data.py).

## Fichiers

| Fichier | Contenu | Re | Source | Statut |
|---|---|---|---|---|
| `ghia1982_u_vertical_centerline.csv` | $u(0.5, y)$, 17 points | 100, 400, 1000, 3200, 5000, 7500, 10000 | Ghia et al. (1982), Table I | transcription secondaire (gist), **non revérifiée sur l'original** |
| `ghia1982_v_horizontal_centerline.csv` | $v(x, 0.5)$, 17 points | idem | Ghia et al. (1982), Table II | idem |
| `erturk2005_u_vertical_centerline.csv` | $u(0.5, y)$, 23 points | 1000, 2500, 5000, 7500, 10000 … 21000 | Erturk et al. (2005), Table 6 | extrait directement du PDF arXiv |
| `erturk2005_v_horizontal_centerline.csv` | $v(x, 0.5)$, 23 points | idem | Erturk et al. (2005), Table 7 | idem |
| `vortices.csv` | $\psi$, $\omega$, $(x, y)$ au centre des tourbillons (principal, BR1, BL1, TL1, BR2…) | 1000 … 21000 | Erturk Tables 1, 4, 5 (+ Ghia et Botella **cités par** Erturk) | voir colonne `note` |
| `bifurcation_literature.csv` | $Re_c$ (et pulsation si connue) | — | plusieurs, voir colonne `status` | **en partie non vérifié** |
| `sources/` | textes bruts d'origine (gists Ghia, pages 27–33 d'Erturk via `pdftotext -raw`) + script de conversion `build_ref.py` | | | traçabilité |

## Erreurs connues dans les tables de Ghia (1982)

Conservées telles quelles dans les CSV, masquées (NaN) par défaut au chargement (`GHIA_SUSPECT`) :

| Table | Re | Point | Valeur imprimée | Pourquoi suspecte |
|---|---|---|---|---|
| I ($u$) | 3200 | $y = 0.4531$ | −0.86636 | incompatible avec les voisins (−0.04272, −0.24427) ; débit $\int u\,dy$ passe de 0.003 à −0.083. Probablement −0.08636. |
| I ($u$) | 10000 | $y = 0.5$ | +0.03111 | signe opposé à toutes les autres colonnes et à Erturk ($u = -0.0287$). Probablement −0.03111. |
| II ($v$) | 400 | $x = 0.9063$ | −0.23827 | rompt la monotonie (−0.22847 → −0.44993) ; signalée douteuse par la transcription d'origine. |

On ne « corrige » pas ces valeurs (pas de source pour la bonne valeur) : on les exclut.

## Ce qui manque encore (à saisir si on accède aux articles originaux)

- **Ghia (1982), Tables III–V** : $\psi$, $\omega$ et positions des tourbillons pour $Re = 100$ et $400$
  (indispensable pour la Partie A !). Non disponibles en source libre vérifiable → **à saisir depuis
  l'article** (bibliothèque de l'école, Q7 de PLAN.md). En attendant, la validation à $Re = 100$ se fait
  sur les profils médians (Tables I–II), qui sont la comparaison standard.
- **Botella & Peyret (1998)** : tableaux de profils à $Re = 1000$ (référence spectrale la plus précise).
  Seules les valeurs du tourbillon principal sont disponibles (via Erturk Table 1).
- **Bifurcation** : vérifier sur les originaux $Re_c$, la pulsation/période du cycle limite, et les
  positions de sondes utilisées (Auteri 2002, Bruneau & Saad 2006, Sahin & Owens 2003).

## Remarques pour l'interprétation

- **Erturk et al. (2005)** obtiennent des solutions *stationnaires* jusqu'à $Re = 21000$ (solveur
  stationnaire, grille 601²). Cela ne contredit pas $Re_c \approx 8000$ : au-delà de $Re_c$ la branche
  stationnaire existe mais est **instable** ; un solveur instationnaire ne peut pas y rester.
  Ces profils au-delà de 8000 ne sont donc **pas** une référence pour nos simulations instationnaires
  de la Partie C (seulement pour la branche stationnaire instable).
- Erturk note aussi qu'à grille grossière (129²) il n'obtient pas de solution stationnaire dès
  $Re = 7500$ : la dispersion des $Re_c$ publiés (7402 → 8031) dépend fortement de la résolution
  — matière directe pour la Partie C (sensibilité aux choix numériques).
- Convention de signe (Erturk) : $\psi < 0$ et $\omega < 0$ dans le tourbillon principal (rotation
  horaire pour un couvercle allant vers $+x$), avec $u = \partial_y\psi$, $v = -\partial_x\psi$,
  $\omega = \partial_x v - \partial_y u$. **Utiliser la même convention dans notre code.**

## Bibliographie

1. U. Ghia, K. N. Ghia, C. T. Shin, *High-Re solutions for incompressible flow using the Navier–Stokes
   equations and a multigrid method*, J. Comput. Phys. **48** (1982) 387–411. doi:10.1016/0021-9991(82)90058-4
2. O. Botella, R. Peyret, *Benchmark spectral results on the lid-driven cavity flow*, Comput. Fluids
   **27** (1998) 421–433. doi:10.1016/S0045-7930(98)00002-4
3. E. Erturk, T. C. Corke, C. Gökçöl, *Numerical solutions of 2-D steady incompressible driven cavity
   flow at high Reynolds numbers*, Int. J. Numer. Meth. Fluids **48** (2005) 747–774. arXiv:cs/0411047
4. F. Auteri, N. Parolini, L. Quartapelle, *Numerical investigation on the stability of singular driven
   cavity flow*, J. Comput. Phys. **183** (2002) 1–25.
5. C.-H. Bruneau, M. Saad, *The 2D lid-driven cavity problem revisited*, Comput. Fluids **35** (2006) 326–348.
6. M. Poliashenko, C. K. Aidun, *A direct method for computation of simple bifurcations*, J. Comput.
   Phys. **121** (1995) 246–260.
7. Y.-F. Peng, Y.-H. Shiau, R. R. Hwang, *Transition in a 2-D lid-driven cavity flow*, Comput. Fluids
   **32** (2003) 337–352.
8. M. Sahin, R. G. Owens, *A novel fully-implicit finite volume method applied to the lid-driven cavity
   problem. Part II: Linear stability analysis*, Int. J. Numer. Meth. Fluids **42** (2003) 79–88.
9. A. Fortin, M. Jardak, J. J. Gervais, R. Pierre, *Localization of Hopf bifurcations in fluid flow
   problems*, Int. J. Numer. Meth. Fluids **24** (1997) 1185–1210.
10. N. Cardoso, P. Bicudo, *Time dependent simulation of the driven lid cavity at high Reynolds number*,
    arXiv:0809.3098 (revue utile des trois « camps » : stationnaire / stabilité / DNS).

Références bibliographiques (volumes, pages) tirées des listes de références d'Erturk (2005) et de
Cardoso & Bicudo (2009), sauf Bruneau & Saad (de mémoire, **à vérifier**).
