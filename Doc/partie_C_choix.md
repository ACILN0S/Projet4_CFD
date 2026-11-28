# Partie C — choix, valeurs retenues et premiers résultats

> **Version périmée (24/11).** Les résultats finaux (tous les calculs, jusqu'à $t = 2000$ et à $N = 256$) sont dans `Doc/partie_C_resultats.md` (25/11).

Document de travail pour l'équipe, état au **24 novembre 2026, 6 h 40**. Il explique ce qui a été décidé pour la Partie C, avec quelles valeurs, pourquoi, et ce que montrent les calculs déjà terminés. Les chiffres viennent du notebook `notebooks/C_unsteady.ipynb` (ré-exécuté le 24/11) et des journaux `results/logs/C_*.log`. Les deux calculs à $N = 256$ tournent encore : leurs résultats ne sont pas dans ce document.

## 1. Ce que demande l'énoncé

L'énoncé (Partie C) rappelle qu'au-delà d'un Reynolds critique $Re_c$ (proche de 8000 dans la littérature), la solution stationnaire de la cavité entraînée devient instable et l'écoulement devient périodique. Il demande :

- de simuler l'écoulement au $Re$ **le plus élevé possible**, en justifiant ce choix ;
- d'étudier la **sensibilité** des résultats, à ce $Re$, aux choix numériques ;
- de **comparer** à la littérature ou à un autre solveur validé.

Les diagnostics restent les mêmes qu'en Parties A et B : l'énergie cinétique $K(t) = \tfrac12 \int (u^2 + v^2)\,dS$ et la puissance injectée par le couvercle $P(t) = \int_0^1 \nu\, \partial_y u|_{y=1}\, U_t\, dx$, avec $U_t(t) = 1 - e^{-t}$ et $\nu = 1/Re$.

Le cadre reste celui de `CLAUDE.md` : strict minimum, uniquement des méthodes du cours, chaque choix justifié par une section du poly. Le plan de la Partie C est l'étape 5 de `PLAN.md`.

## 2. Les dix choix et ce qui a été retenu

La veille, dix choix ont été posés à l'équipe. Le tableau résume les réponses ; les sections suivantes les justifient.

| Choix | Question | Décision |
| --- | --- | --- |
| 1 | Quel Reynolds ? | $Re = 10\,000$ |
| 2 | Plusieurs $Re$ ? | Un seul pour l'instant ; d'autres plus tard peut-être |
| 3 | Quelle grille ? | $N = 128$ d'abord, puis $N = 256$ |
| 4 | Quelles méthodes ? | M1 (projection) et M2 ($\psi$–$\omega$) |
| 5 | Schéma en temps / advection | Euler explicite, advection centrée (+ upwind comme contre-exemple) |
| 6 | Tests de sensibilité numérique | $\Delta t/2$ sur M2 seulement ; tolérance du gradient conjugué non testée |
| 7 | Condition initiale | Départ du repos, comme en Parties A et B |
| 8 | Mesure de la période | Méthode simple des maxima (sans FFT) |
| 9 | Que sauvegarder ? | Tout ce qui est pertinent : historiques, sondes, instantanés |
| 10 | Où calculer ? | Sur le Mac (Apple M1), un processus par cas, en parallèle |

### 2.1 Pourquoi $Re = 10\,000$

- **Au-dessus de $Re_c$ avec de la marge** : les valeurs publiées vont de 7400 à 8050 (`ref_data/README.md`). $Re = 10\,000$ est franchement au-delà, même pour la valeur la plus haute.
- **Des références existent** : Ghia et al. (1982, grille 257²) et Erturk et al. (2005, grille 601²) publient tous deux des profils à $Re = 10\,000$. Aucun $Re$ intermédiaire au-dessus de 8000 n'a des profils chez les deux.
- **Le budget** : la couche limite sous le couvercle a une épaisseur $\delta \sim Re^{-1/2} = 0.01$, soit 1.3 maille à $N = 128$ et 2.6 mailles à $N = 256$. Pour la résoudre à un $Re$ plus haut, il faudrait $N \propto Re^{1/2}$ ; et le pas de temps de stabilité du schéma centré vaut $\Delta t \le 2\nu \propto Re^{-1}$. Le coût d'un calcul croît donc au moins comme $Re \cdot N^2 \propto Re^2$ : $Re = 20\,000$ coûterait 4 fois plus, sur une grille ($N = 256$) qui prend déjà plusieurs jours (§5).

### 2.2 Grilles $N = 128$ puis $N = 256$

- $N = 128$ : la grille la plus fine des Parties A et B, donc des calculs comparables. C'est aussi la grille « à risque » d'Erturk et al. : sur 129², ils trouvent un calcul périodique à $Re = 7500$ alors que la solution est stationnaire sur 257².
- $N = 256$ : la grille de Ghia à $Re = 10\,000$, qui sert à trancher si l'instationnarité vue à $N = 128$ est physique ou numérique.

### 2.3 M1 et M2, Euler explicite, advection centrée

- **M1** (projection, cours §7.3.2) et **M2** ($\psi$–$\omega$, cours §7.3.4, vorticité de paroi par la formule de Thom) : les deux méthodes du projet, inchangées depuis le code de Bastian (commit `ca71e32`).
- **Euler explicite** (§3.1) : le schéma des Parties A et B. Le repli **RK3** (§3.9.4), stable pour l'advection centrée jusqu'à $C \le \sqrt3$ (§5.3.2.1), était prévu dans `PLAN.md` si Euler était trop cher. Il a été **écarté pour l'instant** : à $N = 128$, Euler tient en quelques heures (§5).
- **Advection centrée** : le choix « sage » (ordre 2, sans diffusion numérique). L'**upwind d'ordre 1** est calculé aussi, comme choix « malheureux », pour illustrer l'équation modifiée (§5.4.2) : il ajoute une diffusion numérique qui fait baisser le $Re$ effectif (§7.8).

### 2.4 Condition initiale et durée

- **Départ du repos** ($u = v = 0$, $\psi = \omega = 0$), couvercle en $U_t = 1 - e^{-t}$ : même condition initiale qu'en Parties A et B. Aucune perturbation n'est ajoutée ; si l'écoulement devient instationnaire, c'est que le schéma lui-même amplifie les petits écarts (arrondis, transitoire de démarrage).
- **$t_\text{max} = 300$** unités de temps $L/U$, **sans arrêt sur le résidu** (`residual_tol = 0`) : un régime périodique ne fait jamais tomber le résidu stationnaire, le critère des Parties A et B ne s'applique donc pas. La valeur 300 vient de `PLAN.md` (« simuler $t \approx 300$ ») : en Partie B, il fallait déjà $t = 250$ à $430$ pour atteindre l'état stationnaire à $Re = 5000$. Elle se révèle **trop courte** (§7.1).

### 2.5 Tolérance du gradient conjugué : $10^{-8}$

Même valeur qu'en Parties A et B (ADR 0003 §3) : résidu relatif $\|r_k\|/\|b\| \le 10^{-8}$, départ à chaud depuis la solution du pas précédent. Elle a été vérifiée en Partie A à $N = 128$ : l'écart entre tol $10^{-8}$ et $10^{-10}$ valait $\approx 10^{-10}$, contre $4 \cdot 10^{-3}$ d'écart à Ghia. Elle **n'a pas été re-testée** à $Re = 10\,000$ (choix 6, §4) : chaque calcul centré coûte plusieurs heures, et le test $\Delta t/2$ était plus informatif.

## 3. Le pas de temps

### 3.1 Les trois limites du cours

`cavity/common.py`, fonction `time_step`, applique les limites de stabilité d'Euler explicite avec une marge de 0.8 :

- **Fourier** (diffusion, centré 2D, cours §5.3.3 éq. 5.127) : $Fo = \nu\Delta t/h^2 \le 1/4$, donc $\Delta t \le h^2/(4\nu)$ ;
- **CFL** (advection upwind 2D, cours §5.3.3 éq. 5.131) : $C_x + C_y \le 1$ avec $|u|, |v| \le U = 1$, donc $\Delta t \le h/2$ ;
- **Advection centrée + Euler** (cours §5.3.1.4) : instable seule, stable grâce à la diffusion si $C^2 \le 2Fo$, ce qui donne $\Delta t \le 2\nu/U^2 = 2\nu$.

### 3.2 Valeurs à $Re = 10\,000$ ($\nu = 10^{-4}$), marge 0.8

| $N$ | $h$ | Fourier $0.8\,h^2/4\nu$ | CFL $0.8\,h/2$ | centré $0.8 \times 2\nu$ | $\Delta t$ centré | $\Delta t$ upwind |
| --- | --- | --- | --- | --- | --- | --- |
| 128 | $7.81\cdot10^{-3}$ | $1.22\cdot10^{-1}$ | $3.13\cdot10^{-3}$ | $\mathbf{1.6\cdot10^{-4}}$ | $1.6\cdot10^{-4}$ | $\mathbf{3.125\cdot10^{-3}}$ |
| 256 | $3.91\cdot10^{-3}$ | $3.05\cdot10^{-2}$ | $1.56\cdot10^{-3}$ | $\mathbf{1.6\cdot10^{-4}}$ | $1.6\cdot10^{-4}$ | $1.56\cdot10^{-3}$ |

- **À haut $Re$, la limite de Fourier ne joue plus** : $\nu$ minuscule la rend 1000 fois plus lâche que les autres. C'est l'inverse des Parties A et B (à $Re = 100$, $N = 128$, c'était elle qui fixait $\Delta t$ ; ADR 0003 §1).
- **Schéma centré** : la limite $\Delta t \le 2\nu$ est active, et **elle ne dépend pas de $N$** : $\Delta t = 1.6\cdot10^{-4}$ à $N = 128$ comme à $256$, soit $1.875\cdot10^6$ pas pour $t = 300$.
- **Upwind** : la CFL est active, $\Delta t = 3.125\cdot10^{-3}$ à $N = 128$ : 20 fois moins de pas ($96\,000$).
- **Reynolds de maille** $Uh/\nu = 78$ à $N = 128$ (39 à $N = 256$) : bien au-delà de 2, le seuil au-dessus duquel le schéma centré produit des oscillations nœud à nœud (cours §5.4). C'est le risque principal du schéma centré ici (§7.5).

## 4. Choix 6 : quels tests de sensibilité

Au total, sept cas sont lancés par `scripts/partie_C.py` (un indice désigne toujours le même cas) :

| Cas | Méthode | Advection | $\Delta t$ | $N$ | Rôle |
| --- | --- | --- | --- | --- | --- |
| 0 | M2 | centrée | $\Delta t/2 = 8\cdot10^{-5}$ | 128 | erreur d'Euler sur période et amplitude |
| 1 | M1 | centrée | $1.6\cdot10^{-4}$ | 128 | calcul principal |
| 2 | M2 | centrée | $1.6\cdot10^{-4}$ | 128 | calcul principal |
| 3 | M1 | upwind | $3.125\cdot10^{-3}$ | 128 | choix « malheureux » |
| 4 | M2 | upwind | $3.125\cdot10^{-3}$ | 128 | choix « malheureux » |
| 5 | M2 | centrée | $1.6\cdot10^{-4}$ | 256 | sensibilité au maillage (en cours) |
| 6 | M1 | centrée | $1.6\cdot10^{-4}$ | 256 | sensibilité au maillage (en cours) |

**Pourquoi $\Delta t/2$ sur M2 et pas sur M1.** Le test doit isoler l'erreur d'Euler sur la période et l'amplitude. Pour M2, c'est simple : l'état stationnaire de M2 ne dépend pas de $\Delta t$ (seule la vitesse d'évolution en dépend), donc toute différence entre $\Delta t$ et $\Delta t/2$ est l'erreur en temps. Pour M1, la projection non incrémentale mêle à l'erreur d'Euler une **erreur de projection** proportionnelle à $\Delta t$, qui modifie même l'état stationnaire (docs/adr/0002 §7). Un écart entre $\Delta t$ et $\Delta t/2$ sur M1 serait donc impossible à attribuer à une seule cause. Un seul calcul supplémentaire (7.8 h) suffit sur M2.

**Pourquoi pas la tolérance du gradient conjugué.** Elle a été validée en Partie A (§2.5) et ne dépend pas de $\Delta t$ (résidu relatif) ; un test de plus aurait coûté 2 à 5 h de calcul pour un effet attendu de l'ordre de $10^{-10}$.

## 5. Choix 10 : budget sur l'Apple M1

### 5.1 Un processus par cas, un cœur par processus

- Les noyaux Numba (opérateurs, gradient conjugué) sont des boucles **séquentielles** : un calcul n'utilise qu'un cœur. Paralléliser à l'intérieur d'un calcul (Numba `prange`, MKL multi-thread) n'apporterait presque rien sur une grille de $129^2$ et compliquerait le code.
- Les cas sont **indépendants** : on en lance plusieurs à la fois, un par cœur (`CLAUDE.md` §3, comme en Partie B). L'Apple M1 a **4 cœurs performance** (et 4 cœurs efficacité, bien plus lents) : `scripts/lancer_partie_C.sh` lance 4 cas à la fois (`xargs -P 4`), le 5ᵉ démarre dès qu'un cœur se libère.
- Chaque processus est limité à un thread (`OMP_NUM_THREADS=1`, `NUMBA_NUM_THREADS=1`, `MKL_THREADING_LAYER=SEQUENTIAL`) pour que les 4 calculs ne se disputent pas les cœurs ; `caffeinate -i` empêche la mise en veille du Mac.

### 5.2 Estimation contre réalité ($N = 128$)

| Cas | Pas | ms/pas estimé | Durée estimée | ms/pas réel | it. GC moy. | Durée réelle |
| --- | --- | --- | --- | --- | --- | --- |
| 0 : M2 centré $\Delta t/2$ | $3.75\cdot10^6$ | ~10 | ~11 h | 7.51 | 166 | **7.83 h** |
| 1 : M1 centré | $1.875\cdot10^6$ | ~15 | ~8 h | 3.45 | 42 | **1.80 h** |
| 2 : M2 centré | $1.875\cdot10^6$ | ~10 | ~5–6 h | 8.49 | 186 | **4.42 h** |
| 3 : M1 upwind | $9.6\cdot10^4$ | ~20–30 | < 1 h | 6.32 | 90 | **0.17 h** |
| 4 : M2 upwind | $9.6\cdot10^4$ | ~15–20 | < 1 h | 7.78 | 168 | **0.21 h** |

- Les estimations venaient de bancs de quelques milliers de pas au début du calcul, en plein transitoire, là où le gradient conjugué itère le plus. Ensuite le départ à chaud fait chuter le nombre d'itérations (M1 centré : 42 en moyenne, 21 à la fin), et les calculs ont été **plus rapides que prévu**, surtout M1.
- **M2 coûte plus cher que M1** : son Poisson porte sur $\psi$, avec un second membre $-\omega$ qui change beaucoup d'un pas à l'autre ; il demande ~4 fois plus d'itérations (186 contre 42).
- Pour $\Delta t/2$, un pas coûte un peu moins (166 itérations au lieu de 186) : la solution du pas précédent, point de départ du gradient conjugué, est plus proche.

### 5.3 $N = 256$ : la limite de la machine

Au départ : **93 ms par pas** (518 itérations en moyenne), 4 fois plus de nœuds et plus d'itérations. État le 24/11 à 6 h 40 (journaux) :

| Cas | Avancement | Temps écoulé | It. GC (dernier relevé) | Fin estimée |
| --- | --- | --- | --- | --- |
| 6 : M1, $N = 256$ | $t = 180 / 300$ | 12.4 h | 98 | 24/11 vers midi |
| 5 : M2, $N = 256$ | $t = 80 / 300$ | 12.2 h | 471 | 25/11 dans l'après-midi |

C'est pourquoi on ne monte ni en $Re$ ni en $N$ : un calcul M2 à $N = 256$ prend déjà environ 3 jours pour $t = 300$.

## 6. Ce qui a été ajouté au code (choix 9)

### 6.1 Code de calcul : trois options désactivées par défaut

Le code de Bastian (`ca71e32`) est conservé tel quel. La seule modification de `cavity/` est dans `march_to_steady` (`cavity/common.py`), la boucle en temps commune aux deux méthodes, qui reçoit trois options passées par `run_projection` et `run_psi_omega` :

- `probes` : liste d'indices $(i, j)$ ; $u$ et $v$ y sont enregistrés **à chaque pas** dans des tableaux réservés d'avance (taille `t_max/dt + 1`), pour ne pas ralentir la boucle ;
- `snapshot_times` : instants où l'on copie les champs complets ($u, v, p$ pour M1 ; $\psi, \omega, u, v$ pour M2) ;
- `log_every` : une ligne de journal toutes les 5 unités de temps ($t$, $K$, $P$, itérations du gradient conjugué, minutes écoulées), pour suivre les calculs à distance.

Les trois valent `None` par défaut. **Vérification** : avec les options par défaut, les résultats sont identiques bit à bit à ceux du code d'origine, y compris le cas de la Partie A (M1, $Re = 100$, $N = 32$ : arrêt au même pas, 6223). Les Parties A et B ne sont donc pas touchées.

### 6.2 Scripts

- `scripts/partie_C.py` (~110 lignes) : la liste des 7 cas, les sondes, les instants des instantanés ; trois modes : `liste`, `bench <i>` (3000 pas, coût par pas, rien n'est sauvé) et `<i>` (calcul complet et sauvegarde).
- `scripts/lancer_partie_C.sh` (~20 lignes) : lancement en parallèle, un journal par cas.
- `notebooks/C_unsteady.ipynb` : relit les résultats, mesure période et amplitude, trace les figures. Aucun calcul n'y est lancé.

### 6.3 Ce qui est sauvegardé

Un fichier `results/C_<méthode>_<schéma>_Re10000_N<N>_dt<Δt>.npz` par cas :

- **historiques à chaque pas** : $t$, $K$, $P$, résidu stationnaire, nombre d'itérations du gradient conjugué ($1.875\cdot10^6$ valeurs chacun pour un calcul centré) ;
- **6 sondes**, $u$ et $v$ à chaque pas, aux mêmes points physiques sur les deux grilles :

| Sonde | $(x, y)$ | Rôle |
| --- | --- | --- |
| coin_couvercle_gauche | $(2/128,\ 126/128)$ | à 2 mailles du coin, où le Reynolds de maille est maximal : détecte des zigzags numériques |
| coin_couvercle_droit | $(126/128,\ 126/128)$ | idem, au coin où l'écoulement du couvercle bute sur la paroi |
| centre | $(0.5,\ 0.5)$ | tourbillon principal |
| bas_gauche | $(0.125,\ 0.125)$ | tourbillon secondaire en bas à gauche |
| bas_droit | $(0.875,\ 0.125)$ | tourbillon secondaire en bas à droite |
| haut_gauche | $(0.125,\ 0.875)$ | tourbillon secondaire en haut à gauche |

- **156 instantanés** des champs : tous les 5 unités de temps de $t = 5$ à $280$ (évolution sur tout le calcul), puis tous les 0.2 sur les 20 dernières (environ 18 à 25 images par cycle) ;
- le champ final, les paramètres du cas, la machine, les versions des bibliothèques et la date.

Tailles : 76 et 97 Mo (upwind), 318 et 339 Mo (centré), 594 Mo ($\Delta t/2$), soit environ 1.4 Go pour $N = 128$. Ces fichiers ne sont pas versionnés (`.gitignore`).

## 7. Résultats à $N = 128$

### 7.1 $K(t)$, $P(t)$ et résidu

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/figures/C_K_P_residual.png)

*Figure 1 — $K(t)$, $P(t)$ et résidu stationnaire pour les 5 cas à $N = 128$ (trait plein : centré ; tirets : upwind). La courbe du cas 0 ($\Delta t/2$) est cachée sous celle de M2 centré.*

| Cas | $K(300)$ | $dK/dt$ sur $[250, 300]$ | $P(300)$ | $\lVert R\rVert$ à $t = 300$ |
| --- | --- | --- | --- | --- |
| M2 centré $\Delta t/2$ | $2.341\cdot10^{-2}$ | $1.63\cdot10^{-5}$ | $9.483\cdot10^{-3}$ | 1.77 |
| M1 centré | $4.360\cdot10^{-2}$ | $2.90\cdot10^{-5}$ | $8.319\cdot10^{-3}$ | 0.36 |
| M2 centré | $2.339\cdot10^{-2}$ | $1.64\cdot10^{-5}$ | $9.483\cdot10^{-3}$ | 1.70 |
| M1 upwind | $1.960\cdot10^{-2}$ | $7.0\cdot10^{-8}$ | $7.757\cdot10^{-3}$ | $4.0\cdot10^{-5}$ |
| M2 upwind | $5.930\cdot10^{-3}$ | $-6.4\cdot10^{-9}$ | $7.390\cdot10^{-3}$ | $4.2\cdot10^{-3}$ |

- **$P$ atteint un palier dès $t \approx 50$** : la couche limite sous le couvercle s'établit vite.
- **$K$ croît encore à $t = 300$** pour les calculs centrés : $+3.4\,\%$ (M1) et $+3.5\,\%$ (M2) sur les 50 dernières unités de temps. Le cœur du tourbillon principal se met en rotation très lentement. **À $t = 300$, l'écoulement moyen n'est pas établi.**
- Upwind : $\lVert R\rVert$ décroît exponentiellement et $K$ ne varie plus : écoulement stationnaire.

### 7.2 Choix 8 : mesurer la période et l'amplitude

**Méthode** (fonction `oscillation` du notebook, ~15 lignes, sans transformée de Fourier), sur les 50 dernières unités de temps d'un signal $s(t)$ ($K$ ou une sonde) :

1. retirer la tendance lente par une **parabole** de moindres carrés ;
2. découper le signal en lobes positifs (entre deux passages par zéro) et garder **un seul maximum par lobe** ;
3. **période** $T$ = moyenne des écarts entre maxima successifs (l'écart-type donne l'incertitude) ; **amplitude** d'un cycle = (max − min)/2 sur ce cycle ;
4. **taux** $\sigma$ = pente de $\ln(\text{amplitude})$ en fonction du temps : $\sigma \approx 0$ oscillation entretenue, $\sigma < 0$ amortie, $\sigma > 0$ croissante ;
5. si la variation relative du signal est inférieure à $10^{-5}$, il est déclaré **stationnaire**.

**Pièges rencontrés en la mettant au point :**

- *Une droite ne suffit pas.* La première version retirait une tendance linéaire. Or $K$ croît encore avec une courbure : après une droite, il reste un arc, et un lobe sur deux du signal ne repasse pas sous zéro. Résultat : l'oscillation de $K$ n'était pas détectée, ou la période sortait doublée. La parabole corrige ce défaut.
- *Les harmoniques.* Certaines sondes de M2 ont deux bosses par cycle. En prenant tous les maxima locaux, on mesurait $T \approx 2.5$ ; avec un seul maximum par lobe, on retrouve la fondamentale $T \approx 4.97$, la même que celle de $K$.
- *Le bruit d'arrondi.* Sur les calculs upwind (stationnaires), les variations restantes sont au niveau des arrondis et de la tolérance du gradient conjugué ; sans le seuil $10^{-5}$, elles étaient comptées comme des oscillations.

### 7.3 Période, amplitude et amortissement

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/figures/C_probes.png)

*Figure 2 — $u$ aux 6 sondes, tendance retirée, sur $t \in [250, 300]$, pour les trois calculs centrés.*

| Cas | Signal | Période $T$ | Amplitude | $\sigma$ | Régime |
| --- | --- | --- | --- | --- | --- |
| M2 centré $\Delta t/2$ | $K$ | $4.972 \pm 0.007$ | $7.92\cdot10^{-7}$ | $-1.52\cdot10^{-2}$ | amorti |
| M2 centré $\Delta t/2$ | sonde bas_gauche | $4.970 \pm 0.023$ | $4.47\cdot10^{-4}$ | $-1.47\cdot10^{-2}$ | amorti |
| M1 centré | $K$ | $3.642 \pm 0.016$ | $4.09\cdot10^{-6}$ | $-1.10\cdot10^{-2}$ | amorti |
| M1 centré | sonde coin_couvercle_gauche | $3.643 \pm 0.015$ | $7.64\cdot10^{-3}$ | $-0.89\cdot10^{-2}$ | amorti |
| M2 centré | $K$ | $4.974 \pm 0.007$ | $7.67\cdot10^{-7}$ | $-1.52\cdot10^{-2}$ | amorti |
| M2 centré | sonde bas_gauche | $4.971 \pm 0.024$ | $4.35\cdot10^{-4}$ | $-1.45\cdot10^{-2}$ | amorti |
| M1, M2 upwind | toutes | — | $< 10^{-5}$ relatif | — | stationnaire |

- Dans chaque calcul centré, **$K$ et les 6 sondes oscillent avec la même période** : c'est un seul mode, global.
- **La période dépend de la méthode** : 3.64 pour M1, 4.97 pour M2 (37 % d'écart).
- **Les oscillations sont amorties** : $\sigma$ entre $-0.007$ et $-0.015$ selon la sonde, l'amplitude est divisée par 2 en 50 à 100 unités de temps.

### 7.4 Résultat du test $\Delta t/2$ (cas 0, terminé le 24/11 à 1 h 13)

- **Période** : $4.972 \pm 0.007$ contre $4.974 \pm 0.007$ (écart 0.04 %, dans l'incertitude), sur $K$ comme sur toutes les sondes.
- **Amortissement** : $\sigma = -1.52\cdot10^{-2}$ sur $K$ dans les deux cas.
- **Amplitude** : à 3 % près ($7.92\cdot10^{-7}$ contre $7.67\cdot10^{-7}$ sur $K$).
- **$K$ et $P$** : $K(300)$ à 0.1 % près, $P(300)$ identique à 4 chiffres ; profils (§7.7) identiques à $10^{-3}$ près ; mêmes zigzags (§7.5).

**Conclusion : l'erreur d'Euler explicite est négligeable.** L'écart entre M1 et M2 ne vient pas du pas de temps, mais de la discrétisation en espace (grille, conditions aux parois).

### 7.5 Oscillation physique ou numérique ?

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/figures/C_lid_rows_du.png)

*Figure 3 — Gauche et centre : $u$ sur les rangées $y = 1 - h$, $1 - 2h$, $1 - 4h$, près des coins gauche et droit du couvercle, à $t = 300$. Droite : $u(300) - u(299)$, où l'écoulement bouge.*

- **Zigzags nœud à nœud** sous le couvercle près du coin **droit**, pour M1 et M2 centrés, jusqu'à $\pm 0.2$ sur $u$ sur la première rangée ; M1 en a aussi près du coin gauche. C'est le défaut attendu du schéma centré quand le Reynolds de maille (78) dépasse 2 : **à $N = 128$, la couche limite du couvercle n'est pas résolue**.
- Les zigzags de M2 sont **identiques à $\Delta t/2$** : ils viennent de l'espace, pas du temps. Les profils upwind sont lisses.
- La variation sur une unité de temps couvre toute la couche de cisaillement autour du tourbillon principal et les tourbillons secondaires : le mode est global. Mais il coexiste avec des zigzags numériques : à $N = 128$, **on ne peut pas dire quelle part de l'oscillation est physique**. C'est le rôle des calculs à $N = 256$.

### 7.6 Un cycle de vorticité

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/figures/C_vorticity_cycle_case1.png)

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/figures/C_vorticity_cycle_case2.png)

*Figure 4 — M1 (haut) et M2 (bas) centrés. À gauche, vorticité moyenne sur les 20 dernières unités de temps ; ensuite, fluctuation $\omega - \bar\omega$ à 4 instants espacés d'un quart de période. Échelles bornées (les pics de $\omega$ aux coins du couvercle, en $1/h$, écraseraient tout).*

- $\bar\omega$ : tourbillon principal à cœur de vorticité presque uniforme, entouré d'une fine couche de cisaillement ; stries près du coin haut droit (les zigzags du §7.5).
- Fluctuation : structures organisées qui font le tour de la cavité avec la couche de cisaillement, et structures lisses au centre. Amplitude $\approx 0.24$ pour M1, $\approx 0.04$ pour M2.

### 7.7 Comparaison à Ghia et Erturk

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/figures/C_profiles_ghia_erturk.png)

*Figure 5 — $u(0.5, y)$ et $v(x, 0.5)$ moyennés sur les 20 dernières unités de temps, comparés à Ghia (257²) et Erturk (601²) à $Re = 10\,000$.*

| Cas | max$\lvert\Delta u\rvert$ Ghia | max$\lvert\Delta v\rvert$ Ghia | max$\lvert\Delta u\rvert$ Erturk | max$\lvert\Delta v\rvert$ Erturk |
| --- | --- | --- | --- | --- |
| M2 centré $\Delta t/2$ | 0.144 | 0.159 | 0.156 | 0.198 |
| M1 centré | **0.039** | **0.038** | **0.027** | 0.104 |
| M2 centré | 0.144 | 0.159 | 0.156 | 0.199 |
| M1 upwind | 0.224 | 0.232 | 0.237 | 0.255 |
| M2 upwind | 0.591 | 0.543 | 0.598 | 0.571 |

**Attention** : Ghia et Erturk calculent des solutions **stationnaires** (sans pas de temps physique). Au-delà de $Re_c$, ce sont des solutions stationnaires instables, que notre moyenne temporelle ne reproduit pas forcément. Écart maximal le long de chaque ligne, parois exclues, coquille connue de Ghia à $y = 0.5$ ignorée.

- **M1 centré est le plus proche** : 0.03 à 0.04, sauf $v$ près de la paroi droite comparé à Erturk (0.10).
- **M2 centré s'en écarte nettement** (0.14 à 0.20) : tourbillon principal trop faible. C'est le défaut de la formule de Thom déjà analysé en Partie B (`Doc/partie_B_limites.md` §3) : elle sous-estime toujours $|\omega|$ à la paroi, d'une fraction $h/(3\delta_\omega)$, et $\delta_\omega \propto Re^{-1/2}$, donc à $h$ fixé l'erreur croît comme $\sqrt{Re}$ (−22 % mesuré à $Re = 5000$, $N = 128$). À $Re = 10\,000$, elle est encore plus grande.

### 7.8 L'upwind, choix « malheureux »

L'équation modifiée du schéma upwind (cours §5.4.2.1) fait apparaître une diffusion numérique $D = \tfrac12 U h (1 - C)$, avec $C = U\Delta t / h$ :

- $N = 128$ : $C = 0.40$, $D = 2.34\cdot10^{-3}$, soit **$D/\nu = 23$** ;
- le schéma voit donc un Reynolds effectif $Re_\text{eff} = UL/(\nu + D) \approx$ **410**, très en dessous de $Re_c \approx 8000$.

Comme prévu, M1 et M2 upwind donnent un écoulement **stationnaire** : le schéma supprime la physique qu'on cherche. Leurs profils sont loin des références.

### 7.9 Tableau de sensibilité

| Cas | Régime ($K$) | $T$ ($K$) | Sonde la plus agitée | $T$ | Ampl. | $\bar K$ ($t \in [250, 300]$) | $\bar P$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M2 centré $\Delta t/2$ | amorti | 4.972 | bas_gauche | 4.970 | $4.5\cdot10^{-4}$ | $2.302\cdot10^{-2}$ | $9.483\cdot10^{-3}$ |
| M1 centré | amorti | 3.642 | coin_couvercle_gauche | 3.643 | $7.6\cdot10^{-3}$ | $4.291\cdot10^{-2}$ | $8.309\cdot10^{-3}$ |
| M2 centré | amorti | 4.974 | bas_gauche | 4.971 | $4.3\cdot10^{-4}$ | $2.300\cdot10^{-2}$ | $9.483\cdot10^{-3}$ |
| M1 upwind | stationnaire | — | — | — | — | $1.960\cdot10^{-2}$ | $7.757\cdot10^{-3}$ |
| M2 upwind | stationnaire | — | — | — | — | $5.930\cdot10^{-3}$ | $7.390\cdot10^{-3}$ |

Lecture : à $N = 128$, **le pas de temps ne change rien**, mais **la méthode change tout** (période ; amplitude 5 fois plus grande pour M1 sur $K$, 18 fois sur la sonde la plus agitée ; $K$ moyen d'un facteur 2), et **l'advection** décide du régime lui-même.

## 8. Interprétation honnête

Ce qu'on peut dire, à $t = 300$ et $N = 128$ :

1. **Pas de régime périodique établi.** M1 et M2 centrés oscillent selon un mode global, mais amorti ($\sigma < 0$), avec une période qui dépend de la méthode (3.64 contre 4.97).
2. **L'écoulement moyen n'est pas établi** : $K$ croît encore de 3 % par 50 unités de temps. Tant qu'il évolue, on ne peut conclure ni « périodique » ni « stationnaire » : le mode amorti pourrait aussi bien s'éteindre que se stabiliser une fois l'écoulement moyen établi.
3. **Le pas de temps n'y est pour rien** (test $\Delta t/2$) : les écarts sont spatiaux.
4. **Zigzags numériques** sous le couvercle, pour les deux méthodes : $N = 128$ est trop grossier pour le schéma centré à $Re = 10\,000$ (Reynolds de maille 78, couche limite sur 1.3 maille). C'est cohérent avec Erturk et al. : sur 129², une partie de l'instationnarité est numérique.
5. **M2 est pénalisé par la formule de Thom**, comme en Partie B, et davantage à ce $Re$.
6. **Upwind** : stationnaire à $Re_\text{eff} \approx 410$, exactement ce que prévoit l'équation modifiée.

**Points non expliqués (à ne pas masquer à l'oral) :**

- le cœur de M2 upwind est presque au repos ($K$ 3.3 fois plus faible que M1 upwind, profils plats au centre), alors que les deux schémas ont la même diffusion numérique ;
- le facteur 2 sur $K$ entre M1 et M2 centrés, alors que $P$ ne diffère que de 14 % : la différence est dans le cœur du tourbillon, pas sous le couvercle. Thom l'explique au moins en partie, mais l'ampleur n'est pas chiffrée.

## 9. Décision en attente

L'équipe n'a pas encore tranché ; trois options ont été posées le 23/11. **Aucune n'a été exécutée.**

| Option | Contenu | Coût | Pour | Contre |
| --- | --- | --- | --- | --- |
| **A** | Relancer M1 et M2 centrés à $N = 128$ avec $t_\text{max} = 1000$ | ~6 h (M1), ~15 h (M2), en parallèle | Voir si le mode amorti s'éteint ou s'entretient une fois l'écoulement moyen établi | Toujours à $N = 128$, avec les zigzags numériques |
| **B** | Ajouter une reprise (sauvegarde/redémarrage) pour prolonger les calculs $N = 256$ au-delà de $t = 300$ | du code en plus, puis plusieurs jours | La seule grille qui peut trancher physique/numérique | Exclu par `CLAUDE.md` (« pas de sauvegarde/reprise sauf demande explicite ») ; délai avant le rendu |
| **C** | Garder $t = 300$ et présenter un échec expliqué | rien | Honnête et déjà documenté (§8) ; l'énoncé accepte les échecs expliqués | Pas de conclusion sur le régime |

À noter : les calculs $N = 256$ (cas 5 et 6) finissent de toute façon à $t = 300$ ; leur résultat (zigzags présents ou non, période, amortissement) peut aider à choisir.

## 10. Références

- **Ghia, Ghia & Shin (1982)**, J. Comput. Phys. 48, 387–411 : $\psi$–$\omega$ stationnaire, grille 257² à $Re = 10\,000$ ; profils sur les lignes médianes (`ref_data/ghia1982_*.csv`).
- **Erturk, Corke & Gökçöl (2005)**, Int. J. Numer. Meth. Fluids 48, 747–774 : centré d'ordre 2 + Thom, grilles jusqu'à 601² ; profils à $Re = 10\,000$ (`ref_data/erturk2005_*.csv`). Sur 129², pas de solution stationnaire à $Re = 7500$ (calcul périodique) ; sur 257², stationnaire jusqu'à $Re = 12\,500$.
- **$Re_c$** : 7400 à 8050 selon la méthode et le maillage (Peng et al. 2003 : 7402 ; Poliashenko & Aidun : 7763 ; Fortin et al. 1997 : 7998.5, cités par Erturk et al. ; Auteri et al. 2002 : 8018, à vérifier sur l'article).
- **Période de référence** : aucune n'a encore été relevée dans un article ; on n'en cite pas sans l'avoir lue.
- Documents internes : `PLAN.md` (étape 5), `CLAUDE.md`, `docs/adr/0002`, `0003`, `0004`, `Doc/partie_B_limites.md`, `ref_data/README.md`.
