# Partie C — résultats, choix et interprétation

Document pour les membres du groupe, état au **25 novembre 2026, fin d'après-midi** : tous les calculs sont terminés. Il remplace `Doc/partie_C_choix.md` (24/11), écrit avant les calculs longs et à $N = 256$. Tous les chiffres viennent du notebook `notebooks/C_unsteady.ipynb` (ré-exécuté le 25/11) ; les figures en sont tirées.

## 1. L'essentiel en une page

- **Ce qu'on a fait.** Cavité entraînée à $Re = 10\,000$ (au-delà de $Re_c \approx 8000$), avec nos deux méthodes (M1 projection, M2 $\psi$–$\omega$), Euler explicite, départ du repos ; 11 calculs, environ 95 h de calcul sur le Mac.
- **Le pas de temps n'y est pour rien.** M2 avec $\Delta t/2$ redonne la même période à 0.04 % près : l'erreur d'Euler est négligeable, les écarts viennent de l'espace.
- **L'upwind fait ce que prévoit le cours** : diffusion numérique 23 fois plus grande que $\nu$, Reynolds effectif $\approx 410$, écoulement stationnaire.
- **Il faut des temps très longs.** L'écoulement moyen n'est établi que vers $t \approx 1000$, et le régime final de M1 vers $t \approx 1450$. À $t = 300$ (la durée prévue au départ), toutes les oscillations semblaient amorties.
- **M1, $N = 128$ : un vrai régime périodique** (cycle limite) de période **$T = 2.283$**, amplitude constante de $t = 1450$ à $2000$. C'est le seul régime périodique établi de nos calculs.
- **M2, $N = 128$ : stationnaire** jusqu'à $t = 2000$. Sur la même grille, les deux méthodes ne donnent donc pas le même régime.
- **À $N = 256$, M1 et M2 se rapprochent nettement** (écart sur $K$ : 46 % → 13 % ; erreur de M2 sur les profils divisée par 2.4 à 3.8), mais ces calculs ne vont que jusqu'à $t = 300$ (19 h et 38 h de calcul) : leur régime final n'est pas connu.
- **Conclusion honnête** : nous obtenons une perte de stationnarité et un cycle limite avec M1, mais nous ne pouvons pas affirmer qu'il est physique : il dépend de la méthode à $N = 128$ et sa période dépend encore du maillage.

## 2. L'énoncé et le cadre

L'énoncé (Partie C) : au-delà de $Re_c$, la solution stationnaire devient instable et l'écoulement devient périodique. Il demande de **simuler au $Re$ le plus élevé possible** en le justifiant, d'**étudier la sensibilité aux choix numériques** et de **comparer à la littérature**. Diagnostics : $K(t) = \tfrac12 \int (u^2 + v^2)\,dS$ et $P(t) = \int_0^1 \nu\, \partial_y u|_{y=1}\, U_t\, dx$, couvercle $U_t = 1 - e^{-t}$, $\nu = 1/Re$.

Cadre (`CLAUDE.md`) : strict minimum, uniquement des méthodes du cours (poly *Numerical Methods*), chaque choix cité avec sa section.

## 3. Les choix et leurs raisons

| Choix | Décision | Pourquoi |
| --- | --- | --- |
| Reynolds | $Re = 10\,000$ | Au-dessus de toutes les valeurs publiées de $Re_c$ (7400 à 8050) ; profils de référence chez Ghia (257²) **et** Erturk (601²) |
| Grilles | $N = 128$, puis $256$ | 128 : la grille des Parties A et B ; 256 : celle de Ghia à ce $Re$, pour juger l'effet du maillage |
| Méthodes | M1 (§7.3.2) et M2 (§7.3.4, Thom) | Les deux méthodes du projet, code inchangé |
| Temps | Euler explicite (§3.1) | Celui des Parties A et B ; le test $\Delta t/2$ montre qu'il suffit |
| Advection | centrée ; upwind en contre-exemple | Centré : ordre 2 sans diffusion numérique ; upwind : choix « malheureux » expliqué par l'équation modifiée (§5.4.2) |
| $\Delta t$ | limites de stabilité du cours × 0.8 | Voir ci-dessous |
| Départ | du repos, sans perturbation | Même condition initiale qu'en A et B |
| Durée | 300, puis 1000, puis 2000 | 300 venait de `PLAN.md` ; trop court (§5.3) |
| Gradient conjugué | tol. $10^{-8}$, départ à chaud | Validé en Partie A (écart $10^{-10}$ entre $10^{-8}$ et $10^{-10}$) |
| Période | méthode simple des maxima, sans FFT | Expliquable à l'oral (§4) |

**Le pas de temps** (`common.time_step`). Trois limites du cours (§5.3) : Fourier $\Delta t \le h^2/4\nu$ ; CFL upwind $\Delta t \le h/2$ ; advection centrée + Euler, stable seulement grâce à la diffusion si $C^2 \le 2Fo$, soit $\Delta t \le 2\nu$.

| $N$ | Fourier | CFL | centré $2\nu$ | $\Delta t$ centré | $\Delta t$ upwind | Reynolds de maille $h/\nu$ |
| --- | --- | --- | --- | --- | --- | --- |
| 128 | $1.2\cdot10^{-1}$ | $3.1\cdot10^{-3}$ | $1.6\cdot10^{-4}$ | $\mathbf{1.6\cdot10^{-4}}$ | $\mathbf{3.1\cdot10^{-3}}$ | 78 |
| 256 | $3.1\cdot10^{-2}$ | $1.6\cdot10^{-3}$ | $1.6\cdot10^{-4}$ | $\mathbf{1.6\cdot10^{-4}}$ | $1.6\cdot10^{-3}$ | 39 |

À haut $Re$, la limite de Fourier ne joue plus. Pour le centré, $\Delta t \le 2\nu$ **ne dépend pas de $N$** : $1.875\cdot10^6$ pas pour $t = 300$, à $N = 128$ comme à $256$. Le Reynolds de maille dépasse largement 2 : le schéma centré risque des oscillations nœud à nœud (§5.5).

**Pourquoi pas un $Re$ plus élevé.** La couche limite du couvercle a une épaisseur $\sim Re^{-1/2} = 0.01$ (1.3 maille à $N = 128$). Il faudrait $N \propto Re^{1/2}$, et $\Delta t \propto Re^{-1}$ : le coût croît au moins comme $Re^2$. $Re = 20\,000$ coûterait 4 fois plus, sur une grille qui demande déjà 19 à 38 h pour $t = 300$.

## 4. Les calculs et la mesure du régime

| Cas | Méthode | Advection | $N$ | $t_\text{max}$ | Pas | Durée | it. GC moy. | Rôle |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | M2 | centrée, $\Delta t/2$ | 128 | 300 | $3.75\cdot10^6$ | 7.8 h | 166 | erreur en temps |
| 1 | M1 | centrée | 128 | 300 | $1.88\cdot10^6$ | 1.8 h | 42 | calcul principal |
| 2 | M2 | centrée | 128 | 300 | $1.88\cdot10^6$ | 4.4 h | 186 | calcul principal |
| 3 | M1 | upwind | 128 | 300 | $9.6\cdot10^4$ | 0.2 h | 90 | contre-exemple |
| 4 | M2 | upwind | 128 | 300 | $9.6\cdot10^4$ | 0.2 h | 168 | contre-exemple |
| 5 | M2 | centrée | 256 | 300 | $1.88\cdot10^6$ | **38.0 h** | 440 | maillage |
| 6 | M1 | centrée | 256 | 300 | $1.88\cdot10^6$ | **19.2 h** | 139 | maillage |
| 7 | M1 | centrée | 128 | 1000 | $6.25\cdot10^6$ | 2.9 h | 19 | temps long |
| 8 | M2 | centrée | 128 | 1000 | $6.25\cdot10^6$ | 6.6 h | 83 | temps long |
| 9 | M2 | centrée | 128 | 2000 | $1.25\cdot10^7$ | 6.8 h | 42 | régime final |
| 10 | M1 | centrée | 128 | 2000 | $1.25\cdot10^7$ | 7.4 h | 28 | régime final |

Un processus par cas, un cœur par processus, jusqu'à 4 en parallèle (`scripts/lancer_partie_C.sh`). Les cas 7 à 10 repartent du repos et redonnent **exactement** (écart 0) les cas plus courts sur leur durée commune : les calculs sont reproductibles. Aucun NaN, $|u| \le 1$ partout, gradient conjugué toujours loin de sa limite d'itérations (878 au plus).

**Ce qui est enregistré.** $t$, $K$, $P$, résidu et itérations à chaque pas ; $u$ et $v$ en **6 sondes** à chaque pas (près des deux coins du couvercle, au centre, près des trois autres coins) ; des instantanés des champs (tous les 5, puis tous les 0.2 sur les 20 derniers temps).

**Mesurer le régime** (fonction `oscillation` du notebook), sur une fenêtre de 50 unités de temps : (1) retirer la tendance lente par une parabole ; (2) garder un seul maximum par lobe positif (évite de compter les harmoniques) ; (3) période = écart moyen entre maxima, amplitude = (max − min)/2 par cycle ; (4) taux $\sigma$ = pente de $\ln(\text{amplitude})$ : $\sigma < 0$ amortie, $\sigma > 0$ croissante, $\sigma \approx 0$ **périodique établi** ; (5) variation relative $< 10^{-5}$ : stationnaire. On répète la mesure sur des fenêtres successives pour suivre le régime dans le temps.

## 5. Résultats

### 5.1 Le pas de temps : erreur négligeable

M2 à $\Delta t$ et $\Delta t/2$ (cas 2 et 0) : période $4.974$ contre $4.972$ ($\pm 0.007$), même taux $\sigma = -1.52\cdot10^{-2}$, amplitudes à 3 % près, $K(300)$ à 0.1 %, mêmes profils et mêmes zigzags. **L'erreur d'Euler est négligeable** : tout écart entre méthodes ou maillages est spatial. (Le test est fait sur M2 parce que l'état de M1 dépend de $\Delta t$ par l'erreur de projection, ADR 0002 §7 : on n'aurait pas pu séparer les deux erreurs.)

### 5.2 L'upwind : le choix « malheureux »

Équation modifiée (§5.4.2.1) : diffusion numérique $D = \tfrac12 U h(1 - C)$ avec $C = 0.40$, soit $D = 2.34\cdot10^{-3}$, **$D/\nu = 23$**. Le schéma calcule en fait un écoulement à $Re_\text{eff} = 1/(\nu + D) \approx 410$, très en dessous de $Re_c$ : **stationnaire**, comme prévu. Profils très loin des références (0.22 à 0.60).

### 5.3 Temps long à $N = 128$ : stationnaire pour M2, cycle limite pour M1

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/figures/C_regime_evolution.png)

*Figure 1 — Amplitude (échelle log) et période de l'oscillation de $u$ à la sonde bas gauche, mesurées sur des fenêtres de 50 unités de temps. Traits pleins : $N = 128$ ($t \le 2000$) ; tirets : $N = 256$ ($t \le 300$).*

| Fenêtre finissant à | M2 : $T$ | M2 : amplitude | M1 : $T$ | M1 : amplitude | M1 : $\sigma$ |
| --- | --- | --- | --- | --- | --- |
| 300 | 4.97 | $4.4\cdot10^{-4}$ | 3.64 | $1.7\cdot10^{-3}$ | $-1.1\cdot10^{-2}$ |
| 500 | 4.70 | $1.2\cdot10^{-5}$ | 3.44 | $3.8\cdot10^{-5}$ | $-2.1\cdot10^{-2}$ |
| 700 | — | bruit | 2.288 | $2.6\cdot10^{-5}$ | $+7.9\cdot10^{-3}$ |
| 1000 | — | bruit | 2.284 | $3.4\cdot10^{-4}$ | $+8.9\cdot10^{-3}$ |
| 1300 | — | bruit | 2.282 | $3.8\cdot10^{-3}$ | $+6.0\cdot10^{-3}$ |
| 1500 | — | bruit | 2.282 | $5.8\cdot10^{-3}$ | $+4.6\cdot10^{-4}$ |
| 2000 | — | bruit | **2.283** | $\mathbf{5.8\cdot10^{-3}}$ | $\mathbf{+4\cdot10^{-5}}$ |

- **M2** : l'oscillation de $t = 300$ (période ~5, qui raccourcit en s'amortissant) s'éteint ; après $t \approx 700$, plus rien au-dessus du bruit d'arrondi, jusqu'à $t = 2000$. $K$ ne bouge plus ($dK/dt \approx 8\cdot10^{-10}$).
- **M1** : la première oscillation (période 3.6 → 3.4) s'éteint aussi, mais pas complètement : vers $t \approx 600$, **une deuxième oscillation**, de période 2.28, prend le relais et croît (amplitude doublée en ~80 unités de temps), puis **sature** vers $t \approx 1450$. De 1450 à 2000, amplitude constante et $T = 2.283$ : **cycle limite**, ce qu'on attend après une bifurcation de Hopf (une instabilité qui croît puis sature). La même période se lit sur $K$ et sur les 6 sondes : un seul mode, global.

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/Doc/img_partie_C/probes_case10.png)

*Figure 2 — M1, $N = 128$ : $u$ aux 6 sondes sur $t \in [1950, 2000]$, tendance retirée. Oscillation sinusoïdale d'amplitude constante ($1.9\cdot10^{-2}$ près du coin gauche du couvercle).*

### 5.4 Effet du maillage ($N = 256$, $t \le 300$)

| | M1, $N = 128$ | M1, $N = 256$ | M2, $N = 128$ | M2, $N = 256$ |
| --- | --- | --- | --- | --- |
| $K(300)$ | $4.36\cdot10^{-2}$ | $4.18\cdot10^{-2}$ | $2.34\cdot10^{-2}$ | $3.63\cdot10^{-2}$ |
| $P(300)$ | $8.32\cdot10^{-3}$ | $9.45\cdot10^{-3}$ | $9.48\cdot10^{-3}$ | $1.01\cdot10^{-2}$ |
| Période vers $t = 300$ | 3.64 | ~1.46 | 4.97 | ~2.1 |

- **M1 et M2 se rapprochent** : l'écart sur $K$ passe de 46 % à 13 %. Le désaccord à $N = 128$ venait donc surtout du manque de résolution, en particulier de la formule de Thom (§5.7).
- **Périodes du même ordre** (1.5 et 2.1) mais pas égales, et **pas de régime établi** : chez M2, l'oscillation croît jusqu'à $t \approx 200$ puis s'amortit fortement ; chez M1, la période raccourcit encore (1.80 → 1.58 → 1.46) et $\sigma$ change de signe d'une sonde à l'autre.
- $P$ augmente avec $N$ : c'est attendu, à cause de la singularité aux coins du couvercle ($P$ croît comme $\ln N$, `PLAN.md` §2). On compare donc les $P$ à $N$ fixé.

### 5.5 Physique ou numérique ? Les zigzags sous le couvercle

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/Doc/img_partie_C/lid_rows.png)

*Figure 3 — De haut en bas : M1 $N = 128$ ($t = 2000$), M1 $N = 256$, M2 $N = 256$ ($t = 300$). Gauche et centre : $u$ sur les 1ʳᵉ, 2ᵉ et 4ᵉ rangées sous le couvercle, près des coins. Droite : $u(t) - u(t - 1)$, où l'écoulement bouge.*

- **Zigzags nœud à nœud** près du coin droit du couvercle : le défaut du schéma centré quand le Reynolds de maille dépasse 2 (78 à $N = 128$). Ils sont plus faibles à $N = 256$ (Reynolds de maille 39), surtout pour M2 (une seule rangée, $x \gtrsim 0.93$).
- **Ils ne sont pas la source de l'oscillation** : ils existent aussi dans l'état **stationnaire** de M2, et l'écoulement bouge ailleurs (carte de droite) : le long de la **paroi gauche**, de la paroi du bas et au **coin haut gauche**, pour les deux méthodes et les deux grilles.

### 5.6 La structure de l'oscillation

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/Doc/img_partie_C/vorticity.png)

*Figure 4 — Haut : M1 $N = 128$ en cycle limite ($t \approx 2000$) ; bas : M2 $N = 256$ ($t \approx 300$). À gauche, vorticité moyenne ; ensuite, fluctuation $\omega - \bar\omega$ à 4 instants.*

Tourbillon principal à vorticité presque uniforme, entouré d'une couche de cisaillement. La fluctuation est organisée : **paquets qui descendent le long de la paroi gauche depuis le coin haut gauche**, et bras en spirale vers le cœur (M1). **M2 à $N = 256$ montre le même type de structure** aux mêmes endroits : les deux méthodes voient la même instabilité quand le maillage est assez fin.

### 5.7 $K$, $P$ et comparaison à la littérature

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/figures/C_K_P_residual.png)

*Figure 5 — $K(t)$, $P(t)$ et résidu stationnaire, tous les cas. Le résidu de M1 remonte et se stabilise vers 1.4 (oscillation) ; celui de M2 descend jusqu'à un plancher de $\sim 5\cdot10^{-4}$ fixé par la tolérance du gradient conjugué.*

![](file:///Users/denisvictor522gmail.com/Desktop/projet%204/Projet4_CFD/figures/C_profiles_ghia_erturk.png)

*Figure 6 — $u(0.5, y)$ et $v(x, 0.5)$ moyennés sur les 20 derniers temps, comparés à Ghia (257²) et Erturk (601²).*

| Cas | $\lvert\Delta u\rvert$ Ghia | $\lvert\Delta v\rvert$ Ghia | $\lvert\Delta u\rvert$ Erturk | $\lvert\Delta v\rvert$ Erturk |
| --- | --- | --- | --- | --- |
| M1, $N = 128$, $t = 2000$ | **0.028** | **0.034** | **0.018** | 0.112 |
| M1, $N = 256$, $t = 300$ | 0.044 | 0.040 | 0.033 | **0.029** |
| M2, $N = 128$, $t = 2000$ | 0.134 | 0.150 | 0.152 | 0.189 |
| M2, $N = 256$, $t = 300$ | 0.059 | 0.060 | 0.053 | 0.053 |
| M1 / M2 upwind | 0.22 / 0.59 | 0.23 / 0.54 | 0.24 / 0.60 | 0.26 / 0.57 |

- **M1 est le plus proche des références** (0.02 à 0.04), sauf $v$ près de la paroi droite à $N = 128$ (0.11), corrigé à $N = 256$ (0.03).
- **M2 à $N = 128$ est loin (0.13 à 0.19) bien qu'il soit stationnaire** : son état stationnaire n'est pas celui des références. C'est l'erreur de la formule de Thom vue en Partie B (elle sous-estime $|\omega|$ à la paroi, d'autant plus que $h\sqrt{Re}$ est grand). Elle est divisée par 2.4 à 3.8 à $N = 256$.
- Attention : Ghia et Erturk calculent des solutions **stationnaires** (sans temps physique). Au-delà de $Re_c$, ce sont des solutions stationnaires instables, que notre moyenne temporelle n'a pas à reproduire exactement.

## 6. Interprétation : ce qu'on peut dire, ce qu'on ne peut pas

**Ce qu'on peut affirmer :**

1. Les calculs sont stables, reproductibles, et l'erreur en temps est négligeable.
2. L'upwind tue l'instabilité, comme le prévoit l'équation modifiée.
3. À $Re = 10\,000$, partir du repos impose des temps très longs : $t \approx 1000$ pour l'écoulement moyen, $t \approx 1450$ pour le régime de M1.
4. Avec M1 à $N = 128$, l'écoulement perd sa stationnarité et atteint un **cycle limite** de période $T = 2.283$, d'amplitude constante.
5. L'oscillation vient de la couche de cisaillement le long de la paroi gauche et du coin haut gauche, pas des zigzags du coin droit.

**Ce qu'on ne peut pas affirmer :**

1. **Que ce cycle limite est physique.** À $N = 128$, M2 est stationnaire ; et à $N = 256$, les périodes changent encore (1.5 et 2.1 à $t = 300$). Ce comportement est cohérent avec Erturk et al. : sur 129², ils trouvent un calcul périodique à $Re = 7500$ alors qu'ils obtiennent une solution stationnaire sur 257² jusqu'à $Re = 12\,500$ (une partie de l'instationnarité sur grille grossière est numérique ; nuance : leur solveur stationnaire peut converger vers une solution stationnaire instable).
2. **La valeur de la période physique** : elle n'est pas convergée en maillage, et nous n'avons relevé aucune période de référence dans un article.
3. **Le régime final à $N = 256$** : il faudrait M1 et M2 à $N = 256$ jusqu'à $t \gtrsim 1500$, soit environ 4 à 8 jours de calcul par cas sur cette machine (au rythme mesuré de 19 h et 38 h pour $t = 300$).

**Points non expliqués** (à ne pas masquer) : le cœur de M2 upwind est presque au repos ($K$ 3.3 fois plus faible que M1 upwind) alors que les deux ont la même diffusion numérique ; l'écart entre les périodes de M1 et M2 à $N = 256$.

## 7. Questions probables à l'oral

- **« Pourquoi $Re = 10\,000$ ? »** Au-dessus de toutes les valeurs publiées de $Re_c$, références disponibles chez deux auteurs, et coût en $Re^2$ (§3).
- **« Votre oscillation n'est-elle pas numérique ? »** Trois arguments contre : mode global (une période pour $K$ et toutes les sondes), croissance puis saturation, structure localisée là où est la couche de cisaillement et non aux zigzags. Deux arguments pour : M2 stationnaire sur la même grille, période qui change avec $N$. Nous ne tranchons pas (§6).
- **« Pourquoi pas RK3, qui permet un plus grand $\Delta t$ ? »** Le test $\Delta t/2$ montre que l'erreur en temps est négligeable ; le coût vient du nombre de pas imposé par $\Delta t \le 2\nu$, que RK3 aurait relâché (§5.3.2.1 du cours). C'était le repli prévu dans `PLAN.md` ; il n'a pas été nécessaire à $N = 128$.
- **« Pourquoi M1 et M2 ne donnent-ils pas la même chose ? »** Mêmes équations, mais discrétisations différentes aux parois (pression de Neumann par nœud fantôme pour M1, formule de Thom pour M2). L'écart diminue fortement de $N = 128$ à $256$ : c'est une erreur de discrétisation.
- **« Pourquoi $P$ change-t-il avec $N$ ? »** Le saut de vitesse aux coins du couvercle rend $P$ divergent en continu ; en discret il croît en $\ln N$. On compare les $P$ à $N$ fixé.
- **« Comment mesurez-vous la période ? »** Maxima successifs après retrait d'une tendance parabolique, un maximum par lobe ; taux par la pente de $\ln$(amplitude) (§4). Limite connue : $K$ peut être classé « stationnaire » quand ses variations relatives sont sous $10^{-5}$ alors que les sondes oscillent (cas 7 vers $t = 600$) : le régime se lit sur les sondes.

## 8. Pour refaire les calculs et les figures

- `python scripts/partie_C.py liste` : les 11 cas ; `bash scripts/lancer_partie_C.sh 9 10` : lance des cas en parallèle (journaux `results/logs/C_<i>.log`).
- `notebooks/C_unsteady.ipynb` : relit les fichiers `results/C_*.npz` (0.08 à 2 Go chacun, non versionnés), aucun calcul n'y est lancé ; exécution en moins d'une minute.
- Code de calcul : `cavity/` (fonctions décrites dans `Doc/structure_code.pdf`) ; seule modification pour la Partie C : les options `probes`, `snapshot_times`, `log_every` de `march_to_steady`, désactivées par défaut (résultats des Parties A et B identiques bit à bit).

## 9. Références

- **Ghia, Ghia & Shin (1982)**, J. Comput. Phys. 48, 387–411 : $\psi$–$\omega$ stationnaire, 257² à $Re = 10\,000$ (`ref_data/ghia1982_*.csv`).
- **Erturk, Corke & Gökçöl (2005)**, Int. J. Numer. Meth. Fluids 48, 747–774 : centré d'ordre 2 + Thom, jusqu'à 601² (`ref_data/erturk2005_*.csv`) ; résultats sur 129² et 257² cités au §6.
- **$Re_c$** : 7400 à 8050 selon la méthode et le maillage (Peng et al. 2003 : 7402 ; Poliashenko & Aidun : 7763 ; Fortin et al. 1997 : 7998.5 ; Auteri et al. 2002 : 8018, à vérifier sur l'article) — `ref_data/README.md`.
- **Période de référence** : aucune relevée dans un article ; nous n'en citons pas.
- Internes : `PLAN.md` (étape 5), `CLAUDE.md`, `docs/adr/0002` à `0004`, `Doc/partie_B_limites.md`, `Doc/partie_C_choix.md` (version du 24/11).
