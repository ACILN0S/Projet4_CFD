# Projet 4 — Cavité entraînée à haut Reynolds : plan (version minimale)

> Règle n°1 : **faire le strict minimum demandé par l'énoncé**, avec **uniquement des idées du cours**
> (poly *Numerical Methods*, Vicquelin & Vié, 2025-26). Chaque choix ci-dessous cite le chapitre du cours
> dont il vient. Si une idée n'est pas dans le cours, on ne l'utilise pas.

---

## 1. Le problème

Cavité carrée 2D, Navier–Stokes incompressible (cours §7.2). Trois parois fixes ($u=v=0$), couvercle
$y=L$ à la vitesse $U_t(t) = U(1-e^{-t/\tau})$, $\tau = L/U$. Un seul paramètre : $Re = UL/\nu$.

**Adimensionnement** : $L = U = \rho = 1$ ⇒ $\nu = 1/Re$ et $U_t(t) = 1 - e^{-t}$.

**Diagnostics imposés** (rien d'autre) :

- $K(t) = \tfrac12 \int_\Omega (u^2+v^2)\,dS$
- $P(t) = \int_0^1 \nu\,\partial_y u|_{y=1}\,U_t(t)\,dx$

**Ce que l'énoncé exige** :

- au moins **2 méthodes**, dont une du cours, et montrer l'effet d'un choix « sage » ou « malheureux » ;
- **A** : $Re = 100$ stationnaire, validé (Ghia + l'autre méthode), grille justifiée ;
- **B** : plusieurs $Re < 7500$, topologie, comparaison Ghia, $P(Re)$ + loi d'échelle ;
- **C** : $Re > Re_c \approx 8000$, le plus haut $Re$ abordable (justifié), sensibilité aux choix
  numériques, comparaison littérature. Un échec est acceptable s'il est expliqué.

---

## 2. Les deux méthodes (toutes deux dans le cours §7.3)

Même grille pour les deux : **nœuds colocalisés** $x_i = ih$, $y_j = jh$, $i,j = 0..N$, $h = 1/N$
(comme les workshops et les Projets 1 et 3). Les nœuds de bord sont sur les parois.

|                        | **M1 — Projection (pas fractionnaire)**                                                                                                           | **M2 — Fonction de courant / vorticité**                                                                                                            |
| ---------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Cours                  | §7.3.2 + workshop ch. 15 (cavité)                                                                                                                      | §7.3.4                                                                                                                                                     |
| Inconnues              | $u, v, p$                                                                                                                                              | $\omega, \psi$                                                                                                                                            |
| Un pas de temps        | 1) $u^* = u^n + \Delta t\,[-\text{adv} + \nu\Delta u]$  2) $\Delta p = \frac{1}{\Delta t}\nabla\cdot u^*$  3) $u^{n+1} = u^* - \Delta t\,\nabla p$ | 1) $\omega^{n+1} = \omega^n + \Delta t\,[-\text{adv} + \nu\Delta\omega]$  2) $\Delta\psi = -\omega$  3) $u = \partial_y\psi$, $v = -\partial_x\psi$ |
| Conditions aux limites | $u, v$ imposés aux parois ; $\partial p/\partial n = 0$ partout, couvercle compris (formule décentrée d'ordre 1, §4.2.2.2) ; $p$ à **moyenne nulle** (slides du prof) | $\psi = 0$ aux parois ; $\omega$ de paroi tiré d'un développement de Taylor de $\psi$ (outil du §2.1)                                              |

Briques communes, toutes du cours :

| Brique                                     | Choix                                                                                                               | Cours               |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------- | ------------------- |
| Temps                                      | **Euler explicite**                                                                                           | §3.1, workshop 15  |
| Diffusion $\Delta$                       | centrée ordre 2 (5 points)                                                                                         | §2.3, §4.2.1      |
| Advection                                  | **centrée ordre 2** *ou* **upwind ordre 1** (c'est le « choix sage / malheureux »)                 | §2.2, §5.6.2      |
| Poisson ($p$ ou $\psi$)                | **gradient conjugué sans matrice** (laplacien 5 points appliqué en tranches numpy, départ = solution du pas précédent) ; pour $p$ : **moyenne nulle** retirée à chaque itération (slides du prof) | §4.7.3, workshop 11.3 |
| Pas de temps                               | $\Delta t$ fixé par les conditions de stabilité CFL / Fourier                                                   | §5.3               |
| $\partial_y u$ au couvercle (pour $P$) | différence décentrée d'ordre 2 $(3f_j-4f_{j-1}+f_{j-2})/2h$                                                    | aide-mémoire §2.6 |

**Pas de** grille décalée (MAC), cellules fantômes, DCT/DST, RK3-SSP, QUICK, Numba, Richardson/GCI,
bilan d'énergie, FFT, tests unitaires : rien de tout cela n'est exigé.

### Pourquoi pas de cellules fantômes (contrairement au Projet 3)

Au Projet 3, les fantômes servaient à deux choses : loger des stencils larges ($i \pm 2$, schémas d'ordre 2
et 3) et imposer des conditions aux limites **exactes** par miroir, grâce à la symétrie du tourbillon de
Taylor-Green. Ici, aucune des deux raisons ne tient :

| Brique | Avec fantômes | Gain |
| --- | --- | --- |
| Laplacien 5 points, advection centrée / upwind 1 | stencils à $\pm 1$, nœuds de bord déjà sur les parois | aucun |
| $u, v$ (Dirichlet) | imposés directement sur les nœuds de paroi | aucun |
| $\omega$ de paroi (M2) | $\psi_{N+1} = \psi_{N-1} + 2hU_t$ puis laplacien centré ⇒ **on retombe exactement sur Thom** | aucun |
| $\partial_y u$ au couvercle (pour $P$) | déjà décentré d'ordre 2 | aucun |
| Neumann sur $p$ (M1) | miroir pair $p_{-1} = p_1$ ⇒ Neumann **centrée d'ordre 2** au lieu de l'ordre 1 | **seul gain possible** |
| Coins du couvercle (saut de vitesse) | inchangé | aucun (c'est pourtant la 1ʳᵉ source d'erreur sur $P$) |

Le seul gain (Neumann d'ordre 2 sur $p$) a un coût : la ligne de bord devient $2p_1 - 4p_0 + \dots$, ce qui
casse la symétrie exigée par le gradient conjugué (il faut pondérer les lignes de bord par ½, et ¼ aux coins),
et la moyenne nulle devient une moyenne pondérée (trapèzes). L'erreur d'ordre 1 ne concerne que la rangée
de nœuds collée aux parois : rien ne dit qu'elle limite la précision globale.
**Décision : pas de fantômes**, sauf si l'étude de convergence de l'étape 2 le justifie (voir ci-dessous).

### Singularités aux coins du couvercle (à garder en tête pour le maillage et pour $P(Re)$)

Aux deux coins du haut, la vitesse de paroi saute de $U_t$ à 0. Près du coin, l'écoulement est dominé par la
viscosité (problème du « grattoir » de Taylor, Moffatt — hors poly, référence à vérifier avant de la citer) :
$u$ reste borné mais n'a pas de limite, et $\partial_y u$, $\omega$, $p$ se comportent en $1/r$ ($r$ = distance au coin).

- **$P$ est infini en continu** : $\int dx/x$ diverge en logarithme (de même pour $\varepsilon = \nu\int\omega^2$).
  La grille coupe à $r \approx h$, donc en discret $P \approx P_\text{régulier} + c\,\nu\ln(1/h)$ :
  $P$ **ne converge pas** en maillage, il croît comme $\ln N$. On compare donc toujours les $P$ à $N$ fixé.
- **Biais sur $P(Re)$** : le terme de coin $\propto \nu = Re^{-1}$, le terme de couche limite $\propto Re^{-1/2}$
  ⇒ la pente mesurée sera un peu **plus forte que $-1/2$**, surtout aux $Re$ modérés.
- **Ordre de convergence** : la singularité peut faire baisser la pente de l'écart à Ghia **pour M1 et pour M2**
  (le test M1 vs M2 de l'étape 2 reste valable : le coin touche les deux, la Neumann d'ordre 1 seulement M1).
- **Schéma centré à haut $Re$** : au coin le gradient vaut $\sim U/h$ et le Reynolds de maille $Uh/\nu$ est maximal
  ⇒ risque d'oscillations nœud à nœud, qu'il ne faut pas prendre pour la bifurcation de Hopf.
- **Nœud de coin** : on prend $u = 0$ (il appartient à la paroi fixe). Le stencil 5 points des nœuds intérieurs ne
  l'utilise jamais ; ce choix ne touche que $P$, $\varepsilon$ et les tracés.
- **On ne régularise pas le couvercle** (du type $u = 16x^2(1-x)^2$) : ce serait un autre problème que celui de
  l'énoncé, et on ne pourrait plus comparer à Ghia. On garde la singularité et on explique ses effets.

### Stabilité (tout est dans le cours, c'est aussi la justification du $\Delta t$)

- Diffusion, Euler explicite 2D : $Fo = \nu\Delta t/h^2 \le 1/4$ (§5.3.3).
- Advection **upwind** + Euler : $C_x + C_y \le 1$ (§5.3.3).
- Advection **centrée** + Euler : instable seule (§5.3, ex. 9), stable seulement grâce à la diffusion,
  si $C^2 \le 2\,Fo \le 1$ (§5.3.1.4, équation d'advection-diffusion) ⇒ $\Delta t \lesssim 2\nu = 2/Re$.
  À haut $Re$ c'est très petit : c'est **le coût** de la méthode centrée.
- Upwind : diffusion numérique $D = \tfrac12 U h(1-C)$ (équation modifiée, §5.4.2.1). À $Re = 10^4$,
  $N = 128$ : $D \approx 4\cdot10^{-3} \gg \nu = 10^{-4}$ ⇒ le schéma « voit » un $Re$ effectif ~ 250.
  On s'attend donc à ce que **l'upwind tue l'instationnarité** : c'est notre « choix malheureux »,
  expliqué par le cours.

---

## 3. Plan de travail

### Étape 1 — Code commun + M1 (workshop 15)

- [ ] Grille, $U_t(t)$, laplacien et advection (centrée / upwind) en tranches numpy.
- [ ] Poisson par gradient conjugué sans matrice (Dirichlet pour $\psi$ ; Neumann + moyenne nulle pour $p$).
  Vérification rapide sur une solution connue : Dirichlet, Poisson du workshop 11 ($\phi = (x^2-1)(y^2-1)$) ;
  Neumann, $\phi = \cos(\pi x)\cos(\pi y)$ (dérivée normale nulle aux parois, moyenne nulle).
- [ ] Pas de temps M1, boucle en temps, critère stationnaire $\max|u^{n+1}-u^n| < 10^{-6}$ (même idée
  que le résidu des méthodes itératives, §4.3.2).
- [ ] Diagnostics $K$ (somme des trapèzes, workshop 8) et $P$.

### Étape 2 — Partie A ($Re = 100$)

- [ ] M1 avec $N = 32, 64, 128$ : profils $u(0.5, y)$ et $v(x, 0.5)$ vs Ghia.
- [ ] **Grille justifiée** : écart max à Ghia en fonction de $h$ (log-log, pente ≈ ordre du schéma,
  comme au workshop 9) ; on retient le $N$ où l'écart ne bouge plus.
- [ ] **Signal à surveiller** (après l'étape 3) : si la pente vaut ≈ 1 pour M1 mais ≈ 2 pour M2, le
  premier suspect est la Neumann d'ordre 1 sur $p$ ⇒ alors seulement, passer à la Neumann d'ordre 2 par
  miroir (cf. « Pourquoi pas de cellules fantômes ») et le justifier par cette mesure.
- [ ] $K(t)$, $P(t)$ jusqu'au stationnaire.
- [ ] **Justification du maillage par les coins** (sans calcul en plus : on réutilise $N = 32, 64, 128$) :
  tracer $P_\text{stationnaire}$ en fonction de $\ln N$. Attendu : une droite (croissance en $\ln N$, pas de
  palier) ⇒ $P$ ne sert **pas** à choisir $N$ ; on choisit $N$ sur les profils de Ghia et sur $K$, qui
  convergent, et on fige ce $N$ pour comparer les $P$.
- [ ] Lecture de la pente de convergence : si **M1 et M2** ont une pente faible, suspect = coins ; si seulement
  M1, suspect = Neumann d'ordre 1 (voir ci-dessous).

### Étape 3 — M2 et validation croisée

- [ ] M2 sur la même grille, mêmes briques.
- [ ] $Re = 100$ : M1, M2 et Ghia sur la même figure.

### Étape 4 — Partie B ($Re$ = 400, 1000, 3200, 5000 : colonnes disponibles chez Ghia)

- [ ] $N = 128$ (contrôle à $N = 64$ sur un cas). Lignes de courant (iso-$\psi$) : le tourbillon
  principal se recentre, des tourbillons de coin apparaissent puis grossissent.
- [ ] Profils vs Ghia.
- [ ] $P(Re)$ en log-log et pente par régression linéaire. Attendu : $P \propto Re^{-1/2}$ à grand $Re$
  (couche limite d'épaisseur $\delta \propto Re^{-1/2}$, cf. workshop 12.3 : $\delta_{99} \propto \sqrt{x}$).
- [ ] **Robustesse de la pente** : faire la régression de $P(Re)$ à $N = 64$ **et** à $N = 128$. Si la pente
  bouge peu, elle est fiable ; si elle est plus forte que $-1/2$ ou varie avec $N$, l'expliquer par le terme
  de coin $\propto Re^{-1}\ln N$ (cf. « Singularités aux coins du couvercle »).
- [ ] Remarque à faire dans le rendu : $P$ diverge en continu (coins), donc on compare toujours les $P$
  à $N$ fixé.

### Étape 5 — Partie C ($Re > 8000$)

- [ ] **Budget** : mesurer le temps d'un pas (gradient conjugué compris : noter le nombre moyen
  d'itérations par pas, faible grâce au départ à chaud), en déduire le temps de calcul pour simuler
  $t \approx 300$ avec le $\Delta t$ imposé par la stabilité ⇒ justifie le $Re$ (ex. 8500 ou 10000)
  et le $N$ (ex. 128) retenus. Ordre de grandeur : centré + Euler à $Re = 10^4$ impose
  $\Delta t \le 2\cdot10^{-4}$, soit ~$10^6$ pas. Si c'est trop cher, repli prévu par le cours :
  **RK3** (§3.9.4), stable pour l'advection centrée jusqu'à $C \le \sqrt3$ (§5.3.2.1), sur M2
  (où il s'applique directement à $\omega$).
- [ ] $K(t)$, $P(t)$ : stationnaire ou périodique ? Période mesurée entre deux maxima.
- [ ] **Sensibilité** (au même $Re$) : centré vs upwind, $N = 64$ vs $128$, M1 vs M2.
  Pour chacun : périodique ou non, période, amplitude. L'upwind doit rester stationnaire :
  on l'explique par $D$ (équation modifiée).
- [ ] **Avant de conclure « périodique »** avec le schéma centré : regarder $u$ ou $\omega$ près des coins du
  couvercle. Des oscillations nœud à nœud signalent une instationnarité **numérique**, pas la bifurcation.
- [ ] Tracés de $p$ et $\omega$ : fixer les bornes de l'échelle de couleur (pics en $1/h$ aux coins).
- [ ] Comparaison littérature : $Re_c \approx 8000$ (fourchette 7400–8050, cf. `ref_data/README.md`).

### Étape 6 — Rendu

- [ ] Figures propres, rédaction des choix et des échecs éventuels avec leur explication.

---

## 4. Fichiers

```
Projet4_CFD/
├── cavity/
│   ├── common.py        grille, U_t, opérateurs, Poisson (gradient conjugué), K et P
│   ├── projection.py    méthode M1
│   └── psi_omega.py     méthode M2
├── notebooks/
│   ├── A_Re100.ipynb    (contient aussi la vérification du Poisson)
│   ├── B_sweep_Re.ipynb
│   └── C_unsteady.ipynb
├── ref_data/            2 CSV de Ghia + README
├── docs/pseudocode.md
└── figures/, results/
```

Code en **numpy vectorisé par tranches** (`u[1:-1, 1:-1]`), très commenté, en français.

---

## 5. Questions ouvertes

- ~~Date et format de rendu ?~~ Réponse (slides du prof) : présentation + code à envoyer **avant le
  1ᵉʳ décembre** ; soutenance de **30 min** (10 min de présentation + 20 min de questions), le mercredi
  2 décembre (14h–17h) ou le jeudi 3 décembre (16h30–18h).
- Seul ou en groupe ?
- Peut-on laisser tourner un calcul ~1 h pour la Partie C ?
