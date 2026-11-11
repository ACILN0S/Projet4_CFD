# ADR 0002 — Condition de Neumann sur la pression (M1) : ordre 1, et comment passer à l'ordre 2

- **Statut** : accepté (2026-11-11). Option 2 (nœuds fantômes) appliquée (§6) ; elle rend M1 d'ordre 2 en espace (vérifié à
  $\Delta t$ fixe, §7.3) ; la pente ≈ 1 mesurée avec le $\Delta t$ de stabilité vient de l'erreur de projection, qui dépend de $\Delta t$ (identité
  exacte, §7.3 ; ordre ≈ 0.65 en $\Delta t$). Décision §7.4 : options 1 et 2, réalisées.
- **Mesure à l'origine** : `notebooks/A_Re100.ipynb` §1.c, `figures/A_poisson_convergence.png`,
  `results/A_poisson_verification.csv`.

## 1. Ce qu'on a mesuré

Vérification du Poisson discret sur des solutions connues : on choisit $\phi$, on donne $f = \Delta\phi$ au
solveur, on compare (gradient conjugué arrêté à $10^{-10}$, donc son erreur d'arrêt est négligeable).

| Cas                                                           | Pente observée                   | Lecture                                                                                                                                    |
| ------------------------------------------------------------- | --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| Dirichlet $(x^2-1)(y^2-1)$ sur $[-1,1]^2$ (workshop 11.1) | aucune (erreur $\sim 10^{-12}$) | le laplacien 5 points est **exact** sur ce polynôme ($\partial^4\phi = 0$) : le code est juste, mais ce cas ne mesure pas l'ordre |
| Dirichlet $\sin(\pi x)\sin(\pi y)$                          | **2,00**                    | ordre 2 attendu                                                                                                                            |
| Neumann $\cos(\pi x)\cos(\pi y)$                            | **≈ 1** (0,82 puis 0,92)   | ordre 1, erreur $\approx 2.3h$ à $2.7h$                                                                                               |

## 2. Pourquoi l'ordre 1

- **Deux choses distinctes.** La **discrétisation** (stencil 5 points + condition aux limites) construit le
  système $A\phi = b$ ; le **gradient conjugué** (cours §4.7.3) ne fait que le résoudre, sur tous les nœuds
  intérieurs. La condition aux limites n'entre dans le gradient conjugué qu'à travers $A$ : elle change les
  équations de la première rangée de nœuds intérieurs, dont le stencil touche un nœud de bord.
  L'ordre 1 vient donc de la condition aux limites, **pas** du gradient conjugué.
- **La règle utilisée** : $\phi_0 = \phi_1$ (Neumann décentrée d'ordre 1, cours §4.2.2.2, éq. 4.27-4.28).
  Elle impose $\partial\phi/\partial n = 0$ **entre** les nœuds 0 et 1 : le schéma résout (à l'ordre 2) le problème
  sur le domaine décalé $[h/2,\ 1-h/2]^2$. Estimation 1D : erreur $\approx \tfrac{\pi^2}{4}h \approx 2.5\,h$, du même
  ordre que la mesure.
- **Une erreur née au bord touche tout le volume.** Poisson est elliptique : chaque point dépend de
  tout le domaine et de ses conditions aux limites. Sur la figure, l'erreur / $h$ des trois grilles se
  superpose et forme une pente qui traverse tout le domaine, pas un défaut localisé près des parois.

## 3. Les deux façons de passer à l'ordre 2 (toutes deux dans le cours)

Contrainte commune : le gradient conjugué exige une matrice **symétrique** (définie positive, cours §4.7.3).

### Option 1 — Décentré à 3 points (cours §4.2.2.2, éq. 4.29-4.30)

- Règle : $\phi_0 = (4\phi_1 - \phi_2)/3$ au lieu de $\phi_0 = \phi_1$.
- **Casse la symétrie** : l'équation du nœud 1 reçoit un coefficient sur $\phi_2$, mais pas l'inverse
  ⇒ gradient conjugué non garanti.
- **Casse la compatibilité simple** : « somme du second membre = 0 » ne suffit plus, il faudrait une
  pondération particulière.
- Il faudrait changer de solveur : SOR (cours §4.5.3) accepte une matrice non symétrique, mais coûte
  $O(n^2)$ contre $O(n^{3/2})$ (cours §4.8).

### Option 2 — Nœuds fantômes, miroir pair $\phi_{-1} = \phi_1$ (cf. PLAN.md, « Pourquoi pas de cellules fantômes »)

- Les nœuds de bord deviennent des inconnues ; leur équation est la dérivée normale centrée d'ordre 2 :
  $(2\phi_1 - 4\phi_0 + \phi_{0,j+1} + \phi_{0,j-1})/h^2$.
- Le facteur 2 casse aussi la symétrie, **mais** on la rétablit en pondérant les lignes de bord par ½ et
  celles des coins par ¼ ⇒ **le gradient conjugué reste utilisable**.
- Coût : la moyenne nulle devient une moyenne pondérée (poids des trapèzes), pour le second membre comme
  pour $p$. Quelques lignes de code, aucune idée nouvelle.

**Si on renforce : option 2**, car elle garde le gradient conjugué et ses $O(N)$ itérations.

## 4. Décision

On **garde l'ordre 1** pour l'instant (strict minimum) : un $p$ d'ordre 1 ne fait pas forcément des
vitesses d'ordre 1, et ce sont les vitesses (profils vs Ghia) qui comptent pour la Partie A.

Critère de bascule (étapes 2 et 3, pente de l'écart aux profils de vitesse) :

| M1     | M2     | Suspect                              | Action                                                            |
| ------ | ------ | ------------------------------------ | ----------------------------------------------------------------- |
| ≈ 1   | ≈ 2   | Neumann d'ordre 1 sur $p$          | passer à l'option 2, justifié par la mesure                     |
| faible | faible | singularités des coins du couvercle | garder l'ordre 1, l'expliquer                                     |
| ≈ 2   | ≈ 2   | —                                   | garder l'ordre 1, l'expliquer à l'oral avec la figure du Poisson |

## 5. Mesure de l'étape 3 (2026-11-11, `notebooks/A_Re100.ipynb` §3.c, $Re = 100$, $N = 32, 64, 128$)

Pente de l'écart entre grilles successives $\max|u_N - u_{2N}|$ sur les lignes médianes :

| | M1 (projection) | M2 ($\psi$–$\omega$, Thom) |
| --- | --- | --- |
| $u$ | 1.37 | **2.07** |
| $v$ | 1.42 | **2.02** |

- M2 atteint l'ordre 2 malgré les coins et Thom ⇒ les coins ne limitent pas l'ordre sur les profils.
- M1 seul est en dessous (≈ 1.4, entre 1 et 2) ⇒ **ligne 1 du tableau** : suspect = Neumann d'ordre 1 sur $p$.
- $\max|M1 - M2|$ diminue avec une pente ≈ 1.3–1.4 : les deux méthodes vont vers la même solution, M1 plus lentement.
- **Décision (utilisateur, 2026-11-11)** : appliquer l'option 2 (nœuds fantômes). Voir §6 et §7.

## 6. Option 2 mise en œuvre : ce qu'elle implique (2026-11-11)

Code : `cavity/common.py` (`solve_poisson(..., bc="neumann_ghost")`), `cavity/projection.py` (`divergence_all_nodes`,
`pressure_bc="neumann_ghost"`, défaut de `run_projection`). L'ancienne Neumann d'ordre 1 reste disponible
(`bc="neumann"`, `pressure_bc="neumann"`) pour la comparaison.

### 6.1 Le système de pression

- **Inconnues** : $p$ à **tous** les nœuds, parois et coins compris ($(N+1)^2$ inconnues au lieu de $(N-1)^2$).
- **Équation à chaque nœud** : laplacien 5 points ; aux parois, le voisin manquant est un **fantôme miroir pair**
  $p_{-1,j} = p_{1,j}$, $p_{N+1,j} = p_{N-1,j}$ (et de même en $j$). Cela impose $\partial p/\partial n = 0$ par la dérivée **centrée**
  d'ordre 2 (cours §2.2, éq. 2.9). Ligne de paroi (gauche) : $(2p_{1,j} + p_{0,j+1} + p_{0,j-1} - 4p_{0,j})/h^2$ ; coin :
  $(2p_{1,0} + 2p_{0,1} - 4p_{0,0})/h^2$.
- **Symétrie (exigée par le gradient conjugué)** : le facteur 2 casse la symétrie ; on multiplie chaque ligne par le poids des
  trapèzes $w_i w_j$ (½ sur une paroi, ¼ dans un coin, workshop 8.2). On résout $A p = b$ avec $A = -W\Delta_h$ symétrique
  semi-définie positive et $b = -W\,\text{rhs}$.
- **Existence (compatibilité)** : $A\cdot\mathbf 1 = 0$ ⇒ il faut $\sum w_{ij}\,\text{rhs}_{ij} = 0$, c'est-à-dire $\iint \text{rhs}\,dS = 0$ par les trapèzes :
  la version discrète de la compatibilité continue $\iint \Delta p\,dS = \oint \partial p/\partial n = 0$. On retire donc la **moyenne pondérée**
  du second membre (et non plus la moyenne simple).
- **Unicité** : $p$ à **moyenne pondérée nulle** ($\iint p\,dS = 0$) à chaque itération (slides du prof : « subtract mean(P) », ici
  la moyenne des trapèzes, cohérente avec les inconnues de paroi).
- **Vérification** (`notebooks/A_Re100.ipynb` §1.c, $\phi = \cos\pi x\cos\pi y$) : pente **2.00** sur $N = 16 \to 128$, erreur $\approx 0.82\,h^2$, soit
  ~200 fois plus petite que l'ordre 1 à $N = 64$. Coût : ~40 % d'itérations du gradient conjugué en plus.

### 6.2 Le second membre aux parois (choix validé par l'utilisateur : « option 1 »)

- Les nœuds de paroi étant des inconnues, il faut $\text{rhs} = (\nabla\cdot u^*)/\Delta t$ **aussi sur les parois**. Ce n'est pas un détail :
  la ligne de paroi pèse $\tfrac12 h$ ; une erreur $O(1)$ sur son second membre agit comme une source sur une couche
  d'épaisseur $h$ ⇒ erreur $O(h)$ sur $p$, ce qui ramènerait l'ordre 1.
- **Formule retenue** : chaque dérivée est centrée là où les deux voisins existent, **décentrée d'ordre 2** aux extrémités
  ($(-3f_0 + 4f_1 - f_2)/2h$, workshop 11.2 éq. 11.10, même formule que pour $P$), soit `np.gradient(..., edge_order=2)`.
  Le long d'une paroi : dérivée tangentielle centrée (valeurs de paroi connues), dérivée normale décentrée ; aux coins,
  les deux décentrées. Aux nœuds intérieurs, rien ne change.
- **Écartées** : second membre nul sur les parois (l'écart de $u^*$ à la divergence nulle, $O(\Delta t)$, divisé par $\Delta t$ donne une
  erreur $O(1)$ sur la ligne de paroi) ; fantômes de vitesse (il faudrait choisir la parité de chaque composante).
- **Coins du couvercle** : $u$ saute de $U_t$ à 0 ⇒ la divergence vaut $\sim U_t/h$ aux deux coins et à leurs voisins sur le
  couvercle. Ces nœuds entrent désormais dans l'équation de pression (poids ¼ et ½) : c'est la singularité de coin, acceptée
  (PLAN.md §2) ; la compatibilité est rétablie par la moyenne pondérée. Conséquence visible : $p$ aux coins croît en $1/h$
  ($-1.6, -2.6, -5.5$ au coin gauche pour $N = 32, 64, 128$).
- **Correction de vitesse** : inchangée ($u^{n+1} = u^* - \Delta t\,\nabla p$ aux nœuds intérieurs, gradient centré) ; mais $p$ aux nœuds de
  paroi, utilisé par le gradient des nœuds voisins des parois, est maintenant une vraie inconnue d'ordre 2.

## 7. Résultat sur M1 à $Re = 100$ (2026-11-11) : pente ≈ 1 avec le $\Delta t$ de stabilité, ordre 2 à $\Delta t$ fixe

### 7.1 Mesure

| | M1 Neumann ordre 1 | M1 fantômes | M2 |
| --- | --- | --- | --- |
| $\max|u_{32} - u_{64}|$ | 0.0289 | 0.0056 | 0.0107 |
| $\max|u_{64} - u_{128}|$ | 0.0112 | 0.0028 | 0.0026 |
| pente (u / v) | 1.37 / 1.42 | **1.0 / 0.91** | 2.07 / 2.02 |
| $\max|M1 - M2|$, $N = 32, 64, 128$ (u) | 0.035, 0.015, 0.006 | 0.020, 0.0086, 0.0038 | — |

Les fantômes divisent l'erreur de M1 par ~3 à 5, mais **la pente tombe à ≈ 1**.

### 7.2 Diagnostics (non adoptés, scripts hors projet)

1. **Coins du couvercle ?** Calcul de la divergence avec $u = U_t$ aux deux coins (pas de saut le long du couvercle) : pente
   ≈ 1.1, erreurs plus grandes ⇒ hypothèse **non confirmée**.
2. **Dépendance de l'état stationnaire à $\Delta t$** : on refait $N = 32, 64$ avec $\Delta t/2$ :

   | | $N = 32$ | $N = 64$ |
   | --- | --- | --- |
   | M1 ordre 1 : $\max|\Delta$profils$|$ | $1.2\cdot10^{-2}$ | $4.8\cdot10^{-3}$ |
   | M1 fantômes | $6.0\cdot10^{-3}$ | $2.4\cdot10^{-3}$ |
   | M2 | $1.7\cdot10^{-7}$ | $1.1\cdot10^{-7}$ |

   ⇒ **l'état stationnaire de M1 dépend de $\Delta t$, celui de M2 non** (M2 : erreur d'arrêt seulement). Avec les fantômes, cette
   dépendance est **du même ordre que l'écart entre grilles** : l'étude de grille de M1 mesure en grande partie un effet de $\Delta t$.
   À $N$ fixé, elle est à peu près proportionnelle à $\Delta t$.

### 7.3 Ce qui est démontré, ce qui est réfuté (2026-11-11)

**Démontré (identité exacte).** Notations : $F$ advection + diffusion, $G$ gradient centré aux nœuds intérieurs et $\tilde G p$ ce
gradient prolongé par 0 aux parois (la correction n'y est pas appliquée), $D$ la divergence (`divergence_all_nodes`), $\Delta_h$ le
laplacien à fantômes que résout le gradient conjugué. Un pas : $u^* = u + \Delta t F(u)$ (intérieur), $\Delta_h p = D(u^*)/\Delta t$,
$u^{n+1} = u^* - \Delta t\,\tilde G p$. Au stationnaire $u^{n+1} = u^n = u$ : $F(u) = Gp$ (la quantité de mouvement stationnaire est exacte, sans $\Delta t$),
donc $u^* = u + \Delta t\,\tilde G p$, et l'équation de pression donne
$$D(u) = \Delta t\,(\Delta_h - D\tilde G)\,p .$$
La divergence stationnaire est **proportionnelle à $\Delta t$** et à l'écart entre le laplacien résolu et celui qu'applique la correction.
Vérifié numériquement : $\max|D(u) - \Delta t\,E p| \sim 10^{-4}$ contre $\max|D(u)| \sim 10^{2}$ (erreur d'arrêt). M2 n'a pas d'étape de projection :
son stationnaire ne dépend pas de $\Delta t$ (mesuré : $10^{-7}$). C'est une propriété de la projection **non incrémentale du cours** (§7.3.2).

**Réfuté : l'explication « $\Delta t/h$, donc ordre 1, à cause de $\partial p/\partial n = 0$ non physique ».** Elle prédisait $E p \approx c\,g/h$ sur les parois
($g$ = gradient normal de pression au premier nœud intérieur). Mesure (parois, tiers central) : $E p \approx$ 1.2, 1.3, 1.9 pour
$N = 32, 64, 128$, à peu près constant, alors que $g/h$ double ; $E p/(g/h)$ = 0.70, 0.37, 0.27. Le plus grand $E p$ est aux **coins du
couvercle**. Et l'erreur ne varie pas en $\Delta t/h$ : à $N$ fixé, une division de $\Delta t$ par 2 donne la même variation relative à $N = 32$
et à $64$ ($\approx 1.0\,\Delta t$ fantômes, $1.9\,\Delta t$ ordre 1). (Précision ultérieure, notebook A §2.e : avec trois $\Delta t$ à $N = 64$, l'ordre observé
en $\Delta t$ est **≈ 0.65**, pas 1 ; l'erreur n'est donc pas exactement $C\,\Delta t$.)
**On ne peut donc pas affirmer que la Neumann non physique est la cause** ; on sait seulement que c'est l'écart $\Delta_h - D\tilde G$, qui mêle
le laplacien compact/large de la grille colocalisée, le traitement des parois et les coins. Lequel domine : non déterminé.

**D'où venait la pente ≈ 1 (diagnostic 3, $\Delta t$ fixé à celui de $N = 128$ pour les trois grilles)** :

| | pente $u$ | pente $v$ |
| --- | --- | --- |
| M1 fantômes, $\Delta t$ fixe | **1.91** | 3.53 (écarts ~$5\cdot10^{-4}$, peu fiable) |
| M1 ordre 1, $\Delta t$ fixe | 1.66 | 1.55 |

- Avec $\Delta t$ fixe, la part qui dépend de $\Delta t$ ($E_{\Delta t}$) est la même sur toutes les grilles et disparaît des écarts entre grilles :
  **les fantômes rendent bien M1 d'ordre 2 en espace** ; l'ordre 1 de la Neumann d'origine se voit encore (1.6).
- Avec le $\Delta t$ de stabilité, $\Delta t$ est divisé par 2.56 puis par 4 (CFL à $N = 32$, Fourier ensuite) et $E_{\Delta t}$ change d'une grille à l'autre ;
  comme elle domine les écarts entre grilles, la pente mesurée (≈ 1) est un **artefact de la suite des $\Delta t$**, pas un ordre 1 en $h$.
- **Convergence en $\Delta t$** (notebook A §2.e, $N = 64$, $\Delta t_0$, $\Delta t_0/2$, $\Delta t_0/4$) : écarts $1.8\cdot10^{-3}$ puis $1.1\cdot10^{-3}$, rapport 1.55 ⇒ ordre
  observé en $\Delta t$ **≈ 0.65** (pas 1 : le facteur $(\Delta_h - D\tilde G)p$ dépend lui-même de $\Delta t$ ; non expliqué plus avant). M2 : $10^{-7}$.
- À $N = 128$, l'écart M1 − M2 ($\approx 4\cdot10^{-3}$) est surtout cette part $E_{\Delta t}$ de M1.

### 7.4 Options — **décision (utilisateur, 2026-11-11) : 1 et 2**

Réalisé : M1 garde les fantômes (défaut de `run_projection`) ; le notebook A §2.e chiffre l'erreur de projection (convergence en
$\Delta t$, étude de grille à $\Delta t$ fixe) et le §3.c compare les trois variantes (M1 ordre 1, M1 fantômes, M2). `run_projection` et
`run_psi_omega` acceptent un $\Delta t$ imposé (plus petit que la limite de stabilité) pour ces études.

1. **Garder les fantômes et l'expliquer** : M1 est d'ordre 2 en espace (vérifié à $\Delta t$ fixe), plus une erreur de projection qui dépend
   de $\Delta t$, propre à la méthode du cours, qui domine à $N = 128$ ; M2 sert de référence sans cette erreur. Strict minimum, tout dans le cours.
2. **Chiffrer l'erreur en $\Delta t$ dans le notebook** (convergence en $\Delta t$ à $N$ fixé, et étude de grille à $\Delta t$ fixe comme ci-dessus) :
   c'est ce qui permet de séparer proprement les deux erreurs à l'oral. Reste dans le cours ; coût modéré.
3. **Projection incrémentale** (on garde $-\nabla p^n$ dans la prédiction et on résout le Poisson pour l'incrément) : au stationnaire,
   l'incrément est nul et $D(u^*) = 0$ avec $u^* = u$, donc la part $E_{\Delta t}$ disparaîtrait. **Hors du cours** : seulement si l'utilisateur l'autorise.
   (Non retenue.)
4. Revenir à l'ordre 1 : sans intérêt (moins précis, ordre spatial ≈ 1.6 à $\Delta t$ fixe).
