# ADR 0002 — Condition de Neumann sur la pression (M1) : ordre 1, et comment passer à l'ordre 2

- **Statut** : accepté pour l'instant (2026-11-11, fin de l'étape 1) ; **à revoir** après la comparaison
  M1 / M2 de l'étape 3.
- **Mesure à l'origine** : `notebooks/A_Re100.ipynb` §1.c, `figures/A_poisson_convergence.png`,
  `results/A_poisson_verification.csv`.

## 1. Ce qu'on a mesuré

Vérification du Poisson discret sur des solutions connues : on choisit $\phi$, on donne $f = \Delta\phi$ au
solveur, on compare (gradient conjugué arrêté à $10^{-10}$, donc son erreur d'arrêt est négligeable).

| Cas | Pente observée | Lecture |
| --- | --- | --- |
| Dirichlet $(x^2-1)(y^2-1)$ sur $[-1,1]^2$ (workshop 11.1) | aucune (erreur $\sim 10^{-12}$) | le laplacien 5 points est **exact** sur ce polynôme ($\partial^4\phi = 0$) : le code est juste, mais ce cas ne mesure pas l'ordre |
| Dirichlet $\sin(\pi x)\sin(\pi y)$ | **2,00** | ordre 2 attendu |
| Neumann $\cos(\pi x)\cos(\pi y)$ | **≈ 1** (0,82 puis 0,92) | ordre 1, erreur $\approx 2.3h$ à $2.7h$ |

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

| M1 | M2 | Suspect | Action |
| --- | --- | --- | --- |
| ≈ 1 | ≈ 2 | Neumann d'ordre 1 sur $p$ | passer à l'option 2, justifié par la mesure |
| faible | faible | singularités des coins du couvercle | garder l'ordre 1, l'expliquer |
| ≈ 2 | ≈ 2 | — | garder l'ordre 1, l'expliquer à l'oral avec la figure du Poisson |
