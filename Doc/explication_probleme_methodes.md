# Projet 4 — Le problème et les deux méthodes (explication)

> Fiche de révision des parties 1 et 2 de [PLAN.md](../PLAN.md). Les références §x.y renvoient au poly
> *Numerical Methods* (Vicquelin & Vié, 2025-26).

---

## 1. Le problème

On simule le fluide d'une boîte carrée dont le couvercle glisse. Un seul paramètre, $Re$, fixe le
régime. On suit l'écoulement avec deux grandeurs : $K$ (énergie cinétique) et $P$ (puissance injectée).

### 1.1 La situation physique

La boîte est carrée, de côté $L$, et pleine de fluide. Trois parois restent immobiles : gauche, droite et
fond. La paroi du haut, le **couvercle**, glisse vers la droite. Par frottement visqueux, il entraîne le
fluide juste en dessous, qui entraîne ses voisins. Il se forme un grand **tourbillon** qui tourne dans le
sens horaire.

```
        couvercle  →  →  →   (vitesse U_t(t))
      ┌───────────────────┐
      │    ↗  ────→  ↘    │
      │   ↑            ↓  │
      │   ↑     ⟳      ↓  │   parois fixes : u = v = 0
      │    ↖  ←────  ↙    │
      └───────────────────┘
```

Le couvercle accélère progressivement, avec $\tau = L/U$ :

$$U_t(t) = U\left(1 - e^{-t/\tau}\right)$$

Sa vitesse vaut 0 à $t = 0$ et tend vers $U$ au bout de quelques $\tau$. Ce démarrage doux évite un choc
numérique au premier pas de temps.

### 1.2 Les équations (§7.2)

C'est Navier–Stokes **incompressible** en 2D, avec la vitesse $\mathbf{u} = (u, v)$ et la pression $p$ :

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u} = -\nabla p + \nu\,\Delta\mathbf{u},
\qquad \nabla\cdot\mathbf{u} = 0$$

- $\partial \mathbf{u}/\partial t$ : la variation dans le temps.
- $(\mathbf{u}\cdot\nabla)\mathbf{u}$ : l'**advection**, c'est-à-dire le fluide qui transporte sa propre
  vitesse. C'est le terme non linéaire, et le plus difficile à traiter.
- $-\nabla p$ : la pression qui pousse le fluide.
- $\nu\,\Delta\mathbf{u}$ : la **diffusion** visqueuse, qui lisse les écarts de vitesse.

**Point clé : la pression n'a pas d'équation d'évolution.** Elle prend à chaque instant la valeur qui
garantit $\nabla\cdot\mathbf{u} = 0$. Toute la difficulté numérique vient de là, et les deux méthodes de la
partie 2 sont deux façons de traiter cette contrainte.

### 1.3 Adimensionnement : un seul paramètre, $Re$

On prend $L = U = \rho = 1$. Il reste un seul nombre sans dimension :

$$Re = \frac{UL}{\nu} \;\Rightarrow\; \nu = \frac{1}{Re}, \qquad U_t(t) = 1 - e^{-t}$$

$Re$ compare l'advection, qui crée des structures fines, à la diffusion, qui les efface.

| $Re$ | Régime attendu |
| --- | --- |
| 100 | Stationnaire et lisse ; tourbillon un peu décalé vers le haut à droite |
| 400 à 5000 | Stationnaire ; le tourbillon se recentre, des tourbillons secondaires apparaissent dans les coins du bas |
| $> Re_c \approx 8000$ | La solution stationnaire devient instable et oscille périodiquement (bifurcation de Hopf) : Partie C |

### 1.4 Les deux diagnostics

**Énergie cinétique totale :**

$$K(t) = \frac{1}{2}\int_\Omega \left(u^2 + v^2\right)dS$$

$K$ part de 0 et monte. Ensuite, soit il devient constant (régime stationnaire), soit il oscille (régime
périodique). C'est le signal le plus simple pour voir l'instationnarité.

**Puissance injectée par le couvercle :**

$$P(t) = \int_0^1 \nu \left.\frac{\partial u}{\partial y}\right|_{y=1} U_t(t)\,dx$$

$\nu\,\partial u/\partial y$ est la contrainte de frottement au couvercle (une force par unité de longueur).
On la multiplie par la vitesse du couvercle pour obtenir une puissance, puis on intègre sur sa longueur.

**Lien entre les deux** : $dK/dt = P - \text{dissipation visqueuse}$. Le couvercle apporte de l'énergie et la
viscosité la dissipe. En régime stationnaire, les deux s'équilibrent.

### 1.5 Ce que l'énoncé demande

- **Deux méthodes**, dont au moins une du cours. Il faut aussi montrer qu'un choix numérique peut être
  « sage » ou « malheureux ».
- **Partie A ($Re = 100$)** : valider contre Ghia et al. (1982) et contre l'autre méthode, et justifier la
  grille.
- **Partie B (plusieurs $Re < 7500$)** : décrire la topologie (tourbillons de coin), comparer à Ghia, tracer
  $P(Re)$ et trouver une loi de puissance.
- **Partie C ($Re > 8000$)** : aller le plus haut possible en justifiant la limite, et tester la sensibilité
  aux choix numériques. Un échec est accepté s'il est expliqué.

---

## 2. Les deux méthodes

Les deux méthodes utilisent la même grille et les mêmes briques. Elles ne diffèrent que par le
traitement de la pression : **M1 la calcule, M2 l'élimine**.

### 2.0 La grille commune

Les nœuds sont régulièrement espacés : $x_i = ih$, $y_j = jh$, pour $i, j = 0\dots N$, avec $h = 1/N$.
Cela fait $(N+1)\times(N+1)$ points.

- Les nœuds $i = 0$, $i = N$, $j = 0$ et $j = N$ sont **sur les parois**. Les conditions aux limites
  s'y imposent directement.
- **Colocalisée** veut dire que $u$, $v$ et $p$ sont rangés aux mêmes points. C'est la grille des
  workshops ; on n'utilise pas de grille décalée (MAC), qui n'est pas dans le cours.
- Dans le code, `f[i, j]` correspond à $f(x_i, y_j)$. Les nœuds intérieurs s'écrivent `f[1:-1, 1:-1]`.

---

### 2.1 Méthode M1 : projection, ou pas fractionnaire (§7.3.2, workshop 15)

**Idée** : à chaque pas de temps, on avance d'abord la vitesse **sans la pression**, puis on la
**corrige** pour la rendre à divergence nulle.

1. **Prédiction.** On calcule une vitesse provisoire $\mathbf{u}^*$ en ignorant la pression :
   $$\mathbf{u}^* = \mathbf{u}^n + \Delta t\left[-(\mathbf{u}^n\cdot\nabla)\mathbf{u}^n + \nu\,\Delta\mathbf{u}^n\right]$$
   $\mathbf{u}^*$ tient compte de l'advection et de la viscosité, mais sa divergence n'est **pas nulle** en
   général.
2. **Équation de Poisson pour la pression.** On veut $\mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t\,\nabla p$
   avec $\nabla\cdot\mathbf{u}^{n+1} = 0$. On prend la divergence des deux côtés :
   $$0 = \nabla\cdot\mathbf{u}^* - \Delta t\,\Delta p \quad\Rightarrow\quad \Delta p = \frac{1}{\Delta t}\,\nabla\cdot\mathbf{u}^*$$
   Le second membre est connu, on résout pour obtenir $p$.
3. **Correction.** $\mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t\,\nabla p$. Le gradient de pression retire
   exactement la partie de $\mathbf{u}^*$ qui violait l'incompressibilité. On « projette » $\mathbf{u}^*$
   sur l'ensemble des champs à divergence nulle, d'où le nom.

**Conditions aux limites**

- **Vitesse** : $u = v = 0$ sur les parois fixes ; $u = U_t(t)$ et $v = 0$ sur le couvercle.
- **Pression** : la physique ne fournit pas de condition naturelle. On prend $\partial p/\partial n = 0$,
  discrétisé par une formule décentrée d'ordre 1, par exemple $p_0 = p_1$ (§4.2.2.2).
- **Subtilité** : avec du Neumann partout, $p$ n'est défini **qu'à une constante près**, et la matrice est
  singulière. On **fixe la valeur de $p$ en un point** pour lever l'indétermination.

---

### 2.2 Méthode M2 : fonction de courant et vorticité (§7.3.4)

**Idée** : en 2D, on peut **faire disparaître la pression** et satisfaire l'incompressibilité
automatiquement, en changeant d'inconnues.

**La fonction de courant $\psi$.** On pose

$$u = \frac{\partial\psi}{\partial y}, \qquad v = -\frac{\partial\psi}{\partial x}$$

Alors $\nabla\cdot\mathbf{u} = \partial^2\psi/\partial x\partial y - \partial^2\psi/\partial y\partial x = 0$
**automatiquement**, quel que soit $\psi$. En bonus, les lignes iso-$\psi$ **sont les lignes de courant**,
ce qui donne directement les tracés demandés en Partie B.

**La vorticité** $\omega = \partial v/\partial x - \partial u/\partial y$ mesure la rotation locale du
fluide. On prend le rotationnel de Navier–Stokes. Comme le rotationnel d'un gradient est nul, **le terme
$\nabla p$ disparaît**, et il reste une simple équation d'advection-diffusion :

$$\frac{\partial\omega}{\partial t} + u\,\frac{\partial\omega}{\partial x} + v\,\frac{\partial\omega}{\partial y} = \nu\,\Delta\omega$$

En remplaçant $u$ et $v$ par leur expression en $\psi$ dans la définition de $\omega$, on obtient le lien
entre les deux :

$$\Delta\psi = -\omega$$

**Un pas de temps**

1. On avance la vorticité : $\omega^{n+1} = \omega^n + \Delta t\,[-\text{adv} + \nu\,\Delta\omega^n]$.
2. On résout le Poisson $\Delta\psi = -\omega^{n+1}$.
3. On en déduit $u = \partial\psi/\partial y$ et $v = -\partial\psi/\partial x$ par différences centrées.

**Conditions aux limites (la partie délicate)**

- **Pour $\psi$**, c'est facile : aucun fluide ne traverse les parois, donc $\psi$ est constant sur tout
  le bord, et on prend $\psi = 0$.
- **Pour $\omega$**, on ne connaît pas sa valeur à la paroi. On la reconstruit avec un **développement de
  Taylor** de $\psi$, de la paroi vers le premier nœud intérieur (outil du §2.1). Prenons le couvercle
  ($j = N$) :
  $$\psi_{N-1} = \psi_N - h\,\frac{\partial\psi}{\partial y} + \frac{h^2}{2}\,\frac{\partial^2\psi}{\partial y^2} + \dots$$
  On sait que $\psi_N = 0$ et que $\partial\psi/\partial y = u = U_t$. De plus, $\psi$ est constant le long
  de la paroi, donc $\partial^2\psi/\partial x^2 = 0$ et $\omega_\text{paroi} = -\partial^2\psi/\partial y^2$.
  On isole $\omega$ et on obtient la **formule de Thom** :
  $$\omega_\text{couvercle} = -\frac{2\left(\psi_{N-1} + h\,U_t\right)}{h^2}$$
  Sur les parois fixes, le terme $U_t$ disparaît : $\omega_\text{paroi} = -2\,\psi_\text{voisin}/h^2$.

**Comparaison M1 / M2**

| | M1 (projection) | M2 ($\psi$–$\omega$) |
| --- | --- | --- |
| Pression | calculée à chaque pas | éliminée |
| Incompressibilité | imposée par la projection | automatique |
| Point délicat | condition aux limites sur $p$ | $\omega$ à la paroi (Thom) |
| Généralisable en 3D | oui | non (spécifique au 2D) |

Deux méthodes **construites très différemment** qui donnent le même résultat : c'est une validation croisée
solide.

---

### 2.3 Les briques communes

**Le temps : Euler explicite (§3.1).** $f^{n+1} = f^n + \Delta t \times (\text{second membre évalué en } n)$.
C'est le schéma le plus simple, mais il impose une limite sur $\Delta t$ (voir 2.4).

**La diffusion : schéma centré à 5 points (§2.3).** C'est un schéma d'ordre 2, celui de tous les
workshops :

$$\Delta f \approx \frac{f_{i+1,j} + f_{i-1,j} + f_{i,j+1} + f_{i,j-1} - 4f_{i,j}}{h^2}$$

**L'advection : deux options, notre « choix sage ou malheureux » (§2.2, §5.6.2).**

- **Centrée d'ordre 2** : $u\,\partial f/\partial x \approx u\,(f_{i+1} - f_{i-1})/2h$. Elle est précise et
  n'introduit pas de diffusion parasite, mais elle est délicate du point de vue de la stabilité.
- **Upwind d'ordre 1** : on prend le point situé du côté d'où vient le fluide. Si $u > 0$,
  $u\,\partial f/\partial x \approx u\,(f_i - f_{i-1})/h$. Elle est très robuste, mais introduit une
  **viscosité numérique**, ce qui est grave à haut $Re$ (voir 2.4).

**Le Poisson (pour $p$ ou pour $\psi$) : solveur direct creux (§4.3.1).** On discrétise $\Delta$ avec le
schéma à 5 points, ce qui donne un grand système linéaire $A\,x = b$. La matrice $A$ est creuse : au plus 5
coefficients non nuls par ligne. Surtout, **$A$ ne change jamais**, puisque la grille et les conditions aux
limites restent les mêmes. On la **factorise une seule fois** (décomposition LU, `scipy.sparse.linalg.splu`).
Ensuite, à chaque pas, seul le second membre $b$ change, et la résolution se réduit à deux remontées
triangulaires, ce qui est très rapide.

**$\partial u/\partial y$ au couvercle, pour calculer $P$.** On ne peut pas utiliser une différence
centrée, puisqu'il n'y a pas de nœud au-dessus du couvercle. On prend une formule **décentrée d'ordre 2**
qui n'utilise que les points en dessous :

$$\left.\frac{\partial u}{\partial y}\right|_N \approx \frac{3u_N - 4u_{N-1} + u_{N-2}}{2h}$$

---

### 2.4 La stabilité : comment on choisit $\Delta t$ (§5.3)

On définit deux nombres sans dimension :

- le nombre de **Courant** $C = U\Delta t/h$ : la fraction de maille parcourue par le fluide en un pas ;
- le nombre de **Fourier** $Fo = \nu\Delta t/h^2$ : la même idée, pour la diffusion.

| Cas | Condition de stabilité | Conséquence sur $\Delta t$ |
| --- | --- | --- |
| Diffusion, Euler, 2D (§5.3.3) | $Fo \le 1/4$ | $\Delta t \le h^2/(4\nu)$ |
| Advection upwind, Euler (§5.3.3) | $C_x + C_y \le 1$ | le fluide traverse au plus une maille par pas |
| Advection centrée, Euler, seule (§5.3, ex. 9) | **toujours instable** ($\lvert G\rvert > 1$) | — |
| Advection centrée + diffusion, Euler (§5.3.1.4) | $C^2 \le 2\,Fo \le 1$ | $\Delta t \le 2\nu/U^2 = 2/Re$ |

**Le piège du schéma centré.** On développe $C^2 \le 2Fo$ :

$$\frac{U^2\Delta t^2}{h^2} \le \frac{2\nu\Delta t}{h^2} \quad\Rightarrow\quad \Delta t \le \frac{2\nu}{U^2} = \frac{2}{Re}$$

Les $h^2$ se simplifient : la limite **ne dépend pas de la grille**. À $Re = 10^4$, on a
$\Delta t \le 2\cdot10^{-4}$. Pour simuler jusqu'à $t \approx 300$, il faut environ 1,5 million de pas.
C'est **le coût du schéma centré**, et c'est ce qui limitera le $Re$ atteignable en Partie C.

**Upwind et viscosité numérique : pourquoi ce choix est « malheureux » (§5.4.2.1).** L'**équation
modifiée** montre que le schéma upwind ne résout pas exactement l'équation voulue. Il résout une équation
avec une diffusion supplémentaire :

$$D = \tfrac{1}{2}\,U\,h\,(1 - C)$$

Application à $Re = 10^4$, $N = 128$ : $h \approx 0{,}0078$, donc $D \approx 4\cdot10^{-3}$, alors que
$\nu = 10^{-4}$. La viscosité numérique est **environ 40 fois plus grande que la vraie**. Le schéma
calcule en réalité un écoulement à

$$Re_\text{eff} \approx \frac{1}{\nu + D} \approx 250$$

À $Re \approx 250$, l'écoulement est stationnaire. **L'upwind va donc effacer les oscillations** qu'on
cherche en Partie C. Ce n'est pas un bug : c'est une conséquence directe de l'équation modifiée, et c'est
exactement le « choix malheureux » expliqué par le cours que l'énoncé demande de montrer.

---

## Résumé en 4 phrases

1. On cherche l'écoulement dans une boîte dont le couvercle glisse, avec un seul paramètre, $Re$. On le
   suit grâce à $K$ (énergie) et $P$ (puissance injectée).
2. La difficulté est la pression, qui n'a pas d'équation propre : **M1** la calcule par un Poisson pour
   corriger la vitesse, **M2** l'élimine en passant à $\psi$–$\omega$.
3. Les deux méthodes partagent les mêmes briques : Euler explicite, laplacien à 5 points, Poisson résolu
   par une LU factorisée une seule fois.
4. Le choix de l'advection est le cœur de l'histoire. Le schéma centré est précis mais impose
   $\Delta t \le 2/Re$. L'upwind est robuste mais trop diffusif : il fait « voir » un $Re$ bien plus faible
   et tue les oscillations.
