# ADR 0003 — Pas de temps, critère de stationnarité et tolérance du gradient conjugué (M1, Partie A)

- **Statut** : proposé (2026-11-11, avant l'étape 2), à confirmer par les mesures de l'étape 2.
- **Code concerné** : `cavity/common.py` (`time_step`, `solve_poisson`), `cavity/projection.py`.

## 1. Pourquoi le pas de temps est limité (cours §5.3)

### L'outil : l'analyse de von Neumann (§5.3.1.3)

On injecte une onde $e^{ikx}$ dans le schéma et on calcule le facteur $G$ par lequel elle est multipliée à
chaque pas. Si $|G| > 1$ pour une seule onde, elle grossit sans fin : le calcul explose.

### Diffusion : Euler explicite + centré (§5.3.1.3, éq. 5.27-5.32 ; §5.3.3, éq. 5.127)

- $G = 1 + 2\,Fo\,(\cos kh - 1)$, avec le nombre de Fourier $Fo = \nu\,\Delta t / h^2$.
- Pire cas : l'onde la plus courte que la grille porte (zigzag nœud à nœud, $\cos kh = -1$), $G = 1 - 4Fo$.
- $|G| \le 1$ ⇒ $Fo \le 1/2$ en 1D ; en 2D, les contributions de $x$ et $y$ s'ajoutent ⇒ $Fo \le 1/4$,
  soit $\Delta t \le h^2/(4\nu)$.

### Advection

- Upwind + Euler, 2D : $C_x + C_y \le 1$ (éq. 5.131), avec $|u|, |v| \le 1$ ⇒ $\Delta t \le h/2$.
- Centré + Euler : instable seul (exercice 9) ; stable grâce à la diffusion si $C^2 \le 2Fo \le 1$
  (§5.3.1.4, advection-diffusion) ⇒ $\Delta t \le 2\nu$. (Condition établie en 1D, appliquée telle quelle.)

### Laquelle fixe $\Delta t$ ? Ça dépend de la grille

La limite de diffusion est en $h^2$, celle d'advection en $h$ : en raffinant, la diffusion finit toujours par
l'emporter. Les deux se croisent quand $h^2/(4\nu) = h/2$, soit $h = 2\nu$, c'est-à-dire $N = 1/(2\nu) = Re/2$
($N = 50$ à $Re = 100$).

$Re = 100$ ($\nu = 0.01$), marge 0.8 :

| $N$ | Fourier $h^2/4\nu$ | CFL $h/2$ | centré $2\nu$ | $\Delta t$ retenu | limite active |
| --- | --- | --- | --- | --- | --- |
| 32 | $2.44\cdot10^{-2}$ | $\mathbf{1.56\cdot10^{-2}}$ | $2.00\cdot10^{-2}$ | $1.25\cdot10^{-2}$ | advection (CFL) |
| 64 | $\mathbf{6.10\cdot10^{-3}}$ | $7.81\cdot10^{-3}$ | $2.00\cdot10^{-2}$ | $4.88\cdot10^{-3}$ | diffusion |
| 128 | $\mathbf{1.53\cdot10^{-3}}$ | $3.91\cdot10^{-3}$ | $2.00\cdot10^{-2}$ | $1.22\cdot10^{-3}$ | diffusion |

**À haut $Re$, c'est l'inverse** : $\nu$ minuscule rend la limite de diffusion très lâche ; ce sont la CFL et
surtout $\Delta t \le 2\nu$ (centré) qui coûtent cher. Enjeu de la Partie C.

## 2. Critère de stationnarité

- Pseudo-code : arrêt quand $\max|u^{n+1} - u^n| < 10^{-6}$ et $t > 10$ (rampe du couvercle terminée).
- **Attention** : $\max|u^{n+1}-u^n| = \Delta t\,\max|\partial_t u|$, donc le critère revient à
  $\max|\partial_t u| < 10^{-6}/\Delta t$. Comme $\Delta t$ est ~10 fois plus petit à $N = 128$ qu'à $N = 32$, le
  critère y est ~10 fois **moins sévère**.
- **Décision** : on garde le critère du pseudo-code et on **affiche $\max|\partial_t u|$ à l'arrêt** pour chaque
  grille, afin de vérifier que l'écoulement est bien stationnaire partout.

## 3. Tolérance du gradient conjugué

### Résidu, résidu relatif, tolérance

- Le gradient conjugué est **itératif** : il produit $\phi_1, \phi_2, \dots$ qui se rapprochent de la solution de
  $A\phi = b$. Il faut décider quand s'arrêter.
- On ne connaît pas la solution, donc pas l'erreur ; on mesure de combien l'équation est **violée** : le
  **résidu** $r_k = b - A\phi_k$ (cours §4.3.2 ; workshop 11.1, question 3 : « exit condition controlled by the
  residual »).
- **Résidu relatif** : $\|r_k\| / \|b\|$. Pour la pression, $b = \nabla\cdot u^*/\Delta t$ change d'échelle avec $\Delta t$ et $N$ :
  un seuil relatif garde le même sens d'un cas à l'autre.
- **Tolérance** : on s'arrête dès que $\|r_k\|/\|b\| \le \text{tol}$. Plus petite ⇒ pression plus précise, mais plus
  d'itérations à chaque pas de temps.
- Dans M1 : une pression mal convergée corrige mal $u^*$ ⇒ la vitesse garde une divergence parasite.

### Pourquoi on ne juge pas la tolérance sur $\nabla\cdot u$

Sur grille colocalisée, la correction utilise un gradient **centré** ($p_{i+1} - p_{i-1}$) ; la divergence de ce
gradient est un laplacien **large** (points $i \pm 2$), pas le laplacien compact 5 points que résout le gradient
conjugué. Même avec un Poisson exact, $\nabla\cdot u^{n+1}$ ne s'annule donc pas exactement : il reste une
erreur en $O(h^2)$ sans rapport avec tol.

### Faut-il faire varier tol avec $\Delta t$ et avec $N$ ?

- **Avec $\Delta t$ : non.** $b = \nabla\cdot u^*/\Delta t$ grandit quand $\Delta t$ diminue ; un seuil absolu changerait donc de
  difficulté avec $\Delta t$. Le résidu **relatif** divise par $\|b\|$ : ce changement d'échelle s'annule, la même tol
  garde le même sens quel que soit $\Delta t$.
- **Avec $N$ : oui, en principe**, pour deux raisons qui vont dans le même sens :
  1. *L'erreur à battre diminue.* L'arrêt du gradient conjugué doit rester négligeable devant l'erreur du schéma,
     qui décroît avec $N$ (en $h^2$, plutôt en $h$ pour $p$ avec notre Neumann, cf. ADR 0002). Grille plus fine
     ⇒ Poisson plus précis nécessaire.
  2. *Un même résidu cache une erreur plus grande.* Le résidu mesure la violation de l'équation, pas l'erreur
     sur $p$. Les deux sont reliés par le conditionnement $K(A) = |\lambda|_\max / |\lambda|_\min$ (cours §4.2.3) :
     erreur relative $\le K(A) \times$ résidu relatif. Pour le laplacien, $K(A) = 4N^2/\pi^2$ (éq. 4.46), soit
     $\approx 4\cdot10^2$ à $N = 32$ et $\approx 6.6\cdot10^3$ à $N = 128$. À résidu égal, l'erreur garantie est 16 fois moins
     bonne à $N = 128$ qu'à $N = 32$.
- **Ordre de grandeur avec tol $= 10^{-8}$** : à $N = 128$, la borne vaut $10^{-8} \times 6.6\cdot10^3 \approx 7\cdot10^{-5}$ en relatif,
  du même ordre que $h^2 \approx 6\cdot10^{-5}$. La borne est pessimiste (l'erreur réelle est en général bien plus petite),
  mais rien ne garantit a priori que $10^{-8}$ suffise sur la grille fine : il faut le vérifier **là**.

### Décision

- tol $= 10^{-8}$ fixe dans la boucle en temps (départ à chaud depuis la pression du pas précédent).
- Vérification (workshop 15 : « impact of the convergence of the Poisson solver ») sur la grille **la plus
  exigeante, $N = 128$** : refaire le calcul avec tol $= 10^{-10}$ ; l'écart entre les deux profils doit être très petit
  devant l'écart à Ghia.
  - Si c'est le cas : une tol fixe suffit pour toutes nos grilles, inutile de la faire varier avec $N$.
  - Sinon : passer à tol $= 10^{-10}$ (surcoût modéré : le nombre d'itérations croît comme $\log(1/\text{tol})$, et le
    départ à chaud le limite encore).
