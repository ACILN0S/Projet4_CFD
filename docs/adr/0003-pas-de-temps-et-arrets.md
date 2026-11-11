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
- **Ce qui dépend de $\Delta t$ et ce qui n'en dépend pas** :
  - $\partial_t u$ mesure à quelle vitesse l'écoulement évolue : c'est une propriété de l'écoulement, **pas** du
    pas de temps. On l'estime par $(u^{n+1}-u^n)/\Delta t$ ; la division par $\Delta t$ compense la dépendance.
  - $u^{n+1}-u^n$, lui, dépend de $\Delta t$ (deux fois plus petit si $\Delta t$ l'est).
  - Le critère porte sur $u^{n+1}-u^n$, donc c'est le **seuil** sur $\partial_t u$ qui dépend de $\Delta t$.
- **On connaît d'avance la valeur de $\partial_t u$ à l'arrêt** : $\max|\partial_t u| < 10^{-6}/\Delta t$ par construction, et
  $\Delta t$ est fixé avant le calcul. L'afficher n'apprend donc rien.

  | $N$ ($Re = 100$) | $\Delta t$ | $10^{-6}/\Delta t$ (borne de $\max\lvert\partial_t u\rvert$ à l'arrêt) |
  | --- | --- | --- |
  | 32 | $1.25\cdot10^{-2}$ | $8.0\cdot10^{-5}$ |
  | 64 | $4.88\cdot10^{-3}$ | $2.0\cdot10^{-4}$ |
  | 128 | $1.22\cdot10^{-3}$ | $8.2\cdot10^{-4}$ |

- **Un $\partial_t u$ petit ne prouve pas qu'on est proche du stationnaire.** Près du stationnaire, l'écart décroît
  en $e^{-\lambda t}$ ($\lambda$ : taux d'amortissement du mode le plus lent), donc
  $$|u - u_\infty| \approx \frac{|\partial_t u|}{\lambda}.$$
  Si $\lambda$ est petit, un $\partial_t u$ petit peut cacher un écart à $u_\infty$ important (ex. : $\lambda \sim 0.1$ et
  $\partial_t u = 8\cdot10^{-4}$ ⇒ $|u - u_\infty| \sim 10^{-2}$, de l'ordre de l'écart à Ghia). $\lambda$ n'est pas connu d'avance :
  la borne du tableau ne suffit pas, il faut **vérifier après coup** (voir Décision).
- **$\lambda$ diminue quand $Re$ augmente** : la viscosité amortit de moins en moins vite les perturbations, et
  $\lambda \to 0$ à l'approche de $Re_c \approx 8000$ (bifurcation). Au-delà, il n'y a plus de solution stationnaire :
  l'écoulement est périodique et le critère n'est **jamais** atteint.
- **Pourquoi pas un critère directement sur $\partial_t u$, soit $\max|u^{n+1}-u^n|/\Delta t < 10^{-6}$ ?**
  Il aurait l'avantage d'avoir le même sens sur toutes les grilles, mais :
  1. *Il est beaucoup plus sévère.* Il revient à $\max|u^{n+1}-u^n| < 10^{-6}\,\Delta t$, soit $\approx 1.2\cdot10^{-9}$ à
     $N = 128$ : 1000 fois plus strict que le pseudo-code.
  2. *Il coûte cher.* L'approche du stationnaire est exponentielle ($u - u_\infty \propto e^{-\lambda t}$) : chaque
     facteur 10 de sévérité en plus ajoute une durée fixe $\ln 10/\lambda$ de temps simulé, et à $N = 128$ chaque
     unité de temps demande déjà ~800 pas.
  3. *Il peut ne jamais être atteint.* Le gradient conjugué s'arrête à tol $= 10^{-8}$ (§3) : la pression garde une
     petite erreur à chaque pas, qui fait « bruiter » $u$. Si ce bruit dépasse $10^{-9}$, l'écart $u^{n+1}-u^n$
     plafonne au-dessus du seuil et la boucle ne s'arrête pas.
  - Un critère sur $\partial_t u$ demanderait donc un seuil plus grand (par ex. $10^{-4}$), ce qui revient à choisir
    nous-mêmes une valeur absente du pseudo-code. Et il ne réglerait pas le vrai problème (le $\lambda$ inconnu).

### Décision : une seule règle d'arrêt pour les Parties A, B et C

- **Règle** (celle du pseudo-code, identique pour M1 et M2) : on s'arrête dès que
  - $\max|u^{n+1}-u^n| < 10^{-6}$ et $t > 10$ ⇒ écoulement **stationnaire** ;
  - **ou** $t \ge t_\text{final}$ ⇒ pas de stationnaire atteint.
- **Ce que donne la même règle selon la partie** :
  - A ($Re = 100$) et B ($Re \le 5000$) : stationnaires, l'arrêt se fait par le critère. $t_\text{final}$ ne sert que
    de garde-fou ; l'atteindre est un signal d'alerte à signaler.
  - C ($Re > Re_c$) : schéma centré périodique ⇒ arrêt à $t_\text{final}$ (choisi pour couvrir plusieurs périodes de
    $K(t)$ une fois le régime établi) ; l'upwind, lui, doit rester stationnaire et s'arrêter par le critère.
    On sait **lequel des deux arrêts** a joué : c'est déjà un résultat (stationnaire ou non).
- **Validation du seuil $10^{-6}$, une seule fois, dans le cas le plus défavorable** : critère le plus lâche
  ($N = 128$, plus petit $\Delta t$) et amortissement le plus lent (plus haut $Re$ stationnaire, $Re = 5000$).
  On prolonge ce calcul jusqu'à $\max|u^{n+1}-u^n| < 10^{-7}$ et on compare les profils de $u$ au centre :
  - écart très petit devant l'écart à Ghia ⇒ $10^{-6}$ suffit pour **tous** les cas stationnaires de A et B
    (ils sont moins défavorables) ;
  - sinon ⇒ passer à $10^{-7}$ partout.
  Au début de la Partie A, on peut faire le même test à $Re = 100$, $N = 128$ (peu coûteux) pour avancer sans attendre B.

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
