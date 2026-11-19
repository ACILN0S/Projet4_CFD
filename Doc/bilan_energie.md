# Bilan d'énergie de la cavité : $dK/dt = P - \varepsilon$ (Partie B, 2026-11-15)

> Contrôle demandé en Partie B : au stationnaire, la puissance injectée $P$ doit être dissipée par la viscosité,
> $P = \varepsilon$. Ce qui est **démontré** est marqué *(démontré)* ; ce qui repose sur un argument d'échelle ou sur la seule
> observation est marqué comme tel. Calculs faits sur les champs finaux de la Partie B (`results/B_*_tref30.npz`).
> Notations : $\nu = 1/Re$, $h = 1/N$, sommes sur les nœuds intérieurs $i, j = 1..N-1$ sauf mention.

---

## 1. Les observations (résumé)

1. **Le bilan discret n'est pas fermé** : au stationnaire ($dK/dt \approx 0$ à mieux que 1 %), l'écart $(P - \varepsilon)/P$ vaut
   de $-36$ % à $+15$ % selon la méthode, $Re$ et $N$ ; **M1** a $\varepsilon > P$ (source d'énergie numérique), **M2** a $\varepsilon < P$
   (dissipation numérique en plus de la viscosité).
2. **L'écart dépend du maillage et tend vers 0 quand on raffine** : c'est une erreur de discrétisation, pas un défaut de
   convergence en temps. M1, $Re = 5000$ : $-36$ %, $-17$ %, $-8.6$ % à $N = 64, 128, 256$ (ordre $\approx 1$) ; M2 : $+15$ %, $+14$ %, $+5.4$ %.
3. **Retirer les coins du couvercle ne ferme pas le bilan.** Ni en les retirant de $\varepsilon$ seul (on enlève de la dissipation
   sans enlever l'injection correspondante), ni en retirant **la même zone** de $P$ et de $\varepsilon$ : près des coins, $P$ et $\varepsilon$
   portent la même part (5–6 % chacun à $Re = 5000$), qui se compense. **Il n'y a pas de zone coupable dans $P$ ni dans $\varepsilon$.**
4. **Mais pour M1, le défaut du schéma lui-même est localisé** (§4) : il vient surtout d'un terme de **travail de la pression aux
   parois**, nul en continu, dont 97 à 100 % est porté par le couvercle à moins de $0.05$ des deux coins. Ce terme n'appartient ni à
   $P$ ni à $\varepsilon$ : c'est pourquoi retirer les coins de $P$ et de $\varepsilon$ ne le supprime pas.

## 2. Le bilan continu

### 2.1 $dK/dt = P - \varepsilon$ *(démontré, écoulement régulier)*

On multiplie Navier–Stokes, $\partial_t \mathbf u + (\mathbf u\cdot\nabla)\mathbf u = -\nabla p + \nu\Delta\mathbf u$, par $\mathbf u$ et on intègre sur $\Omega$ ($\nabla\cdot\mathbf u = 0$) :

- $\int_\Omega \mathbf u\cdot\partial_t\mathbf u = dK/dt$ ;
- advection : $\mathbf u\cdot(\mathbf u\cdot\nabla)\mathbf u = \nabla\cdot(\tfrac12|\mathbf u|^2\mathbf u)$, donc $\int = \oint \tfrac12|\mathbf u|^2\,\mathbf u\cdot\mathbf n = 0$ ($\mathbf u\cdot\mathbf n = 0$ sur toutes les parois) ;
- pression : $\mathbf u\cdot\nabla p = \nabla\cdot(p\mathbf u)$, donc $\int = \oint p\,\mathbf u\cdot\mathbf n = 0$ ;
- viscosité : $\int \mathbf u\cdot\Delta\mathbf u = -\int|\nabla\mathbf u|^2 + \oint \mathbf u\cdot\partial_n\mathbf u$. Seul le couvercle a $\mathbf u \ne 0$ : $\mathbf u = (U_t, 0)$, $\partial_n = \partial_y$,
  d'où $\nu\oint = \int_0^1 \nu\,U_t\,\partial_y u\,dx = P$.

Donc $\boxed{dK/dt = P - \varepsilon}$, avec $\varepsilon = \nu\int_\Omega|\nabla\mathbf u|^2\,dS$.

### 2.2 $\int|\nabla\mathbf u|^2 = \int\omega^2$ dans la cavité *(démontré, écoulement régulier)*

En 2D incompressible, $|\nabla\mathbf u|^2 - \omega^2 = u_x^2 + v_y^2 + 2u_yv_x = -2(u_xv_y - u_yv_x)$ (on a utilisé $u_x = -v_y$).
Le jacobien $J = u_xv_y - u_yv_x$ est une divergence : $J = \partial_x(u\,v_y) - \partial_y(u\,v_x)$. Donc $\int_\Omega J = \oint u\,(v_y n_x - v_x n_y)\,ds = \pm\oint u\,\partial_s v\,ds$
($\partial_s$ : dérivée le long de la paroi). Sur les parois latérales et le fond, $u = 0$ ; sur le couvercle, $v = 0$ le long de la paroi, donc $\partial_s v = 0$.
Le terme de bord est nul : $\varepsilon = \nu\int|\nabla\mathbf u|^2 = \nu\int\omega^2$. (En discret, les deux diffèrent de 4 à 9 % : erreur de discrétisation.)

### 2.3 Les coins du couvercle : $P$ et $\varepsilon$ sont infinis, leur différence ne l'est pas *(argument d'échelle)*

Près d'un coin (distance $r$), l'écoulement est dominé par la viscosité (Stokes) et auto-semblable : $\mathbf u \approx U\,\mathbf g(\theta)$
(problème du « grattoir » de Taylor ; référence à vérifier, cf. PLAN §2). Alors $|\nabla\mathbf u| \sim U/r$ et $p \sim \nu U/r$.
- $\varepsilon$ sur le quart de disque $\delta < r < R$ : $\nu\int (U/r)^2\, r\,dr\,d\theta \propto \nu U^2\ln(R/\delta)$ : **divergence logarithmique** ;
- $P$ : $\partial_y u \sim U/x$ sur le couvercle, $\int_\delta \nu U\,dx/x \propto \nu U^2\ln(1/\delta)$ : **même divergence**. Sur une grille, $\delta \sim h$ : $P$ et $\varepsilon$
  croissent en $\ln N$ (vérifié en Partie A pour $P$).
- **Bilan sur le domaine privé des deux quarts de disque** $\Omega_\delta$ : $dK_\delta/dt = P_\delta - \varepsilon_\delta + F_\delta$, où $F_\delta$ est le flux d'énergie à travers
  les deux arcs (advection, travail de la pression, diffusion). Avec la forme auto-semblable, le terme visqueux $\nu\,\mathbf u\cdot\partial_r\mathbf u$ est nul
  ($\mathbf u$ ne dépend pas de $r$), l'advection est en $O(U^3\delta)$, mais le travail de la pression vaut $\int p\,u_r\,\delta\,d\theta = O(\nu U^2)$ :
  **fini et non nul quand $\delta \to 0$**. Conséquence : retirer la même zone de $P$ et de $\varepsilon$ laisse un écart $O(\nu U^2)$, soit $O(\nu/P) = O(Re^{-1/2})$
  en relatif, donc **plus visible aux bas $Re$**. C'est ce qu'on observe : à $Re = 100$, M2 (bilan fermé à 0.1 %) passe à $-2$ à $-4$ % quand on retire
  2 à 5 nœuds par coin.

## 3. Le bilan discret : ce qui est démontré

### 3.1 Le temps ne joue aucun rôle au stationnaire *(démontré)*

Au stationnaire, $u^{n+1} = u^n$ : les équations discrètes stationnaires sont satisfaites exactement, quel que soit l'intégrateur en
temps. L'écart du bilan au stationnaire est donc une **erreur d'espace** (formules aux différences, conditions aux limites, intégrales).
(Pour M1, l'état stationnaire dépend quand même de $\Delta t$ à travers la projection : ADR 0002 §7.3.)

### 3.2 M1 : une identité d'énergie discrète exacte *(démontré, vérifié à $10^{-16}$)*

À chaque pas, le résidu stationnaire de M1 vaut exactement $R = -A_h(\mathbf u) + \nu L_h\mathbf u - G_h p$ aux nœuds intérieurs
($A_h$ : advection centrée, $L_h$ : laplacien 5 points, $G_h$ : gradient centré). On le multiplie par $h^2\mathbf u$ et on somme.
**Sommation par parties** (analogue discret de l'intégration par parties), en 1D avec $u_0, u_N$ donnés :
$$\sum_{i=1}^{N-1} u_i(u_{i+1} - 2u_i + u_{i-1}) = -\sum_{i=0}^{N-1}(u_{i+1}-u_i)^2 + u_N(u_N - u_{N-1}) - u_0(u_1 - u_0).$$
En 2D, les termes de bord ne subsistent qu'au couvercle (seule paroi où $\mathbf u \ne 0$). On obtient, **exactement** :
$$P_h - \varepsilon_h = \mathcal A + \mathcal G + h^2\sum\mathbf u\cdot R,$$
avec $P_h = \nu\sum_{i=1}^{N-1} U\,(U - u_{i,N-1})$ (puissance au couvercle, $\partial_y u$ par différence d'ordre 1), $\varepsilon_h = \nu\sum_{\text{arêtes}}(\Delta u)^2 + (\Delta v)^2$
(dissipation sur les arêtes de la grille), $\mathcal A = h^2\sum\mathbf u\cdot A_h(\mathbf u)$, $\mathcal G = h^2\sum\mathbf u\cdot G_h p$. En continu, $\mathcal A$ et $\mathcal G$ sont nuls (§2.1) ;
en discret, ils ne le sont pas : **ce sont eux, et eux seuls, qui empêchent le bilan de se fermer.**

Puis, par une seconde sommation par parties, $\mathcal G = -h^2\sum p\,\mathrm{div}_h\mathbf u + \mathcal B$, avec le **travail discret de la pression aux parois**
$$\mathcal B = \tfrac h2\sum_j\big(u_{N-1,j}\,p_{N,j} - u_{1,j}\,p_{0,j}\big) + \tfrac h2\sum_i\big(v_{i,N-1}\,p_{i,N} - v_{i,1}\,p_{i,0}\big),$$
la pression **sur la paroi** multipliée par la vitesse normale **au premier nœud intérieur**. En continu, ce terme est $\oint p\,\mathbf u\cdot\mathbf n = 0$.

Vérification sur les champs finaux (en % de $P$ ; les deux identités sont satisfaites à $10^{-16}$ près, $h^2\sum\mathbf u\cdot R$ est au plus $10^{-4}$ %) :

| $Re$, $N$ | $(P-\varepsilon)/P$ mesuré | $(P_h-\varepsilon_h)/P$ | $\mathcal A/P$ | $-h^2\sum p\,\mathrm{div}_h\mathbf u\,/P$ | $\mathcal B/P$ |
| --- | --- | --- | --- | --- | --- |
| 1000, 64 | $-14.3$ | $-19.5$ | $-2.5$ | $-0.9$ | $-16.0$ |
| 1000, 128 | $-9.1$ | $-13.3$ | $-1.0$ | $-1.8$ | $-10.5$ |
| 5000, 64 | $-35.8$ | $-45.0$ | $-16.9$ | $-0.04$ | $-28.1$ |
| 5000, 128 | $-17.0$ | $-23.2$ | $-6.1$ | $-0.08$ | $-17.0$ |
| 5000, 256 | $-8.6$ | $-12.9$ | $-1.7$ | $-0.4$ | $-10.8$ |

- **Le terme dominant est $\mathcal B$** ; la divergence discrète non nulle de M1 (grille colocalisée) ne pèse presque rien dans le bilan.
- $\mathcal A$ (l'advection centrée sous forme non conservative ne conserve pas exactement l'énergie) décroît vite : $-17 \to -6 \to -1.7$ % (ordre $\approx 1.5$–$1.8$).
- $\mathcal B$ décroît lentement : $-28 \to -17 \to -10.8$ % (ordre $\approx 0.7$).
- La différence entre le bilan « mesuré » ($P$ d'ordre 2 au couvercle, $\varepsilon$ par $\omega$ et les trapèzes) et le bilan discret exact vient des formules de mesure
  ($P - P_h$ et $\varepsilon - \varepsilon_h$ valent chacun 8 à 28 % de $P$ et se compensent en partie).

### 3.3 Pourquoi $\mathcal B$ est porté par les coins du couvercle *(démontré hors des coins ; observé aux coins)*

Sur une paroi où la vitesse tangentielle est constante le long de la paroi (toutes nos parois, couvercle compris : $u = U_t$), la vitesse
normale vérifie $u_n = 0$ et $\partial_n u_n = -\partial_s u_s = 0$ (incompressibilité). Par Taylor (cours §2.1), au premier nœud intérieur,
$u_n(h) = \tfrac{h^2}{2}\partial_n^2 u_n + O(h^3)$ : là où l'écoulement est régulier, la contribution de $\mathcal B$ est $O(h^2)$. Aux deux coins du couvercle,
la vitesse tangentielle saute de $U_t$ à 0 et $p \sim 1/r$ : l'argument tombe, et c'est là que $\mathcal B$ se concentre. Observé : 97 à 100 % de $\mathcal B$ vient
du couvercle à moins de $0.05$ d'un coin ($Re = 1000$ et $5000$, $N = 64$ à $256$). D'où l'ordre faible de $\mathcal B$ ($\approx 0.7$), imposé par la singularité.

### 3.4 Pourquoi on ne peut pas « localiser » l'erreur avec les densités $\mathbf u\cdot\nabla p$ ou $\mathbf u\cdot(\mathbf u\cdot\nabla)\mathbf u$

Même en continu, ces densités ne sont pas nulles localement (ce sont des divergences : seule leur **intégrale** s'annule). Leurs parts par
zone se compensent fortement (pour $\mathcal A$ : $-644$ % le long des parois, $+747$ % à l'intérieur, à $Re = 1000$). Seules des quantités nulles en
continu **en chaque point** ($\mathrm{div}_h\mathbf u$) ou **en chaque point de paroi** ($u_n$ au premier nœud, dans $\mathcal B$) mesurent une erreur locale.

### 3.5 M2 : pas d'identité équivalente *(démontré)*

M2 ne discrétise pas l'équation de quantité de mouvement mais celle de $\omega$ : multiplier son schéma par $\mathbf u$ ne donne pas de bilan de $K$, et
l'écart de M2 ($+5$ à $+15$ % à haut $Re$) ne peut pas être décomposé de la même façon. Fait démontré en passant : la vitesse de M2,
$u = \delta_y\psi$, $v = -\delta_x\psi$ (différences centrées), a une divergence centrée **exactement nulle à tous les nœuds intérieurs** : les différences
centrées commutent, et aux nœuds voisins des parois la divergence n'utilise que $u$ sur les parois latérales et $v$ sur le couvercle et le fond, dont
les valeurs imposées (0) sont exactement celles que donnerait $\psi = 0$ ; la vitesse du couvercle $U$ n'y intervient pas (vérifié : $\max|\mathrm{div}_h\mathbf u| = 7\cdot10^{-15}$).
Le défaut de M2 ne vient donc pas de la divergence, mais des parois (formule de Thom, d'ordre 1 localement) et de l'advection, sans qu'on puisse le
chiffrer terme à terme.

## 4. Conclusion

- Le bilan $dK/dt = P - \varepsilon$ est exact en continu ; en discret, il ne se ferme qu'**à l'erreur de maillage près**, et cette erreur décroît avec $N$
  (ordre $\approx 1$ pour M1). Le critère d'arrêt n'y est pour rien ($dK/dt \approx 0$ à mieux que 1 %).
- **Dans $P$ et $\varepsilon$, aucune zone n'est coupable** : leurs parts près des coins se compensent, et les retirer n'améliore pas le bilan (il reste un
  terme de flux $O(\nu U^2)$ à travers la frontière de la zone retirée, §2.3).
- **Pour M1, l'origine du défaut est identifiée et démontrée** : le travail discret de la pression aux parois ($\mathcal B$), qui en continu est nul, et qui est
  concentré aux coins du couvercle ; l'advection centrée y ajoute un défaut qui disparaît vite quand on raffine.
- Pour M2, l'écart décroît nettement seulement au-delà de $N = 128$, comme ses profils (notebook B §4).
