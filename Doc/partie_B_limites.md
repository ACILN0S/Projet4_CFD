# Partie B — difficultés rencontrées, diagnostics et limites

Ce document accompagne le notebook B (`notebooks/B_sweep_Re.ipynb`), qui ne garde que les résultats. Il rassemble ce qu'on a dû
comprendre ou corriger en chemin, avec les chiffres qui le justifient. Le bilan d'énergie (dK/dt = P − ε) est traité à part :
[bilan_energie.md](bilan_energie.md).

Réglages communs : grille colocalisée, Euler explicite, advection centrée, gradient conjugué (tol 1e-8), Δt = limites de stabilité
du cours. Les calculs ont été faits sur Ruche ; les scripts de diagnostic (calibration, Woods, Δt) et la version optimisée du code
qui les a fait tourner sont archivés dans `Projet4_CFD_version_optimisee_2026-11-19.zip` (dossier parent). Ils redonnent les mêmes
résultats que le code livré à 2·10⁻⁹ près (vérifié sur M1 et M2, Re = 400, N = 64).

---

## 1. Critère d'arrêt : à calibrer à haut Re

**Problème.** L'ancien critère (‖R‖ < 10⁻² ‖R(10)‖, choisi en Partie A à Re = 100) arrêtait beaucoup trop tôt à haut Re : à
Re = 5000, M1 s'arrêtait à t = 98 au lieu de t ≈ 430, avec une erreur d'arrêt de 245 % de l'écart à Ghia. Le grand tourbillon
finit de s'accélérer très lentement (il faut ~190 unités de temps pour diviser le résidu par 10).

**Méthode.** Un calcul long par cas à Re = 5000 (le plus défavorable), jusqu'à t = 1000, en gardant l'historique de ‖R‖ et les
profils médians. L'état final sert de référence ; pour chaque couple (t_ref, ε), on retrouve l'instant où le critère aurait arrêté
le calcul et l'écart des profils à cet instant. Règle : erreur d'arrêt ≤ 10 % de l'écart à Ghia.

| Cas | t_ref = 10, ε = 1e-2 (ancien) | t_ref = 10, ε = 1e-3 | **t_ref = 30, ε = 1e-3** | t_ref = 10, ε = 1e-4 |
|---|---|---|---|---|
| M1, N = 128 | t = 98 : 245 % | t = 247 : 40 % | **t = 427 : 4.6 %** | t = 438 : 4.0 % |
| M1, N = 64  | t = 151 : 65 % | t = 314 : 9.3 % | **t = 365 : 5.1 %** | t = 508 : 0.9 % |
| M2, N = 128 (R sur ω) | t = 156 : 15 % | t = 229 : 5.8 % | **t = 266 : 3.7 %** | t = 341 : 1.5 % |

**Décision** : t_ref = 30, ε = 10⁻³ pour tous les Re de la Partie B (ADR 0003 §2). Avec ce couple, K varie de moins de 1.5 % et
P de moins de 10⁻⁴ en relatif après l'arrêt.

**Résidu de M2 mesuré sur ω.** Mesuré d'abord sur u, v, le résidu de M2 valait **exactement 0** à certains pas : le gradient
conjugué, parti du ψ du pas précédent, était déjà sous sa tolérance et faisait 0 itération, donc ψ, u, v ne bougeaient pas alors
que ω avait évolué. 13 calculs M2 sur 23 s'étaient arrêtés à tort sur un tel pas. On mesure donc le résidu de l'équation que M2
fait évoluer, $R = (\omega^{n+1} - \omega^n)/\Delta t$, le pendant exact du $R = (u^{n+1} - u^n)/\Delta t$ de M1. Vérifié avec le même couple (tableau,
dernière ligne). Tous les calculs M2 ont été refaits ; ils redonnent les anciennes valeurs à 3·10⁻⁴ près.

**Limite qui reste.** À Re = 100, l'écoulement est stationnaire avant t_ref = 30 ; ‖R(30)‖ est déjà au niveau du bruit du gradient
conjugué et ne peut plus être divisé par 1000. M1 à N = 128 va donc jusqu'au garde-fou t = 1000 (résultat juste, plus convergé que
nécessaire). C'est la limite d'un critère relatif ; un seuil absolu en plus serait un réglage de plus à justifier, pour un gain nul.

---

## 2. M2 à N = 64 : une fausse solution stationnaire pour Re ≥ 3200

**Constat.**

| Re | M2, N = 64 : K / ψ_min | M2, N = 128 : K / ψ_min | M1, N = 64 : K |
|---|---|---|---|
| 1000 | 0.0363 / −0.1052 | 0.0423 / −0.1155 | 0.0439 |
| 2500 | 0.0266 / −0.0879 | 0.0403 / −0.1121 | 0.0455 |
| 3200 | **0.0080 / −0.0306** | 0.0388 / −0.1097 | 0.0471 |
| 5000 | **0.0057 / −0.0212** | 0.0345 / −0.1026 | 0.0535 |

M2 à N = 64 se dégrade dès Re ≈ 700, puis s'effondre entre Re = 2500 et 3200 : le tourbillon principal perd 75 % de son
intensité (Ghia : ψ_min ≈ −0.12). Le calcul est pourtant bien stationnaire : poussé jusqu'à une réduction du résidu de 10⁷
(calcul long du §1), il reste sur cette solution.

**Pourquoi « fausse ».** On ne peut pas démontrer que le problème continu n'a pas cette solution. Mais elle n'est pas convergée en
maillage (N = 128 et 256 donnent un tourbillon 4 à 5 fois plus intense), M1 sur la même grille ne la donne pas, et elle contredit
Ghia (1982) et Erturk (2005).

**Ce que ce n'est pas** (test M2, N = 64, schéma centré contre upwind) :

| Re | Péclet de maille Uh/ν | centré : ψ_min | upwind : ψ_min |
|---|---|---|---|
| 1500 | 23 | −0.0996 | −0.0755 |
| 3200 | 50 | −0.0306 | −0.0458 |
| 5000 | 78 | −0.0212 | −0.0347 |

- **Pas de la diffusion numérique** : l'upwind ajoute une viscosité numérique de plusieurs fois ν. Une solution trop diffusée
  ressemblerait à celle d'un Re plus faible, donc à un tourbillon intense ; on observe l'inverse, avec les deux schémas.
- **Pas le schéma d'advection** : deux schémas très différents échouent de la même façon.
- **Pas seulement les oscillations du schéma centré** : M1 a des oscillations nœud à nœud comparables (amplitude 0.24 contre 0.30
  à Re = 3200) et reste correct.

Ce qui distingue M1 de M2 à grille et schéma identiques, c'est la condition de paroi : §3.

---

## 3. Pourquoi M2 se dégrade plus vite que M1 : la formule de Thom

### 3.1 Erreur de la formule de Thom (démonstration)

M1 impose u = U exactement au couvercle. M2 n'a pas de condition physique sur ω ; il la reconstruit à chaque pas par la formule de
Thom. Notations : n = distance à la paroi vers l'intérieur, ω = ∂v/∂x − ∂u/∂y, donc Δψ = −ω. Le long de la paroi, ψ est constant
et la vitesse de paroi uniforme, donc ∂ψ/∂n = ∓u_s, ∂²ψ/∂n² = −ω_w, ∂³ψ/∂n³ = −∂ω/∂n|_w. Taylor :

$$\psi_1 = \psi_w + h\,\partial_n\psi + \frac{h^2}{2}\,\partial_n^2\psi + \frac{h^3}{6}\,\partial_n^3\psi + O(h^4).$$

Thom ne garde que les deux premiers termes, ω_Thom = −2(ψ₁ − ψ_w − h ∂ψ/∂n)/h², d'où

$$\omega_{\text{Thom}} = \omega_w + \frac{h}{3}\,\partial_n\omega\big|_w + O(h^2).$$

- L'erreur locale est **d'ordre 1**. Dans la couche limite, |ω| décroît en s'éloignant de la paroi : Thom **sous-estime toujours
  |ω_w|**, d'une fraction h/(3δ_ω), où δ_ω est l'épaisseur de la couche limite.
- δ_ω ∝ Re^(−1/2), donc **à h fixé, l'erreur croît comme √Re**.
- Thom donne une solution globalement d'ordre 2 (Erturk et al. 2005, §3 ; Huang et Wetton 1996), mais c'est un résultat
  asymptotique qui suppose h ≪ δ_ω. Ici h/δ_ω ≈ 0.3 à Re = 3200, N = 128.

**Vérification** : formule de Thom appliquée au profil de M1 sous le couvercle (qui n'utilise pas Thom), moyenne sur x ∈ [0.1, 0.9].

| Re | profil M1 | pas H | erreur de Thom sur ω_w |
|---|---|---|---|
| 100 | N = 128 | 1/64 | −2.3 % |
| 1000 | N = 128 | 1/64 | −7.3 % |
| 3200 | N = 128 | 1/64 | −16.6 % |
| 5000 | N = 128 | 1/64 | −22.1 % |
| 3200 | N = 256 | 1/128 | −5.2 % |
| 5000 | N = 256 | 1/128 | −7.5 % |

Signe prévu dans tous les cas ; erreur/√Re presque constante (0.21 à 0.31) ; diviser H par 2 divise l'erreur par 2.4 à 2.6.

### 3.2 Pourquoi la vitesse sous le couvercle est trop faible (démonstration)

Au couvercle, ψ(1) − ψ(1−h) = h ū₁, où ū₁ est la vitesse moyenne dans la première maille. Thom s'écrit donc exactement
$$\bar u_1 = U - \frac{h}{2}\,|\omega_w| ,$$
c'est la vitesse moyenne d'un profil **linéaire**, la tangente à la paroi. Le vrai profil vérifie
$$\bar u_1 = U - \frac{h}{2}\,|\omega_w| + R, \qquad R \approx \frac{h^2}{6}\,\partial_n^2 u > 0,$$
 parce que le profil est **convexe** sous le couvercle. Le signe vient de la
physique : à la paroi, u = U et v = 0, et l'équation de quantité de mouvement stationnaire se réduit à ν ∂²u/∂y² = ∂p/∂x ; le
couvercle pousse le fluide vers le coin haut droit, où la pression monte, donc ∂p/∂x > 0.

**Conclusion** : pour une vorticité de paroi donnée, Thom impose dans la première maille une vitesse trop faible, d'une fraction
≈ h/(3δ_ω), qui croît comme √Re.

**Mesures** (u(0.5, y = 0.953), juste sous la couche limite ; profils complets :
[figures/B_lid_boundary_layer.png](../figures/B_lid_boundary_layer.png)) :

| Cas | Ghia | M1 | M2 |
|---|---|---|---|
| Re = 3200, N = 128 | 0.461 | 0.453 | **0.409** (−11 %) |
| Re = 3200, N = 256 | 0.461 | 0.466 | 0.448 |
| Re = 5000, N = 128 | 0.460 | 0.466 | **0.396** (−14 %) |
| Re = 5000, N = 256 | 0.460 | 0.479 | 0.456 |
| Re = 3200, N = 64 | 0.461 | 0.450 | **−0.133** (fausse solution) |

- Le déficit de M2 au premier nœud a le signe et la tendance prévus ; il vaut 1.5 à 2 fois la prédiction (h²/6) ∂²u/∂n².
- Le saut de vitesse à travers la couche limite est donc plus grand chez M2, d'où une vorticité de paroi plus forte que celle de M1
  (rapport 1.08 à Re = 3200, N = 128), et moins de fluide entraîné : le tourbillon principal est trop faible.
- Raffiner corrige le défaut : à Re = 3200, le déficit de M2 (par rapport à M1, N = 256) passe de −0.057 à N = 128 à −0.018 à N = 256.
- La fausse solution de N = 64 est le cas extrême : le fluide sous le couvercle est presque à l'arrêt.

**Ce qui n'est pas démontré** : le passage du déficit local (première maille) au déficit, plus grand, de tout le jet sous le
couvercle. C'est la réponse globale de l'écoulement discret.

**Test de confirmation : M2 avec la formule de Woods** (ω_w = 3/2 ω_Thom − ½ ω₁, erreur locale d'ordre 2 ; hors méthode retenue,
`scripts/test_woods.py` dans l'archive). Appliquée aux profils de M1, elle se trompe de −5 % à +1 % sur ω_w, contre −2 à −22 % pour
Thom. Calculs complets (Ruche, 2026-11-19) :

| Cas | u(0.5, 0.953) Thom → Woods (Ghia ≈ 0.46) | écart max à Ghia Thom → Woods | ψ_min Thom → Woods |
|---|---|---|---|
| Re = 1000, N = 64 | 0.412 → 0.424 | 0.062 → 0.049 | −0.105 → −0.108 |
| Re = 1000, N = 128 | 0.457 → 0.459 | 0.014 → 0.011 | −0.116 → −0.116 |
| Re = 3200, N = 64 | **−0.133 → 0.355** | **0.59 → 0.12** | **−0.031 → −0.099** |
| Re = 5000, N = 64 | fausse solution → **diverge** (t = 24) | | |

- **Woods supprime la fausse solution** de N = 64 à Re = 3200 et réduit l'erreur dans tous les cas : la formule de Thom est bien
  une cause du défaut de M2, et la cause de la bascule vers la fausse solution.
- **Mais Woods ne ramène pas M2 au niveau de M1** (Re = 3200, N = 64 : u = 0.355 contre 0.450 pour M1). M2 a donc d'autres erreurs
  près de la paroi quand la couche limite tient en une ou deux mailles (vitesses tirées de ψ par différences centrées, transport de ω
  dans une couche limite très raide) ; on ne les a pas isolées.
- **Woods est moins stable** : il diverge à Re = 5000, N = 64, comme le signale la littérature pour les schémas explicites. C'est une
  raison de garder Thom dans la méthode retenue.
- Woods à N = 128, Re = 3200 et 5000 : en cours (à compléter).

**Pourquoi M2 n'a pas « raison » et M1 tort.** Le signe de l'erreur de Thom ne repose que sur la convexité du vrai profil ; à
Re = 5000, entre N = 128 et 256, la vitesse sous le couvercle de M2 bouge 5 fois plus que celle de M1 et rejoint M1 ; par rapport à
Erturk, M1 est 3 à 5 fois plus précis que M2 à grille égale (§4).

---

## 4. Ghia a sa propre erreur ; comparaison à Erturk

**Ghia n'est pas exact à haut Re.** Ghia et al. (1982) n'ont pas d'étude de convergence chiffrée (129² jusqu'à Re = 3200, 257² au-delà).
Erturk et al. (2005) ont calculé sur 401², 513², 601² avec extrapolation de Richardson. Tourbillon principal (`ref_data/primary_vortex.csv`) :

| Re | ψ_min Ghia | ψ_min Erturk | écart |
|---|---|---|---|
| 1000 | −0.11793 | −0.11894 | 0.9 % (Botella et Peyret, spectral : −0.118937) |
| 5000 | −0.11897 | −0.12223 | **2.7 %** |

**Conséquence : M1 semblait ne pas converger.** Par rapport à Ghia, M1 s'éloigne en raffinant (Re = 5000 : écart max 0.012 à
N = 128, 0.019 à N = 256). Par rapport aux profils d'Erturk (Tables 6 et 7, extraits dans `ref_data/erturk2005_*.csv`), il converge
régulièrement :

| Re = 5000, max\|Δu\| / max\|Δv\| à Erturk | N = 64 | N = 128 | N = 256 |
|---|---|---|---|
| M1 | 0.057 / 0.26 | 0.014 / 0.013 | **0.005 / 0.005** |
| M2 | 0.61 / 0.58 (faux) | 0.084 / 0.11 | 0.019 / 0.026 |

Sur ψ_min, l'écart de M1 à Erturk change de signe mais diminue (+5.8 %, −2.0 %, +0.6 % pour N = 64, 128, 256) : convergence
oscillante. Le très bon accord de M1 avec Ghia à N = 128 était en partie une compensation d'erreurs.

**Re = 3200** (Erturk ne publie pas ce Re) : ψ_min de M1 vaut −0.121, −0.121, −0.124 pour N = 64, 128, 256 ; celui de M2,
−0.031 (faux), −0.110, −0.119 (Richardson sur N = 128, 256 : −0.121). Le temps n'est pas en cause (à l'arrêt, K ne bouge plus que de
0.1 à 0.3 %). **Explication la plus probable, déjà établie en Partie A** (ADR 0002 §7.3) : la solution stationnaire de M1 dépend
de Δt (projection non incrémentale du cours : au stationnaire, la divergence vaut Δt (Δ_h − D G̃) p, ordre ≈ 0.65 en Δt), celle de M2
non. Or à haut Re, Δt = 0.8 × 2ν est **le même pour toutes les grilles** : cette part de l'erreur de M1 ne diminue pas quand on
raffine, d'où une convergence en maillage non monotone. Test de confirmation : M1 à Re = 3200 avec Δt/2 et Δt/4.
Premier résultat (N = 64) : **la solution de M1 dépend bien de Δt à haut Re** : K = 0.0471, 0.0476, 0.0479 pour Δt, Δt/2, Δt/4
(+1.0 % puis +0.7 %), et le tourbillon s'intensifie quand Δt diminue. N = 128 et 256 : en cours (à compléter).

---

## 5. P et les coins du couvercle

P dépend de N à cause de la singularité des coins du couvercle (saut de vitesse) : terme en ν ln N (Partie A). On ne compare donc
les P qu'à N fixé. Dans P(Re), ce terme donne une correction en Re⁻¹ qui rend la pente globale plus forte que −1/2 (≈ −0.7) ; les
pentes locales tendent vers −1/2 (M2 : −0.54 à partir de Re = 2500 ; M1 : −0.62, terme de coin plus grand).

---

## 6. Ce qui reste ouvert

- Résultat du test Woods (§3.2) et du test sur Δt de M1 (§4).
- Pour les profils à Re = 3200, pas de référence plus fine que Ghia.
