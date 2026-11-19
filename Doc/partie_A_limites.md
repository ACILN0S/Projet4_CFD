# Partie A — choix intermédiaires, difficultés et vérifications complémentaires

Ce document accompagne le notebook A (`notebooks/A_Re100.ipynb`), qui ne garde que les résultats. Les chiffres ci-dessous viennent de
l'ancienne version longue du notebook (archivée dans `Projet4_CFD_version_optimisee_2026-11-19.zip`, dossier parent), calculée avec le
premier critère d'arrêt de la Partie A (ε = 10⁻², t_ref = 10). Décisions détaillées : ADR 0002 (pression), ADR 0003 (Δt, arrêt, tolérance).

---

## 1. Condition de Neumann sur la pression : l'ordre 1 abandonné pour des nœuds fantômes

**Première version.** Neumann décentrée d'ordre 1 (cours §4.2.2.2) : « bord = voisin intérieur », φ₀ = φ₁. Test sur
φ = cos(πx)cos(πy) : **ordre 1** (pentes 0.82 puis 0.92), erreur maximale ≈ 2.3 h à 2.7 h.

**Pourquoi.** φ₀ = φ₁ impose ∂φ/∂n = 0 **entre** les nœuds 0 et 1 : le schéma résout (à l'ordre 2) le problème de Neumann sur le domaine
décalé [h/2, 1 − h/2]². En 1D, l'écart à la vraie solution est une pente linéaire d'amplitude ≈ (π²/2) h (x − ½), soit une erreur
maximale ≈ (π²/4) h ≈ 2.5 h : c'est ce qu'on mesure. Les courbes « erreur / h » des trois grilles se superposent : l'erreur est globale
et régulière, pas limitée au bord.

**Version retenue.** Nœuds fantômes miroirs φ₋₁ = φ₁ (dérivée normale centrée d'ordre 2). Les nœuds de paroi deviennent des inconnues ;
les lignes de paroi sont pondérées par les poids des trapèzes (½ sur une paroi, ¼ dans un coin) pour garder la matrice symétrique, ce
qu'exige le gradient conjugué ; compatibilité et moyenne nulle deviennent des moyennes pondérées (∬ φ dS = 0). Résultat : **ordre 2.00**,
erreur ≈ 0.82 h², environ 200 fois plus petite à N = 64, pour ~40 % d'itérations en plus.

**Comment on l'a trouvé.** La comparaison M1 / M2 : avec la Neumann d'ordre 1, M1 convergeait à l'ordre ≈ 1.4 et M2 à l'ordre 2. Une
pente faible pour M1 seul désignait ce qui est propre à M1, la pression (ADR 0002 §5).

**Piège de lecture.** À N = 128, la version d'ordre 1 était **plus proche de Ghia** (0.004 / 0.004) que la version à fantômes (0.006 / 0.012).
Elle n'était pas plus juste : plus loin de M2 et de pente moins bonne, elle traversait les valeurs de Ghia en convergeant. On juge donc
sur l'écart entre grilles et sur l'écart à l'autre méthode, pas sur Ghia seul.

---

## 2. L'erreur de projection de M1 dépend de Δt

**Constat.** Avec le Δt de stabilité, l'écart entre grilles successives de M1 converge à l'ordre **≈ 1.0**, alors que M2 converge à l'ordre 2.

**Explication (identité exacte, ADR 0002 §7.3).** Un pas de M1 : $u^* = u + \Delta t\,F(u)$, $\Delta_h p = D(u^*)/\Delta t$, $u^{n+1} = u^* - \Delta t\,\tilde G p$, où $\Delta_h$ est le
laplacien que résout le gradient conjugué et D G̃ celui qu'applique la correction (gradient centré, non appliqué aux parois). Au
stationnaire, F(u) = G p (sans Δt) et

$$D(u) = \Delta t\,\bigl(\Delta_h - D\tilde G\bigr)\,p .$$

Sur grille colocalisée, les deux laplaciens ne coïncident pas : la solution stationnaire de M1 garde une divergence proportionnelle à Δt et
**dépend de Δt**. M2 n'a pas d'étape de projection : sa solution stationnaire n'en dépend pas. C'est une propriété de la projection non
incrémentale du cours (§7.3.2).

**Mesures (N = 64).**

| | Δt → Δt/2 | Δt/2 → Δt/4 |
|---|---|---|
| M1, écart max des profils | 1.8·10⁻³ | 1.1·10⁻³ |
| M2, écart max des profils | 1·10⁻⁷ (erreur d'arrêt) | – |

Ordre observé en Δt ≈ 0.65 (rapport 1.55 et non 2). **À Δt fixe** (celui de N = 128 imposé aux trois grilles), M1 converge à l'ordre
**1.91** en espace : la Neumann à fantômes l'a bien rendu d'ordre 2.

**Conséquence.** Avec le Δt de stabilité, Δt change avec N (CFL à N = 32, puis Fourier) : la part d'erreur en Δt change d'une grille à
l'autre et domine les écarts entre grilles. Ce n'est pas un ordre 1 en h. Dans la Partie B, à haut Re, Δt = 0.8 × 2ν est **le même pour
toutes les grilles** : cette erreur ne diminue plus quand on raffine (voir `partie_B_limites.md` §4).

**Explication réfutée.** « Erreur en Δt/h à cause de ∂p/∂n = 0 non physique » : elle prévoyait (Δ_h − D G̃) p ≈ c g/h aux parois. Mesuré :
≈ constant (1.2, 1.3, 1.9 pour N = 32, 64, 128) alors que g/h double ; le plus grand écart est aux coins du couvercle.

---

## 3. Critère d'arrêt

- **Résidu stationnaire** $R = (f^{n+1} - f^n)/\Delta t$, norme L2 aux nœuds intérieurs (u, v pour M1 ; ω pour M2). Avec Euler explicite, R est
  exactement le membre de droite discret : il mesure de combien le champ viole les équations stationnaires, et ne dépend pas de Δt. Le
  fluide partant du repos, le résidu initial est nul : la référence est prise à t_ref.
- **À Re = 100** : passer de ε = 10⁻² à 10⁻³ changeait les profils de 2.5·10⁻⁵, soit 0.4 % de l'écart à Ghia. Le résidu décroît en
  $e^{-\lambda t}$ avec **λ ≈ 0.54, le même pour M1, M2 et les trois grilles** : c'est une propriété de l'écoulement (mode le plus lent), et les
  deux méthodes s'arrêtaient au même instant (t ≈ 18.5).
- **Critère final** : ε = 10⁻³, t_ref = 30, choisi au cas le plus défavorable de la Partie B (Re = 5000) et appliqué partout, y compris à
  Re = 100. Conséquence à Re = 100 : l'écoulement est stationnaire avant t = 30, ‖R(30)‖ est déjà au niveau du bruit du gradient conjugué,
  et M1 à N = 128 va jusqu'au garde-fou t = 1000 (résultat juste, plus convergé que nécessaire).
- **Premier pas de M2** : ‖R‖ = 0, car la vorticité de paroi (Thom) n'est calculée qu'à la fin du pas ; M2 « voit » le couvercle avec un
  pas de retard, sans effet visible dès t ≈ 0.1.

---

## 4. Tolérance du gradient conjugué

tol = 10⁻⁸ (résidu relatif) contre 10⁻¹⁰ à N = 128 : les profils changent de ~10⁻¹⁰, sept ordres de grandeur sous l'écart à Ghia, et P
ne change pas. tol = 10⁻⁸ suffit pour toutes les grilles ; 10⁻¹⁰ doublerait le nombre d'itérations pour rien (ADR 0003 §3).

---

## 5. P et les coins du couvercle

Aux deux coins du couvercle, la vitesse de paroi saute de U à 0 : ∂u/∂y ~ 1/r, donc ∫ ν ∂u/∂y dx diverge comme ln en continu. La grille
coupe à r ≈ h : $P \approx P_{\text{régulier}} + c\,\nu \ln N$. Mesuré (M1) : P = 0.154, 0.182, 0.212 pour N = 32, 64, 128, accroissements 0.028 puis 0.030,
c ≈ 4.2 ν. P ne converge pas en maillage et ne sert pas à choisir N ; on compare toujours les P à N et méthode fixés. K, lui, converge.

---

## 6. Erreur propre de Ghia

L'écart de M2 à Ghia **ne diminue plus** entre N = 64 et 128 (0.003 → 0.004 en u), alors que M2 converge à l'ordre 2 et que son N = 128 est
à ~10⁻³ de sa limite. On mesure alors l'erreur de Ghia (129², ordre 2), de l'ordre de 4 à 8·10⁻³, plus grande sur le minimum de v.

---

## 7. Coût

Numba (recommandé par le prof) : noyaux compilés, boucles explicites reprenant les formules du cours. Vérifié contre l'ancienne version
numpy pure : mêmes champs à 2·10⁻¹⁰ près, 3 à 5 fois plus rapide. À Re = 100 : ~1 s (N = 32), ~15 s (N = 64), quelques minutes (N = 128) ;
chaque doublement de N coûte ~16 fois plus (4 fois plus de nœuds, Δt ÷ 4 par la condition de Fourier).
