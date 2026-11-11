# Données de référence — le strict minimum

Articles complets dans [`doc_ref/`](../doc_ref/) : Ghia et al. (1982) et Erturk et al. (2005, `0411047v3.pdf`).

## Profils de Ghia et al. (1982)

Deux fichiers, une colonne par Re (100, 400, 1000, 3200, 5000, 7500, 10000) :

- `ghia1982_u_vertical_centerline.csv` : u(x = 0.5, y), Table I de l'article
- `ghia1982_v_horizontal_centerline.csv` : v(x, y = 0.5), Table II de l'article

Lecture en une ligne :

```python
data = np.loadtxt("ref_data/ghia1982_u_vertical_centerline.csv", delimiter=",")
y, u_re100 = data[:, 0], data[:, 1]      # colonnes : y, Re100, Re400, Re1000, ...
```

On s'en sert pour valider la Partie A (Re = 100) et pour comparer les profils de la Partie B.

### Ce que Ghia a fait (vérifié sur l'article, §3–4)

- **Formulation** ψ–ω (comme notre M2), **calcul stationnaire** (le pas de temps ne sert qu'à relaxer) :
  pas de K(t), pas de P.
- **Ordre 2** : diffusion centrée ; advection upwind d'ordre 1 + terme d'ordre 2 en correction différée
  (ordre 2 à convergence) ; ω de paroi par un nœud fantôme (ψ_y d'ordre 3, ω d'ordre 2, éq. 3.4–3.5).
- **Grilles uniformes** : **129×129 pour Re ≤ 3200**, **257×257 pour Re = 5000, 7500, 10000**.
- **Pas d'étude de convergence chiffrée** : seul contrôle, les lignes de courant à Re = 400 sur 129 et 257
  (visuel). Constats : le pas de la grille fine est le paramètre dominant ; ψ et ω au centre des tourbillons
  sont bien plus sensibles au maillage que les profils de vitesse.

**Conséquence pour nous** : Ghia a sa propre erreur de discrétisation. Quand notre écart à Ghia arrête de
diminuer, ce palier peut être **le sien**. On complète donc par l'écart entre nos grilles N et 2N.

### Vérification des CSV

Les CSV viennent d'une copie en ligne (gists ivan-pi). Ils ont été recoupés avec le texte extrait de
l'article (`doc_ref/`) : les 119 valeurs de u sont retrouvées ; pour v, deux écarts :
- v, Re = 1000, x = 0.9063 : la copie donnait −0.51500, l'article donne **−0.51550** ⇒ **corrigé** dans le CSV ;
- v, Re = 400, x = 0.5 : CSV 0.05186, l'extraction lit 0.05188 (confusion 6/8 possible) ⇒ laissé, **à
  vérifier à l'œil** sur le PDF (écart de 2·10⁻⁵, sans effet sur nos comparaisons).

**Coquilles connues dans l'article** (on ignore simplement ces points) :
- u, Re = 3200, y = 0.4531 : −0.86636 (sans doute −0.08636)
- u, Re = 10000, y = 0.5 : +0.03111 (sans doute −0.03111)
- v, Re = 400, x = 0.9063 : −0.23827 (douteuse)

## Erturk et al. (2005) : la vraie étude de maillage

- ψ–ω stationnaire, **centré d'ordre 2 + Thom** aux parois (comme notre M2), résidus < 10⁻¹⁰.
  Ils rappellent (§3) que Thom, d'ordre 1 localement, **préserve l'ordre 2 global** : utile pour l'oral.
- Grilles 401, 513, 601 + **extrapolation de Richardson répétée** (h⁴ puis h⁶ en théorie) sur ψ et ω au
  centre du tourbillon principal (Tables 2–3). Pas de K, pas de P.
- **Résultat utile pour la Partie C** : sur 129×129, pas de solution stationnaire à Re = 7500 (calcul
  périodique) ; sur 257×257, stationnaire jusqu'à Re = 12 500. Une partie de l'instationnarité sur grille
  grossière est donc **numérique**. Nuance : un solveur stationnaire peut converger vers une solution
  stationnaire **instable**, ce qui ne contredit pas Re_c ≈ 8000.
- Profils u, v sur les lignes médianes en 601×601 (Tables 6–7) : référence plus précise que Ghia pour
  Re ≥ 1000, si besoin (non extraits ici).

## Re critique (Partie C)

L'énoncé donne Re_c ≈ 8000. Les valeurs publiées vont d'environ **7400 à 8050** selon la méthode et le
maillage, et cette dispersion vaut la peine d'être commentée. Deux points d'appui suffisent :

- Peng et al. (2003) : Re_c ≈ 7402 ; Fortin et al. (1997) : Re_c ≈ 7998.5 ; Poliashenko & Aidun : ≈ 7763
  (valeurs citées par Erturk et al. 2005, §4, vérifié dans `doc_ref/`)
- Auteri et al. (2002), J. Comput. Phys. 183:1–25 : Re_c ≈ 8018 (**à vérifier** sur l'article)

## Références

1. U. Ghia, K. N. Ghia, C. T. Shin, J. Comput. Phys. 48 (1982) 387–411.
2. E. Erturk, T. C. Corke, C. Gökçöl, Int. J. Numer. Meth. Fluids 48 (2005) 747–774 (arXiv:cs/0411047).
3. F. Auteri, N. Parolini, L. Quartapelle, J. Comput. Phys. 183 (2002) 1–25.
