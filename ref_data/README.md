# Données de référence — le strict minimum

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

**Coquilles connues dans l'article** (on ignore simplement ces points) :
- u, Re = 3200, y = 0.4531 : −0.86636 (sans doute −0.08636)
- u, Re = 10000, y = 0.5 : +0.03111 (sans doute −0.03111)
- v, Re = 400, x = 0.9063 : −0.23827 (douteuse)

## Re critique (Partie C)

L'énoncé donne Re_c ≈ 8000. Les valeurs publiées vont d'environ **7400 à 8050** selon la méthode et le
maillage, et cette dispersion vaut la peine d'être commentée. Deux points d'appui suffisent :

- Peng et al. (2003) : Re_c ≈ 7402 ; Fortin et al. (1997) : Re_c ≈ 7998.5
  (valeurs citées par Erturk et al. 2005, arXiv:cs/0411047)
- Auteri et al. (2002), J. Comput. Phys. 183:1–25 : Re_c ≈ 8018 (**à vérifier** sur l'article)

## Références

1. U. Ghia, K. N. Ghia, C. T. Shin, J. Comput. Phys. 48 (1982) 387–411.
2. E. Erturk, T. C. Corke, C. Gökçöl, Int. J. Numer. Meth. Fluids 48 (2005) 747–774 (arXiv:cs/0411047).
3. F. Auteri, N. Parolini, L. Quartapelle, J. Comput. Phys. 183 (2002) 1–25.
