# Source du PDF d'oral

`oral_A_B.html` produit `../oral_parties_A_B.pdf` (mise en page reprise du poly *Numerical Methods* : Times, formules KaTeX,
cadres « À dire », encadrés de questions). Formules chargées depuis cdn.jsdelivr.net (connexion nécessaire). Régénérer, depuis ce dossier :

    "C:\Program Files\Google\Chrome\Application\chrome.exe" --headless=new --no-pdf-header-footer --allow-file-access-from-files --virtual-time-budget=20000 --print-to-pdf=..\oral_parties_A_B.pdf oral_A_B.html

Les figures sont lues dans `../../figures/` : exécuter d'abord les notebooks A et B.

## Documents de limites

`md_to_pdf.py` convertit un `.md` de `Doc/` en PDF dans la même mise en page (page de garde, table des matières, tableaux,
formules `$...$` et `$$...$$` en KaTeX). Il faut le paquet `markdown` (`pip install markdown`). Depuis ce dossier :

    python md_to_pdf.py ../partie_A_limites.md "Partie A : choix intermédiaires et vérifications complémentaires"
    python md_to_pdf.py ../partie_B_limites.md "Partie B : difficultés rencontrées, diagnostics et limites"
