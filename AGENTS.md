# AGENTS.md

Les conventions pour les assistants IA (style, Numba, docstrings en français, tests, workflow)
sont dans [CLAUDE.md](CLAUDE.md), qui fait foi. Ce fichier n'existe que pour les outils qui
cherchent `AGENTS.md` ; ne pas y dupliquer de règles.

Ordre de priorité : énoncé ([Doc/FinalProject_LidDrivenCavity_HighRe.pdf](Doc/FinalProject_LidDrivenCavity_HighRe.pdf))
> [PLAN.md](PLAN.md) > [CLAUDE.md](CLAUDE.md).

Commande de test : `python -m pytest tests -q` dans l'environnement conda `CFD`
(sous Windows, si `conda` n'est pas dans le PATH : `~/miniconda3/envs/CFD/python.exe -m pytest tests -q`).
