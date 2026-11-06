"""Chargement des données de référence de la littérature (dossier ``ref_data/``).

Les fichiers CSV commencent par des lignes de commentaire ``#`` (source, tableau, remarques),
puis une ligne d'en-tête, par exemple ``y,Re100,Re400,...`` : la première colonne est la
coordonnée le long de la ligne médiane, les suivantes une colonne par nombre de Reynolds.

Toutes les grandeurs sont adimensionnées avec L = U = 1 (comme dans notre code).
"""

from pathlib import Path

import numpy as np

REF_DIR = Path(__file__).resolve().parent.parent / "ref_data"

# Profils disponibles : nom court -> fichier CSV
PROFILE_FILES = {
    "ghia_u": "ghia1982_u_vertical_centerline.csv",    # u(0.5, y), Ghia Table I
    "ghia_v": "ghia1982_v_horizontal_centerline.csv",  # v(x, 0.5), Ghia Table II
    "erturk_u": "erturk2005_u_vertical_centerline.csv",    # u(0.5, y), Erturk Table 6
    "erturk_v": "erturk2005_v_horizontal_centerline.csv",  # v(x, 0.5), Erturk Table 7
}

# Valeurs de Ghia (1982) manifestement erronées dans la table publiée : (profil, Re, coordonnée).
# Elles sont conservées telles quelles dans les CSV et masquées (NaN) au chargement par défaut.
#  - Table I, Re = 3200, y = 0.4531 : -0.86636 (|u| incompatible avec les voisins -0.04272 et
#    -0.24427 ; probablement -0.08636, coquille d'un chiffre).
#  - Table I, Re = 10000, y = 0.5 : +0.03111 (signe opposé à toutes les autres colonnes et à
#    Erturk 2005, u = -0.0287 à Re = 10000 ; probablement -0.03111).
#  - Table II, Re = 400, x = 0.9063 : -0.23827 (rompt la monotonie entre -0.22847 et -0.44993 ;
#    signalée comme douteuse par la transcription d'origine).
GHIA_SUSPECT = (
    ("ghia_u", 3200, 0.4531),
    ("ghia_u", 10000, 0.5000),
    ("ghia_v", 400, 0.9063),
)


def _read_csv(path):
    """Lit un CSV à en-tête simple, en ignorant les lignes de commentaire ``#``."""
    lines = [ln.strip() for ln in Path(path).read_text(encoding="utf-8").splitlines()]
    lines = [ln for ln in lines if ln and not ln.startswith("#")]
    header = lines[0].split(",")
    data = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    return header, data


def load_profile(name, mask_suspect=True):
    """Charge un profil de vitesse le long d'une ligne médiane.

    Paramètres
    ----------
    name : str
        Clé de ``PROFILE_FILES`` (``"ghia_u"``, ``"ghia_v"``, ``"erturk_u"``, ``"erturk_v"``).
    mask_suspect : bool
        Si vrai, remplace par NaN les valeurs listées dans ``GHIA_SUSPECT``.

    Retour
    ------
    dict
        ``{Re: (coord, values)}`` avec ``coord`` trié par ordre croissant (tableaux numpy).
    """
    if name not in PROFILE_FILES:
        raise KeyError(f"Profil inconnu '{name}'. Choix possibles : {sorted(PROFILE_FILES)}")
    header, data = _read_csv(REF_DIR / PROFILE_FILES[name])
    order = np.argsort(data[:, 0])
    coord = data[order, 0]

    profiles = {}
    for col, label in enumerate(header[1:], start=1):
        reynolds = int(label.removeprefix("Re"))
        values = data[order, col].copy()
        if mask_suspect:
            for prof, re_bad, c_bad in GHIA_SUSPECT:
                if prof == name and re_bad == reynolds:
                    values[np.isclose(coord, c_bad)] = np.nan
        profiles[reynolds] = (coord.copy(), values)
    return profiles


def load_vortices():
    """Charge ``vortices.csv`` : liste de dicts (source, Re, vortex, psi, omega, x, y, grid, note).

    Les champs numériques absents (position non publiée) valent ``None``.
    """
    lines = (REF_DIR / "vortices.csv").read_text(encoding="utf-8").splitlines()
    lines = [ln for ln in lines if ln and not ln.startswith("#")]
    keys = lines[0].split(",")
    numeric = {"psi", "omega", "x", "y"}
    records = []
    for ln in lines[1:]:
        raw = dict(zip(keys, ln.split(",")))
        rec = {k: (float(v) if v else None) if k in numeric else v for k, v in raw.items()}
        rec["Re"] = int(raw["Re"])
        records.append(rec)
    return records
