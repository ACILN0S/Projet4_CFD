"""Convertit les tableaux extraits (gists Ghia, texte pdftotext d'Erturk) en CSV ref_data/."""
import re
import sys
from pathlib import Path

DL = Path(sys.argv[1])
OUT = Path(sys.argv[2])

GHIA = ("Ghia U., Ghia K.N., Shin C.T. (1982), High-Re solutions for incompressible flow using the "
        "Navier-Stokes equations and a multigrid method, J. Comput. Phys. 48(3):387-411, "
        "doi:10.1016/0021-9991(82)90058-4")
ERTURK = ("Erturk E., Corke T.C., Gokcol C. (2005), Numerical solutions of 2-D steady incompressible "
          "driven cavity flow at high Reynolds numbers, Int. J. Numer. Meth. Fluids 48(7):747-774, "
          "arXiv:cs/0411047")


def numeric_rows(text, ncols):
    rows = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) != ncols:
            continue
        try:
            rows.append([float(p) for p in parts])
        except ValueError:
            continue
    return rows


def write_wide(path, header_lines, coord, res, rows, fmt_coord="{:.4f}", fmt_val="{:.5f}"):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        for h in header_lines:
            f.write(f"# {h}\n")
        f.write(",".join([coord] + [f"Re{r}" for r in res]) + "\n")
        for r in rows:
            f.write(",".join([fmt_coord.format(r[0])] + [fmt_val.format(v) for v in r[1:]]) + "\n")


ghia_re = [100, 400, 1000, 3200, 5000, 7500, 10000]
u = numeric_rows((DL / "ghiau.txt").read_text(), 8)
v = numeric_rows((DL / "ghiav.txt").read_text(), 8)
assert len(u) == 17 and len(v) == 17, (len(u), len(v))

common_ghia = [
    f"Source : {GHIA}",
    "Grille de Ghia : 129x129 (Re <= 3200) et 257x257 (Re >= 5000) ; methode psi-omega, multigrille.",
    "Transcription : gists de I. Pisarev (ivan-pi) 3e9326d18a366ffe6a8e5bfda6353219 / caa6c6737d36a9140fbcf2ea59c78b3c,",
    "  recoupee a la main avec la table publiee ; valeurs recopiees TELLES QU'IMPRIMEES (coquilles incluses).",
    "Adimensionnement : L = U = 1. Valeurs suspectes : voir ref_data/README.md et cavity/ref_data.py (GHIA_SUSPECT).",
]
write_wide(OUT / "ghia1982_u_vertical_centerline.csv",
           ["Ghia (1982) TABLE I : u(x = 0.5, y) le long de la ligne verticale mediane."] + common_ghia,
           "y", ghia_re, u)
write_wide(OUT / "ghia1982_v_horizontal_centerline.csv",
           ["Ghia (1982) TABLE II : v(x, y = 0.5) le long de la ligne horizontale mediane."] + common_ghia,
           "x", ghia_re, v)

erturk_re = [1000, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 21000]
common_erturk = [
    f"Source : {ERTURK}",
    "Grille 601x601, psi-omega differences finies ordre 2 (solution stationnaire).",
    "Extraction : pdftotext -raw de l'arXiv cs/0411047 (v1), pages 32-33, aucune retouche.",
]
eu = numeric_rows((DL / "p32.txt").read_text(), 11)
ev = numeric_rows((DL / "p33.txt").read_text(), 11)
assert len(eu) == 23 and len(ev) == 23, (len(eu), len(ev))
write_wide(OUT / "erturk2005_u_vertical_centerline.csv",
           ["Erturk et al. (2005) Table 6 : u(x = 0.5, y)."] + common_erturk, "y", erturk_re, eu,
           fmt_coord="{:.3f}", fmt_val="{:.4f}")
write_wide(OUT / "erturk2005_v_horizontal_centerline.csv",
           ["Erturk et al. (2005) Table 7 : v(x, y = 0.5)."] + common_erturk, "x", erturk_re, ev,
           fmt_coord="{:.3f}", fmt_val="{:.4f}")

# --- Tourbillons (Erturk Table 5) ---------------------------------------------------------
lines = (DL / "p31.txt").read_text().splitlines()
vortex_rows = []
label_re = re.compile(r"^(Primary|BR1|BL1|BR2|BL2|TL1|BR3|BL3|TL2)\s+(.*)$")


def tok(s):
    return s.split()


i = 0
while i < len(lines):
    m = label_re.match(lines[i].strip())
    if not m:
        i += 1
        continue
    name = m.group(1)
    if name == "Primary":
        psi = tok(m.group(2))
        omg = tok(lines[i + 1].replace("Vortex", ""))
        xy = lines[i + 2]
        name = "primary"
    else:
        psi = tok(lines[i - 1])
        omg = tok(m.group(2))
        xy = lines[i + 1]
    xys = re.findall(r"\(([\d.]+) , ([\d.]+)\)|(?<![\d.E])-(?![\d.])", xy.replace("(x , y)", ""))
    assert len(psi) == len(omg) == len(xys) == 10, (name, len(psi), len(omg), len(xys))
    for re_, p, o, (x, y) in zip(erturk_re, psi, omg, xys):
        if p == "-":
            continue
        vortex_rows.append(("Erturk2005", re_, name, float(p), float(o), float(x), float(y), "601x601", "Table 5"))
    i += 1

# Re = 1000, comparaison (Erturk Table 1, signes tels qu'imprimes : psi > 0 en valeur absolue)
t1 = [
    ("Erturk2005_extrap", 1000, "primary", -0.118942, -2.067213, None, None, "extrapole h6", "Erturk Table 1"),
    ("BotellaPeyret1998", 1000, "primary", -0.1189366, -2.067753, 0.5308, 0.5652, "Chebyshev N=160",
     "cite via Erturk Table 1 (non verifie sur l'original)"),
    ("Ghia1982", 1000, "primary", -0.117929, -2.04968, 0.5313, 0.5625, "129x129",
     "cite via Erturk Table 1 (non verifie sur l'original)"),
]
# Erturk Table 4 : valeurs de Ghia au centre du tourbillon principal
t4 = [("Ghia1982", r, "primary", p, o, None, None, "257x257", "cite via Erturk Table 4 (non verifie sur l'original)")
      for r, p, o in [(5000, -0.118966, -1.86016), (7500, -0.119976, -1.87987), (10000, -0.119731, -1.88082)]]

with open(OUT / "vortices.csv", "w", encoding="utf-8", newline="\n") as f:
    f.write("# Proprietes des tourbillons : psi et omega au centre, position (x, y).\n")
    f.write("# Convention de signe (Erturk) : psi < 0 et omega < 0 dans le tourbillon principal (rotation horaire),\n")
    f.write("#   psi > 0 dans les tourbillons secondaires (BR1, BL1, TL1), < 0 dans les tertiaires (BR2, BL2, TL2).\n")
    f.write("# Noms : B/T = bottom/top, L/R = left/right, 1/2/3 = secondaire/tertiaire/quaternaire.\n")
    f.write(f"# Sources : {ERTURK} (Tables 1, 4, 5) ; Botella & Peyret (1998) Comput. Fluids 27:421-433 et\n")
    f.write("#   Ghia (1982) CITES PAR Erturk (colonne note) : a reverifier sur les articles originaux si possible.\n")
    f.write("# Signes de Table 1 (imprimes en valeur absolue par Erturk) ramenes a la convention ci-dessus.\n")
    f.write("source,Re,vortex,psi,omega,x,y,grid,note\n")
    for r in vortex_rows + t1 + t4:
        vals = ["" if v is None else str(v) for v in r]
        f.write(",".join(vals) + "\n")

print("OK", len(vortex_rows), "tourbillons Erturk")
