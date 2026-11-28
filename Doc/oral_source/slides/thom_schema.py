"""Schéma : sous le couvercle, Thom remplace le vrai profil (convexe) par sa tangente."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RED, BORD, M2, LIL = "#B51127", "#63003C", "#FF7F0E", "#967DAB"
plt.rcParams.update({"font.family": "Arial", "font.size": 15})
fig, ax = plt.subplots(figsize=(6.2, 4.6), dpi=200)
s = np.linspace(0, 1.35, 200)          # s = n / h, distance au couvercle en mailles
true = 1 - 0.95 * s + 0.24 * s**2      # profil convexe
tang = 1 - 0.95 * s                    # tangente à la paroi (Thom)
y = -s                                  # couvercle en haut

for yy, lab in ((0, "couvercle (nœud 0)"), (-1, "nœud 1  (n = h)")):
    ax.axhline(yy, color="0.75", lw=1, ls=":", zorder=0)
    ax.text(0.03 if yy == 0 else 0.72, yy + 0.04, lab, color="0.35", fontsize=12)
ax.plot([0, 1.15], [0, 0], color=RED, lw=7, solid_capstyle="butt")
ax.annotate("", xy=(1.17, 0.13), xytext=(0.88, 0.13),
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.5, mutation_scale=20))
ax.text(0.70, 0.11, "U", color=RED, fontsize=17, fontweight="bold", va="center")

m = s <= 1
ax.fill_betweenx(y[m], tang[m], true[m], color=LIL, alpha=0.45, lw=0)
ax.plot(true, y, color=BORD, lw=3, label="vrai profil (convexe)")
ax.plot(tang, y, color=M2, lw=3, ls="--", label="Thom : profil linéaire")
ax.annotate("R ≈ (h²/6) ∂²u/∂n² > 0", xy=(0.36, -0.75), xytext=(0.62, -0.62),
            color=BORD, fontsize=14,
            arrowprops=dict(arrowstyle="-", color=BORD, lw=1.2))
ax.set_xlim(-0.02, 1.25); ax.set_ylim(-1.5, 0.3)
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_xlabel("u", fontsize=15); ax.set_ylabel("y", fontsize=15, rotation=0, labelpad=12)
ax.legend(loc="lower right", frameon=False, fontsize=13)
fig.tight_layout()
fig.savefig(sys.argv[1], transparent=True)
