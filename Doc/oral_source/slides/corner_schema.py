"""Schéma : au coin du couvercle, u saute de U à 0 ; les isovaleurs de u partent toutes du coin,
donc ∂u/∂y ~ 1/r, et la grille coupe la singularité à r ≈ h."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

RED, BORD, LIL = "#B51127", "#63003C", "#967DAB"
plt.rcParams.update({"font.family": "Arial", "font.size": 15})
fig, ax = plt.subplots(figsize=(5.4, 4.4), dpi=200)
cmap = LinearSegmentedColormap.from_list("cs", [BORD, LIL, RED])

# isovaleurs de u : rayons issus du coin (0, 1), de u = 0 (paroi) à u = 1 (couvercle)
for k, frac in enumerate(np.linspace(0.08, 0.92, 9)):
    th = -np.pi / 2 * (1 - frac)          # angle mesuré depuis le couvercle
    r = np.linspace(0, 1.25, 50)
    ax.plot(r * np.cos(th), 1 + r * np.sin(th), color=cmap(frac), lw=2)
ax.plot([0, 1.4], [1, 1], color=RED, lw=8, solid_capstyle="butt")
ax.plot([0, 0], [1, -0.4], color=BORD, lw=8, solid_capstyle="butt")
ax.annotate("", xy=(1.35, 1.14), xytext=(0.85, 1.14),
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.5, mutation_scale=20))
ax.text(0.30, 1.14, "u = U", color=RED, fontsize=16, fontweight="bold", va="center")
ax.text(-0.08, 0.05, "u = 0", color=BORD, fontsize=16, fontweight="bold", rotation=90, ha="right")
circ = plt.Circle((0, 1), 0.22, fill=False, ls="--", lw=1.8, color="0.25")
ax.add_patch(circ)
ax.text(-0.2, 1.22, "r ≈ h", fontsize=15, color="0.2", ha="right")
ax.text(0.68, -0.42, "isovaleurs de u :\n∂u/∂y ~ 1/r", fontsize=14, color=BORD)
ax.set_xlim(-0.6, 1.5); ax.set_ylim(-0.5, 1.3)
ax.set_aspect("equal"); ax.axis("off")
fig.tight_layout()
fig.savefig(sys.argv[1], transparent=True)
