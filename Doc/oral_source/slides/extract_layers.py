"""Découpe l'image des anneaux (slide de titre du template) en calques par couleur,
chacun recentré sur son centre de rotation, pour les animer en rotation dans PowerPoint."""
import json, sys
import numpy as np
from PIL import Image
from scipy import ndimage

SRC, OUT = sys.argv[1], sys.argv[2]
im = np.asarray(Image.open(SRC).convert("RGBA")).astype(np.int32)
H, W = im.shape[:2]
rgb, a = im[..., :3], im[..., 3]

# disque blanc : centre et rayon
white = (a > 200) & (rgb.min(-1) > 245)
lab, n = ndimage.label(white)
big = np.argmax(np.bincount(lab.ravel())[1:]) + 1
ys, xs = np.nonzero(lab == big)
# ajustement d'un cercle (moindres carrés) sur le bord du disque blanc rempli
filled = ndimage.binary_fill_holes(lab == big)
edge = filled & ~ndimage.binary_erosion(filled)
ys, xs = np.nonzero(edge)
def fit(xs, ys):
    A_ = np.c_[2 * xs, 2 * ys, np.ones_like(xs)]
    sol, *_ = np.linalg.lstsq(A_, (xs**2 + ys**2).astype(float), rcond=None)
    return sol[0], sol[1], np.sqrt(sol[2] + sol[0]**2 + sol[1]**2)
# ajustement itératif : on écarte les points des encoches (arc lilas superposé)
for tol in (200, 60, 15, 4):
    cx, cy, R = fit(xs, ys)
    keep = np.abs(np.hypot(xs - cx, ys - cy) - R) < tol
    xs, ys = xs[keep], ys[keep]
cx, cy, R = fit(xs, ys)
res = np.hypot(xs - cx, ys - cy) - R
print("fit residual pct", np.percentile(np.abs(res), [50, 90, 99]))
print("disk", cx, cy, R)

yy, xx = np.mgrid[0:H, 0:W]
rr = np.hypot(xx - cx, yy - cy)
pal = {"lilac": (150, 125, 171), "dark": (99, 0, 60), "teal": (0, 172, 169), "pink": (245, 181, 210)}
names = list(pal)
d = np.stack([np.abs(rgb - np.array(pal[k])).sum(-1) for k in names], -1)
cls = np.argmin(d, -1)
vis = a > 0
inside = rr <= R + 6

layers = {}
for i, k in enumerate(names):
    layers[k] = vis & ~inside & (cls == i)
# l'arc lilas recouvre en partie le disque : on garde ses pixels francs à l'intérieur (hors emblème)
box = (xx > 2250) & (xx < 2700) & (yy > 770) & (yy < 1185)
lil_strict = vis & (np.abs(rgb - np.array(pal["lilac"])).sum(-1) < 90) & ~box
layers["lilac"] = layers["lilac"] | (inside & lil_strict)
layers["disk"] = vis & inside

# séparation de l'anneau lilas en ses deux arcs
lab, n = ndimage.label(ndimage.binary_dilation(layers["lilac"], iterations=2))
sizes = np.bincount(lab.ravel())[1:]
order = np.argsort(sizes)[::-1]
lil = layers.pop("lilac")
layers["lilacA"] = lil & (lab == order[0] + 1)
layers["lilacB"] = lil & (lab == order[1] + 1)
print("lilac comps", sizes[order[:3]])

# emblème (deux « c ») : boîte repérée sur l'image, centres = centres des pastilles
emb_zone = (xx > 2250) & (xx < 2700) & (yy > 770) & (yy < 1185)
# rouge et lilas mélangés au blanc : R - G vaut 163 t (rouge) contre 25 t (lilas)
colored = (a > 0) & (rgb.min(-1) < 245)
red = colored & (rgb[..., 0] > rgb[..., 2])  # rouge : R > B ; lilas : B > R
lilc = colored & ~red
EMB_C = {"embA": (2397, 913), "embB": (2555, 1043)}
embl = {}
notwhite = ~((rgb.min(-1) > 250) & (a > 200))
for k, m in (("embA", lilc & emb_zone), ("embB", red & emb_zone)):
    embl[k] = m
meta = {"W": W, "H": H, "cx": cx, "cy": cy, "R": R, "layers": {}}
SC = 2  # réduction de résolution

def save(name, mask, ccx, ccy, fill_white_mask=None):
    rgba = im.copy()
    rgba[..., 3] = np.where(mask, rgba[..., 3], 0)
    if fill_white_mask is not None:
        rgba[fill_white_mask, :3] = 255
        rgba[fill_white_mask, 3] = 255
    ys, xs = np.nonzero(rgba[..., 3] > 0)
    half = int(np.ceil(max(np.abs(xs - ccx).max(), np.abs(ys - ccy).max()))) + 4
    x0, y0 = int(round(ccx)) - half, int(round(ccy)) - half
    canvas = np.zeros((2 * half, 2 * half, 4), np.uint8)
    sx0, sy0 = max(x0, 0), max(y0, 0)
    sx1, sy1 = min(x0 + 2 * half, W), min(y0 + 2 * half, H)
    canvas[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = rgba[sy0:sy1, sx0:sx1]
    img = Image.fromarray(canvas).resize((2 * half // SC, 2 * half // SC), Image.LANCZOS)
    img.save(f"{OUT}/{name}.png", optimize=True)
    meta["layers"][name] = {"x0": x0, "y0": y0, "size": 2 * half}
    print(name, x0, y0, 2 * half)

for k in ("lilacA", "lilacB", "dark", "teal", "pink"):
    save(k, layers[k], cx, cy)
emb_all = embl["embA"] | embl["embB"]
# disque : cercle blanc complet + texte du logo (le lilas qui le recouvre est retiré)
lil_in = inside & (a > 0) & (rgb[..., 2] > rgb[..., 0] + 5) & (rgb.min(-1) < 245) & ~emb_zone
disk_mask = (rr <= R) | (filled & (rr <= R + 60))  # + zones blanches d'origine sous l'arc lilas
# hors du logo (texte + emblème), le disque est blanc pur : supprime le filet rouge du bord
txt = disk_mask & (rr < R - 40) & (np.abs(rgb - np.array(pal["dark"])).sum(-1) < 150) & (a > 0)
ty, tx = np.nonzero(txt | emb_all)
logo_box = (xx >= tx.min() - 20) & (xx <= tx.max() + 20) & (yy >= ty.min() - 20) & (yy <= ty.max() + 20)
print("logo box", tx.min(), tx.max(), ty.min(), ty.max())
fillw = (emb_all & inside) | lil_in | (disk_mask & (a < 255)) | (disk_mask & ~logo_box)
save("disk", disk_mask, cx, cy, fill_white_mask=fillw)
for k, m in embl.items():
    save(k, m, *EMB_C[k])
json.dump(meta, open(f"{OUT}/meta.json", "w"), indent=1)
