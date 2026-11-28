"""Construit la présentation des Parties A et B sur le template CentraleSupélec.

Usage : python build_deck.py <template.pptx> <dossier figures> <sortie.pptx>
Pré-requis (même dossier) : layers/ (extract_layers.py), img/ (thom_schema.py, corner_schema.py).
"""
import json
import sys
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.chart.data import CategoryChartData
from pptx.dml.color import RGBColor
from pptx.enum.chart import XL_CHART_TYPE, XL_LABEL_POSITION, XL_MARKER_STYLE, XL_LEGEND_POSITION
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Emu, Inches, Pt

import anim
from lxml import etree

TEMPLATE, FIGS, OUT = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
HERE = Path(__file__).parent
IMG, LAYERS = HERE / "img", HERE / "layers"

# --- palette du template + couleurs des figures (M1 bleu, M2 orange, comme matplotlib) -------------
RED, BORD, LILAC, TEAL, PINK = "B51127", "63003C", "967DAB", "00ACA9", "F5B5D2"
TINT, GREY, DARK, WHITE = "F3EFF6", "595959", "262626", "FFFFFF"
M1C, M2C = "1F77B4", "FF7F0E"
FONT = "Arial"


def rgb(h):
    return RGBColor.from_string(h)


SUP, SUB = {"base": 30000}, {"base": -25000}


def sup(t, **kw):
    return (t, {**SUP, **kw})


def sub(t, **kw):
    return (t, {**SUB, **kw})


# --- découpe des figures multi-panneaux --------------------------------------------------------
def crops():
    def crop(name, boxes, tag):
        im = Image.open(FIGS / f"{name}.png")
        for k, b in enumerate(boxes):
            im.crop(b).save(IMG / f"{tag}_{k}.png")
    crop("A_flow_K_P", [(40, 0, 600, 600), (680, 0, 1240, 600), (1262, 0, 1905, 600), (1905, 0, 2550, 600)], "flowA")
    crop("A_profiles_Ghia", [(0, 0, 825, 690), (828, 0, 1650, 690)], "profA")
    crop("A_poisson_convergence", [(0, 0, 825, 630)], "poisson")
    crop("A_grid_convergence", [(52, 0, 825, 630)], "gridA")
    crop("A_M1_P_lnN", [(0, 0, 825, 600)], "PlnN")
    crop("B_P_Re", [(0, 0, 975, 720), (975, 0, 1950, 720)], "PRe")
    cols = [0, 451, 896, 1342, 1787, 2233, 2700]
    crop("B_streamlines", [(cols[i], 56, cols[i + 1], 500) for i in range(6)], "strM1")
    crop("B_streamlines", [(cols[i], 502, cols[i + 1], 960) for i in range(6)], "strM2")
    pc = [0, 606, 1206, 1800, 2394, 3000]
    crop("B_profiles_Ghia", [(pc[i], 62, pc[i + 1], 1200) for i in range(5)], "profB")


# --- primitives de dessin ----------------------------------------------------------------------
def text(slide, x, y, w, h, paras, size=16, color=DARK, bold=False, align="l", anchor="t",
         italic=False, name=None, spacing=4):
    """paras : liste de paragraphes ; un paragraphe = str ou liste de (texte, {options de run})."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    if name:
        tb.name = name
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Inches(0.04)
    tf.margin_top = tf.margin_bottom = Inches(0.02)
    tf.vertical_anchor = {"t": MSO_ANCHOR.TOP, "m": MSO_ANCHOR.MIDDLE, "b": MSO_ANCHOR.BOTTOM}[anchor]
    _fill_tf(tf, paras, size, color, bold, align, italic, spacing)
    return tb


def _fill_tf(tf, paras, size, color, bold, align, italic=False, spacing=4):
    if isinstance(paras, str):
        paras = [paras]
    for i, p in enumerate(paras):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        para.alignment = {"l": PP_ALIGN.LEFT, "c": PP_ALIGN.CENTER, "r": PP_ALIGN.RIGHT}[align]
        para.space_after = Pt(spacing)
        runs = [(p, {})] if isinstance(p, str) else p
        for t, o in runs:
            r = para.add_run()
            r.text = t
            f = r.font
            f.name = FONT
            f.size = Pt(o.get("size", size))
            f.bold = o.get("bold", bold)
            f.italic = o.get("italic", italic)
            f.color.rgb = rgb(o.get("color", color))
            if "base" in o:
                f._rPr.set("baseline", str(o["base"]))


def box(slide, x, y, w, h, fill=TINT, line=None, lw=1.5, shape=MSO_SHAPE.ROUNDED_RECTANGLE,
        radius=0.12, paras=None, size=14, color=DARK, bold=False, align="c", anchor="m", name=None,
        dash=False, margin=0.1):
    sh = slide.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    if name:
        sh.name = name
    if shape == MSO_SHAPE.ROUNDED_RECTANGLE:
        sh.adjustments[0] = radius
    if fill:
        sh.fill.solid()
        sh.fill.fore_color.rgb = rgb(fill)
    else:
        sh.fill.background()
    if line:
        sh.line.color.rgb = rgb(line)
        sh.line.width = Pt(lw)
        if dash:
            sh.line.dash_style = 4  # tiret
    else:
        sh.line.fill.background()
    sh.shadow.inherit = False
    tf = sh.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Inches(margin)
    tf.margin_top = tf.margin_bottom = Inches(0.04)
    tf.vertical_anchor = {"t": MSO_ANCHOR.TOP, "m": MSO_ANCHOR.MIDDLE, "b": MSO_ANCHOR.BOTTOM}[anchor]
    if paras:
        _fill_tf(tf, paras, size, color, bold, align)
    return sh


def badge(slide, x, y, d, label, fill=RED, color=WHITE, size=16):
    """Pastille ronde : le motif visuel du deck (écho aux cercles du logo)."""
    return box(slide, x, y, d, d, fill=fill, shape=MSO_SHAPE.OVAL, paras=[label], size=size,
               color=color, bold=True, margin=0)


def pic(slide, path, x, y, w=None, h=None, name=None):
    p = slide.shapes.add_picture(str(path), Inches(x), Inches(y),
                                 Inches(w) if w else None, Inches(h) if h else None)
    if name:
        p.name = name
    return p


def arrow(slide, x, y, w, h, color=LILAC):
    sh = slide.shapes.add_shape(MSO_SHAPE.CHEVRON, Inches(x), Inches(y), Inches(w), Inches(h))
    sh.fill.solid()
    sh.fill.fore_color.rgb = rgb(color)
    sh.line.fill.background()
    sh.shadow.inherit = False
    return sh


def line(slide, x1, y1, x2, y2, color=BORD, w=2):
    c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1), Inches(y1), Inches(x2), Inches(y2))
    c.line.color.rgb = rgb(color)
    c.line.width = Pt(w)
    return c


def card(slide, x, y, w, h, head, body, badge_label=None, badge_fill=RED, head_color=BORD,
         fill=TINT, size=14, head_size=16):
    """Carte : fond teinté, pastille optionnelle, titre gras, corps. Renvoie les formes à animer."""
    bg = box(slide, x, y, w, h, fill=fill)
    shapes = [bg]
    tx = x + 0.2
    if badge_label:
        shapes.append(badge(slide, x + 0.18, y + 0.18, 0.5, badge_label, fill=badge_fill, size=13))
        tx = x + 0.82
    paras = [[(head, {"bold": True, "color": head_color, "size": head_size})]]
    paras += [p for p in ([body] if isinstance(body, str) else body)]
    shapes.append(text(slide, tx, y + 0.14, x + w - tx - 0.15, h - 0.2, paras, size=size, color=DARK))
    return shapes


def ids(*shapes):
    return [s.shape_id for s in shapes]


def notes(slide, txt):
    slide.notes_slide.notes_text_frame.text = txt


def chart_style(chart, title):
    chart.has_title = True
    tf = chart.chart_title.text_frame
    tf.text = ""
    _fill_tf(tf, [title] if isinstance(title, str) else [title], 14, BORD, True, "c", spacing=0)
    chart.font.name = FONT
    chart.font.size = Pt(12)
    chart.font.color.rgb = rgb(GREY)


def quiet_axes(chart, fmt=None, gridlines=True):
    va = chart.value_axis
    va.has_major_gridlines = gridlines
    if gridlines:
        va.major_gridlines.format.line.color.rgb = rgb("D9D9D9")
        va.major_gridlines.format.line.width = Pt(0.75)
    va.format.line.fill.background()
    va.tick_labels.font.size = Pt(12)
    if fmt:
        va.tick_labels.number_format = fmt
        va.tick_labels.number_format_is_linked = False
    ca = chart.category_axis
    ca.format.line.color.rgb = rgb("A6A6A6")
    ca.tick_labels.font.size = Pt(12)
    ca.has_major_gridlines = False


# --- construction ------------------------------------------------------------------------------
crops()
prs = Presentation(TEMPLATE)
# on vide le template (6 diapositives d'exemple)
sldIdLst = prs.slides._sldIdLst
for sldId in list(sldIdLst):
    prs.part.drop_rel(sldId.rId)
    sldIdLst.remove(sldId)
L = {l.name: l for l in prs.slide_layouts}
LAY_TITLE, LAY_SECTION, LAY_CONTENT = L["Diapositive de titre 2"], L["Titre de section"], L["Titre seul"]
meta = json.load(open(LAYERS / "meta.json"))
# image des anneaux telle que placée dans la disposition « Diapositive de titre 2 »
K = 5007429 / (meta["W"] * (1 - 0.28279))  # EMU par pixel de l'image source
OX, OY = 7184571, 67081


def content_slide(title):
    s = prs.slides.add_slide(LAY_CONTENT)
    tf = s.shapes.title.text_frame
    tf.text = ""
    runs = [(title, {})] if isinstance(title, str) else title
    para = tf.paragraphs[0]
    for t, o in runs:
        r = para.add_run()
        r.text = t
        r.font.name = FONT
        r.font.bold = True
        if "base" in o:
            r.font._rPr.set("baseline", str(o["base"]))
    return s


def ring_slide(title, subtitle, note):
    """Diapositive rouge du template, anneaux reconstitués en calques qui tournent à l'arrivée."""
    s = prs.slides.add_slide(LAY_TITLE)
    s._element.set("showMasterSp", "0")  # masque l'image fixe de la disposition
    bg = box(s, 0, 0, 13.333, 7.5, fill=RED, shape=MSO_SHAPE.RECTANGLE, name="Fond rouge")
    layer = {}
    for name in ("dark", "pink", "teal", "disk", "lilacA", "lilacB", "embA", "embB"):
        m = meta["layers"][name]
        size = Emu(int(m["size"] * K))
        p = s.shapes.add_picture(str(LAYERS / f"{name}.png"), Emu(int(OX + m["x0"] * K)),
                                 Emu(int(OY + m["y0"] * K)), size, size)
        p.name = f"Anneau {name}"
        layer[name] = p
    # placeholders du titre au premier plan
    tree = s.shapes._spTree
    for ph in list(s.placeholders):
        tree.remove(ph._element)
        tree.append(ph._element)
    tree.remove(bg._element)
    tree.insert(2, bg._element)
    t = s.shapes.title
    t.text_frame.text = ""
    _fill_tf(t.text_frame, title, 48, WHITE, True, "l")
    sub = [ph for ph in s.placeholders if ph.placeholder_format.idx == 1][0]
    sub.text_frame.text = ""
    _fill_tf(sub.text_frame, subtitle, 20, WHITE, False, "l", spacing=6)
    g = [
        (layer["dark"].shape_id, "spin_in", 0, 1700, {"angle": -150, "bld": "pic"}),
        (layer["lilacA"].shape_id, "spin_in", 150, 1600, {"angle": 200, "bld": "pic"}),
        (layer["lilacB"].shape_id, "spin_in", 300, 1700, {"angle": -260, "bld": "pic"}),
        (layer["teal"].shape_id, "spin_in", 250, 1900, {"angle": 130, "bld": "pic"}),
        (layer["pink"].shape_id, "spin_in", 400, 1800, {"angle": -170, "bld": "pic"}),
        (layer["disk"].shape_id, "zoom", 600, 800, {"bld": "pic"}),
        (layer["embA"].shape_id, "spin_in", 1250, 1100, {"angle": -360, "bld": "pic"}),
        (layer["embB"].shape_id, "spin_in", 1250, 1100, {"angle": 360, "bld": "pic"}),
        (t.shape_id, "fly_l", 1500, 800),
        (sub.shape_id, "rise", 2100, 700),
    ]
    anim.apply(s, [{"trigger": "auto", "items": g}], transition="fade")
    notes(s, note)
    return s


def section_slide(title, subtitle, steps, note):
    s = prs.slides.add_slide(LAY_SECTION)
    s.shapes.title.text = ""
    _fill_tf(s.shapes.title.text_frame, title, 54, RED, True, "l")
    body = [ph for ph in s.placeholders if ph.placeholder_format.idx == 1][0]
    body.text_frame.text = ""
    _fill_tf(body.text_frame, subtitle, 22, BORD, False, "l")
    items = [(s.shapes.title.shape_id, "fly_l", 0, 700), (body.shape_id, "rise", 500, 600)]
    n = len(steps)
    w = 11.6 / max(n, 1)
    if not subtitle:          # intercalaire réduit au titre : pas de zone de sous-titre
        body._element.getparent().remove(body._element)
        items = items[:1]
    for i, st in enumerate(steps):
        x = 0.9 + i * w
        b = badge(s, x, 4.95, 0.62, str(i + 1), fill=[RED, LILAC, BORD, TEAL, PINK][i % 5])
        t = text(s, x + 0.75, 4.92, w - 0.85, 0.75, st, size=15, color=DARK, anchor="m")
        items += [(b.shape_id, "zoom", 900 + 200 * i, 450), (t.shape_id, "fade", 1000 + 200 * i, 450)]
    anim.apply(s, [{"trigger": "auto", "items": items}], transition="push_u")
    notes(s, note)
    return s


# =============================================================================================
# 1. Titre
ring_slide(["Cavité entraînée", "à haut Reynolds"],
           ["Projet 4 de CFD · Parties A et B",
            [("Projection (M1) et ψ–ω (M2), de Re = 100 à 5000", {"size": 16})]],
           "Présentation : deux solveurs du cours, M1 (projection) et M2 (fonction de courant–vorticité). "
           "À Re = 100 ils convergent vers la même solution, celle de Ghia, ce qui valide les deux. "
           "Jusqu'à Re = 5000, ils décrivent l'évolution de l'écoulement et donnent la loi P ∝ Re^-1/2 ; "
           "M2 se dégrade plus vite que M1 quand Re augmente, et on sait expliquer pourquoi.")

# =============================================================================================
# 2. Le problème
s = content_slide("Le problème : une cavité carrée entraînée par son couvercle")
sq = box(s, 0.95, 2.1, 4.1, 4.1, fill=WHITE, line=BORD, lw=3, shape=MSO_SHAPE.RECTANGLE)
lid = box(s, 0.95, 1.92, 4.1, 0.2, fill=RED, shape=MSO_SHAPE.RECTANGLE)
arr = s.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW, Inches(1.25), Inches(1.35), Inches(1.7), Inches(0.42))
arr.fill.solid(); arr.fill.fore_color.rgb = rgb(RED); arr.line.fill.background(); arr.shadow.inherit = False
ut = text(s, 3.1, 1.3, 2.9, 0.5, [[("U", {}), sub("t"), ("(t) = 1 − e", {}), sup("−t")]], bold=True,
          size=20, color=RED, anchor="m")
vortex = s.shapes.add_shape(MSO_SHAPE.CIRCULAR_ARROW, Inches(1.85), Inches(2.75), Inches(2.5), Inches(2.5))
vortex.fill.solid(); vortex.fill.fore_color.rgb = rgb(LILAC); vortex.line.fill.background()
vortex.shadow.inherit = False
def small_vortex(x, y, d, gap_deg):
    """Tourbillon secondaire : flèche circulaire presque fermée (≈ 290°), retournée (flipH) pour tourner dans le
    sens trigonométrique, inverse du principal. gap_deg = direction de l'ouverture à l'écran (degrés, sens
    horaire depuis +x) : on l'oriente vers l'intérieur de la cavité, symétrique d'un coin à l'autre."""
    v = s.shapes.add_shape(MSO_SHAPE.CIRCULAR_ARROW, Inches(x), Inches(y), Inches(d), Inches(d))
    v.fill.solid(); v.fill.fore_color.rgb = rgb(PINK); v.line.fill.background(); v.shadow.inherit = False
    A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"
    v._element.spPr.find(A + "xfrm").set("flipH", "1")
    gc = (180 - gap_deg) % 360                      # centre de l'ouverture avant retournement
    end, start = (gc - 35) % 360, (gc + 35) % 360   # arc horaire de start à end, ouverture de 70°
    av = v._element.spPr.find(A + "prstGeom").find(A + "avLst")
    for g in list(av):
        av.remove(g)
    for name, val in (("adj1", 16000), ("adj2", 1142319), ("adj3", int(end * 60000)),
                      ("adj4", int(start * 60000)), ("adj5", 18000)):
        etree.SubElement(av, A + "gd", name=name, fmla=f"val {val}")
    return v


c1 = small_vortex(1.08, 5.42, 0.68, 315)   # coin bas gauche : ouverture vers le haut à droite
c2 = small_vortex(4.17, 5.25, 0.8, 225)    # coin bas droit : ouverture vers le haut à gauche
walls = text(s, 1.8, 5.75, 2.35, 0.4, "parois fixes : u = v = 0", size=12, color=GREY, align="c")
dim = text(s, 0.95, 6.25, 4.1, 0.4, "L = 1", size=14, color=GREY, align="c")
cards = []
cards.append(card(s, 6.0, 1.45, 6.75, 1.45, "Navier–Stokes incompressible 2D",
                  ["Adimensionné : L = U = ρ = 1, ν = 1/Re", "Fluide au repos à t = 0, couvercle en rampe"],
                  badge_label="NS", badge_fill=BORD))
cards.append(card(s, 6.0, 3.1, 6.75, 1.45, "K(t) : énergie cinétique",
                  [[("K(t) = ½ ∬ (u² + v²) dS", {"size": 18, "color": DARK})],
                   "intensité de l'écoulement dans la cavité"], badge_label="K", badge_fill=LILAC))
cards.append(card(s, 6.0, 4.75, 6.75, 1.45, "P(t) : puissance injectée par le couvercle",
                  [[("P(t) = ∫₀¹ ν ∂u/∂y|", {"size": 18}), sub("y=1", size=18), (" U", {"size": 18}),
                    sub("t", size=18), (" dx", {"size": 18})],
                   "travail du couvercle qui entretient le mouvement"], badge_label="P", badge_fill=TEAL))
g0 = [(sq.shape_id, "fade", 0, 500), (lid.shape_id, "wipe_l", 300, 600), (arr.shape_id, "fly_l", 600, 700),
      (ut.shape_id, "fade", 900, 500), (vortex.shape_id, "zoom", 1000, 700),
      (vortex.shape_id, "spin", 1700, 1600), (c1.shape_id, "spin_in", 1900, 1300, {"angle": 360}),   # rotation 360 → 0 : sens trigo
      (c2.shape_id, "spin_in", 2050, 1300, {"angle": 360}),
      (walls.shape_id, "fade", 2100, 400), (dim.shape_id, "fade", 2100, 400)]
g1 = [(sid, "rise", 250 * i + 60 * j, 500) for i, cd in enumerate(cards) for j, sid in enumerate(ids(*cd))]
anim.apply(s, [{"trigger": "auto", "items": g0}, {"trigger": "click", "items": g1}])
notes(s, "Cavité carrée, couvercle qui démarre en rampe U_t = 1 − e^-t. Adimensionnement L = U = ρ = 1, ν = 1/Re. "
         "Les deux diagnostics de l'énoncé : K, l'énergie cinétique, et P, la puissance que le couvercle injecte. "
         "[clic] cartes.")

# =============================================================================================
# 3. Les deux méthodes
s = content_slide("Deux méthodes du cours sur une base commune")
base_lbl = text(s, 0.6, 1.25, 5, 0.35, "Base commune (mêmes briques pour M1 et M2)", size=14, color=GREY,
                bold=True)
chips = ["Grille colocalisée\n(N+1)² nœuds", "Euler explicite\n(cours §3.1)", "Schémas centrés\nordre 2",
         "Poisson : gradient\nconjugué (§4.7.3)", "Δt : limites de\nstabilité (§5.3)"]
chip_shapes = []
for i, c in enumerate(chips):
    chip_shapes.append(box(s, 0.6 + i * 2.46, 1.65, 2.32, 0.85, fill=TINT, paras=c.split("\n"), size=14,
                           color=BORD, bold=False))


def lane(y, color, label, head, steps):
    b = badge(s, 0.6, y + 0.45, 0.85, label, fill=color, size=18)
    h = text(s, 1.65, y, 11, 0.4, head, size=16, color=color, bold=True)
    n = len(steps)
    gap = 0.42
    w = (11.05 - gap * (n - 1)) / n
    out = []
    for i, (t1, t2) in enumerate(steps):
        x = 1.65 + i * (w + gap)
        bx = box(s, x, y + 0.42, w, 0.95, fill=WHITE, line=color, lw=2,
                 paras=[[(t1, {"bold": True, "size": 14, "color": DARK})], [(t2, {"size": 13, "color": GREY})]])
        out.append(bx)
        if i < n - 1:
            out.append(arrow(s, x + w + 0.1, y + 0.72, 0.22, 0.35, color=color))
    return b, h, out


m1b, m1h, m1s = lane(2.75, M1C, "M1", "Projection : u, v, p   (cours §7.3.2, workshop 15)",
                     [("Prédiction sans pression", "u* = uⁿ + Δt F(uⁿ)"),
                      ("Poisson pour p (Neumann)", "Δp = ∇·u* / Δt"),
                      ("Correction", "uⁿ⁺¹ = u* − Δt ∇p")])
m2b, m2h, m2s = lane(4.25, M2C, "M2", "Fonction de courant – vorticité ψ–ω   (cours §7.3.4)",
                     [("Transport de ω", "ωⁿ⁺¹ = ωⁿ + Δt G(ωⁿ)"),
                      ("Poisson pour ψ", "Δψ = −ω (Dirichlet)"),
                      ("Vitesses", "u = ∂ψ/∂y, v = −∂ψ/∂x"),
                      ("ω de paroi", "formule de Thom")])
msg = text(s, 0.6, 5.85, 12.1, 0.5,
           [[("Même grille, mêmes schémas, même Poisson ⇒ ", {}),
             ("leurs écarts viennent de la méthode elle-même.", {"bold": True})]],
           size=17, color=BORD, align="c", anchor="m")
anim.apply(s, [
    {"trigger": "click", "items": [(base_lbl.shape_id, "fade", 0, 400)] +
        [(c.shape_id, "rise", 120 * i, 450) for i, c in enumerate(chip_shapes)]},
    {"trigger": "click", "items": [(m1b.shape_id, "zoom", 0, 450), (m1h.shape_id, "fade", 150, 400)] +
        [(x.shape_id, "wipe_l", 300 + 220 * i, 400) for i, x in enumerate(m1s)]},
    {"trigger": "click", "items": [(m2b.shape_id, "zoom", 0, 450), (m2h.shape_id, "fade", 150, 400)] +
        [(x.shape_id, "wipe_l", 300 + 200 * i, 400) for i, x in enumerate(m2s)]},
    {"trigger": "click", "items": [(msg.shape_id, "rise", 0, 600)]},
])
notes(s, "Deux méthodes du cours qui traitent l'incompressibilité de deux façons : M1 avec une pression, "
         "M2 en l'éliminant par la fonction de courant. Elles partagent tout le reste (grille colocalisée, Euler "
         "explicite, schémas centrés d'ordre 2, gradient conjugué sans matrice, départ à chaud, tolérance 1e-8, Δt "
         "des limites de stabilité). Leurs écarts viennent donc de la méthode elle-même : on peut valider l'une par "
         "l'autre, puis comparer un choix « sage » et un choix « malheureux ».")

# =============================================================================================
# 4. Section A
section_slide("Partie A", "Re = 100 : vérifier chaque brique, puis valider la solution",
              ["Poisson vérifié", "Écoulement stationnaire", "Validation croisée", "Convergence en maillage"],
              "Partie A : Re = 100. On vérifie d'abord la brique la plus délicate, puis on valide la solution.")

# =============================================================================================
# 5. Poisson
s = content_slide("A · Le solveur de Poisson est vérifié avant tout")
fig = pic(s, IMG / "poisson_0.png", 0.45, 1.3, h=4.85)
cap = text(s, 0.45, 6.2, 6.3, 0.4, "Δφ = −2π²φ, solutions connues sin·sin (Dirichlet) et cos·cos (Neumann)",
           size=12, color=GREY)
st1 = text(s, 7.1, 1.3, 5.6, 0.9, [[("ordre 2", {"size": 44, "bold": True, "color": RED})]], anchor="b")
st2 = text(s, 7.1, 2.2, 5.6, 1.0, ["erreur ∝ h² en Dirichlet et en Neumann",
                                   "itérations du gradient conjugué ∝ N : 33 → 68 → 134"], size=15, color=DARK)
pb = card(s, 7.1, 3.35, 5.65, 2.95, "Problème rencontré : Neumann d'ordre 1",
          ["« bord = voisin » φ₀ = φ₁ place la paroi en h/2 → ordre 1 seulement",
           [("Correction : nœuds fantômes miroirs φ₋₁ = φ₁ → ", {}),
            ("ordre 2, erreur ÷ 200 à N = 64", {"bold": True, "color": RED})]],
          badge_label="!", badge_fill=RED, size=14)
# mini-schéma du nœud fantôme
wall = line(s, 9.6, 5.3, 9.6, 5.82, color=BORD, w=3)
gh = box(s, 8.5, 5.47, 0.3, 0.3, fill=WHITE, line=LILAC, lw=2, shape=MSO_SHAPE.OVAL, dash=True)
n0 = box(s, 9.45, 5.47, 0.3, 0.3, fill=BORD, shape=MSO_SHAPE.OVAL)
n1 = box(s, 10.4, 5.47, 0.3, 0.3, fill=BORD, shape=MSO_SHAPE.OVAL)
lb_ = [text(s, xc - 0.5, 5.8, 1.0, 0.3, t_, size=11, color=GREY, align="c")
       for xc, t_ in ((8.65, "fantôme"), (9.6, "paroi"), (10.55, "intérieur"))]
mirror = text(s, 10.95, 5.42, 1.7, 0.4, [[("φ", {}), sub("−1"), (" = φ", {}), sub("1")]], size=16, color=LILAC,
              bold=True)
anim.apply(s, [
    {"trigger": "auto", "items": [(fig.shape_id, "fade", 0, 700, {"bld": "pic"}), (cap.shape_id, "fade", 300, 500)]},
    {"trigger": "click", "items": [(st1.shape_id, "zoom", 0, 500), (st2.shape_id, "rise", 300, 500)]},
    {"trigger": "click", "items": [(x, "rise", 80 * i, 500) for i, x in enumerate(ids(*pb))] +
        [(wall.shape_id, "fade", 600, 300), (n0.shape_id, "zoom", 700, 300), (n1.shape_id, "zoom", 800, 300),
         (gh.shape_id, "zoom", 1100, 400), (mirror.shape_id, "wipe_l", 1400, 500)] +
        [(x.shape_id, "fade", 900, 400) for x in lb_]},
])
notes(s, "Avant Navier–Stokes, on vérifie le Poisson sur des solutions connues, avec les deux conditions utilisées : "
         "Dirichlet pour ψ, Neumann pour la pression. Ordre 2 exactement, et un nombre d'itérations qui double quand N "
         "double, comme le prévoit le cours. Neumann : compatibilité (on retire la moyenne du second membre), unicité "
         "(moyenne nulle). Problème : la première version « bord = voisin » n'était que d'ordre 1 car elle place la paroi "
         "en h/2. Nœuds fantômes miroirs, pondérés par les poids des trapèzes pour garder la symétrie (exigée par le "
         "gradient conjugué) : ordre 2, erreur 200 fois plus petite à N = 64. On l'a trouvé en comparant M1 et M2.")

# =============================================================================================
# 6. Écoulement Re = 100
s = content_slide("A · L'écoulement stationnaire à Re = 100")
hh = 2.95
xs = 0.45
panels = []
for k in range(4):
    im = Image.open(IMG / f"flowA_{k}.png")
    w = hh * im.width / im.height
    panels.append(pic(s, IMG / f"flowA_{k}.png", xs, 1.3, h=hh))
    xs += w + 0.12
cds = [card(s, 0.45 + i * 4.18, 4.6, 3.98, 1.6, h_, b_, badge_label=bl, badge_fill=bf)
       for i, (h_, b_, bl, bf) in enumerate([
           ("Un tourbillon principal", "décalé vers le coin haut droit, plus deux petits tourbillons de coin",
            "1", BORD),
           ("Stationnaire vers t ≈ 8", "K et P suivent la rampe du couvercle : simple amortissement visqueux",
            "2", LILAC),
           ("Arrêt du calcul", ["résidu stationnaire divisé par 1000 :",
                                [("R = (f", {}), sup("n+1"), (" − f", {}), sup("n"), (")/Δt", {})]],
            "3", TEAL)])]
anim.apply(s, [
    {"trigger": "auto", "items": [(panels[0].shape_id, "fade", 0, 600, {"bld": "pic"}),
                                  (panels[1].shape_id, "fade", 250, 600, {"bld": "pic"}),
                                  (panels[2].shape_id, "wipe_l", 700, 1300, {"bld": "pic"}),
                                  (panels[3].shape_id, "wipe_l", 1100, 1300, {"bld": "pic"})]},
    {"trigger": "click", "items": [(x, "rise", 250 * i + 60 * j, 500) for i, cd in enumerate(cds)
                                   for j, x in enumerate(ids(*cd))]},
])
notes(s, "À Re = 100 : tourbillon principal décalé vers le haut à droite, dans le sens de l'entraînement, et deux "
         "petits tourbillons de coin en bas. K et P suivent la rampe du couvercle puis se stabilisent vers t ≈ 8 : "
         "le transitoire est un simple amortissement visqueux. On arrête quand le résidu stationnaire (de combien le "
         "champ viole les équations stationnaires) est divisé par mille.")

# =============================================================================================
# 7. Validation
s = content_slide("A · Validation : M1, M2 et Ghia se superposent")
pa = pic(s, IMG / "profA_0.png", 0.4, 1.3, h=3.45)
pb_ = pic(s, IMG / "profA_1.png", 4.6, 1.3, h=3.45)
cd = CategoryChartData()
cd.categories = ["N = 32", "N = 64", "N = 128"]
cd.add_series("max|Δu| M1 − M2", (0.020, 0.009, 0.004))
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(8.95), Inches(1.3), Inches(3.85), Inches(3.45), cd)
ch = gf.chart
chart_style(ch, "Écart max M1 − M2 (u)")
ch.has_legend = False
pl = ch.plots[0]
pl.gap_width = 70
pl.has_data_labels = True
pl.data_labels.number_format = "0.000"
pl.data_labels.number_format_is_linked = False
pl.data_labels.position = XL_LABEL_POSITION.OUTSIDE_END
pl.data_labels.font.size = Pt(13)
pl.data_labels.font.bold = True
pl.data_labels.font.color.rgb = rgb(BORD)
pl.series[0].format.fill.solid()
pl.series[0].format.fill.fore_color.rgb = rgb(LILAC)
pl.series[0].points[2].format.fill.solid()
pl.series[0].points[2].format.fill.fore_color.rgb = rgb(RED)
quiet_axes(ch, "0.000")
ch.value_axis.maximum_scale = 0.024
ch.value_axis.visible = False
ch.value_axis.has_major_gridlines = False
cA = card(s, 0.4, 5.0, 6.05, 1.45, "Validation croisée",
          "Deux méthodes différentes convergent vers la même solution : celle de Ghia (écart ≈ 10⁻²)",
          badge_label="✓", badge_fill=TEAL)
cB = card(s, 6.75, 5.0, 6.05, 1.45, "Ghia a sa propre erreur",
          "M2 − Ghia cesse de baisser à N = 128 (0,004) alors que M2 converge encore",
          badge_label="!", badge_fill=BORD)
anim.apply(s, [
    {"trigger": "auto", "items": [(pa.shape_id, "wipe_l", 0, 900, {"bld": "pic"}),
                                  (pb_.shape_id, "wipe_l", 400, 900, {"bld": "pic"})]},
    {"trigger": "click", "items": [(gf.shape_id, "wipe_u", 0, 800, {"bld": "chart"})]},
    {"trigger": "click", "items": [(x, "rise", 70 * j, 500) for j, x in enumerate(ids(*cA))]},
    {"trigger": "click", "items": [(x, "rise", 70 * j, 500) for j, x in enumerate(ids(*cB))]},
])
notes(s, "L'énoncé demande de valider par la littérature ou par un autre solveur : on fait les deux. Nos deux méthodes "
         "et Ghia se superposent, écart max ≈ 1e-2. Surtout, l'écart M1 − M2 diminue quand on raffine, jusqu'à 0.004 : "
         "deux méthodes différentes convergent vers la même solution, celle de Ghia. L'écart de M2 à Ghia cesse de "
         "diminuer à N = 128 alors que M2 converge : on atteint l'erreur propre de Ghia (129², ordre 2).")

# =============================================================================================
# 8. Convergence en maillage
s = content_slide("A · Convergence en maillage : N = 128 suffit")
gfig = pic(s, IMG / "gridA_0.png", 0.45, 1.3, h=4.6)
gcap = text(s, 0.45, 5.95, 5.9, 0.5, "Écart max entre deux grilles successives (N et 2N) : ne dépend d'aucune référence",
            size=12, color=GREY)
k2 = card(s, 6.75, 1.3, 6.0, 1.15, "M2 : ordre 2", "à ≈ 10⁻³ de sa limite à N = 128",
          badge_label="M2", badge_fill=M2C)
k1 = card(s, 6.75, 2.6, 6.0, 1.15, "M1 : ordre ≈ 1 seulement ?", "à Δt fixé : ordre 1,91 en espace",
          badge_label="M1", badge_fill=M1C)
expl = box(s, 6.75, 3.9, 6.0, 1.45, fill=WHITE, line=BORD, lw=1.5, align="l", anchor="m",
           paras=[[("Au stationnaire, la projection du cours vérifie", {"size": 14, "color": GREY})],
                  [("∇ₕ·u = Δt (Δₕ − D G̃) p", {"size": 19, "bold": True, "color": BORD})],
                  [("⇒ la solution de M1 dépend de Δt, qui change avec N", {"size": 14, "color": DARK})]])
concl = box(s, 6.75, 5.55, 6.0, 0.8, fill=RED, paras=[[("N = 128 : M1 et M2 à 0,004 l'un de l'autre", {})]],
            size=17, color=WHITE, bold=True)
anim.apply(s, [
    {"trigger": "auto", "items": [(gfig.shape_id, "fade", 0, 700, {"bld": "pic"}), (gcap.shape_id, "fade", 300, 500)]},
    {"trigger": "click", "items": [(x, "rise", 70 * j, 500) for j, x in enumerate(ids(*k2))]},
    {"trigger": "click", "items": [(x, "rise", 70 * j, 500) for j, x in enumerate(ids(*k1))] +
        [(expl.shape_id, "wipe_d", 500, 700)]},
    {"trigger": "click", "items": [(concl.shape_id, "zoom", 0, 500)]},
])
notes(s, "Pour justifier la grille sans dépendre de Ghia, on mesure l'écart entre deux grilles successives. M2 converge "
         "à l'ordre 2 : à N = 128 il est à ~1e-3 de sa limite. M1 ne montre que l'ordre 1 ici. Ce n'est pas un défaut en "
         "espace : à Δt fixé, M1 est d'ordre 1.91. Au stationnaire, D(u) = Δt (Δ_h − D G̃) p : le laplacien résolu par "
         "le Poisson et celui qu'applique la correction diffèrent sur grille colocalisée, donc la solution de M1 dépend "
         "de Δt, qui change avec N (CFL puis Fourier). M2, sans projection, n'en dépend pas. N = 128 suffit.")

# =============================================================================================
# 9. P et coins
s = content_slide("A · Problème : P ne converge pas en maillage")
csch = pic(s, IMG / "corner_schema.png", 0.4, 1.25, h=4.0)
plot = pic(s, IMG / "PlnN_0.png", 6.6, 1.3, h=3.9)
ban = box(s, 0.45, 5.4, 12.4, 1.1, fill=TINT, align="l", anchor="m", margin=0.3,
          paras=[[("Coins du couvercle : u saute de U à 0 ⇒ ∫ ν ∂u/∂y dx diverge en ln. ", {"color": DARK}),
                  ("En discret : P ≈ P", {"bold": True, "color": RED}), sub("régulier", bold=True, color=RED),
                  (" + c ν ln N", {"bold": True, "color": RED}), (" (c ≈ 4,2)", {})],
                 [("⇒ on ne compare les P qu'à N fixé ; K, lui, converge.", {"bold": True, "color": BORD})]],
          size=16)
anim.apply(s, [
    {"trigger": "auto", "items": [(csch.shape_id, "zoom", 0, 800, {"bld": "pic"})]},
    {"trigger": "click", "items": [(plot.shape_id, "wipe_l", 0, 900, {"bld": "pic"})]},
    {"trigger": "click", "items": [(ban.shape_id, "rise", 0, 600)]},
])
notes(s, "Point important pour la suite : K converge avec la grille, P non. Aux deux coins du couvercle la vitesse "
         "saute de 1 à 0 : ∂u/∂y ~ 1/r et l'intégrale diverge comme ln en continu. La grille coupe à r ≈ h, d'où "
         "P ≈ P_régulier + c ν ln N. Mesuré (M1) : P = 0.154, 0.182, 0.212 pour N = 32, 64, 128. On ne compare donc les P "
         "qu'à grille et méthode fixées.")

# =============================================================================================
# 10. Section B
section_slide("Partie B", "Re de 100 à 5000 : balayer, comparer, comprendre",
              ["Balayage et critère d'arrêt", "Évolution de l'écoulement", "Littérature et défaut de M2",
               "Loi d'échelle P(Re)"],
              "Partie B : balayage en Reynolds jusqu'à 5000.")

# =============================================================================================
# 11. Balayage
s = content_slide("B · Un balayage en Reynolds sur le calculateur Ruche")
stats = [("11", "valeurs de Re,\nde 100 à 5000"), ("48", "calculs : M1 et M2,\nN = 64, 128 (et 256)"),
         ("≈ 10⁶", "pas de temps à Re = 5000\n(Δt ≤ 2ν, schéma centré)"), ("22 h", "pour le cas le plus long\n(M2, N = 256)")]
st_shapes = []
for i, (num, lab) in enumerate(stats):
    x = 0.5 + i * 3.1
    a_ = text(s, x, 1.4, 2.9, 1.0, [[(num, {"size": 48, "bold": True, "color": RED})]], align="c", anchor="b")
    b_ = text(s, x, 2.45, 2.9, 0.9, lab.split("\n"), size=15, color=DARK, align="c", spacing=0)
    st_shapes.append((a_, b_))
import math
RES = [100, 200, 400, 700, 1000, 1500, 2000, 2500, 3200, 4000, 5000]
GHIA = {100, 400, 1000, 3200, 5000}
x0, x1 = 1.0, 12.2


def xpos(re):
    return x0 + (math.log10(re) - 2) / (math.log10(5000) - 2) * (x1 - x0)


axis = line(s, x0 - 0.3, 4.75, x1 + 0.3, 4.75, color=BORD, w=2.5)
dots = []
for re in RES:
    xc = xpos(re)
    d = box(s, xc - 0.16, 4.59, 0.32, 0.32, fill=RED if re in GHIA else LILAC, shape=MSO_SHAPE.OVAL)
    lbl = text(s, xc - 0.45, 5.0, 0.9, 0.35, str(re), size=13, color=DARK, align="c")
    dots.append((d, lbl))
br = box(s, xpos(3200) - 0.3, 3.9, xpos(5000) - xpos(3200) + 0.6, 0.45, fill=TINT,
         paras=["+ N = 256"], size=13, color=BORD, bold=True)
leg = text(s, 1.0, 5.55, 11.2, 0.45,
           [[("●", {"color": RED}), (" valeurs publiées par Ghia et al. (1982)     ", {}),
             ("●", {"color": LILAC}), (" valeurs intermédiaires     échelle logarithmique", {})]],
           size=13, color=GREY, align="c")
anim.apply(s, [
    {"trigger": "auto", "items": [(x, k, 250 * i + 150 * j, 500) for i, pr in enumerate(st_shapes)
                                  for j, (x, k) in enumerate(zip(ids(*pr), ("zoom", "fade")))]},
    {"trigger": "click", "items": [(axis.shape_id, "wipe_l", 0, 800)] +
        [(x, "zoom" if j == 0 else "fade", 300 + 110 * i, 350) for i, pr in enumerate(dots)
         for j, x in enumerate(ids(*pr))] + [(br.shape_id, "rise", 1700, 500), (leg.shape_id, "fade", 1700, 500)]},
])
notes(s, "11 valeurs de Re de 100 à 5000 (celles de Ghia plus des valeurs intermédiaires), M1 et M2, N = 64 et 128, "
         "plus N = 256 à Re = 3200 et 5000 : 48 calculs sur Ruche. À haut Re, Δt ≤ 2ν (schéma centré + Euler, "
         "C² ≤ 2 Fo) ne dépend pas de N : environ 1e6 pas à Re = 5000, jusqu'à 22 h de calcul pour un cas.")

# =============================================================================================
# 12. Critère d'arrêt
s = content_slide("B · Problème : le critère d'arrêt de la Partie A s'arrête trop tôt")
cd = CategoryChartData()
cd.categories = ["Partie A\n‖R(10)‖ ÷ 10²", "‖R(10)‖ ÷ 10³", "retenu\n‖R(30)‖ ÷ 10³", "‖R(10)‖ ÷ 10⁴"]
cd.add_series("erreur d'arrêt", (245, 40, 4.6, 4.0))
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(0.4), Inches(1.3), Inches(6.9), Inches(4.7), cd)
ch = gf.chart
chart_style(ch, "Erreur d'arrêt, en % de l'écart à Ghia")
ccap = text(s, 0.4, 6.05, 6.9, 0.4, "M1, N = 128, Re = 5000 · arrêt quand ‖R‖ est divisé par 10², 10³ ou 10⁴ "
            "depuis t = 10 ou 30", size=12, color=GREY, align="c")
ch.has_legend = False
pl = ch.plots[0]
pl.gap_width = 60
pl.has_data_labels = True
pl.data_labels.number_format = '0.0" %"'
pl.data_labels.number_format_is_linked = False
pl.data_labels.position = XL_LABEL_POSITION.OUTSIDE_END
pl.data_labels.font.size = Pt(14)
pl.data_labels.font.bold = True
pl.data_labels.font.color.rgb = rgb(BORD)
ser = pl.series[0]
ser.format.fill.solid()
ser.format.fill.fore_color.rgb = rgb(LILAC)
for i, col in ((0, RED), (2, TEAL)):
    ser.points[i].format.fill.solid()
    ser.points[i].format.fill.fore_color.rgb = rgb(col)
quiet_axes(ch, '0" %"')
ch.value_axis.maximum_scale = 280
ch.value_axis.minimum_scale = 0
ch.value_axis.major_unit = 50
p1 = card(s, 7.65, 1.3, 5.1, 1.3, "Constat", "à Re = 5000, M1 s'arrêtait à t = 98 au lieu de ≈ 430",
          badge_label="!", badge_fill=RED)
p2 = card(s, 7.65, 2.75, 5.1, 1.55, "Calibration",
          "calcul long jusqu'à t = 1000 comme référence ; règle : erreur d'arrêt ≤ 10 % de l'écart à Ghia",
          badge_label="?", badge_fill=LILAC)
dec = box(s, 7.65, 4.45, 5.1, 0.75, fill=TEAL, paras=[[("Retenu partout : ‖R‖ < 10", {}), sup("−3"), (" ‖R(t = 30)‖", {})]],
          size=16, color=WHITE, bold=True)
p3 = card(s, 7.65, 5.35, 5.1, 1.05, "Résidu de M2 mesuré sur ω",
          "sur u, v : 0 itération de CG ⇒ 13 calculs sur 23 arrêtés à tort", badge_label="M2", badge_fill=M2C,
          size=13, head_size=14)
anim.apply(s, [
    {"trigger": "auto", "items": [(gf.shape_id, "wipe_u", 0, 1000, {"bld": "chart"}),
                                  (ccap.shape_id, "fade", 600, 400)]},
    {"trigger": "click", "items": [(x, "rise", 70 * j, 500) for j, x in enumerate(ids(*p1))]},
    {"trigger": "click", "items": [(x, "rise", 70 * j, 500) for j, x in enumerate(ids(*p2))] +
        [(dec.shape_id, "zoom", 600, 500)]},
    {"trigger": "click", "items": [(x, "rise", 70 * j, 500) for j, x in enumerate(ids(*p3))]},
])
notes(s, "À haut Re, l'approche du stationnaire est beaucoup plus lente : le critère de la Partie A arrêtait le calcul "
         "quatre fois trop tôt à Re = 5000 (t = 98 au lieu de ≈ 430, erreur 245 % de l'écart à Ghia). Recalibré au cas "
         "le plus défavorable, par un calcul long jusqu'à t = 1000 qui sert de référence : résidu divisé par mille depuis "
         "t = 30, erreur d'arrêt < 5 % de l'écart à Ghia. Même critère pour tout le projet. Pour M2, le résidu est mesuré "
         "sur ω, la variable qu'il fait évoluer : sur u, v, le CG parti du ψ précédent faisait parfois 0 itération.")

# =============================================================================================
# 13. Évolution de l'écoulement
s = content_slide("B · Quand Re augmente, l'écoulement se réorganise")
row_lbl1 = text(s, 0.3, 2.0, 0.5, 0.6, "M1", size=15, color=M1C, bold=True)
row_lbl2 = text(s, 0.3, 4.0, 0.5, 0.6, "M2", size=15, color=M2C, bold=True)
cols_ = []
for k in range(6):
    x = 0.85 + k * 2.05
    a_ = pic(s, IMG / f"strM1_{k}.png", x, 1.3, w=1.95)
    b_ = pic(s, IMG / f"strM2_{k}.png", x, 3.28, w=1.95)
    cols_.append((a_, b_))
chips = [("Le centre descend", "(0,62 ; 0,73) → (0,52 ; 0,54)"),
         ("Tourbillon haut-gauche", "apparaît entre Re = 1000 et 1500"),
         ("M1 et M2 : même topologie", "centres à 1–2 mailles de Ghia")]
chip_sh = [card(s, 0.45 + i * 4.18, 5.45, 3.98, 1.0, h_, b_, size=13, head_size=15) for i, (h_, b_) in enumerate(chips)]
anim.apply(s, [
    {"trigger": "auto", "items": [(row_lbl1.shape_id, "fade", 0, 300), (row_lbl2.shape_id, "fade", 0, 300)] +
        [(x, "fade", 350 * k, 500, {"bld": "pic"}) for k, pr in enumerate(cols_) for x in ids(*pr)]},
    {"trigger": "click", "items": [(x, "rise", 250 * i + 60 * j, 500) for i, cd_ in enumerate(chip_sh)
                                   for j, x in enumerate(ids(*cd_))]},
])
notes(s, "N = 128 ; en noir le tourbillon principal, en rouge les secondaires. Le tourbillon principal descend vers le "
         "centre, de (0.62, 0.73) à Re = 100 à (0.52, 0.54) à Re = 5000, et s'arrondit. Les deux tourbillons de coin du "
         "bas grossissent, un troisième apparaît en haut à gauche entre Re = 1000 et 1500, ce qui précise la fourchette "
         "de Ghia (entre 1000 et 3200). Même topologie pour les deux méthodes ; centres à 1–2 mailles de ceux de Ghia.")

# =============================================================================================
# 14. Comparaison à la littérature
s = content_slide("B · Comparaison à la littérature : Ghia et Erturk")
pp = []
for i, k in enumerate((2, 3, 4)):
    pp.append(pic(s, IMG / f"profB_{k}.png", 0.35 + i * 2.72, 1.25, h=4.95))
pcap = text(s, 0.35, 6.25, 8.2, 0.4, "Profils médians u(0,5 ; y) et v(x ; 0,5), N = 128 (N = 256 à Re = 5000) · "
            "○ Ghia (1982) · △ Erturk (2005)", size=11, color=GREY)
q1 = card(s, 8.75, 1.3, 4.05, 2.45, "M2 se dégrade plus vite",
          [[("3×", {"size": 36, "bold": True, "color": RED}), ("  plus loin de Ghia", {"size": 15})],
           "que M1 à Re = 3200 (0,060 contre 0,020)"], badge_label="M2", badge_fill=M2C)
q2 = card(s, 8.75, 3.95, 4.05, 2.45, "Ghia n'est pas exact",
          [[("2,7 %", {"size": 36, "bold": True, "color": RED}), ("  sur ψ", {"size": 15}), sub("min", size=15)],
           "à Re = 5000 (Erturk, 601² + Richardson) ; M1 − Erturk : 0,005"], badge_label="!", badge_fill=BORD)
anim.apply(s, [
    {"trigger": "auto", "items": [(p.shape_id, "wipe_d", 300 * i, 800, {"bld": "pic"}) for i, p in enumerate(pp)] +
        [(pcap.shape_id, "fade", 900, 400)]},
    {"trigger": "click", "items": [(x, "rise", 70 * j, 500) for j, x in enumerate(ids(*q1))]},
    {"trigger": "click", "items": [(x, "rise", 70 * j, 500) for j, x in enumerate(ids(*q2))]},
])
notes(s, "On compare sur la grille de Ghia : 129² jusqu'à Re = 3200, 257² à Re = 5000. Les deux méthodes reproduisent "
         "Ghia à 1e-2 près. Deux messages. M2 se dégrade plus vite que M1 : à Re = 3200 il est trois fois plus loin "
         "(0.060 contre 0.020) et il lui faut la grille fine pour le rejoindre. Ensuite, Ghia n'est pas exact : Erturk "
         "et al. (601², Richardson) montrent que Ghia sous-estime le tourbillon de 2.7 % à Re = 5000 ; par rapport à "
         "Erturk, M1 est à 0.005 près.")

# =============================================================================================
# 15. Thom
s = content_slide("B · Pourquoi M2 décroche : la formule de Thom")
sch = pic(s, IMG / "thom_schema.png", 0.35, 1.2, h=3.95)
form = box(s, 0.45, 5.25, 5.6, 1.2, fill=TINT, anchor="m",
           paras=[[("ω", {}), ("Thom", {"size": 12}), (" = ω", {}), ("w", {"size": 12}),
                   (" + (h/3) ∂ₙω|", {}), ("w", {"size": 12}), (" + O(h²)", {})],
                  [("erreur d'ordre 1, qui sous-estime toujours |ω", {"size": 13, "color": GREY}),
                   sub("w", size=13, color=GREY), ("|", {"size": 13, "color": GREY})]],
           size=19, color=BORD, bold=True)
cd = CategoryChartData()
cd.categories = ["Re = 100", "Re = 1000", "Re = 3200", "Re = 5000"]
cd.add_series("erreur", (2.3, 7.3, 16.6, 22.1))
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(6.5), Inches(1.3), Inches(6.3), Inches(3.3), cd)
ch = gf.chart
chart_style(ch, [("Sous-estimation de |ω", {}), sub("w"), ("| par Thom", {})])
tcap = text(s, 6.5, 4.62, 6.3, 0.3, "formule appliquée aux profils de M1 (qui n'utilise pas Thom), h = 1/64",
            size=11, color=GREY, align="c")
ch.has_legend = False
pl = ch.plots[0]
pl.gap_width = 60
pl.has_data_labels = True
pl.data_labels.number_format = '"−"0.0" %"'
pl.data_labels.number_format_is_linked = False
pl.data_labels.position = XL_LABEL_POSITION.OUTSIDE_END
pl.data_labels.font.size = Pt(14)
pl.data_labels.font.bold = True
pl.data_labels.font.color.rgb = rgb(BORD)
pl.series[0].format.fill.solid()
pl.series[0].format.fill.fore_color.rgb = rgb(M2C)
quiet_axes(ch, '0" %"')
ch.value_axis.maximum_scale = 26
ch.value_axis.visible = False
ch.value_axis.has_major_gridlines = False
why = card(s, 6.5, 5.0, 6.3, 1.45, "Couche limite δ ∝ 1/√Re ⇒ erreur ∝ h √Re",
           "sous le couvercle, M2 impose une vitesse trop faible : u = 0,41 contre 0,45 (M1) et 0,46 (Ghia) "
           "à Re = 3200, N = 128", badge_label="√", badge_fill=RED, size=13)
anim.apply(s, [
    {"trigger": "auto", "items": [(sch.shape_id, "fade", 0, 800, {"bld": "pic"})]},
    {"trigger": "click", "items": [(form.shape_id, "rise", 0, 600)]},
    {"trigger": "click", "items": [(gf.shape_id, "wipe_u", 0, 900, {"bld": "chart"}),
                                   (tcap.shape_id, "fade", 500, 400)]},
    {"trigger": "click", "items": [(x, "rise", 70 * j, 500) for j, x in enumerate(ids(*why))]},
])
notes(s, "La seule différence de fond est au couvercle : M1 impose u = U, M2 reconstruit ω de paroi par Thom. Taylor : "
         "ω_Thom = ω_w + h/3 ∂nω + O(h²). Erreur d'ordre 1 qui sous-estime toujours |ω_w| (|ω| décroît en s'éloignant de "
         "la paroi), d'une fraction h/(3δ) avec δ ~ Re^-1/2 : proportionnelle à h√Re. Mesurée sur les profils de M1 : "
         "−2 % à Re = 100, −22 % à Re = 5000. Thom s'écrit ū1 = U − h/2 |ω_w| : la moyenne d'un profil linéaire ; le vrai "
         "profil est convexe (ν ∂²u/∂y² = ∂p/∂x > 0), donc Thom impose une vitesse trop faible. L'ordre 2 global de "
         "Thom suppose h ≪ δ ; ici h/δ ≈ 0.3.")

# =============================================================================================
# 16. Fausse solution
s = content_slide("B · Cas extrême : une fausse solution stationnaire (M2, N = 64)")
cd = CategoryChartData()
cd.categories = ["1000", "2500", "3200", "5000"]
cd.add_series("M2, N = 128", (0.1155, 0.1121, 0.1097, 0.1026))
cd.add_series("M2, N = 64", (0.1052, 0.0879, 0.0306, 0.0212))
gf = s.shapes.add_chart(XL_CHART_TYPE.LINE_MARKERS, Inches(0.4), Inches(1.3), Inches(6.6), Inches(5.1), cd)
ch = gf.chart
chart_style(ch, [("Intensité du tourbillon principal |ψ", {}), sub("min"), ("|", {})])
ch.has_legend = True
ch.legend.position = XL_LEGEND_POSITION.BOTTOM
ch.legend.include_in_layout = False
ch.legend.font.size = Pt(13)
for sr, col, mk in zip(ch.plots[0].series, (M2C, RED), (XL_MARKER_STYLE.CIRCLE, XL_MARKER_STYLE.SQUARE)):
    sr.format.line.color.rgb = rgb(col)
    sr.format.line.width = Pt(3)
    sr.smooth = False
    sr.marker.style = mk
    sr.marker.size = 10
    sr.marker.format.fill.solid()
    sr.marker.format.fill.fore_color.rgb = rgb(col)
    sr.marker.format.line.color.rgb = rgb(col)
quiet_axes(ch, "0.00")
ch.value_axis.minimum_scale = 0
ch.value_axis.maximum_scale = 0.14
ch.category_axis.has_title = True
ch.category_axis.axis_title.text_frame.text = "Re"
ch.category_axis.axis_title.text_frame.paragraphs[0].runs[0].font.size = Pt(12)
top = card(s, 7.35, 1.3, 5.45, 1.5, "Le tourbillon perd 75 % pour Re ≥ 3200",
           "pourtant le calcul est bien stationnaire (résidu ÷ 10⁷) ; ni M1 ni N = 128 ne la donnent",
           badge_label="M2", badge_fill=M2C, size=13)
rows = [("✗", BORD, "Diffusion numérique ?", "non : l'upwind échoue aussi, et diffuser renforcerait le tourbillon"),
        ("✗", BORD, "Oscillations du centré ?", "non : M1 en a autant et reste juste"),
        ("✓", TEAL, "Condition de paroi (Thom)", [("oui : avec Woods (ordre 2 local), ψ", {}), sub("min"),
                                            (" passe de −0,031 à −0,099", {})])]
row_sh = []
for i, (sym, col, h_, b_) in enumerate(rows):
    y = 3.0 + i * 0.92
    bd = badge(s, 7.45, y + 0.12, 0.55, sym, fill=col, size=18)
    tx = text(s, 8.15, y, 4.65, 0.85, [[(h_, {"bold": True, "color": BORD})], b_], size=13, color=DARK, spacing=1)
    row_sh.append((bd, tx))
lim = text(s, 7.35, 5.85, 5.45, 0.6, "Woods ne corrige que ~½ du défaut et diverge à Re = 5000 : on garde Thom, "
           "qui exige N ≳ 2√Re", size=12, color=GREY, italic=True)
anim.apply(s, [
    {"trigger": "auto", "items": [(gf.shape_id, "wipe_l", 0, 1300, {"bld": "chart"})]},
    {"trigger": "click", "items": [(x, "rise", 70 * j, 500) for j, x in enumerate(ids(*top))]},
    {"trigger": "click", "items": [(x, "zoom" if j == 0 else "fade", 250 * i + 100 * j, 450)
                                   for i, pr in enumerate(row_sh[:2]) for j, x in enumerate(ids(*pr))]},
    {"trigger": "click", "items": [(row_sh[2][0].shape_id, "zoom", 0, 500), (row_sh[2][1].shape_id, "fade", 100, 450),
                                   (lim.shape_id, "fade", 700, 500)]},
])
notes(s, "Pour Re ≥ 3200, M2 à N = 64 converge vers un état stationnaire où le tourbillon a perdu 75 % de son intensité "
         "(ψ_min = −0.03 au lieu de −0.12), fluide presque à l'arrêt sous le couvercle. Vraiment stationnaire (résidu "
         "divisé par 1e7). Fausse : pas convergée en maillage, absente chez M1 sur la même grille, contredit Ghia et "
         "Erturk. Ce n'est pas de la diffusion numérique (l'upwind échoue aussi ; une solution trop diffusée aurait un "
         "tourbillon intense). Test de Woods (ω_w = 3/2 ω_Thom − 1/2 ω_1, ordre 2 local, hors méthode retenue) : supprime "
         "la fausse solution, donc Thom est en cause. Mais Woods n'explique qu'environ la moitié du défaut et diverge "
         "à Re = 5000, N = 64 : on garde Thom. M2 + Thom demande N ≳ 2√Re.")

# =============================================================================================
# 17. P(Re)
s = content_slide([("B · P(Re) : la loi d'échelle en Re", {}), sup("−1/2")])
pr0 = pic(s, IMG / "PRe_0.png", 0.45, 1.25, h=3.95)
pr1 = pic(s, IMG / "PRe_1.png", 6.95, 1.25, h=3.95)
steps = [("Couche limite", [("δ ~ Re", {}), sup("−1/2")]),
         ("Puissance", [("P ~ ν ∂u/∂y ~ ν/δ ~ Re", {}), sup("−1/2")]),
         ("Pente mesurée ≈ −0,7", [("coins : terme en ν ln N ∝ Re", {}), sup("−1"),
                                   (" ; pentes locales → −0,54", {})])]
st_sh = []
for i, (h_, b_) in enumerate(steps):
    w = [3.2, 4.0, 4.4][i]
    x = [0.45, 4.1, 8.55][i]
    st_sh.append(box(s, x, 5.4, w, 1.05, fill=[TINT, TINT, RED][i], anchor="m",
                     paras=[[(h_, {"bold": True, "size": 15, "color": [BORD, BORD, WHITE][i]})],
                            [(t_, {**o_, "size": 13, "color": [DARK, DARK, WHITE][i]}) for t_, o_ in b_]]))
    if i < 2:
        st_sh.append(arrow(s, x + w + 0.09, 5.75, 0.22, 0.35, color=LILAC))
anim.apply(s, [
    {"trigger": "auto", "items": [(pr0.shape_id, "fade", 0, 700, {"bld": "pic"})]},
    {"trigger": "click", "items": [(pr1.shape_id, "wipe_l", 0, 900, {"bld": "pic"})]},
    {"trigger": "click", "items": [(x.shape_id, "wipe_l", 250 * i, 450) for i, x in enumerate(st_sh)]},
])
notes(s, "Le couvercle entraîne le fluide à travers une couche limite d'épaisseur δ ~ Re^-1/2, donc ∂u/∂y ~ 1/δ et "
         "P ~ ν/δ ~ Re^-1/2. La pente globale mesurée est plus forte (≈ −0.7) parce que les coins du couvercle ajoutent un "
         "terme en ν ln N, c'est-à-dire en Re^-1. Les pentes locales tendent vers −1/2 (−0.54 pour M2 au-delà de 2500), "
         "et le modèle à deux termes P√Re = a + c Re^-1/2 est une droite (figure de droite).")

# =============================================================================================
# 18. Bilan
s = content_slide("Bilan des Parties A et B")
cols3 = [("✓", TEAL, "Vérifié et validé", ["Poisson d'ordre 2", "M1 ≈ M2 ≈ Ghia à Re = 100", "N = 128 suffit"]),
         ("↗", LILAC, "Jusqu'à Re = 5000", ["topologie qui évolue", "tourbillon haut-gauche entre 1000 et 1500",
                                             "P ∝ 1/√Re aux grands Re"]),
         ("⚖", RED, "Sage ou malheureux", ["M1 : le plus précis à grille égale",
                                                     "M2 + Thom : exige N ≳ 2√Re"])]
col_sh = []
for i, (sym, colr, head, items) in enumerate(cols3):
    x = 0.5 + i * 4.18
    bg_ = box(s, x, 1.45, 3.95, 3.85, fill=TINT)
    bd = badge(s, x + 1.53, 1.7, 0.9, sym, fill=colr, size=26)
    hd = text(s, x + 0.2, 2.75, 3.55, 0.5, head, size=18, color=BORD, bold=True, align="c")
    bl = text(s, x + 0.3, 3.35, 3.35, 1.9, [[("• " + it, {})] for it in items], size=15, color=DARK, spacing=6)
    col_sh.append((bg_, bd, hd, bl))
nxt = box(s, 0.5, 5.55, 12.3, 0.9, fill=BORD, anchor="m",
          paras=[[("Partie C  ·  Re > Re", {"bold": True}), sub("c", bold=True), (" ≈ 8000 : ", {"bold": True}),
                  ("distinguer une instationnarité physique d'un artefact de grille", {})]],
          size=17, color=WHITE)
anim.apply(s, [
    {"trigger": "auto", "items": [(x, ["fade", "zoom", "fade", "fade"][j], 350 * i + 120 * j, 500)
                                  for i, c in enumerate(col_sh) for j, x in enumerate(ids(*c))]},
    {"trigger": "click", "items": [(nxt.shape_id, "wipe_l", 0, 800)]},
])
notes(s, "Deux solveurs du cours, vérifiés brique par brique et validés l'un par l'autre à Re = 100. Jusqu'à Re = 5000, "
         "ils décrivent l'évolution de la topologie et donnent P ∝ Re^-1/2. M1 est le plus précis à grille égale ; M2 "
         "demande une grille plus fine à haut Re, à cause de sa condition de paroi. Pour la Partie C, au-delà de "
         "Re ≈ 8000 : distinguer une instationnarité physique d'un artefact de grille (Erturk : périodique à 7500 sur "
         "129², stationnaire sur 257²).")

# =============================================================================================
# 19. Section C (à compléter)
section_slide("Partie C", "", [], "Partie C : à venir.")

# =============================================================================================
# 20. Fin
ring_slide(["Merci !"], ["Questions ?"], "Questions : voir la section « À savoir pour les questions » de "
           "Doc/oral_parties_A_B.pdf.")

SLDNUM = ('<p:sp xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" '
          'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"><p:nvSpPr><p:cNvPr id="{id}" '
          'name="Numéro de diapositive"/><p:cNvSpPr><a:spLocks noGrp="1"/></p:cNvSpPr><p:nvPr><p:ph type="sldNum" '
          'sz="quarter" idx="12"/></p:nvPr></p:nvSpPr><p:spPr/><p:txBody><a:bodyPr/><a:lstStyle/><a:p>'
          '<a:fld id="{{52338FEE-DEA1-0B41-B2EF-9F30B23FA466}}" type="slidenum"><a:rPr lang="fr-FR"/><a:t>{n}</a:t>'
          '</a:fld></a:p></p:txBody></p:sp>')
for n, sl in enumerate(prs.slides, start=1):
    if sl.slide_layout.name in ("Titre seul", "Titre de section"):
        new_id = max(sh.shape_id for sh in sl.shapes) + 1
        sl.shapes._spTree.append(etree.fromstring(SLDNUM.format(id=new_id, n=n)))

prs.save(OUT)
print("écrit :", OUT, "—", len(prs.slides), "diapositives")
