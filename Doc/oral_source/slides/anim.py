"""Animations et transitions PowerPoint écrites directement en XML (p:timing, p:transition).

Une diapositive reçoit une liste de groupes : le premier peut démarrer seul à l'arrivée
sur la diapositive (trigger="auto"), les suivants démarrent au clic. Dans un groupe, chaque
effet a son propre délai (ms) : tout est en « avec la précédente » + délai.
"""
from lxml import etree

NS_P = "http://schemas.openxmlformats.org/presentationml/2006/main"
NS_A = "http://schemas.openxmlformats.org/drawingml/2006/main"
NSDECL = f'xmlns:p="{NS_P}" xmlns:a="{NS_A}"'


class _Ids:
    def __init__(self):
        self.n = 2

    def next(self):
        self.n += 1
        return self.n


def _tgt(spid):
    return f'<p:tgtEl><p:spTgt spid="{spid}"/></p:tgtEl>'


def _val(v):
    if isinstance(v, (int, float)):
        return f'<p:val><p:fltVal val="{v}"/></p:val>'
    return f'<p:val><p:strVal val="{v}"/></p:val>'


def _set_visible(ids, spid):
    return (f'<p:set><p:cBhvr><p:cTn id="{ids.next()}" dur="1" fill="hold"><p:stCondLst>'
            f'<p:cond delay="0"/></p:stCondLst></p:cTn>{_tgt(spid)}<p:attrNameLst>'
            f'<p:attrName>style.visibility</p:attrName></p:attrNameLst></p:cBhvr>'
            f'<p:to><p:strVal val="visible"/></p:to></p:set>')


def _filter(ids, spid, flt, dur):
    return (f'<p:animEffect transition="in" filter="{flt}"><p:cBhvr><p:cTn id="{ids.next()}" '
            f'dur="{dur}"/>{_tgt(spid)}</p:cBhvr></p:animEffect>')


def _anim(ids, spid, attr, v0, v1, dur, decel=0):
    dec = f' decel="{decel}"' if decel else ""
    return (f'<p:anim calcmode="lin" valueType="num"><p:cBhvr additive="base"><p:cTn id="{ids.next()}" '
            f'dur="{dur}" fill="hold"{dec}/>{_tgt(spid)}<p:attrNameLst><p:attrName>{attr}</p:attrName>'
            f'</p:attrNameLst></p:cBhvr><p:tavLst><p:tav tm="0">{_val(v0)}</p:tav>'
            f'<p:tav tm="100000">{_val(v1)}</p:tav></p:tavLst></p:anim>')


def _rot(ids, spid, by_deg, dur):
    return (f'<p:animRot by="{int(by_deg * 60000)}"><p:cBhvr><p:cTn id="{ids.next()}" dur="{dur}" '
            f'fill="hold" accel="30000" decel="30000"/>{_tgt(spid)}<p:attrNameLst><p:attrName>r</p:attrName>'
            f'</p:attrNameLst></p:cBhvr></p:animRot>')


# effet -> (presetID, presetClass, presetSubtype, comportements)
def _behaviours(ids, kind, spid, dur, opt):
    if kind == "fade":
        return 10, "entr", 0, _set_visible(ids, spid) + _filter(ids, spid, "fade", dur)
    if kind.startswith("wipe"):
        side = {"wipe_l": ("left", 8), "wipe_r": ("right", 2), "wipe_u": ("up", 4), "wipe_d": ("down", 1)}[kind]
        return 22, "entr", side[1], _set_visible(ids, spid) + _filter(ids, spid, f"wipe({side[0]})", dur)
    if kind == "rise":  # apparition + léger glissement vers le haut
        return (42, "entr", 0, _set_visible(ids, spid) + _filter(ids, spid, "fade", dur)
                + _anim(ids, spid, "ppt_y", "#ppt_y+0.06", "#ppt_y", dur, decel=100000))
    if kind == "zoom":
        return (53, "entr", 16, _set_visible(ids, spid) + _filter(ids, spid, "fade", dur)
                + _anim(ids, spid, "ppt_w", 0, "#ppt_w", dur, decel=60000)
                + _anim(ids, spid, "ppt_h", 0, "#ppt_h", dur, decel=60000))
    if kind == "spin_in":  # apparition en tournant autour de son centre jusqu'à sa position
        ang = opt.get("angle", -180)
        return (49, "entr", 0, _set_visible(ids, spid) + _filter(ids, spid, "fade", int(dur * 0.45))
                + _anim(ids, spid, "style.rotation", ang, 0, dur, decel=90000))
    if kind == "fly_l":
        return (2, "entr", 8, _set_visible(ids, spid)
                + _anim(ids, spid, "ppt_x", "0-#ppt_w/2", "#ppt_x", dur, decel=100000)
                + _anim(ids, spid, "ppt_y", "#ppt_y", "#ppt_y", dur))
    if kind == "spin":  # accentuation : un tour complet
        return 8, "emph", 0, _rot(ids, spid, opt.get("angle", 360), dur)
    raise ValueError(kind)


def timing_xml(groups):
    """groups : [{"trigger": "auto"|"click", "items": [(spid, kind, delay_ms, dur_ms, opt), ...]}]"""
    ids = _Ids()
    pars, bld = [], {}
    for gi, g in enumerate(groups):
        auto = g["trigger"] == "auto"
        assert not auto or gi == 0, "seul le premier groupe peut être automatique"
        outer_id, mid_id = ids.next(), ids.next()
        cond = ('<p:cond delay="indefinite"/><p:cond evt="onBegin" delay="0"><p:tn val="2"/></p:cond>'
                if auto else '<p:cond delay="indefinite"/>')
        effects = []
        for k, it in enumerate(g["items"]):
            spid, kind, delay, dur = it[:4]
            opt = it[4] if len(it) > 4 else {}
            eff_id = ids.next()
            pid, pcls, psub, beh = _behaviours(ids, kind, spid, dur, opt)
            node = ("afterEffect" if auto else "clickEffect") if k == 0 else "withEffect"
            effects.append(
                f'<p:par><p:cTn id="{eff_id}" presetID="{pid}" presetClass="{pcls}" presetSubtype="{psub}" '
                f'fill="hold" grpId="0" nodeType="{node}"><p:stCondLst><p:cond delay="{delay}"/></p:stCondLst>'
                f'<p:childTnLst>{beh}</p:childTnLst></p:cTn></p:par>')
            if pcls == "entr":
                bld.setdefault(spid, opt.get("bld", "sp"))
        pars.append(
            f'<p:par><p:cTn id="{outer_id}" fill="hold"><p:stCondLst>{cond}</p:stCondLst><p:childTnLst>'
            f'<p:par><p:cTn id="{mid_id}" fill="hold"><p:stCondLst><p:cond delay="0"/></p:stCondLst>'
            f'<p:childTnLst>{"".join(effects)}</p:childTnLst></p:cTn></p:par></p:childTnLst></p:cTn></p:par>')
    bld_xml = ""
    for spid, typ in bld.items():
        if typ == "sp":
            bld_xml += f'<p:bldP spid="{spid}" grpId="0" animBg="1"/>'
        elif typ == "chart":
            bld_xml += f'<p:bldGraphic spid="{spid}" grpId="0"><p:bldAsOne/></p:bldGraphic>'
    xml = (f'<p:timing {NSDECL}><p:tnLst><p:par><p:cTn id="1" dur="indefinite" restart="never" nodeType="tmRoot">'
           f'<p:childTnLst><p:seq concurrent="1" nextAc="seek"><p:cTn id="2" dur="indefinite" nodeType="mainSeq">'
           f'<p:childTnLst>{"".join(pars)}</p:childTnLst></p:cTn><p:prevCondLst><p:cond evt="onPrev" delay="0">'
           f'<p:tgtEl><p:sldTgt/></p:tgtEl></p:cond></p:prevCondLst><p:nextCondLst><p:cond evt="onNext" delay="0">'
           f'<p:tgtEl><p:sldTgt/></p:tgtEl></p:cond></p:nextCondLst></p:seq></p:childTnLst></p:cTn></p:par>'
           f'</p:tnLst>{"<p:bldLst>" + bld_xml + "</p:bldLst>" if bld_xml else ""}</p:timing>')
    return etree.fromstring(xml)


TRANSITIONS = {
    "fade": '<p:transition {ns} spd="med"><p:fade/></p:transition>',
    "push_u": '<p:transition {ns} spd="slow"><p:push dir="u"/></p:transition>',
    "wipe_r": '<p:transition {ns} spd="slow"><p:wipe dir="r"/></p:transition>',
    "circle": '<p:transition {ns} spd="slow"><p:circle/></p:transition>',
}


def apply(slide, groups=None, transition="fade"):
    """Ajoute transition puis animations, dans l'ordre exigé par le schéma (après clrMapOvr)."""
    sld = slide._element
    for tag in ("transition", "timing"):
        for old in sld.findall(f"{{{NS_P}}}{tag}"):
            sld.remove(old)
    anchor = sld.find(f"{{{NS_P}}}clrMapOvr")
    idx = list(sld).index(anchor) + 1
    if transition:
        sld.insert(idx, etree.fromstring(TRANSITIONS[transition].format(ns=NSDECL)))
        idx += 1
    if groups:
        sld.insert(idx, timing_xml(groups))
