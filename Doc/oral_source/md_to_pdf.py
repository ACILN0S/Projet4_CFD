"""Convertit un document Markdown de Doc/ en PDF, avec la mise en page du poly « Numerical Methods ».

Usage (depuis Doc/oral_source/) : python md_to_pdf.py ../partie_A_limites.md "Sous-titre de la page de garde"
Sortie : le PDF à côté du .md (même nom). Le HTML intermédiaire est écrit dans ce dossier.
Dépendances : paquet Python `markdown` (pip install markdown), Google Chrome (impression en PDF),
connexion internet (KaTeX, chargé depuis cdn.jsdelivr.net, pour les formules $...$ et $$...$$).
"""
import html
import re
import subprocess
import sys
from pathlib import Path

import markdown

HERE = Path(__file__).resolve().parent
CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"

CSS = r"""
  @page { size: A4; margin: 24mm 27mm 22mm 27mm;
    @top-left { content: "__HEADER__"; font-family: "Times New Roman", serif; font-style: italic; font-size: 8.5pt; letter-spacing: .3px; }
    @top-right { content: counter(page); font-family: "Times New Roman", serif; font-size: 10pt; } }
  @page :first { @top-left { content: none; } @top-right { content: none; } }
  body { font-family: "Times New Roman", Times, serif; font-size: 11pt; line-height: 1.36; color: #000; background: #fff; margin: 0;
         text-align: justify; hyphens: auto; }
  p { margin: 0 0 6pt 0; }
  ul, ol { margin: 2pt 0 7pt 0; padding-left: 18pt; }
  li { margin: 1.5pt 0; }
  code { font-family: "Courier New", monospace; font-size: 9.5pt; }
  a { color: #000; text-decoration: none; }
  .titlepage { break-after: page; text-align: center; padding-top: 18mm; }
  .titlepage img { width: 62mm; }
  .titlepage hr { border: 0; border-top: 1.1pt solid #000; margin: 0; }
  .titlepage .t { font-size: 20pt; font-weight: bold; margin: 11pt 0; }
  .titlepage .s { margin-top: 12mm; font-size: 12pt; }
  .titlepage .auth { margin-top: 26mm; font-size: 12pt; }
  .titlepage .year { margin-top: 60mm; font-size: 12pt; }
  .chapter { padding-top: 22mm; }
  .chapter h1 { font-size: 24pt; font-weight: bold; margin: 0 0 22pt 0; text-align: left; line-height: 1.15; }
  .minitoc { margin: 0 0 16pt 0; font-size: 9.5pt; }
  .minitoc .h { font-weight: bold; font-size: 11pt; border-bottom: .6pt solid #000; padding-bottom: 2pt; margin-bottom: 3pt; }
  .minitoc ol { list-style: none; margin: 0; padding: 0 0 4pt 30pt; border-bottom: .6pt solid #000; }
  .minitoc li { font-weight: bold; margin: 1pt 0; }
  h2 { font-size: 14pt; font-weight: bold; margin: 18pt 0 7pt 0; break-after: avoid; text-align: left; }
  h3 { font-size: 11.5pt; font-weight: bold; margin: 12pt 0 5pt 0; break-after: avoid; text-align: left; }
  h2 .n, h3 .n { display: inline-block; min-width: 30pt; }
  hr { display: none; }
  table { border-collapse: collapse; margin: 8pt auto 10pt auto; font-size: 9.6pt; break-inside: avoid; }
  thead th { border-top: 1.2pt solid #000; border-bottom: .6pt solid #000; padding: 2pt 7pt; font-weight: normal; }
  tbody td { padding: 1.5pt 7pt; text-align: center; }
  tbody tr:last-child td { border-bottom: 1.2pt solid #000; }
  .katex { font-size: 1.04em; }
  .katex-display { margin: 8pt 0; }
"""


def protect_math(text):
    # Les formules $...$ / $$...$$ ne doivent pas passer par Markdown (souligné, barres obliques) : on les met de côté
    store = []

    def keep(m):
        store.append(m.group(0))
        return f"\x00MATH{len(store) - 1}\x00"
    text = re.sub(r"\$\$.+?\$\$", keep, text, flags=re.S)
    text = re.sub(r"(?<!\$)\$[^$\n]+?\$(?!\$)", keep, text)
    return text, store


def number_heading(m, level):
    # « 1. Titre » ou « 3.1 Titre » : le numéro du Markdown est gardé (les renvois « §3.2 » du texte en dépendent)
    title = m.group(1)
    num = re.match(r"^(\d+(?:\.\d+)*)\.?\s+(.*)$", title)
    if num:
        return f'<h{level}><span class="n">{num.group(1)}</span>{num.group(2)}</h{level}>'
    return f"<h{level}>{title}</h{level}>"


def build(md_path, subtitle):
    src = Path(md_path).read_text(encoding="utf-8")
    title = re.search(r"^# (.+)$", src, flags=re.M).group(1).strip()
    src = re.sub(r"^# .+$\n", "", src, count=1, flags=re.M)
    sections = re.findall(r"^## (.+)$", src, flags=re.M)
    body, store = protect_math(src)
    body = markdown.markdown(body, extensions=["tables"])
    body = re.sub(r"<h2>(.+?)</h2>", lambda m: number_heading(m, 2), body)
    body = re.sub(r"<h3>(.+?)</h3>", lambda m: number_heading(m, 3), body)
    body = re.sub("\x00MATH(\\d+)\x00", lambda m: html.escape(store[int(m.group(1))], quote=False), body)
    toc = "".join(f"<li>{html.escape(s)}</li>" for s in sections)
    header = title.split(" — ")[0].upper() + " — DOCUMENT D'ACCOMPAGNEMENT"
    page = f"""<!doctype html>
<html lang="fr"><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/contrib/auto-render.min.js"
  onload="renderMathInElement(document.body,{{delimiters:[{{left:'$$',right:'$$',display:true}},{{left:'$',right:'$',display:false}}],throwOnError:false}});"></script>
<style>{CSS.replace("__HEADER__", header)}</style></head>
<body>
<div class="titlepage"><img src="logo.png" alt="CentraleSupélec"><div style="height:16mm"></div>
<hr><div class="t">Projet 4 — Cavité entraînée à haut Reynolds</div><hr>
<div class="s">{html.escape(subtitle)}</div><div class="auth">Équipe du Projet 4</div><div class="year">2026-2027</div></div>
<div class="chapter"><h1>{html.escape(title)}</h1>
<div class="minitoc"><div class="h">Contenu</div><ol>{toc}</ol></div></div>
{body}
</body></html>"""
    html_path = HERE / (Path(md_path).stem + ".html")
    html_path.write_text(page, encoding="utf-8")
    pdf_path = Path(md_path).with_suffix(".pdf").resolve()
    subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", "--allow-file-access-from-files",
                    "--virtual-time-budget=20000", f"--print-to-pdf={pdf_path}", html_path.as_uri()],
                   check=True, capture_output=True)
    print("écrit :", pdf_path)


if __name__ == "__main__":
    build(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "Document d'accompagnement")
