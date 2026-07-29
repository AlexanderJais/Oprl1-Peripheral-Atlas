"""Assemble the manuscript into one rendered PDF.

Introduction, Results with Figures 1 and 2 and their legends, Methods,
References, and Table S1 as a supplement. Every section is read from the
markdown file that holds it, so the PDF has no copy of the text to drift from
the source; the only thing this script decides is order, page furniture and
where the figures fall.

Figures are placed after the Results section that first cites them, each on its
own page with its legend beneath it, which is where a reader wants them and is
close enough to a journal's own placement to read as a manuscript.

Run from the repository root: python3 paper/build_pdf.py
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPER = ROOT / "paper"
OUT = PAPER / "Oprl1_peripheral_atlas.pdf"

import markdown as md
from weasyprint import CSS, HTML

CSS_TEXT = """
@page {
  size: A4;
  margin: 20mm 18mm 18mm 18mm;
  @bottom-center {
    content: counter(page);
    font-family: "Nimbus Sans", Helvetica, Arial, sans-serif;
    font-size: 8pt; color: #555;
  }
}
@page landscape { size: A4 landscape; margin: 14mm; }

html { font-size: 10.5pt; }
body {
  font-family: "Nimbus Roman", "Times New Roman", Times, serif;
  line-height: 1.55; text-align: justify; hyphens: auto; color: #111;
}
h1, h2, h3 { font-family: "Nimbus Sans", Helvetica, Arial, sans-serif; }
h1.title {
  font-size: 17pt; line-height: 1.25; text-align: left; margin: 0 0 2mm 0;
  hyphens: none;
}
h2 {
  font-size: 12pt; margin: 7mm 0 2mm 0; text-align: left; hyphens: none;
  break-after: avoid;
}
h3 {
  font-size: 10.5pt; margin: 5mm 0 1.5mm 0; text-align: left; hyphens: none;
  break-after: avoid;
}
p { margin: 0 0 2.4mm 0; }
em { font-style: italic; }
code {
  font-family: "Nimbus Mono PS", "DejaVu Sans Mono", monospace;
  font-size: 0.88em; background: #f2f2f2; padding: 0 0.15em;
}
a { color: inherit; text-decoration: none; }

.byline {
  font-family: "Nimbus Sans", Helvetica, Arial, sans-serif;
  font-size: 9pt; color: #555; margin: 0 0 8mm 0; text-align: left;
}

/* References: hanging indent, and never justified into rivers. */
#references p {
  text-align: left; hyphens: none;
  padding-left: 7mm; text-indent: -7mm; margin-bottom: 1.8mm;
}

figure { break-before: page; break-inside: avoid; margin: 0; }
figure img { width: 100%; height: auto; }
figcaption {
  font-size: 8.6pt; line-height: 1.42; text-align: justify; margin-top: 4mm;
}
figcaption .lead {
  font-family: "Nimbus Sans", Helvetica, Arial, sans-serif;
  font-weight: bold; font-size: 9pt;
}

/* Table S1 is 27 columns; it gets landscape pages and small type. */
.supplement { page: landscape; break-before: page; }
.supplement h2 { font-size: 11pt; }
.supplement table { border-collapse: collapse; width: 100%; font-size: 5.6pt; }
.supplement th, .supplement td {
  border: 0.4pt solid #bbb; padding: 0.8pt 1.6pt; text-align: left;
  font-family: "Nimbus Sans", Helvetica, Arial, sans-serif; line-height: 1.15;
}
.supplement th { background: #eee; font-weight: bold; }
.supplement p, .supplement li { font-size: 8pt; text-align: left; }
.supplement ul { margin: 2mm 0 0 0; padding-left: 5mm; }
/* Each table, and the list of reasons, opens its own landscape page rather
   than starting three rows from the foot of the one before it. */
.supplement h3 { break-before: page; font-size: 9.5pt; margin: 0 0 2.5mm 0; }
.supplement h3.first { break-before: avoid; }
.supplement figure { break-before: page; }
.supplement figcaption { font-size: 8pt; }
.supplement table { break-inside: avoid; }
"""


def render(text, strip_first_heading=True):
    """Markdown to HTML, optionally without the file's own leading h1."""
    if strip_first_heading:
        text = re.sub(r"\A#\s+.*?\n", "", text, count=1)
    return md.markdown(text, extensions=["tables", "attr_list"])


def read(name):
    return (PAPER / f"{name}.md").read_text()


def demote(text):
    """A file's own ## subheadings sit under the section heading, not beside it."""
    return re.sub(r"^## ", "### ", text, flags=re.M)


def split_legends(text):
    """The legend file, as {figure number: (title, body html)}."""
    out = {}
    for block in re.split(r"^## ", text, flags=re.M)[1:]:
        head, _, body = block.partition("\n")
        n = re.match(r"Figure (\w+)\.", head)
        if n:
            out[n.group(1)] = (head.strip(), render(body, False))
    return out


def figure_block(number, legends, image):
    """One figure on its own page, legend beneath."""
    title, body = legends[number]
    lead, _, rest = title.partition(". ")
    rest = md.markdown(rest)[3:-4] if rest else ""
    return (
        f'<figure id="figure{number}">'
        f'<img src="{image.as_uri()}" alt="Figure {number}">'
        f'<figcaption><span class="lead">{lead}. {rest}</span> {body}</figcaption>'
        f"</figure>"
    )


def results_with_figures(legends):
    """Results, with each figure and legend after the section that cites it."""
    text = read("results")
    parts = re.split(r"^## ", text, flags=re.M)
    head = render(parts[0])
    blocks, placed = [], set()
    for i, block in enumerate(parts[1:], start=1):
        title, _, body = block.partition("\n")
        blocks.append(f"<h3>{title.strip()}</h3>" + render(body, False))
        if str(i) in legends:
            blocks.append(figure_block(str(i), legends, PAPER / "figures"
                                       / f"Figure{i}.png"))
            placed.add(str(i))
    # A numbered figure whose Results section is not yet written still belongs
    # in the manuscript; it follows the section before it rather than vanishing.
    for n in sorted(k for k in legends if k.isdigit() and k not in placed):
        blocks.append(figure_block(n, legends, PAPER / "figures" / f"Figure{n}.png"))
    return head + "".join(blocks)


def supplement(legends):
    """The supplemental figures and Table S1, on landscape pages."""
    figs = "".join(
        figure_block(n, legends, PAPER / "figures" / f"Figure{n}.png")
        for n in ("S1", "S2", "S3") if n in legends)
    table = (PAPER / "tables" / "TableS1.md").read_text()
    # Drop the rendered file's own title and its pointer back to the sources.
    table = re.sub(r"\A# Table S1\n+\*.*?\*\n", "", table, flags=re.S)
    legend = re.sub(r"\A# Table legends\n+", "", read("table_legends"))
    legend = demote(legend)
    html = render(legend, False) + render(demote(table), False)
    # The legend's own heading follows the section title and must not be pushed
    # onto a page of its own by the rule that separates the tables.
    html = html.replace("<h3>", '<h3 class="first">', 1)
    # The list of reasons belongs under table (b), not on a page of its own.
    html = html.replace("<h3>Why an ambient", '<h3 class="first">Why an ambient')
    return ('<section class="supplement">'
            "<h2>Supplemental information</h2>"
            + figs + html + "</section>")


def main() -> int:
    intro = read("introduction")
    title = re.match(r"#\s+(.*)", intro).group(1)
    legends = split_legends(read("figure_legends"))
    missing = [n for n in ("1", "2") if n not in legends]
    if missing:
        raise SystemExit(f"no legend for Figure {', '.join(missing)}")

    body = "".join([
        f'<h1 class="title">{md.markdown(title)[3:-4]}</h1>',
        '<p class="byline">Manuscript draft. Figures and tables are generated '
        "from the analysis in this repository.</p>",
        "<h2>Introduction</h2>",
        render(re.sub(r"^## Introduction\n", "", intro.split("\n", 1)[1],
                      flags=re.M), False),
        "<h2>Results</h2>", results_with_figures(legends),
        "<h2>Discussion</h2>", render(demote(read("discussion"))),
        "<h2>Methods</h2>", render(demote(read("methods"))),
        '<h2 id="references">References</h2>',
        f'<div id="references">{render(read("references"))}</div>',
        supplement(legends),
    ])

    html = f"<html><head><meta charset='utf-8'><title>{title}</title></head>" \
           f"<body>{body}</body></html>"
    HTML(string=html, base_url=str(PAPER)).write_pdf(
        OUT, stylesheets=[CSS(string=CSS_TEXT)])
    print(f"  [pdf ] paper/{OUT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
