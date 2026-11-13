#!/usr/bin/env python3
"""Put the swab pilot figures into one slide deck (.pptx), in the layout of the TEQUILA-seq progress decks: white
16:9 slides, bold centred titles, one figure per slide, short bullets, footnotes with the criteria.

Every number in a title or bullet is read from the result tables when the deck is built (read counts, read filter,
kit comparison, isoform and splicing-event tables), so the text always matches the figures. A slide is left out
when its figure does not exist.
"""
import argparse
import csv
import os
import re
from collections import defaultdict

import numpy as np
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.parts.image import Image
from pptx.util import Emu, Inches, Pt

FONT = "Arial"
BLACK = RGBColor(0, 0, 0)
GREY = RGBColor(0x40, 0x40, 0x40)
RED = RGBColor(0xC0, 0x00, 0x00)
SW, SH = 13.333, 7.5

BACKGROUND = {
    "STXBP1": {
        "title": "Developmental and epileptic encephalopathy due to STXBP1 haploinsufficiency",
        "left": ["STXBP1 (MUNC18-1) binds syntaxin-1 and other members of the SNARE complex at the presynaptic "
                 "membrane",
                 "Facilitates synaptic vesicle fusion and is critical for neurotransmitter release"],
        "right": ["De novo mutations in STXBP1 cause early infantile epileptic encephalopathy",
                  "Clinical features include seizures, intellectual disability, speech issues, movement disorders, "
                  "hypotonia, spasticity, feeding difficulties and behaviour issues",
                  "Estimated incidence of at least 1:30,000; one of the most common genes for epileptic "
                  "encephalopathies and related neurodevelopmental disorders"],
        "refs": "Saitsu H, et al. Nature Genetics. 2008;40:782-788.  Suri M, et al. Molecular Genetics and Genomic "
                "Medicine. 2017;5(5):495-507.",
    },
    "SYNGAP1": {
        "title": "Developmental and epileptic encephalopathy due to SYNGAP1 haploinsufficiency",
        "left": ["SYNGAP1 suppresses Ras-ERK1/2 signalling, which is needed for membrane insertion of AMPA receptors",
                 "Regulates synaptic plasticity and is essential for proper cognitive development"],
        "right": ["Truncating mutations in SYNGAP1 (leading to haploinsufficiency) cause intellectual disability, "
                  "epilepsy, autism and sensory processing defects",
                  "De novo SYNGAP1 mutations were identified in nonsyndromic intellectual disability and in epileptic "
                  "encephalopathies"],
        "refs": "Hamdan FF, et al. N Engl J Med. 2009;360:599-605.  Carvill GL, et al. Nature Genetics. "
                "2013;45:825-830.  Weldon M, et al. J Neurodev Disord. 2018;10(1):6.",
    },
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True, help="comma-separated")
    p.add_argument("--figures-dir", required=True, help="results/figures (steps 09 and 10)")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--read-counts", required=True)
    p.add_argument("--run-label", required=True, help="flowcell / run name for the sequencing slide")
    p.add_argument("--public-figures", help="figures of scripts/public_data_figures.py")
    p.add_argument("--rin", help="optional TSV with sample_id, RIN of this pilot's swab RNA")
    p.add_argument("--full-figures", help="figures of scripts/full_path_figures.py (all reads aligned)")
    p.add_argument("--expression-figures", help="figures of scripts/expression_qc.py (11_expression_qc.sh)")
    p.add_argument("--variant-figures", help="figures of scripts/variant_coverage.py")
    p.add_argument("--rin-figures", help="figures of scripts/rin_figures.py")
    p.add_argument("--sex-figures", help="figures of scripts/sex_check.py")
    p.add_argument("--out", required=True)
    return p.parse_args()


def panel_size():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "panel_genes.txt")
    return sum(1 for l in open(path) if l.strip() and not l.startswith("#")) if os.path.exists(path) else "NA"


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fmt_range(values, digits=0, pct=False):
    vals = sorted(values)
    if not vals:
        return "NA"
    f = "{:,.%df}" % digits
    lo, hi = f.format(vals[0]), f.format(vals[-1])
    return (lo if lo == hi else "{}-{}".format(lo, hi)) + ("%" if pct else "")


# ------------------------------------------------------------------------------------------- slide helpers

class Deck:
    def __init__(self, genes):
        self.prs = Presentation()
        self.prs.slide_width, self.prs.slide_height = Inches(SW), Inches(SH)
        self.blank = self.prs.slide_layouts[6]
        self.italic = re.compile(r"\b({})\b".format("|".join(map(re.escape, genes))))

    def runs(self, paragraph, text, size, bold=False, color=BLACK, italic=False):
        for k, part in enumerate(self.italic.split(text)):
            if not part:
                continue
            r = paragraph.add_run()
            r.text = part
            f = r.font
            f.name, f.size, f.bold, f.color.rgb = FONT, Pt(size), bold, color
            f.italic = italic or k % 2 == 1

    def slide(self, title, size=30):
        s = self.prs.slides.add_slide(self.blank)
        if title:
            tb = s.shapes.add_textbox(Inches(0.4), Inches(0.2), Inches(SW - 0.8), Inches(1.25))
            tf = tb.text_frame
            tf.word_wrap, tf.vertical_anchor = True, MSO_ANCHOR.MIDDLE
            p = tf.paragraphs[0]
            p.alignment = PP_ALIGN.CENTER
            self.runs(p, title, size, bold=True)
        return s

    def divider(self, text, sub=None):
        s = self.prs.slides.add_slide(self.blank)
        tb = s.shapes.add_textbox(Inches(1), Inches(2.6), Inches(SW - 2), Inches(2.0))
        tf = tb.text_frame
        tf.word_wrap, tf.vertical_anchor = True, MSO_ANCHOR.MIDDLE
        p = tf.paragraphs[0]
        p.alignment = PP_ALIGN.CENTER
        self.runs(p, text, 36, bold=True)
        if sub:
            p = tf.add_paragraph()
            p.alignment = PP_ALIGN.CENTER
            p.space_before = Pt(14)
            self.runs(p, sub, 18, color=GREY)
        return s

    def caption(self, s, text, x, y, w, size=16):
        tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(0.4))
        p = tb.text_frame.paragraphs[0]
        p.alignment = PP_ALIGN.CENTER
        self.runs(p, text, size, bold=True)

    def picture(self, s, path, x, y, w, h):
        px, py = Image.from_file(path).size
        scale = min(w / px, h / py)
        pw, ph = px * scale, py * scale
        s.shapes.add_picture(path, Inches(x + (w - pw) / 2), Inches(y + (h - ph) / 2), Inches(pw), Inches(ph))

    def bullets(self, s, items, x, y, w, h, size=16, color=BLACK):
        """items: text, or (level, text)."""
        tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
        tf = tb.text_frame
        tf.word_wrap = True
        for i, it in enumerate(items):
            level, text = it if isinstance(it, tuple) else (0, it)
            p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
            p.space_after = Pt(8)
            pPr = p._p.get_or_add_pPr()
            pPr.set("marL", str(int(Inches(0.3 + 0.4 * level))))
            pPr.set("indent", str(-int(Inches(0.25))))
            bu = pPr.makeelement(qn("a:buChar"), {"char": "•"})
            pPr.append(bu)
            self.runs(p, text, size - 2 * level, color=color)
        return tb

    def note(self, s, text, size=12, y=None, color=BLACK):
        tb = s.shapes.add_textbox(Inches(0.3), Inches(y if y is not None else SH - 0.75), Inches(SW - 0.6),
                                  Inches(0.65))
        tf = tb.text_frame
        tf.word_wrap, tf.vertical_anchor = True, MSO_ANCHOR.BOTTOM
        self.runs(tf.paragraphs[0], text, size, color=color)

    def table(self, s, rows, x, y, w, col_w, size=12):
        t = s.shapes.add_table(len(rows), len(rows[0]), Inches(x), Inches(y), Inches(w),
                               Inches(0.38 * len(rows))).table
        for j, cw in enumerate(col_w):
            t.columns[j].width = Inches(cw)
        for i, row in enumerate(rows):
            for j, val in enumerate(row):
                c = t.cell(i, j)
                c.fill.solid()
                c.fill.fore_color.rgb = RGBColor(0xDC, 0xE6, 0xF2) if i == 0 else RGBColor(0xFF, 0xFF, 0xFF)
                c.margin_left = c.margin_right = Inches(0.06)
                tf = c.text_frame
                tf.paragraphs[0].alignment = PP_ALIGN.CENTER
                self.runs(tf.paragraphs[0], str(val), size, bold=(i == 0))
        return t


# ------------------------------------------------------------------------------------------- content

def main():
    args = parse_args()
    genes = args.genes.split(",")
    fig = args.figures_dir
    deck_dir = os.path.join(fig, "deck")
    D = Deck(genes)
    sheet = read_tsv(args.sample_sheet)
    total = defaultdict(int)
    for r in read_tsv(args.read_counts):
        total[r["sample_id"]] += int(r["n_reads"])
    samples = [r["sample_id"] for r in sheet]
    filt = {(r["gene"], r["sample_id"]): r for r in read_tsv(os.path.join(fig, "read_filter.tsv"))}
    kit = {(r["measure"], r["gene"]): r for r in read_tsv(os.path.join(fig, "kit_comparison.tsv"))}
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    kit_text = {r["condition"]: r["treatment"] for r in sheet}

    # 1. title
    s = D.slide(None)
    tb = s.shapes.add_textbox(Inches(1), Inches(2.3), Inches(SW - 2), Inches(2.2))
    tf = tb.text_frame
    tf.word_wrap, tf.vertical_anchor = True, MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    D.runs(p, "Targeted long-read RNA sequencing of epilepsy genes in buccal swab RNA", 40, bold=True)
    p = tf.add_paragraph()
    p.alignment = PP_ALIGN.CENTER
    p.space_before = Pt(18)
    D.runs(p, "TEQUILA-seq pilot: {} in {} donors, {} swab kits".format(
        " and ".join(genes), len({r["donor"] for r in sheet}), len(kits)), 20, color=GREY)

    # 2. study design
    s = D.slide("Study design: the same donors sampled with two swab and RNA extraction kits")
    rows = [["Sample", "Donor", "Swab + RNA extraction", "Barcode"]] + [
        [r["sample_id"].replace("_", " "), r["donor"], r["treatment"], r["barcode"].replace("barcode", "NB")]
        for r in sheet]
    D.table(s, rows, 0.5, 1.7, 7.4, [1.3, 0.8, 4.3, 1.0], size=12)
    D.bullets(s, ["Questions for low-RIN swab RNA:",
                  (1, "Are {} detected, and at what depth?".format(" and ".join(genes))),
                  (1, "How much of each transcript do the reads cover? Are any reads full length?"),
                  (1, "Which isoforms and splicing events can be measured, and do the kits differ?"),
                  "Kits compared within donor (paired, n = {})".format(len({r["donor"] for r in sheet}))],
              8.3, 1.7, 4.6, 4.5, size=16)

    if args.public_figures and os.path.isdir(args.public_figures):
        public_slides(D, args, genes, samples)

    # 3. sequencing (reads per sample)
    path = os.path.join(deck_dir, "reads_per_sample.png")
    if os.path.exists(path):
        s = D.slide("Targeted long-read RNA sequencing of {} in buccal swabs".format(" and ".join(genes)))
        D.picture(s, path, 2.2, 1.45, 8.9, 3.85)
        per = [total[x] / 1e6 for x in samples]
        D.bullets(s, [
            "Prepared a multiplexed TEQUILA-seq library ({}-gene epilepsy and neurodevelopmental disorder probe "
            "panel; genes of interest {}) from buccal swab RNA of {} donors collected with {} kits ({}), based on "
            "ONT sequencing kit SQK-NBD114.24".format(
                panel_size(), " and ".join(genes), len({r["donor"] for r in sheet}), len(kits),
                "; ".join("{}: {}".format(k, kit_text[k]) for k in kits)),
            "Sequenced on {}: {:.1f} M reads in total, {} M per sample".format(
                args.run_label, sum(per), fmt_range(per, 1))], 0.4, 5.3, SW - 0.8, 2.0, size=14)

    if args.full_figures and os.path.isdir(args.full_figures):
        full_path_slides(D, args, genes, samples, filt)

    # 4. detection: reads of the two genes per swab
    path = os.path.join(deck_dir, "gene_reads.png")
    if os.path.exists(path) and filt:
        detected = {g: sum(1 for x in samples if int(filt.get((g, x), {}).get("kept", 0)) > 0) for g in genes}
        all_det = all(detected[g] == len(samples) for g in genes)
        s = D.slide("{} detected in {}".format(
            " and ".join(genes) + (" are" if len(genes) > 1 else " is"),
            "all {} swabs".format(len(samples)) if all_det else "; ".join(
                "{} in {} of {} swabs".format(g, detected[g], len(samples)) for g in genes)))
        D.picture(s, path, 0.6, 1.4, SW - 1.2, 4.1)
        kit_of = {r["sample_id"]: r["condition"] for r in sheet}
        items = []
        for g in genes:
            parts = []
            for k in kits:
                xs = [x for x in samples if kit_of[x] == k and (g, x) in filt]
                parts.append("{} {} reads ({} per million)".format(
                    k, fmt_range([int(filt[(g, x)]["kept"]) for x in xs]),
                    fmt_range([int(filt[(g, x)]["kept"]) / total[x] * 1e6 for x in xs])))
            items.append("{}: {}".format(g, "; ".join(parts)))
        D.bullets(s, items, 0.6, 5.6, SW - 1.2, 1.3, size=15)

    # 5. fragmentation (Figure 2)
    path = os.path.join(fig, "figure2.png")
    if os.path.exists(path):
        fl = {}
        for g in genes:
            k = kit.get(("full_length_reads", g))
            fl[g] = sum(int(float(v)) for part in k["per_donor"].split(";") for v in part.split(":")[1].split("/")) \
                if k else None
        title = "Swab RNA is fragmented: " + ", ".join(
            "{:,} full-length {} reads".format(fl[g], g) for g in genes if fl[g] is not None) + " in total" \
            if all(fl[g] is not None for g in genes) else "Fragmentation of {} transcripts".format(" and ".join(genes))
        s = D.slide(title)
        D.picture(s, path, 0.3, 1.5, 8.6, 3.6)
        items = []
        for g in genes:
            k = kit.get(("median_transcript_coverage", g))
            if k:
                cov = [float(v) * 100 for part in k["per_donor"].split(";") for v in part.split(":")[1].split("/")]
                items.append("{}: reads cover a median {} of the canonical transcript".format(
                    g, fmt_range(cov, 0, True)))
        items.append("Full length = the read reaches the first and the last exon of the canonical (MANE) transcript")
        D.bullets(s, items, 8.8, 1.8, 4.2, 4.8, size=15)

    rin_slide(D, args, genes)
    if args.variant_figures and os.path.isdir(args.variant_figures):
        variant_slides(D, args, genes, samples)

    for g in genes:
        bg = BACKGROUND.get(g)
        model = os.path.join(deck_dir, "{}_gene_model.png".format(g))
        if bg:
            s = D.slide(bg["title"])
            if os.path.exists(model):
                D.picture(s, model, 0.4, 1.5, 6.2, 1.6)
            D.bullets(s, bg["left"], 0.5, 3.4, 6.2, 3.0, size=16)
            D.bullets(s, bg["right"], 6.9, 1.6, 6.0, 4.9, size=16)
            D.note(s, bg["refs"], size=12)
        add_gene_slides(D, g, fig, deck_dir, samples, sheet)

    # summary
    s = D.slide("Summary")
    items = []
    if filt:
        for g in genes:
            kept = [int(filt[(g, x)]["kept"]) for x in samples if (g, x) in filt]
            items.append("{}: detected in {} of {} swabs ({} reads per swab)".format(
                g, sum(1 for v in kept if v > 0), len(samples), fmt_range(kept)))
    for g in genes:
        drawn = read_tsv(os.path.join(deck_dir, "{}_events_drawn.tsv".format(g)))
        if drawn:
            items.append("{}: {} splicing event(s) with >= 20 reads, see the event slides".format(g, len(drawn)))
    items += ["Isoform proportions from fragmented swab RNA should be read with care: partial reads are assigned to "
              "isoforms they are compatible with",
              "Open question: was DNase treatment used for each kit? Some gene reads may come from genomic DNA"]
    D.bullets(s, items, 0.8, 1.6, SW - 1.6, 5.4, size=18)

    # appendix: the earlier isoform CPM figures
    cpm_figs = []
    for g in genes:
        for name in sorted(os.listdir(os.path.join(fig, g))) if os.path.isdir(os.path.join(fig, g)) else []:
            m = re.match(r"{}_(\w+)_isoform_abundance(_by_donor)?\.png$".format(g), name)
            if m:
                cpm_figs.append((g, name, m))
    if cpm_figs:
        D.divider("Appendix: isoform abundance (CPM)")
    for g, name, m in cpm_figs:
        s = D.slide("{} isoform abundance (CPM), {}{}".format(
            g, m.group(1).replace("_", " "), ", by donor" if m.group(2) else ""), size=26)
        D.picture(s, os.path.join(fig, g, name), 0.5, 1.45, SW - 1.0, 5.5)
        D.note(s, "Isoform abundance = isoform reads per million reads sequenced; the 4 most abundant isoforms "
                  "are shown, all others are summed as Other", size=11)

    results = os.path.dirname(os.path.abspath(fig))
    report_figs = {"structure_heatmap": "isoform structures and proportions",
                   "stacked_proportions": "isoform proportions per sample",
                   "isoform_by_donor": "isoform proportions by donor",
                   "gene_reads": "reads assigned to the gene"}
    blocks = []
    for g in genes:
        gres = os.path.join(results, g)
        for method in sorted(m for m in os.listdir(gres) if os.path.isdir(os.path.join(gres, m))) \
                if os.path.isdir(gres) else []:
            for key, what in report_figs.items():
                q = os.path.join(gres, method, "{}_{}_{}.png".format(g, method, key))
                if os.path.exists(q):
                    blocks.append(("{} isoform report ({}): {}".format(g, method.replace("_", " "), what), q))
        q = os.path.join(gres, "method_comparison.png")
        if os.path.exists(q):
            blocks.append(("{} isoform proportions compared between quantification methods".format(g), q))
    frag = os.path.join(results, "fragmentation", "filtered")
    frag_figs = [("{} fragmentation".format(g), os.path.join(frag, "{}_fragmentation.png".format(g)))
                 for g in genes] + [("Detection and full-length reads, all genes and swabs",
                                     os.path.join(frag, "overview.png"))]
    frag_figs = [(t, q) for t, q in frag_figs if os.path.exists(q)]
    if blocks or frag_figs:
        D.divider("Appendix: isoform and fragmentation reports")
    for t, q in blocks + frag_figs:
        s = D.slide(t, size=26)
        D.picture(s, q, 0.5, 1.45, SW - 1.0, 5.5)

    D.prs.save(args.out)
    print("{} slides -> {}".format(len(D.prs.slides), args.out))


def full_path_slides(D, args, genes, samples, fast_filt):
    ff = args.full_figures
    vals = defaultdict(dict)
    for r in read_tsv(os.path.join(ff, "full_path_figures.tsv")):
        vals[(r["measure"], r["gene"])][r["sample_id"]] = r["value"]
    p = os.path.join(ff, "capture_qc.png")
    if os.path.exists(p):
        mr = [float(v) for v in vals[("mapping_rate_percent", "")].values()]
        ml = [float(v) for v in vals[("median_read_length", "")].values()]
        s = D.slide("{} of reads map to the genome; median read length {} nt".format(fmt_range(mr, 0, True),
                                                                                  fmt_range(ml)))
        D.picture(s, p, 1.5, 1.45, SW - 3.0, 4.6)
        D.note(s, "All reads aligned to the genome with minimap2 (full path)", size=11)
    expression_slides(D, args, genes, "capture")
    p = os.path.join(ff, "top_genes.png")
    if os.path.exists(p):
        ranked = [r for r in read_tsv(os.path.join(ff, "full_path_figures.tsv"))
                  if r["measure"] == "top_gene_percent_of_mapped"]
        names = [r["gene"] for r in ranked]
        s = D.slide("Genes collecting the most reads across all swabs")
        D.picture(s, p, 0.4, 1.45, 7.6, 5.2)
        targets = {r["gene"]: float(r["value"]) for r in read_tsv(os.path.join(ff, "full_path_figures.tsv"))
                   if r["measure"] == "target_percent_of_mapped"}
        pooled = [float(r["value"]) for r in read_tsv(os.path.join(ff, "full_path_figures.tsv"))
                  if r["measure"] == "intergenic_percent_of_mapped_pooled"]
        items = []
        if pooled:
            items.append("{:.0f}% of mapped reads fall outside annotated genes".format(pooled[0]))
        if ranked:
            items.append("Most reads in one gene: {} ({:.1f}% of mapped reads)".format(names[0],
                                                                                    float(ranked[0]["value"])))
        for g in genes:
            if g in names:
                items.append("{}: rank {} ({:.3f}% of mapped reads)".format(
                    g, names.index(g) + 1, float(ranked[names.index(g)]["value"])))
            elif g in targets:
                items.append("{}: {:.3f}% of mapped reads, not among the top {} genes".format(g, targets[g],
                                                                                             len(names)))
        items.append("A read is assigned to the gene containing its alignment midpoint; ENSG IDs are unnamed genes")
        D.bullets(s, items, 8.3, 1.7, 4.7, 5.0, size=15)
    p = os.path.join(ff, "bait_recovery.png")
    if os.path.exists(p):
        rec = {}
        for g in genes:
            pairs = [tuple(map(int, v.split("/"))) for v in vals[("kept_full_vs_fast", g)].values()]
            rec[g] = [b / a * 100 for a, b in pairs if a > 0]
        s = D.slide("k-mer baiting recovers " + "; ".join("{} of {} gene reads".format(fmt_range(rec[g], 0, True), g)
                                                          for g in genes if rec[g]) + " found by aligning all reads",
                    size=26)
        D.picture(s, p, 0.6, 1.45, 6.6, 5.3)
        items = ["Full path: every read aligned to the genome",
                 "Fast path: only reads carrying 31-mers unique to the gene, aligned afterwards",
                 "Both count the gene's reads the same way, so the fast-path numbers in this deck are complete"]
        D.bullets(s, items, 8.8, 1.7, 4.3, 5.0, size=15)
    expression_slides(D, args, genes, "samples")
    sex_slide(D, args)


def rin_slide(D, args, genes):
    rf = args.rin_figures
    p = os.path.join(rf, "rin_integrity.png") if rf else ""
    if not p or not os.path.exists(p):
        return
    rows = read_tsv(os.path.join(rf, "rin_integrity.tsv"))
    rin = {r["sample_id"]: float(r["value"]) for r in rows if r["measure"] == "RIN"}
    top = max(rin, key=rin.get)
    s = D.slide("RIN of the swabs and the read-level integrity of {}".format(" and ".join(genes)), size=28)
    D.picture(s, p, 0.4, 1.45, SW - 0.8, 3.8)
    items = ["RIN {} (mean {:.1f}); highest: {} ({:.1f})".format(fmt_range(list(rin.values()), 1),
                                                                 sum(rin.values()) / len(rin), top.replace("_", " "),
                                                                 rin[top])]
    for r in rows:
        if r["measure"].startswith("spearman_rho_median_tx_coverage"):
            items.append("{}: Spearman rho between RIN and median transcript coverage of a read = {}".format(
                r["gene"], r["value"]))
    items.append("n = 8 swabs; correlations are exploratory")
    D.bullets(s, items, 0.5, 5.35, SW - 1.0, 1.6, size=14)


def sex_slide(D, args):
    sf = args.sex_figures
    p = os.path.join(sf, "sex_check.png") if sf else ""
    if not p or not os.path.exists(p):
        return
    rows = read_tsv(os.path.join(sf, "sex_check.tsv"))
    donors = {r["donor"]: r["donor_swabs_agree"] == "True" for r in rows}
    unclear = sorted({r["donor"] for r in rows if r["call"] == "unclear"})
    s = D.slide("Sample identity check: sex-marker expression agrees between the two swabs of {} of {} donors{}".format(
        sum(donors.values()), len(donors),
        "; donor {} unclear in both".format(", ".join(unclear)) if unclear else ""), size=26)
    D.picture(s, p, 0.8, 1.5, 6.0, 5.0)
    D.bullets(s, ["XIST against Y-linked genes (RPS4Y1, DDX3Y, KDM5D, UTY, EIF1AY, USP9Y, ZFY), counts per million "
                  "from all reads",
                  "The two swabs of a donor should fall together: a split would point to a sample swap"],
              7.2, 1.8, 5.7, 4.5, size=15)


def variant_slides(D, args, genes, samples):
    vf = args.variant_figures
    rows = read_tsv(os.path.join(vf, "variant_coverage.tsv"))
    if not rows:
        return
    def covered(r, t):
        return max(int(r["{}_reads".format(s)]) for s in samples if "{}_reads".format(s) in r) >= t
    p = os.path.join(vf, "variant_coverage_summary.png")
    if os.path.exists(p):
        parts = []
        for g in genes:
            vs = [r for r in rows if r["gene"] == g]
            if vs:
                parts.append("{} {} of {}".format(g, sum(covered(r, 20) for r in vs), len(vs)))
        s = D.slide("Known patient variants covered by >= 20 reads in at least one swab: " + ", ".join(parts),
                    size=26)
        D.picture(s, p, 0.5, 1.45, SW - 1.0, 4.3)
        n_match = sum(r["reference_check"] == "match" for r in rows)
        n_mis = sum(r["reference_check"].startswith("MISMATCH") for r in rows)
        D.bullets(s, [
            "Distinct c. variants of {} from the lab's patient list; copy-number changes, exon deletions and other "
            "genes left out".format(" and ".join(genes)),
            "c. positions on the MANE Select transcripts; reference base confirmed for {} variants, {} mismatches"
            .format(n_match, n_mis),
            "Splice-site variants count reads with a base at the site and reads splicing across its intron"],
            0.5, 5.8, SW - 1.0, 1.4, size=14)
    for g in genes:
        p = os.path.join(vf, "variant_coverage_{}.png".format(g))
        vs = [r for r in rows if r["gene"] == g]
        if not os.path.exists(p) or not vs:
            continue
        five = [r for r in vs if not covered(r, 1)]
        s = D.slide("{}: {} of {} patient variants covered by >= 20 reads, {} not covered at all".format(
            g, sum(covered(r, 20) for r in vs), len(vs), len(five)), size=26)
        D.picture(s, p, 0.2, 1.35, SW - 0.4, 5.6)
        D.note(s, "Reads per variant and swab (gene reads); variants ordered 5' to 3' along the transcript; top: summed "
                  "coverage of the 8 swabs with the variant positions", size=11)


def expression_slides(D, args, genes, part):
    ef = args.expression_figures
    if not ef or not os.path.isdir(ef):
        return
    vals = defaultdict(dict)
    for r in read_tsv(os.path.join(ef, "expression_qc.tsv")):
        vals[r["measure"]][r["sample_or_pair"]] = r["value"]
    n_panel = sum(1 for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "panel_genes.txt"))
                  if l.strip() and not l.startswith("#"))
    p = os.path.join(ef, "panel_capture.png")
    if part == "capture" and os.path.exists(p):
        on = [float(v) for v in vals["on_target_percent"].values()]
        s = D.slide("On-target rate on the {}-gene panel: {} of all reads".format(n_panel, fmt_range(on, 1, True)))
        D.picture(s, p, 0.4, 1.45, SW - 0.8, 4.0)
        items = ["On-target = reads overlapping any of the {} probe-panel genes, as % of all reads sequenced".format(
            n_panel)]
        for part, what in (("Outside_annotated_genes", "outside annotated genes"), ("Unmapped", "unmapped")):
            v = [float(x) for x in vals.get("reads_percent_" + part, {}).values()]
            if v:
                items.append("Reads {}: {}".format(what, fmt_range(v, 0, True)))
        D.bullets(s, items, 0.5, 5.5, SW - 1.0, 1.4, size=14)
    for name, title, note in (
            ("pca", "Samples cluster by {} in gene expression", "PCA of log2(CPM + 1); genes with CPM >= 1 in at "
             "least 2 swabs; all reads aligned (full path)"),
            ("correlation", "Gene expression correlation between swabs", "Pearson correlation of log2(CPM + 1)"),
            ("kit_scatter", "Panel gene expression, Zymo against Genotek, per donor",
             "Each point is a panel gene; log2(CPM + 1); dashed line: equal expression")):
        p = os.path.join(ef, name + ".png")
        if part != "samples" or not os.path.exists(p):
            continue
        if name == "pca":
            title = title.format(pca_grouping(vals))
        s = D.slide(title)
        D.picture(s, p, 0.4, 1.45, SW - 0.8, 5.2)
        D.note(s, note, size=11)


def pca_grouping(vals):
    """'kit' or 'donor' or 'neither': which grouping separates the swabs along PC1 (all genes)."""
    pts = {k: float(v.split(",")[0]) for k, v in vals.get("PC1_PC2_all", {}).items()}
    if len(pts) < 4:
        return "kit or donor"
    by_kit, by_donor = defaultdict(list), defaultdict(list)
    for sid, x in pts.items():
        kit, donor = sid.split("_", 1) if "_" in sid else (sid, "")
        by_kit[kit].append(x)
        by_donor[donor].append(x)
    def between(groups):
        allv = np.array([x for g in groups.values() for x in g])
        ss_tot = ((allv - allv.mean()) ** 2).sum() or 1
        return sum(len(g) * (np.mean(g) - allv.mean()) ** 2 for g in groups.values()) / ss_tot
    k, d = between(by_kit), between(by_donor)
    return "kit" if k > d else "donor"


def public_slides(D, args, genes, samples):
    pub = args.public_figures
    mart = {}
    raw = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "public_data", "martorella2023_summary.tsv")
    for r in csv.DictReader((l for l in open(raw) if not l.startswith("#")), delimiter="\t"):
        mart[(r["measure"], r["sample_type"])] = r
    rins = []
    if args.rin and os.path.exists(args.rin):
        rins = [float(r["RIN"]) for r in read_tsv(args.rin) if r.get("RIN") and r["sample_id"] in samples]
    ours = ("This pilot: RIN {} across {} swabs, mean {:.1f}; all well below RIN 7, like the published buccal "
            "swabs (mean 2.5)".format(fmt_range(rins, 1), len(rins), sum(rins) / len(rins)) if rins
            else "This pilot: swab RNA RIN about 1-3 (per-sample values to be added)")
    cite = ("Martorella M, et al. Evaluation of noninvasive biospecimens for transcriptome studies. BMC Genomics. "
            "2023;24:790.")
    p = os.path.join(pub, "public_martorella2023.png")
    if os.path.exists(p):
        b_rin = mart[("mean_RIN", "Buccal swab")]["value"]
        q = mart[("illumina_preseq_qc", "Buccal swab")]
        if rins:
            s = D.slide("Our swab RNA is as degraded as reported for buccal swabs: RIN {} (mean {:.1f}) here, mean {} "
                        "in Martorella et al. 2023".format(fmt_range(rins, 1), sum(rins) / len(rins), b_rin), size=26)
        else:
            s = D.slide("Buccal swab RNA is highly degraded in published data too: mean RIN {}, and {} of {} swabs "
                        "failed standard library QC".format(b_rin, q["n_fail"], q["n_total"]), size=26)
        D.picture(s, p, 0.4, 1.45, SW - 0.8, 3.6)
        D.bullets(s, [
            "{} healthy donors sampled 3-4 times: buccal swab, saliva, hair follicles and urine; short-read RNA-seq"
            .format(19),
            "Buccal swabs captured a median {}% of OMIM disease genes (> 0.1 TPM), against {}% for hair and {}% for "
            "urine".format(mart[("omim_capture_percent", "Buccal swab")]["value"],
                           mart[("omim_capture_percent", "Hair follicles")]["value"],
                           mart[("omim_capture_percent", "Urine cell pellet")]["value"]),
            ours], 0.5, 5.05, SW - 1.0, 1.6, size=14)
        D.note(s, cite + " Values redrawn from the paper's text.", size=11)

    expr = read_tsv(os.path.join(pub, "public_expression.tsv"))
    val = {(r["source"], r["tissue"], r["gene"]): float(r["value"]) for r in expr}
    p = os.path.join(pub, "public_accessible.png")
    if os.path.exists(p) and expr:
        def g(t, gene, src="GTEx v10"):
            return val.get((src, t, gene))
        sal = [g("Minor Salivary Gland", x) for x in genes]
        s = D.slide("{} {} expressed at low levels in oral tissue: {} TPM in minor salivary gland".format(
            " and ".join(genes), "are" if len(genes) > 1 else "is", " and ".join("{:.1f}".format(v) for v in sal if v is not None)), size=26)
        D.picture(s, p, 0.3, 1.45, 7.9, 4.6)
        items = ["Brain cortex: {} TPM".format(" and ".join("{:.0f}".format(g("Brain - Cortex", x) or 0)
                                                            for x in genes)),
                 "Whole blood: {} TPM".format(" and ".join("{:.1f}".format(g("Whole Blood", x) or 0) for x in genes)),
                 "No public buccal swab or saliva dataset reports either gene, and no long-read RNA study of these "
                 "sample types was found",
                 "TEQUILA-seq reads per million (targeted capture) are not comparable to transcriptome-wide TPM"]
        D.bullets(s, items, 8.4, 1.6, 4.6, 5.0, size=15)
        D.note(s, "GTEx v10 median TPM per tissue (GTEx Portal API, retrieved 2026-11-12)", size=11)
    for gene in genes:
        p = os.path.join(pub, "public_gtex_{}.png".format(gene))
        if os.path.exists(p):
            s = D.slide("{} expression across GTEx v10 tissues".format(gene))
            D.picture(s, p, 0.3, 1.45, SW - 0.6, 5.0)
            D.note(s, "GTEx v10 short-read RNA-seq, TPM per sample; boxes sorted by median; whiskers 1.5 x IQR; "
                      "accessible tissues in colour", size=11)
    if expr:
        s = D.slide("Public expression of {} in oral and accessible tissues".format(" and ".join(genes)))
        rows = [["Source", "Tissue"] + genes]
        for src, tissues in (("GTEx v10 (median TPM)", ["Minor Salivary Gland", "Esophagus - Mucosa", "Whole Blood",
                                                         "Cultured fibroblasts", "Brain - Cortex"]),
                             ("HPA (nTPM)", ["salivary gland", "tongue", "esophagus", "cerebral cortex"])):
            key = "GTEx v10" if src.startswith("GTEx") else "Human Protein Atlas"
            for t in tissues:
                vals = [val.get((key, t, x)) for x in genes]
                if any(v is not None for v in vals):
                    rows.append([src, t[0].upper() + t[1:]] + ["NA" if v is None else "{:.1f}".format(v)
                                                               for v in vals])
        D.table(s, rows, 0.5, 1.6, 6.6, [2.2, 2.4] + [1.0] * len(genes), size=12)
        D.bullets(s, [
            "Clinical RNA-seq for rare-disease diagnosis (genome-wide, not these two genes) has used blood, "
            "fibroblasts, lymphoblastoid cells and muscle, not buccal swabs (Cummings 2017; Fresard 2019; Aicher "
            "2020; Murdock 2021)",
            "Buccal and saliva RNA are degraded and partly non-human, so fewer genes are captured (Martorella 2023)",
            "Detection of these two genes in buccal swab RNA has not, to our knowledge, been reported"],
            7.5, 1.6, 5.4, 5.0, size=15)
        D.note(s, "GTEx Portal API v2 and Human Protein Atlas search API, retrieved 2026-11-12. Cummings BB, et al. Sci "
                  "Transl Med 2017; Fresard L, et al. Nat Med 2019; Aicher JK, et al. Genet Med 2020; Murdock DR, et "
                  "al. J Clin Invest 2021; " + cite, size=10)


def add_gene_slides(D, g, fig, deck_dir, samples, sheet):
    gdir = os.path.join(fig, g)
    p = os.path.join(gdir, "{}_coverage_lines.png".format(g))
    if os.path.exists(p):
        s = D.slide("Read coverage along {} in each swab".format(g))
        D.picture(s, p, 2.0, 1.45, 9.3, 5.4)
        D.note(s, "Reads of the gene summed per position; introns shortened. Coverage adds up many short reads, so it "
                  "does not show fragmentation by itself: see the full-length reads and read-length slide", size=11)
    p = os.path.join(gdir, "{}_coverage_tracks.png".format(g))
    if os.path.exists(p):
        s = D.slide("{} read coverage per swab".format(g))
        D.picture(s, p, 2.0, 1.45, 9.3, 5.4)
        D.note(s, "Reads of the gene, genomic coordinates; each track is scaled to its own maximum [0 - max]",
               size=11)

    iso = read_tsv(os.path.join(deck_dir, "{}_isoforms.tsv".format(g)))
    crit = ("*Only displaying transcripts detected with an isoform fraction >= 5% (minimum of 20 reads in "
            "denominator) and read support >= 3 in at least one sample")
    p = os.path.join(deck_dir, "{}_isoform_fraction.png".format(g))
    if os.path.exists(p) and iso:
        method = iso[0]["method"].replace("_", " ")
        ok = [x for x in samples if iso[0].get("{}_fraction".format(x), "NA") != "NA"]
        canon = [r for r in iso if r["category"] == "canonical"]
        if canon and ok:
            fr = [float(canon[0]["{}_fraction".format(x)]) * 100 for x in ok]
            title = "The canonical {} isoform accounts for {} of {} isoform reads".format(g, fmt_range(fr, 0, True), g)
        elif not ok:
            title = "Too few {} reads were assigned to isoforms (fewer than 20 in every swab)".format(g)
        else:
            title = "{} isoform usage in buccal swabs".format(g)
        s = D.slide(title, size=28)
        D.picture(s, p, 0.3, 1.45, SW - 0.6, 5.3)
        D.note(s, crit + "; isoforms quantified by {}".format(method), size=11)
    p = os.path.join(deck_dir, "{}_isoform_category.png".format(g))
    if os.path.exists(p) and iso:
        ok = [x for x in samples if iso[0].get("{}_fraction".format(x), "NA") != "NA"]
        ptc = [r for r in iso if r["category"] == "Premature termination codon"]
        if ptc and ok:
            fr = [sum(float(r["{}_fraction".format(x)]) for r in ptc) * 100 for x in ok]
            title = "{} isoforms with a premature termination codon: {} of isoform reads".format(
                g, fmt_range(fr, 0, True))
        else:
            title = "{} isoforms grouped by open reading frame".format(g)
        s = D.slide(title, size=28)
        D.picture(s, p, 0.3, 1.45, SW - 0.6, 5.3)
        D.note(s, "ORF translated from the canonical start codon; premature termination codon = stop codon > 50 nt "
                  "upstream of the last exon-exon junction; incomplete ORF = no canonical start codon or no stop "
                  "codon. " + crit[1:], size=11)

    drawn = read_tsv(os.path.join(deck_dir, "{}_events_drawn.tsv".format(g)))
    levels = defaultdict(dict)
    for r in read_tsv(os.path.join(deck_dir, "{}_splicing_events.tsv".format(g))):
        levels[r["event_id"]][r["sample_id"]] = r
    for d in drawn:
        p = os.path.join(deck_dir, d["figure"] + ".png")
        if not os.path.exists(p):
            continue
        lv = [float(r["level_percent"]) for r in levels[d["event_id"]].values() if r["level_percent"] != "NA"]
        one = next(iter(levels[d["event_id"]].values()))
        is_ce = d["type"] == "Cassette exon"
        what = "{} {}{}".format(g, "" if one["in_gencode"] == "True" else "novel ",
                                "cassette exon" if is_ce else d["type"].lower())
        title = "{}: {} {} in {} swab(s) with >= 20 reads".format(
            what[0].upper() + what[1:] if not what.startswith(g) else what,
            "inclusion" if is_ce else "usage", fmt_range(lv, 0, True), len(lv)) if lv else what
        s = D.slide(title, size=28)
        D.picture(s, p, 0.3, 1.4, SW - 0.6, 5.2)
        D.note(s, "{}:{:,}-{:,}. *{} levels were estimated as the number of {} reads divided by the number of {} "
                  "reads (minimum of 20 reads in denominator)".format(
                      one["chrom"], int(one["start"]), int(one["end"]),
                      "Exon inclusion" if is_ce else "Splice site usage",
                      "inclusion junction" if is_ce else "alternative junction",
                      "inclusion or skipping junction" if is_ce else "alternative or canonical junction"), size=11)

    p = os.path.join(gdir, "{}_figure.png".format(g))
    if os.path.exists(p):
        s = D.slide("{} sashimi plots, reads and coverage per swab".format(g))
        D.picture(s, p, 0.4, 1.4, 6.6, 5.7)
        D.bullets(s, ["a  Sashimi plot per swab: coverage and junction arcs with their read counts (arcs below the "
                      "per-swab floor are hidden)",
                      "b  IGV-style reads, 25 per swab; full-length reads darker",
                      "c  Gene-body coverage along the canonical transcript, 5' to 3'",
                      "d  Reads on canonical exons; the number above each bar is the count of "
                      "full-length reads"], 7.3, 1.7, 5.6, 5.0, size=15)
        D.note(s, "Reads of the gene; canonical (MANE) gene model; introns shortened", size=11)
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    by = {(r["donor"], r["condition"]): r["sample_id"] for r in sheet}
    for d in dict.fromkeys(r["donor"] for r in sheet):
        for kind, what in (("all_reads", "reads of the gene"),):
            pics = [(k, os.path.join(gdir, "igv", "{}_{}_{}.png".format(g, by[(d, k)], kind)))
                    for k in kits if (d, k) in by]
            pics = [(k, q) for k, q in pics if os.path.exists(q)]
            if not pics:
                continue
            s = D.slide("{} IGV view, donor {}".format(g, d), size=28)
            w = (SW - 0.6 - 0.2 * (len(pics) - 1)) / len(pics)
            for i, (k, q) in enumerate(pics):
                x = 0.3 + i * (w + 0.2)
                D.caption(s, k, x, 1.4, w)
                D.picture(s, q, x, 1.85, w, 4.9)
            D.note(s, "Reads of {}; only {} annotated".format(g, g), size=11)


if __name__ == "__main__":
    main()
