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
    p.add_argument("--out", required=True)
    return p.parse_args()


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
            "Prepared a multiplexed TEQUILA-seq library ({} target gene{}) from buccal swab RNA of {} donors "
            "collected with {} kits ({}), based on ONT sequencing kit SQK-NBD114.24".format(
                len(genes), "s" if len(genes) > 1 else "", len({r["donor"] for r in sheet}), len(kits),
                "; ".join("{}: {}".format(k, kit_text[k]) for k in kits)),
            "Sequenced on {}: {:.1f} M reads in total, {} M per sample".format(
                args.run_label, sum(per), fmt_range(per, 1))], 0.4, 5.3, SW - 0.8, 2.0, size=14)

    if args.full_figures and os.path.isdir(args.full_figures):
        full_path_slides(D, args, genes, samples, filt)

    # 4. detection (Figure 1)
    path = os.path.join(fig, "figure1.png")
    if os.path.exists(path) and filt:
        detected = {g: sum(1 for x in samples if int(filt.get((g, x), {}).get("kept", 0)) > 0) for g in genes}
        all_det = all(detected[g] == len(samples) for g in genes)
        s = D.slide("{} detected in {}".format(
            " and ".join(genes) + (" are" if len(genes) > 1 else " is"),
            "all {} swabs".format(len(samples)) if all_det else "; ".join(
                "{} in {} of {} swabs".format(g, detected[g], len(samples)) for g in genes)))
        D.picture(s, path, 0.4, 1.4, SW - 0.8, 3.85)
        items = []
        for g in genes:
            kept = [int(filt[(g, x)]["kept"]) for x in samples if (g, x) in filt]
            rpm = [int(filt[(g, x)]["kept"]) / total[x] * 1e6 for x in samples if (g, x) in filt]
            spl = [int(filt[(g, x)]["spliced_kept"]) / max(int(filt[(g, x)]["kept"]), 1) * 100
                   for x in samples if (g, x) in filt]
            items.append("{}: {} gene reads per sample ({} per million); {} spliced".format(
                g, fmt_range(kept), fmt_range(rpm), fmt_range(spl, 0, True)))
            k = kit.get(("reads_per_million", g))
            if k and len(kits) == 2:
                items.append((1, "{} vs {}: mean log2 ratio {}, higher in {} of {} donors (paired t p = {})".format(
                    kits[1], kits[0], k["mean_log2_ratio_{}_vs_{}".format(kits[1], kits[0])],
                    k["donors_higher_in_{}".format(kits[1])], k["n_donors"], k["paired_t_p"])))
        D.bullets(s, items, 0.6, 5.3, SW - 1.2, 1.6, size=14)
        D.note(s, "Gene reads: primary alignments with >= 50% of aligned bases on the gene's exons, more than on any "
                  "other gene, and (spliced reads) on the gene's strand", size=11)

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
        D.picture(s, path, 0.3, 1.45, 8.3, 5.3)
        items = []
        for g in genes:
            k = kit.get(("median_transcript_coverage", g))
            if k:
                cov = [float(v) * 100 for part in k["per_donor"].split(";") for v in part.split(":")[1].split("/")]
                items.append("{}: reads cover a median {} of the canonical transcript".format(
                    g, fmt_range(cov, 0, True)))
        items.append("Full length = the read reaches the first and the last exon of the canonical (MANE) transcript")
        D.bullets(s, items, 8.8, 1.8, 4.2, 4.8, size=15)

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
            spl = [int(filt[(g, x)]["spliced_kept"]) for x in samples if (g, x) in filt]
            items.append("{}: detected in {} of {} swabs ({} gene reads, {} spliced reads per sample)".format(
                g, sum(1 for v in kept if v > 0), len(samples), fmt_range(kept), fmt_range(spl)))
    for g in genes:
        drawn = read_tsv(os.path.join(deck_dir, "{}_events_drawn.tsv".format(g)))
        if drawn:
            items.append("{}: {} splicing event(s) with >= 20 reads, see the event slides".format(g, len(drawn)))
    items += ["Isoform proportions from fragmented swab RNA should be read with care: partial reads are assigned to "
              "isoforms they are compatible with",
              "Open question: was DNase treatment used for each kit? Unspliced reads may include genomic DNA"]
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
    frag_figs = [("{} fragmentation (gene-filtered reads)".format(g), os.path.join(frag, "{}_fragmentation.png".format(g)))
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
        on = {g: [float(v) for v in vals[("on_target_percent", g)].values()] for g in genes}
        s = D.slide("On-target rate: " + "; ".join("{} {} of all reads".format(g, fmt_range(on[g], 3, True))
                                                     for g in genes if on[g]), size=28)
        D.picture(s, p, 0.4, 1.45, SW - 0.8, 3.6)
        mr = [float(v) for v in vals[("mapping_rate_percent", "")].values()]
        ml = [float(v) for v in vals[("median_read_length", "")].values()]
        inter = vals[("intergenic_percent_of_mapped", "")]
        kit_of = {r["sample_id"]: r["condition"] for r in read_tsv(args.sample_sheet)}
        by_kit = defaultdict(list)
        for x, v in inter.items():
            by_kit[kit_of.get(x, "")].append(float(v))
        items = ["All reads aligned to the genome (minimap2, full path): {} mapped; median read length {} nt".format(
                     fmt_range(mr, 0, True), fmt_range(ml)),
                 "On-target = reads overlapping the gene's annotated span, as % of all reads sequenced"]
        if inter:
            items.append("Reads outside annotated genes: " + "; ".join(
                "{} {}".format(k, fmt_range(v, 0, True)) for k, v in by_kit.items()) + " of mapped reads")
        D.bullets(s, items, 0.5, 5.15, SW - 1.0, 1.6, size=14)
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
        D.picture(s, p, 0.4, 1.45, 8.2, 5.2)
        items = []
        for g in genes:
            pairs = [tuple(map(int, v.split("/"))) for v in vals[("spliced_kept_full_vs_fast", g)].values()]
            sp = [b / a * 100 for a, b in pairs if a > 0]
            if sp:
                items.append("{}: spliced gene reads recovered {}".format(g, fmt_range(sp, 0, True)))
        items.append("Both paths use the same gene read filter; fast path = reads carrying gene-unique 31-mers, "
                     "aligned afterwards")
        D.bullets(s, items, 8.8, 1.7, 4.3, 5.0, size=15)


def public_slides(D, args, genes, samples):
    pub = args.public_figures
    mart = {}
    raw = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "public_data", "martorella2023_summary.tsv")
    for r in csv.DictReader((l for l in open(raw) if not l.startswith("#")), delimiter="\t"):
        mart[(r["measure"], r["sample_type"])] = r
    rins = []
    if args.rin and os.path.exists(args.rin):
        rins = [float(r["RIN"]) for r in read_tsv(args.rin) if r.get("RIN") and r["sample_id"] in samples]
    ours = ("This pilot: RIN {} across {} swabs".format(fmt_range(rins, 1), len(rins)) if rins
            else "This pilot: swab RNA RIN about 1-3 (per-sample values to be added)")
    cite = ("Martorella M, et al. Evaluation of noninvasive biospecimens for transcriptome studies. BMC Genomics. "
            "2023;24:790.")
    p = os.path.join(pub, "public_martorella2023.png")
    if os.path.exists(p):
        b_rin = mart[("mean_RIN", "Buccal swab")]["value"]
        q = mart[("illumina_preseq_qc", "Buccal swab")]
        s = D.slide("Buccal swab RNA is highly degraded in published data too: mean RIN {}, and {} of {} swabs failed "
                    "standard library QC".format(b_rin, q["n_fail"], q["n_total"]), size=26)
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
                                                         "Cells - Cultured fibroblasts", "Brain - Cortex"]),
                             ("HPA (nTPM)", ["salivary gland", "tongue", "esophagus", "cerebral cortex"])):
            key = "GTEx v10" if src.startswith("GTEx") else "Human Protein Atlas"
            for t in tissues:
                vals = [val.get((key, t, x)) for x in genes]
                if any(v is not None for v in vals):
                    rows.append([src, t[0].upper() + t[1:]] + ["NA" if v is None else "{:.1f}".format(v)
                                                               for v in vals])
        D.table(s, rows, 0.5, 1.6, 6.6, [2.2, 2.4] + [1.0] * len(genes), size=12)
        D.bullets(s, [
            "RNA diagnostics in clinically accessible tissues has used blood, fibroblasts, lymphoblastoid cells and "
            "muscle, not buccal swabs (Cummings 2017; Fresard 2019; Aicher 2020; Murdock 2021)",
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
        s = D.slide("Read coverage along {}: all reads and spliced reads".format(g))
        D.picture(s, p, 2.0, 1.45, 9.3, 5.4)
        D.note(s, "Gene-filtered reads; introns shortened; spliced reads = reads with at least one intron", size=11)
    p = os.path.join(gdir, "{}_coverage_tracks.png".format(g))
    if os.path.exists(p):
        s = D.slide("{} read coverage per swab".format(g))
        D.picture(s, p, 2.0, 1.45, 9.3, 5.4)
        D.note(s, "Gene-filtered reads, genomic coordinates; each track is scaled to its own maximum [0 - max]",
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
                  "reads (minimum of 20 reads in denominator); gene-filtered reads".format(
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
                      "d  Reads on canonical exons, spliced and unspliced; the number above each bar is the count of "
                      "full-length reads"], 7.3, 1.7, 5.6, 5.0, size=15)
        D.note(s, "Gene-filtered reads; canonical (MANE) gene model; introns shortened", size=11)
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    by = {(r["donor"], r["condition"]): r["sample_id"] for r in sheet}
    for d in dict.fromkeys(r["donor"] for r in sheet):
        for kind, what in (("all_reads", "all gene-filtered reads"), ("spliced_reads", "spliced reads only")):
            pics = [(k, os.path.join(gdir, "igv", "{}_{}_{}.png".format(g, by[(d, k)], kind)))
                    for k in kits if (d, k) in by]
            pics = [(k, q) for k, q in pics if os.path.exists(q)]
            if not pics:
                continue
            s = D.slide("{} IGV view, donor {}: {}".format(g, d, what), size=28)
            w = (SW - 0.6 - 0.2 * (len(pics) - 1)) / len(pics)
            for i, (k, q) in enumerate(pics):
                x = 0.3 + i * (w + 0.2)
                D.caption(s, k, x, 1.4, w)
                D.picture(s, q, x, 1.85, w, 4.9)
            D.note(s, "Gene-filtered reads; only {} annotated".format(g), size=11)


if __name__ == "__main__":
    main()
