#!/usr/bin/env python3
"""How well the sample reads cover known patient variants of the target genes.

Input: a variant table (tab- or |-separated, header with a gene column and a nucleotide (c.) column; protein and
pathogenicity columns are used when present). c. positions are read on the gene's canonical (MANE Select)
transcript from the GTF. Rows that are not c. variants of a target gene (copy-number changes, exon deletions, other
genes) are listed in variants_skipped.tsv with the reason. The table holds patient information: keep it outside the
repository.

For every variant and sample (gene reads, i.e. the gene-filtered BAMs):
  * exonic variants: reads with an aligned base at the variant (for a multi-base variant: at its first and last
    exonic base)
  * intronic / splice-site variants (c.N+k, c.N-k): reads with an aligned base at the site, plus reads that splice
    across the intron holding the site (those show the junction the variant would affect); both are counted
The reference base of every substitution (and the deleted / duplicated bases when given) is checked against the
genome, which catches a c. numbering on a different transcript.

Outputs (--outdir): variant_coverage.tsv, variants_skipped.tsv, variant_coverage_<GENE>.{pdf,png} (coverage of all
samples with the variant positions, and reads per variant and sample), variant_coverage_summary.{pdf,png}.
"""
import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import ConnectionPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import IntronShrink, canonical_model  # noqa: E402

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.linewidth": 0.8, "pdf.fonttype": 42, "savefig.facecolor": "white"})
CONSEQ_COLORS = {"Missense": "#3A8DC4", "Nonsense": "#C00000", "Frameshift": "#EC7B02", "In-frame indel": "#57B45B",
                 "Splice site / intronic": "#8E44AD", "Synonymous": "#999999", "Other": "#555555"}
DEPTH_BINS = [0, 1, 10, 20, 10 ** 9]
DEPTH_COLORS = ["#F2F2F2", "#FED976", "#FD8D3C", "#BD0026"]
COMP = str.maketrans("ACGTN", "TGCAN")
POS = r"[-*]?\d+(?:[+-]\d+)?"
HGVS = re.compile(r"c\.(" + POS + r"(?:_" + POS + r")?)((?:[ACGT]>[ACGT])|delins[ACGT0-9]*|del[ACGT0-9]*(?:ins[ACGT0-9]+)?|"
                  r"dup[ACGT0-9]*|ins[ACGT0-9]+)?")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--variants", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--bam-pattern", required=True, help="gene reads: path with {gene} and {sample}")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


# ------------------------------------------------------------------------------------------------ input

def read_variants(path, genes):
    text = open(path, encoding="utf-8", errors="replace").read().splitlines()
    sep = "|" if "|" in text[0] and "\t" not in text[0] else "\t"
    header = [h.strip().lower() for h in text[0].split(sep)]

    def col(*keys):
        for i, h in enumerate(header):
            if any(k in h for k in keys):
                return i
        return None
    ig, ic, ip, ipath = col("gene"), col("nucleotide", "c."), col("amino", "protein"), col("pathogen")
    kept, skipped = [], []
    for line in text[1:]:
        f = [x.strip() for x in line.split(sep)]
        if not any(f):
            continue
        f += [""] * (len(header) - len(f))
        raw_gene = f[ig] if ig is not None else ""
        m = re.match(r"([A-Za-z0-9]+)", raw_gene)
        gene = m.group(1).upper() if m else ""
        blob = " ".join(x for x in (f[ic] if ic is not None else "", raw_gene) if x)
        squeezed = re.sub(r"[\s()]", "", blob)
        hit = HGVS.search(squeezed)
        protein = f[ip] if ip is not None and f[ip] else (re.search(r"p\.[^\s,]+", blob).group(0)
                                                          if re.search(r"p\.[^\s,]+", blob) else "")
        rec = {"gene": gene, "given": blob.strip(), "protein": protein.strip(),
               "pathogenicity": f[ipath] if ipath is not None else ""}
        if gene not in genes:
            skipped.append((rec, "not a target gene" if gene else "no gene"))
        elif not hit:
            skipped.append((rec, "not a c. variant (copy-number change, exon deletion or breakpoints)"))
        else:
            rec["hgvs"] = "c." + hit.group(1) + (hit.group(2) or "")
            rec["pos"], rec["change"] = hit.group(1), hit.group(2) or ""
            kept.append(rec)
    return kept, skipped


# ------------------------------------------------------------------------------------------------ c. -> genome

class CdnaMap:
    def __init__(self, model):
        self.strand = model["strand"]
        exons = model["exons"] if self.strand == "+" else model["exons"][::-1]
        self.tx = []  # genomic position of every transcript base, 5' -> 3'
        for s, e in exons:
            rng = range(s, e + 1) if self.strand == "+" else range(e, s - 1, -1)
            self.tx.extend(rng)
        cds = sorted(model["cds"])
        start = cds[0][0] if self.strand == "+" else cds[-1][1]
        self.c1 = self.tx.index(start)
        self.cds_len = sum(e - s + 1 for s, e in cds)
        self.exonic = set(self.tx)
        self.dir = 1 if self.strand == "+" else -1

    def base(self, spec):
        """c. position spec ('123', '-45', '*10', '123+5', '124-2') -> (genomic position, intronic offset)."""
        m = re.fullmatch(r"([-*]?)(\d+)(?:([+-])(\d+))?", spec)
        kind, n = m.group(1), int(m.group(2))
        if kind == "-":
            i = self.c1 - n
        elif kind == "*":
            i = self.c1 + self.cds_len + 2 + n  # after the stop codon
        else:
            i = self.c1 + n - 1
        if not 0 <= i < len(self.tx):
            raise ValueError("outside the transcript")
        g = self.tx[i]
        off = 0
        if m.group(3):
            off = int(m.group(4)) * (1 if m.group(3) == "+" else -1)
        return g + off * self.dir, off, g


def locate(rec, cmap):
    pos, change = rec["pos"], rec["change"]
    parts = pos.split("_")
    # "c.1022-1030delins" / "c.2955-2958del": a range written with a dash
    m = re.fullmatch(r"(\d+)-(\d+)", pos)
    if m and change[:3] in ("del", "dup", "ins") and int(m.group(2)) > int(m.group(1)):
        parts = [m.group(1), m.group(2)]
        rec["note"] = "range written with '-' read as {}_{}".format(*parts)
    ends = [cmap.base(p) for p in parts]
    gpos = sorted(g for g, _, _ in ends)
    rec["g_start"], rec["g_end"] = gpos[0], gpos[-1]
    rec["intronic"] = any(off != 0 for _, off, _ in ends)
    rec["anchor"] = [anchor for _, _, anchor in ends]  # nearest exonic base of each end
    return rec


def ref_check(rec, fasta, chrom, strand):
    ch = rec["change"]
    if ">" in ch:
        want = ch[0]
        got = fasta.fetch(chrom, rec["g_start"] - 1, rec["g_start"]).upper()
        got = got.translate(COMP) if strand == "-" else got
        return "match" if got == want else "MISMATCH (genome {})".format(got)
    m = re.match(r"(?:del|dup)([ACGT]+)", ch)
    if m and not rec["intronic"]:
        want = m.group(1)
        seq = fasta.fetch(chrom, rec["g_start"] - 1, rec["g_end"]).upper()
        seq = seq.translate(COMP)[::-1] if strand == "-" else seq
        if len(want) == len(seq):
            return "match" if seq == want else "MISMATCH (genome {})".format(seq)
    return "not checked"


def consequence(rec):
    """variant type from the protein change, else from the nucleotide change."""
    p = re.sub(r"^P\.", "p.", rec["protein"].replace(" ", ""))
    if rec["intronic"] or "p.?" in p or "p.(?)" in p:
        return "Splice site / intronic"
    if "fs" in p.lower():
        return "Frameshift"
    if "=" in p:
        return "Synonymous"
    if re.search(r"\*|Ter|X\d*\)?\.?$", p):
        return "Nonsense"
    if re.search(r"del|dup|ins", p):
        return "In-frame indel"
    if re.search(r"p\.\(?[A-Za-z]{1,3}\d+[A-Za-z]{1,3}\)?\.?$", p):
        return "Missense"
    ch = rec["change"]
    if ch.startswith(("del", "dup", "ins")):
        n = len(re.sub(r"[^ACGT]", "", ch.replace("del", "", 1)))
        return "Frameshift" if n % 3 else "In-frame indel"
    return "Other"


# ------------------------------------------------------------------------------------------------ reads

def depth(bam_path, chrom, rec):
    """(reads with aligned bases at the variant, reads splicing across the intron holding an intronic site)."""
    if not os.path.exists(bam_path):
        return 0, 0
    lo, hi = rec["g_start"], rec["g_end"]
    if rec["intronic"]:
        sites = {lo, hi}
    else:
        sites = {min(rec["anchor"]), max(rec["anchor"])}
    n_base = n_junc = 0
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(chrom, lo - 1, hi):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            covered = set()
            for s, e in r.get_blocks():
                for x in sites:
                    if s < x <= e:
                        covered.add(x)
            if covered == sites or (rec["intronic"] and covered):
                n_base += 1
                continue
            if rec["intronic"]:
                pos = r.reference_start
                for op, ln in r.cigartuples or ():
                    if op == 3 and pos < lo and hi <= pos + ln:
                        n_junc += 1
                        break
                    if op in (0, 2, 3, 7, 8):
                        pos += ln
    return n_base, n_junc


# ------------------------------------------------------------------------------------------------ figures

def gene_figure(gene, model, variants, samples, kits_of, cov_sum, glo, out):
    if not variants:
        return
    shrink = IntronShrink(model["exons"], intron_fraction=0.35)
    order = sorted(variants, key=lambda v: v["g_start"] * (1 if model["strand"] == "+" else -1))
    W, H = 13.0, 6.4
    fig = plt.figure(figsize=(W, H))
    # coverage of all samples with the variant positions
    ax = fig.add_axes([0.7 / W, 4.55 / H, (W - 1.0) / W, 1.5 / H])
    lo, hi = model["exons"][0][0], model["exons"][-1][1]
    pos = np.arange(lo, hi + 1)
    x = shrink(pos + 0.5)
    bins = 900
    idx = np.minimum((x * bins).astype(int), bins - 1)
    y = np.bincount(idx, weights=cov_sum[lo - glo: hi - glo + 1], minlength=bins) / np.maximum(
        np.bincount(idx, minlength=bins), 1)
    ax.fill_between(np.linspace(0, 1, bins + 1), np.append(y, y[-1]), step="post", color="#708090", lw=0)
    for v in order:
        ax.axvline(float(shrink(v["g_start"])), color=CONSEQ_COLORS[v["consequence"]], lw=0.8, alpha=0.9)
    ax.set_xlim(-0.005, 1.005)
    ax.set_ylabel("Reads, all samples")
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.set_xticks([])
    ax.set_title("{}: read coverage of the 8 samples (summed) and patient variant positions".format(gene),
                 loc="left", fontsize=10, fontstyle="normal")
    # gene model
    axm = fig.add_axes([0.7 / W, 4.2 / H, (W - 1.0) / W, 0.3 / H])
    axm.set_xlim(-0.005, 1.005)
    axm.set_ylim(-1, 1)
    axm.axis("off")
    axm.plot([shrink(lo), shrink(hi + 1)], [0, 0], color="black", lw=0.8)
    for s, e in model["exons"]:
        axm.add_patch(Rectangle((shrink(s), -0.6), max(shrink(e + 1) - shrink(s), 0.0015), 1.2, color="black"))
    # heatmap: reads per variant and sample
    axh = fig.add_axes([0.7 / W, 1.55 / H, (W - 1.0) / W, 1.55 / H])
    mat = np.array([[v["depth"][s] for v in order] for s in samples], float)
    cmap = ListedColormap(DEPTH_COLORS)
    axh.imshow(mat, aspect="auto", cmap=cmap, norm=BoundaryNorm(DEPTH_BINS, cmap.N), interpolation="nearest")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            axh.text(j, i, "{:.0f}".format(mat[i, j]), ha="center", va="center", fontsize=5,
                     color="white" if mat[i, j] >= 20 else "black")
    axh.set_yticks(range(len(samples)))
    axh.set_yticklabels([s.replace("_", " ") for s in samples], fontsize=8)
    axh.set_xticks(range(len(order)))
    axh.set_xticklabels([v["hgvs"] for v in order], rotation=90, fontsize=6)
    for lbl, v in zip(axh.get_xticklabels(), order):
        lbl.set_color(CONSEQ_COLORS[v["consequence"]])
    axh.tick_params(length=0)
    n = len(order)
    for j, v in enumerate(order):
        fig.add_artist(ConnectionPatch((float(shrink(v["g_start"])), -0.6), (j, -0.5), "data", "data", axesA=axm,
                                       axesB=axh, color=CONSEQ_COLORS[v["consequence"]], lw=0.4, alpha=0.6))
    # legends
    handles = [Rectangle((0, 0), 1, 1, color=c) for c in CONSEQ_COLORS.values()]
    used = [k for k in CONSEQ_COLORS if any(v["consequence"] == k for v in order)]
    fig.legend([handles[list(CONSEQ_COLORS).index(k)] for k in used], used, loc="lower left",
               bbox_to_anchor=(0.7 / W, 0.0), ncol=len(used), frameon=False, fontsize=8, title="Variant type",
               title_fontsize=8)
    fig.legend([Rectangle((0, 0), 1, 1, color=c, ec="#999999", lw=0.4) for c in DEPTH_COLORS],
               ["0", "1-9", "10-19", ">= 20"], loc="lower right", bbox_to_anchor=(1 - 0.3 / W, 0.0), ncol=4,
               frameon=False, fontsize=8, title="Reads covering the variant", title_fontsize=8)
    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png", dpi=250)
    plt.close(fig)


def summary_figure(genes, variants, samples, kits_of, out):
    fig, axes = plt.subplots(1, len(genes), figsize=(5.0 * len(genes), 3.4), squeeze=False)
    fig.subplots_adjust(left=0.07, right=0.95, top=0.88, bottom=0.32, wspace=0.25)
    thresholds = [(1, "#FED976"), (10, "#FD8D3C"), (20, "#BD0026")]
    for ax, g in zip(axes[0], genes):
        vs = [v for v in variants if v["gene"] == g]
        x = np.arange(len(samples))
        w = 0.8 / len(thresholds)
        for k, (t, c) in enumerate(thresholds):
            ax.bar(x - 0.4 + w * (k + 0.5), [sum(v["depth"][s] >= t for v in vs) for s in samples], width=w,
                   color=c, label=">= {} reads".format(t))
        ax.axhline(len(vs), color="#555555", lw=0.8, ls="--")
        ax.text(len(samples) - 0.5, len(vs), "{} variants".format(len(vs)), ha="right", va="bottom", fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels([s.replace("_", " ") for s in samples], rotation=-30, ha="left", rotation_mode="anchor")
        ax.set_ylabel("Patient variants covered")
        ax.set_ylim(0, len(vs) * 1.15 + 1)
        ax.set_title(g, fontstyle="italic", fontsize=11)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=8, loc="lower center", ncol=len(thresholds), bbox_to_anchor=(0.5, 0.0))
    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png", dpi=300)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ main

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    genes = args.genes.split(",")
    with open(args.sample_sheet) as fh:
        sheet = list(csv.DictReader(fh, delimiter="\t"))
    conds = list(dict.fromkeys(r["condition"] for r in sheet))
    samples = [r["sample_id"] for c in conds for r in sheet if r["condition"] == c]
    kits_of = {r["sample_id"]: r["condition"] for r in sheet}
    fasta = pysam.FastaFile(args.fasta)
    kept, skipped = read_variants(args.variants, set(genes))

    # one entry per distinct variant; n_listed counts how often it appears in the table
    by_key = {}
    for rec in kept:
        key = (rec["gene"], rec["hgvs"])
        if key in by_key:
            by_key[key]["n_listed"] += 1
        else:
            rec["n_listed"] = 1
            by_key[key] = rec
    variants, models = [], {}
    for (g, _), rec in by_key.items():
        if g not in models:
            models[g] = canonical_model(args.gtf, g)
        m = models[g]
        try:
            locate(rec, CdnaMap(m))
        except (ValueError, AttributeError) as e:
            skipped.append((rec, "c. position not on the canonical transcript ({})".format(e)))
            continue
        rec["chrom"] = m["chrom"]
        rec["ref_check"] = ref_check(rec, fasta, m["chrom"], m["strand"])
        rec["consequence"] = consequence(rec)
        rec["depth"], rec["base"], rec["junction"] = {}, {}, {}
        for s in samples:
            b, j = depth(args.bam_pattern.format(gene=g, sample=s), m["chrom"], rec)
            rec["base"][s], rec["junction"][s], rec["depth"][s] = b, j, b + j
        variants.append(rec)

    with open(os.path.join(args.outdir, "variant_coverage.tsv"), "w") as fh:
        fh.write("\t".join(["gene", "variant", "as_listed", "protein", "pathogenicity", "times_listed", "type",
                            "location", "chrom", "start", "end", "reference_check", "note"] +
                           ["{}_reads".format(s) for s in samples] +
                           ["{}_junction_reads".format(s) for s in samples] +
                           ["all_samples_reads", "samples_with_20_reads"]) + "\n")
        for v in sorted(variants, key=lambda v: (v["gene"], v["g_start"])):
            fh.write("\t".join(map(str, [
                v["gene"], v["hgvs"], v["given"], v["protein"], v["pathogenicity"], v["n_listed"], v["consequence"],
                "intronic/splice" if v["intronic"] else "exonic", v["chrom"], v["g_start"], v["g_end"],
                v["ref_check"], v.get("note", "")] + [v["depth"][s] for s in samples] +
                [v["junction"][s] for s in samples] +
                [sum(v["depth"].values()), sum(1 for s in samples if v["depth"][s] >= 20)])) + "\n")
    with open(os.path.join(args.outdir, "variants_skipped.tsv"), "w") as fh:
        fh.write("gene\tas_listed\treason\n")
        for rec, why in skipped:
            fh.write("{}\t{}\t{}\n".format(rec["gene"], rec["given"], why))

    for g in genes:
        if g not in models:
            continue
        m = models[g]
        lo, hi = m["exons"][0][0], m["exons"][-1][1]
        cov = np.zeros(hi - lo + 2)
        for s in samples:
            path = args.bam_pattern.format(gene=g, sample=s)
            if not os.path.exists(path):
                continue
            with pysam.AlignmentFile(path) as bam:
                for r in bam.fetch(m["chrom"], lo - 1, hi):
                    if r.is_unmapped or r.is_secondary or r.is_supplementary:
                        continue
                    for a, b in r.get_blocks():
                        a, b = max(a + 1, lo), min(b, hi)
                        if a <= b:
                            cov[a - lo] += 1
                            cov[b - lo + 1] -= 1
        cov = np.cumsum(cov)[:-1]
        gene_figure(g, m, [v for v in variants if v["gene"] == g], samples, kits_of, cov, lo,
                    os.path.join(args.outdir, "variant_coverage_{}".format(g)))
    summary_figure([g for g in genes if g in models], variants, samples, kits_of,
                   os.path.join(args.outdir, "variant_coverage_summary"))
    for g in genes:
        vs = [v for v in variants if v["gene"] == g]
        if vs:
            print("{}: variants covered by >= 20 reads per sample: {}".format(g, ", ".join(
                "{} {}".format(sm, sum(v["depth"][sm] >= 20 for v in vs)) for sm in samples)))
    mism = [v for v in variants if v["ref_check"].startswith("MISMATCH")]
    for g in genes:
        vs = [v for v in variants if v["gene"] == g]
        if vs:
            print("{}: {} distinct variants; covered by >= 1 / >= 10 / >= 20 reads in at least one sample: {} / {} / {}".format(
                g, len(vs), *[sum(1 for v in vs if max(v["depth"].values()) >= t) for t in (1, 10, 20)]))
    print("skipped rows: {} (variants_skipped.tsv); reference check: {} match, {} mismatch, {} not checked".format(
        len(skipped), sum(v["ref_check"] == "match" for v in variants), len(mism),
        sum(v["ref_check"] == "not checked" for v in variants)))
    for v in mism:
        print("  reference mismatch: {} {} ({})".format(v["gene"], v["hgvs"], v["ref_check"]))


if __name__ == "__main__":
    main()
