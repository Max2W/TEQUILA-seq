#!/usr/bin/env python3
"""Do accessible tissues use the same isoforms of the panel genes as brain? (public data; GTEx v8 transcript medians)

For each panel gene and tissue, isoform usage = each transcript's median TPM divided by the gene's summed transcript
TPM (public_data/gtex_v8_transcript_medians.tsv.gz, from scripts/fetch_gtex_isoforms.py). Each tissue is compared
with brain cortex, for genes with a summed transcript TPM >= --min-tpm in both:
  overlap          sum over transcripts of min(usage in cortex, usage in tissue): 1 = the same isoform mix, 0 = none
                   shared
  same_dominant    the most used transcript is the same in both
  mane_share       usage of the MANE Select transcript in the tissue (and in cortex)
Brain cerebellum against cortex is the within-brain reference. Tissues stand in for accessible specimens: lung for
nasal epithelium (no public transcript-level nasal data), esophagus mucosa for buccal, skin for hair follicles,
minor salivary gland for saliva; fibroblasts, EBV lymphocytes and whole blood as such.

Outputs (--outdir): isoform_similarity.tsv (gene x tissue), isoform_usage_targets.tsv (usage of every transcript of
the --genes per tissue, with transcript names from --gtf), isoform_similarity.{pdf,png} and
isoform_targets.{pdf,png}, isoform_similarity_figures.tsv (summary numbers).
"""
import argparse
import csv
import gzip
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "savefig.facecolor": "white"})
REF = "Brain_Cortex"
TISSUES = [  # id, label, stands for, colour
    ("Brain_Cerebellum", "Brain cerebellum", "within-brain reference", "#555555"),
    ("Lung", "Lung", "airway: nasal", "#1B7837"),
    ("Cells_Cultured_fibroblasts", "Fibroblasts", "skin biopsy", "#F4A261"),
    ("Minor_Salivary_Gland", "Salivary gland", "saliva", "#E8957B"),
    ("Skin_Not_Sun_Exposed_Suprapubic", "Skin", "hair follicle", "#D6604D"),
    ("Esophagus_Mucosa", "Esophagus mucosa", "squamous: buccal", "#708090"),
    ("Cells_EBV-transformed_lymphocytes", "EBV lymphocytes", "from blood", "#F6C28B"),
    ("Whole_Blood", "Whole blood", "blood draw", "#C0392B"),
]
HEAT = LinearSegmentedColormap.from_list("deck", ["#FFFFFF", "#F3E79B", "#FAB27B", "#EB7F86", "#CE6693", "#A059A0",
                                                  "#5C53A5"])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--medians", required=True, help="public_data/gtex_v8_transcript_medians.tsv.gz")
    p.add_argument("--gtf", required=True, help="GENCODE GTF (plain or .gz): MANE Select and transcript names")
    p.add_argument("--genes", required=True, help="highlighted genes, comma-separated")
    p.add_argument("--min-tpm", type=float, default=5.0)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def transcripts(gtf, genes):
    """unversioned transcript id -> (name, is MANE Select) for the given gene names."""
    out = {}
    op = gzip.open if gtf.endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "transcript":
                continue
            a = f[8]
            gname = a.split('gene_name "', 1)[1].split('"', 1)[0]
            if gname not in genes:
                continue
            tid = a.split('transcript_id "', 1)[1].split('"', 1)[0].split(".")[0]
            tname = a.split('transcript_name "', 1)[1].split('"', 1)[0] if 'transcript_name "' in a else tid
            out[tid] = (tname, 'tag "MANE_Select"' in a)
    return out


def classic(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def save(fig, path):
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=300)
    plt.close(fig)


def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    targets = a.genes.split(",")
    use = defaultdict(lambda: defaultdict(dict))   # gene -> tissue -> transcript -> TPM
    with gzip.open(a.medians, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            use[r["gene"]][r["tissue"]][r["transcript_id"].split(".")[0]] = float(r["median_tpm"])
    tx = transcripts(a.gtf, set(use))
    rows, out = [], []
    for g, by_t in use.items():
        ref = by_t.get(REF, {})
        rt = sum(ref.values())
        mane = [t for t in ref if tx.get(t, ("", False))[1]]
        for tid, lab, _, _ in TISSUES:
            q = by_t.get(tid, {})
            qt = sum(q.values())
            row = {"gene": g, "tissue": tid, "cortex_tpm": round(rt, 2), "tissue_tpm": round(qt, 2),
                   "overlap": "", "same_dominant": "", "mane_share_cortex": "", "mane_share_tissue": ""}
            if rt >= a.min_tpm and qt >= a.min_tpm:
                ts = set(ref) | set(q)
                pu = {t: ref.get(t, 0) / rt for t in ts}
                qu = {t: q.get(t, 0) / qt for t in ts}
                row["overlap"] = round(sum(min(pu[t], qu[t]) for t in ts), 3)
                row["same_dominant"] = max(pu, key=pu.get) == max(qu, key=qu.get)
                if mane:
                    row["mane_share_cortex"] = round(pu[mane[0]], 3)
                    row["mane_share_tissue"] = round(qu[mane[0]], 3)
            rows.append(row)
    cols = list(rows[0])
    with open(os.path.join(a.outdir, "isoform_similarity.tsv"), "w") as fh:
        fh.write("# Isoform usage of each panel gene in a tissue against brain cortex (GTEx v8 transcript medians); "
                 "genes with summed transcript TPM >= {:g} in both. Made by scripts/isoform_similarity.py.\n".format(
                     a.min_tpm))
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(os.path.join(a.outdir, "isoform_usage_targets.tsv"), "w") as fh:
        fh.write("gene\ttranscript_id\ttranscript_name\tmane_select\ttissue\tmedian_tpm\tusage\n")
        for g in targets:
            for tid in [REF] + [t[0] for t in TISSUES]:
                q = use[g].get(tid, {})
                tot = sum(q.values())
                for t, v in sorted(q.items(), key=lambda kv: -kv[1]):
                    name, mane = tx.get(t, (t, False))
                    fh.write("{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.4f}\n".format(g, t, name, mane, tid, v,
                                                                       v / tot if tot else 0))

    # ---- figure 1: overlap per tissue and the share with the same dominant isoform
    by = defaultdict(list)
    dom = defaultdict(list)
    for r in rows:
        if r["overlap"] != "":
            by[r["tissue"]].append(r["overlap"])
            dom[r["tissue"]].append(r["same_dominant"])
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2), gridspec_kw={"width_ratios": [1.5, 1, 1]})
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.3, top=0.86, wspace=0.35)
    ax = axes[0]
    data = [by[t[0]] for t in TISSUES]
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color="black"))
    for patch, t in zip(bp["boxes"], TISSUES):
        patch.set_facecolor(t[3])
        patch.set_alpha(0.8)
    for i, d in enumerate(data):
        ax.text(i + 1, 1.03, "n={}".format(len(d)), ha="center", fontsize=7, color="#555555")
        out.append(("median_overlap", TISSUES[i][0], round(float(np.median(d)), 3) if d else ""))
        out.append(("genes_compared", TISSUES[i][0], len(d)))
    ax.set_xticks(range(1, len(TISSUES) + 1))
    ax.set_xticklabels(["{}\n({})".format(t[1], t[2]) for t in TISSUES], rotation=40, ha="right",
                       rotation_mode="anchor", fontsize=7.5)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Isoform overlap with brain cortex")
    ax.set_title("a  Isoform usage similarity, per panel gene", loc="left", fontsize=10, fontweight="bold")
    classic(ax)
    ax = axes[1]
    v = [np.mean(dom[t[0]]) * 100 if dom[t[0]] else 0 for t in TISSUES]
    ax.bar(range(len(TISSUES)), v, color=[t[3] for t in TISSUES], edgecolor="black", lw=0.3)
    for i, x in enumerate(v):
        ax.text(i, x + 1, "{:.0f}%".format(x), ha="center", fontsize=7.5)
        out.append(("percent_same_dominant", TISSUES[i][0], round(x, 1)))
    ax.set_xticks(range(len(TISSUES)))
    ax.set_xticklabels([t[1] for t in TISSUES], rotation=40, ha="right", rotation_mode="anchor", fontsize=7.5)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of genes")
    ax.set_title("b  Same dominant isoform as cortex", loc="left", fontsize=10, fontweight="bold")
    classic(ax)
    ax = axes[2]
    lo = []
    for t in TISSUES:
        d = np.array(by[t[0]])
        lo.append((d < 0.5).mean() * 100 if len(d) else 0)
    ax.bar(range(len(TISSUES)), lo, color=[t[3] for t in TISSUES], edgecolor="black", lw=0.3)
    for i, x in enumerate(lo):
        ax.text(i, x + 1, "{:.0f}%".format(x), ha="center", fontsize=7.5)
        out.append(("percent_overlap_below_0.5", TISSUES[i][0], round(x, 1)))
    ax.set_xticks(range(len(TISSUES)))
    ax.set_xticklabels([t[1] for t in TISSUES], rotation=40, ha="right", rotation_mode="anchor", fontsize=7.5)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of genes")
    ax.set_title("c  Expressed, but overlap < 0.5", loc="left", fontsize=10, fontweight="bold")
    classic(ax)
    fig.text(0.06, 0.01, "GTEx v8 transcript medians; genes with summed transcript TPM >= {:g} in cortex and the "
                         "tissue. Overlap = sum of min(usage in cortex, usage in tissue).".format(a.min_tpm),
             fontsize=7.5, color="#555555")
    save(fig, os.path.join(a.outdir, "isoform_similarity"))

    # ---- figure 2: target genes, overlap per tissue; isoform usage of the first two genes
    tg = [g for g in targets if g in use]
    look = {(r["gene"], r["tissue"]): r for r in rows}
    M = np.full((len(tg), len(TISSUES)), np.nan)
    for i, g in enumerate(tg):
        for j, t in enumerate(TISSUES):
            r = look.get((g, t[0]))
            if r and r["overlap"] != "":
                M[i, j] = r["overlap"]
            out.append(("target_overlap", "{}|{}".format(g, t[0]), "" if np.isnan(M[i, j]) else M[i, j]))
    stack = [g for g in tg[:2]]
    fig = plt.figure(figsize=(13.0, 5.6))
    ax = fig.add_axes([0.07, 0.38, 0.27, 0.54])
    im = ax.imshow(np.nan_to_num(M, nan=0), aspect="auto", cmap=HEAT, vmin=0, vmax=1)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, "-" if np.isnan(M[i, j]) else "{:.2f}".format(M[i, j]), ha="center", va="center",
                    fontsize=7, color="white" if not np.isnan(M[i, j]) and M[i, j] > 0.7 else "black")
    ax.set_xticks(range(len(TISSUES)))
    ax.set_xticklabels([t[1] for t in TISSUES], rotation=45, ha="right", rotation_mode="anchor", fontsize=7.5)
    ax.set_yticks(range(len(tg)))
    ax.set_yticklabels(tg, fontstyle="italic")
    ax.set_title("a  Isoform overlap with brain cortex", loc="left", fontsize=10, fontweight="bold")
    pal = ["#2166AC", "#67A9CF", "#D1E5F0", "#F4A582", "#BDBDBD"]
    cols = [(REF, "Brain cortex")] + [(t[0], t[1]) for t in TISSUES]
    for k, g in enumerate(stack):
        ax = fig.add_axes([0.42 + k * 0.29, 0.38, 0.24, 0.54])
        ref = use[g].get(REF, {})
        top = [t for t, _ in sorted(ref.items(), key=lambda kv: -kv[1])[:4]]
        for j, (tid, lab) in enumerate(cols):
            q = use[g].get(tid, {})
            tot = sum(q.values())
            b = 0
            for n, t in enumerate(top + ["other"]):
                v = (q.get(t, 0) if t != "other" else sum(x for tt, x in q.items() if tt not in top)) / tot \
                    if tot else 0
                ax.bar(j, v, bottom=b, color=pal[n], edgecolor="white", lw=0.3, width=0.8,
                       alpha=1.0 if tot >= a.min_tpm else 0.3,
                       label=(tx.get(t, (t, False))[0] + (" (MANE)" if tx.get(t, ("", False))[1] else "")
                              if t != "other" else "Other") if j == 0 else None)
                b += v
            ax.text(j, 1.02, "{:.0f}".format(tot), ha="center", fontsize=6.5, color="#555555")
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([c[1] for c in cols], rotation=45, ha="right", rotation_mode="anchor", fontsize=7.5)
        ax.set_ylim(0, 1.1)
        ax.set_ylabel("Isoform usage" if k == 0 else "")
        ax.set_title("{}  {} isoforms (TPM above bars)".format("bc"[k], g), loc="left", fontsize=10,
                     fontweight="bold", fontstyle="normal")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.33), ncol=3, frameon=False, fontsize=7,
                  handlelength=1.2, columnspacing=0.8)
        classic(ax)
    fig.text(0.07, 0.008, "GTEx v8 (GENCODE v26 transcripts, RSEM on short reads). -: below {:g} TPM in cortex or the "
                         "tissue (faded bars). Isoforms: the four most used in cortex.".format(a.min_tpm), fontsize=7.5,
             color="#555555")
    save(fig, os.path.join(a.outdir, "isoform_targets"))

    with open(os.path.join(a.outdir, "isoform_similarity_figures.tsv"), "w") as fh:
        fh.write("measure\tkey\tvalue\n")
        for r in out:
            fh.write("{}\t{}\t{}\n".format(*r))
    for r in out:
        if r[0] in ("median_overlap", "percent_same_dominant", "percent_overlap_below_0.5", "genes_compared"):
            print("  ", r)


if __name__ == "__main__":
    main()
