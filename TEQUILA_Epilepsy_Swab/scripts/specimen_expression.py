#!/usr/bin/env python3
"""Expression of the probe-panel genes measured in real accessible specimens, from open public RNA-seq (GEO).

Each dataset is reduced to its healthy (or control) samples, turned into TPM the same way (counts divided by the
length of the gene's merged GENCODE exons), scaled over protein-coding genes only (so that small RNAs such as 7SL in
total-RNA libraries do not shrink every other gene), and summarised as the median TPM per gene. Datasets (files from https://ftp.ncbi.nlm.nih.gov/geo/series/<GSEnnn>/<GSE>/suppl/ and the GEO series
matrix for the sample annotation):
  buccal     GSE330022  buccal brushings (soft toothbrush, both cheeks), raw counts; group "Non-user" (n = 24).
                        George J, ... Besaratinia A, Front Oncol 2026 (PMID 42306798)
  nasal_kids GSE118761  nasal epithelial brushings of children, raw counts; nasal, "Non-atopic non-asthmatic".
                        Kicic A, ... AusREC, J Allergy Clin Immunol 2020 (PMID 32113981)
  hair       GSE224010  plucked scalp hair follicles, raw counts; Hair1-5 (n = 5).
                        Takeuchi S, ... Akiyama M, J Dermatol Sci 2023 (PMID 37355462)
  blood      GSE112087  PAXgene whole blood, raw counts (Ensembl 70), two lanes per sample summed; "Healthy" (n = 29).
                        Oon S, ... Wilson NJ, Clin Transl Immunology 2019 (PMID 31890206)
  urine      GSE235813  cultured urine-derived progenitor renal cells, TPM; urine samples (n = 37). Donors are
                        kidney-transplant recipients, not healthy. Cruzado JM, ... Berdasco M, Am J Transplant 2025
                        (PMID 39521058)
Saliva: no open gene-level matrix with annotated healthy samples was usable; it stays represented by GTEx minor
salivary gland. Output (committed under public_data/): specimen_panel_tpm.tsv (panel genes x specimens, median TPM)
and specimen_datasets.tsv (dataset, samples used).
"""
import argparse
import gzip
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd

ALIAS = {"AFG2A": ["SPATA5"], "AARS1": ["AARS"], "EPRS1": ["EPRS"], "QARS1": ["QARS"], "H1-4": ["HIST1H1E"],
         "H3-3A": ["H3F3A"], "NEXMIF": ["KIAA2022"], "LINS1": ["LINS"], "PLPBP": ["PROSC"], "KMT5B": ["SUV420H1"],
         "NSD2": ["WHSC1"]}


CODING = set()


def gencode(gtf):
    """merged exon length per gene name and per unversioned gene id; id -> name."""
    iv, id2name = defaultdict(list), {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "exon":
                continue
            gid = f[8].split('gene_id "', 1)[1].split('"', 1)[0].split(".")[0]
            name = f[8].split('gene_name "', 1)[1].split('"', 1)[0]
            id2name[gid] = name
            if 'gene_type "protein_coding"' in f[8]:
                CODING.add(name)
            iv[gid].append((int(f[3]) - 1, int(f[4])))
    by_id = {}
    for gid, v in iv.items():
        v.sort()
        tot, cs, ce = 0, -1, -1
        for s, e in v:
            if s > ce:
                tot += max(0, ce - cs)
                cs, ce = s, e
            else:
                ce = max(ce, e)
        by_id[gid] = tot + max(0, ce - cs)
    by_name = {}
    for gid, n in by_id.items():
        nm = id2name[gid]
        by_name[nm] = max(by_name.get(nm, 0), n)
    return by_id, by_name, id2name


def series_annotation(path):
    """GEO series matrix -> {field: [values per sample]} (characteristics keyed by their 'key:' prefix)."""
    ann = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            f = line.rstrip("\n").split("\t")
            vals = [x.strip('"') for x in f[1:]]
            key = f[0][8:]
            if key == "characteristics_ch1" and vals and ":" in vals[0]:
                key = vals[0].split(":", 1)[0].strip()
                vals = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
            ann.setdefault(key, vals)
    return ann


def to_tpm(counts, lengths):
    """counts: DataFrame genes x samples (index = gene symbol); lengths: Series by symbol."""
    lens = lengths.reindex(counts.index)
    ok = lens.notna() & (lens > 0) & counts.index.isin(CODING)
    rpk = counts[ok].div(lens[ok] / 1000, axis=0)
    return rpk.div(rpk.sum(axis=0), axis=1) * 1e6


def collapse(df):
    df = df[df.index.notna()]
    return df.groupby(level=0).sum()


def load(name, d, lens_name, lens_id, id2name):
    if name == "buccal":
        ann = series_annotation(os.path.join(d, "GSE330022_series_matrix.txt.gz"))
        keep = [t for t, g in zip(ann["title"], ann["group"]) if g == "Non-user"]
        x = pd.read_excel(os.path.join(d, "GSE330022_RawCounts_AllSamples.xlsx"))
        x = x.set_index("gene_name")[keep].apply(pd.to_numeric, errors="coerce").fillna(0)
        return to_tpm(collapse(x), lens_name), keep
    if name == "nasal_kids":
        ann = series_annotation(os.path.join(d, "GSE118761_series_matrix.txt.gz"))
        ids = [v.replace("Sample name:", "").strip() for v in ann["description"]]
        keep = [i for i, t, dx in zip(ids, ann["tissue"], ann["diagnosis"])
                if t == "Nasal epithelium" and dx == "Non-atopic non-asthmatic"]
        x = pd.read_csv(os.path.join(d, "GSE118761_genecounts.csv.gz"))
        sym = x.columns[-2]
        x = x.set_index(sym)[keep]
        return to_tpm(collapse(x), lens_name), keep
    if name == "hair":
        x = pd.read_excel(os.path.join(d, "GSE224010_DESeq2_all_result.xlsx"))
        x = x.set_index(x.columns[0])
        keep = [c for c in x.columns if re.match(r"Hair\d+\.counts$", str(c))]
        return to_tpm(collapse(x[keep]), lens_name), keep
    if name == "blood":
        ann = series_annotation(os.path.join(d, "GSE112087_series_matrix.txt.gz"))
        lanes = [t.split(" ", 1)[1] for t, ph in zip(ann["title"], ann["phenotype"]) if ph == "Healthy"]
        x = pd.read_csv(os.path.join(d, "GSE112087_counts-matrix-EnsembIDs-GRCh37.p10.txt.gz"), sep="\t", index_col=0)
        x = x[[c for c in lanes if c in x.columns]]
        x = x.T.groupby(lambda c: c.rsplit("_", 3)[0]).sum().T      # sum the two lanes of a sample
        keep = list(x.columns)
        x.index = [id2name.get(i.split(".")[0]) for i in x.index]
        return to_tpm(collapse(x), lens_name), keep
    if name == "urine":
        x = pd.read_csv(os.path.join(d, "GSE235813_tpm_genes.gencV35.csv.gz"), index_col=0)
        keep = [c for c in x.columns if re.fullmatch(r"\d{3}", str(c))]
        x = collapse(x[keep])
        x = x[x.index.isin(CODING)]
        return x.div(x.sum(axis=0), axis=1) * 1e6, keep
    raise ValueError(name)


LABELS = {"buccal": "Buccal brushing", "nasal_kids": "Nasal brushing, children", "hair": "Plucked hair follicles",
          "blood": "Whole blood (PAXgene)", "urine": "Urine-derived cells (cultured)"}
SOURCE = {"buccal": "GSE330022", "nasal_kids": "GSE118761", "hair": "GSE224010", "blood": "GSE112087",
          "urine": "GSE235813"}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", required=True, help="folder with the GEO files named as above (*_series_matrix.txt.gz)")
    p.add_argument("--gtf", required=True, help="GENCODE annotation (gz)")
    p.add_argument("--panel", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    lens_id, lens_name, id2name = gencode(a.gtf)
    lens_name = pd.Series(lens_name)
    panel = [l.strip() for l in open(a.panel) if l.strip() and not l.startswith("#")]
    med, info = {}, []
    for name in LABELS:
        tpm, keep = load(name, a.data_dir, lens_name, lens_id, id2name)
        m = tpm.median(axis=1)
        med[name] = m
        info.append((name, LABELS[name], SOURCE[name], len(keep), int((m.reindex(panel).notna()).sum())))
        print("{}: {} samples, {} genes, {} panel genes".format(name, len(keep), tpm.shape[0], info[-1][-1]))
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "specimen_panel_tpm.tsv"), "w") as fh:
        fh.write("# Median TPM of each panel gene in open public RNA-seq of accessible specimens (healthy or control "
                 "samples; see specimen_datasets.tsv). Made by scripts/specimen_expression.py.\n")
        fh.write("gene\t" + "\t".join(LABELS[n] for n in LABELS) + "\n")
        for g in panel:
            vals = []
            for n in LABELS:
                v = next((med[n][x] for x in [g] + ALIAS.get(g, []) if x in med[n].index), None)
                vals.append("NA" if v is None or pd.isna(v) else "{:.2f}".format(v))
            fh.write(g + "\t" + "\t".join(vals) + "\n")
    with open(os.path.join(a.outdir, "specimen_datasets.tsv"), "w") as fh:
        fh.write("dataset\tspecimen\tgeo\tsamples\tpanel_genes_found\n")
        for r in info:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
