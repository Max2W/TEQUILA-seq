#!/usr/bin/env python3
"""Download GTEx v8 transcript- and junction-level medians for the panel genes (GTEx Portal API v2), for the
isoform-similarity and variant-testability analyses. Run once with internet; outputs are committed under public_data/.

  gtex_v8_transcript_medians.tsv.gz   gene, gencode id, transcript id, tissue, median TPM (all panel genes)
  raw/gtex_v8_junctions_<GENE>.tsv    junction (chrom_start_end, 1-based intron), tissue, median read count
                                      (--junction-genes only)
Tissues: brain cortex and cerebellum, and the GTEx tissues used for accessible specimens.
"""
import argparse
import csv
import gzip
import json
import os
import time
import urllib.parse
import urllib.request

API = "https://gtexportal.org/api/v2/"
TISSUES = ["Brain_Cortex", "Brain_Cerebellum", "Lung", "Esophagus_Mucosa", "Minor_Salivary_Gland",
           "Skin_Not_Sun_Exposed_Suprapubic", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes",
           "Whole_Blood"]
OLD = {"AFG2A": "SPATA5", "AARS1": "AARS", "EPRS1": "EPRS", "QARS1": "QARS", "H1-4": "HIST1H1E", "H3-3A": "H3F3A",
       "NEXMIF": "KIAA2022", "LINS1": "LINS", "PLPBP": "PROSC", "KMT5B": "SUV420H1", "NSD2": "WHSC1"}


def get(endpoint, params):
    url = API + endpoint + "?" + urllib.parse.urlencode(params, doseq=True)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            if attempt == 4:
                raise
            time.sleep(3 * (attempt + 1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel", required=True)
    p.add_argument("--junction-genes", default="", help="comma-separated genes for the junction table")
    p.add_argument("--outdir", required=True, help="public_data/")
    a = p.parse_args()
    panel = [l.strip() for l in open(a.panel) if l.strip() and not l.startswith("#")]
    query = {OLD.get(g, g): g for g in panel}
    ids = {}
    names = list(query)
    for i in range(0, len(names), 50):
        d = get("reference/gene", {"geneId": names[i:i + 50], "gencodeVersion": "v26",
                                   "genomeBuild": "GRCh38/hg38", "itemsPerPage": 500})
        for x in d["data"]:
            if x["geneSymbol"] in query:
                ids[query[x["geneSymbol"]]] = x["gencodeId"]
    print("{} of {} panel genes found in GENCODE v26".format(len(ids), len(panel)))
    rows = []
    genes = [g for g in panel if g in ids]
    for i in range(0, len(genes), 20):
        chunk = genes[i:i + 20]
        d = get("expression/medianTranscriptExpression", {"gencodeId": [ids[g] for g in chunk],
                                                          "datasetId": "gtex_v8", "tissueSiteDetailId": TISSUES,
                                                          "itemsPerPage": 100000})
        back = {ids[g]: g for g in chunk}
        for x in d["data"]:
            rows.append((back.get(x["gencodeId"], x["geneSymbol"]), x["gencodeId"], x["transcriptId"],
                         x["tissueSiteDetailId"], "{:.4g}".format(x["median"])))
        print("  transcripts: {} / {} genes".format(min(i + 20, len(genes)), len(genes)), flush=True)
    with gzip.open(os.path.join(a.outdir, "gtex_v8_transcript_medians.tsv.gz"), "wt") as fh:
        fh.write("gene\tgencode_id\ttranscript_id\ttissue\tmedian_tpm\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    for g in [x for x in a.junction_genes.split(",") if x]:
        d = get("expression/medianJunctionExpression", {"gencodeId": ids[g], "datasetId": "gtex_v8",
                                                        "tissueSiteDetailId": TISSUES, "itemsPerPage": 100000})
        with open(os.path.join(a.outdir, "raw", "gtex_v8_junctions_{}.tsv".format(g)), "w") as fh:
            fh.write("# GTEx v8 median junction read counts per tissue (Portal API v2, medianJunctionExpression), "
                     "{} {}\njunction\ttissue\tmedian_reads\n".format(g, ids[g]))
            for x in d["data"]:
                fh.write("{}\t{}\t{:.4g}\n".format(x["junctionId"], x["tissueSiteDetailId"], x["median"]))
        print("  junctions: {} ({} rows)".format(g, len(d["data"])))


if __name__ == "__main__":
    main()
