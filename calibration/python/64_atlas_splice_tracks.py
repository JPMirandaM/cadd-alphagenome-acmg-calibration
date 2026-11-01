"""Paper 2, amendment A-P2-4 (clinical outcomes): Atlas splicing track scores for the expert-panel RNA-study
candidates (splicing-stratum current VUS with at least PP3 Moderate; list written by script 63).

For each variant and each scorer (SPLICE_SITE_USAGE, SPLICE_JUNCTIONS) reports the track with the largest absolute
raw score: raw score, track name, biosample, ontology CURIE and gene, plus the number of tracks returned. The API
key is read only from ALPHAGENOME_API_KEY. Outputs are under AlphaGenome non-commercial terms and are not versioned.

Usage: uv run --python 3.12 --with alphagenome python -I 64_atlas_splice_tracks.py IN.tsv OUT.tsv [--workers 8]
"""
import argparse, csv, os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from alphagenome.atlas import atlas
from alphagenome.data import genome

SCORERS = ["SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
p = argparse.ArgumentParser()
p.add_argument("input"); p.add_argument("output"); p.add_argument("--workers", type=int, default=8)
a = p.parse_args()
with open(a.input, newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
client = atlas.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)


def col(df, names, i):
    for n in names:
        if n in df.columns:
            return str(df[n].iloc[i])
    return ""


def fetch(r):
    out = {}
    for attempt in range(3):
        try:
            s = client.query_variant(genome.Variant.from_str(r["variant"]), requested_scorers=SCORERS)
            for sc in SCORERS:
                ad = s.get(sc)
                if ad is None or ad.X is None or ad.X.size == 0:
                    out[sc] = ["", "", "", "", "", 0]
                    continue
                x = np.asarray(ad.X, dtype=float)
                i, j = np.unravel_index(np.nanargmax(np.abs(x)), x.shape)
                out[sc] = [f"{x[i, j]:.4g}", col(ad.var, ["name", "track_name"], j), col(ad.var, ["biosample_name"], j),
                           col(ad.var, ["ontology_curie"], j), col(ad.obs, ["gene_name", "gene_id"], i), int(x.shape[1])]
            return r, out, ""
        except Exception as e:  # recorded per variant, never fatal
            err = f"{type(e).__name__}: {str(e)[:200]}"
    return r, None, err


with ThreadPoolExecutor(max_workers=a.workers) as ex:
    res = list(ex.map(fetch, rows))
with open(a.output, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    hdr = list(rows[0].keys()) if rows else ["variant"]
    w.writerow(hdr + [f"{sc}_{c}" for sc in SCORERS for c in ["max_abs_raw", "track", "biosample", "ontology", "gene", "n_tracks"]] + ["error"])
    for r, out, err in res:
        vals = [v for sc in SCORERS for v in (out[sc] if out else [""] * 6)]
        w.writerow([r[h] for h in hdr] + vals + [err])
print("done", len(res), "variants;", sum(1 for _, o, _ in res if o is None), "errors")
