"""Paper 2: retrieves AVI scores and the 18 SHAP feature attributions from the AlphaGenome Atlas API.

Input: a TSV with a column `variant` (GRCh38, 1-based, chr:pos:ref>alt). Output: gzipped TSV with
avi_raw, avi_cdf, the 18 attributions (fi_<name>) and an error column, one row per variant, plus a
one-line timing summary appended to <output>.timing.tsv. The API key is read only from the
environment variable ALPHAGENOME_API_KEY. Outputs are under AlphaGenome non-commercial terms and are
not versioned (see .gitignore).

Usage: uv run --python 3.12 --with alphagenome python -I 51_atlas_attributions.py IN.tsv OUT.tsv.gz
       [--workers 8] [--first N] [--skip N] [--label NAME]
"""
import argparse, csv, gzip, os, sys, time
from concurrent.futures import ThreadPoolExecutor

from alphagenome.atlas import atlas
from alphagenome.data import genome

p = argparse.ArgumentParser()
p.add_argument("input"); p.add_argument("output")
p.add_argument("--workers", type=int, default=8)
p.add_argument("--skip", type=int, default=0)
p.add_argument("--first", type=int, default=0)
p.add_argument("--label", default="")
a = p.parse_args()

with open(a.input, newline="") as f:
    variants = [r["variant"] for r in csv.DictReader(f, delimiter="\t")]
variants = variants[a.skip:]
if a.first:
    variants = variants[: a.first]

client = atlas.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)
SCORERS = ["AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE"]
names = None


def fetch(v):
    try:
        s = client.query_variant(genome.Variant.from_str(v), requested_scorers=SCORERS)
        avi, fi = s["AVI_SCORE"], s["AVI_SCORE_FEATURE_IMPORTANCE"]
        returned = str(avi.obs["variant"].iloc[0]) if "variant" in avi.obs else v
        if returned.replace("chr", "") != v.replace("chr", ""):
            return v, None, None, None, f"returned variant {returned}"
        return v, float(avi.X[0, 0]), float(avi.layers["quantiles"][0, 0]), (list(fi.var["name"]), [float(x) for x in fi.X[0]]), ""
    except Exception as e:  # recorded per variant, never fatal
        return v, None, None, None, f"{type(e).__name__}: {str(e)[:200]}"


t0 = time.time()
with ThreadPoolExecutor(max_workers=a.workers) as ex:
    rows = list(ex.map(fetch, variants))
elapsed = time.time() - t0

for r in rows:
    if r[3] is not None:
        names = r[3][0]; break
names = names or []
mode = "at" if os.path.exists(a.output) else "wt"
with gzip.open(a.output, mode, newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    if mode == "wt":
        w.writerow(["variant", "avi_raw", "avi_cdf"] + [f"fi_{n}" for n in names] + ["error"])
    for v, raw, cdf, fi, err in rows:
        vals = fi[1] if fi else [""] * len(names)
        w.writerow([v, "" if raw is None else raw, "" if cdf is None else cdf] + vals + [err])

n_err = sum(1 for r in rows if r[4])
timing = a.output + ".timing.tsv"
new = not os.path.exists(timing)
with open(timing, "a") as f:
    if new:
        f.write("label\tworkers\tvariants\terrors\tseconds\tvariants_per_second\n")
    f.write(f"{a.label}\t{a.workers}\t{len(rows)}\t{n_err}\t{elapsed:.1f}\t{len(rows) / elapsed:.2f}\n")
print(f"{a.label}: {len(rows)} variants, {n_err} errors, {elapsed:.1f} s, {len(rows) / elapsed:.2f} variants/s with {a.workers} workers")
if n_err:
    print("first errors:", [r[4] for r in rows if r[4]][:3])
