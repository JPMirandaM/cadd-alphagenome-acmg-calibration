"""Paper 2: re-queries the variants whose Atlas lookup failed in chunks written by 55_atlas_batch.py.

For each chunk_*.tsv.gz in the given output directories, rows with a non-empty `error` other than a returned-
variant mismatch are queried again with a small pool and backoff on RESOURCE_EXHAUSTED (quota); recovered rows
replace the failed ones and the chunk is rewritten atomically. Rows that still fail keep their error.
A summary per directory is appended to <dir>/repair.tsv. The API key is read only from ALPHAGENOME_API_KEY.

Usage: uv run --python 3.12 --with alphagenome python -I 56_atlas_repair.py DIR [DIR ...] [--workers 8]
"""
import argparse, csv, glob, gzip, os, time
from concurrent.futures import ThreadPoolExecutor

from alphagenome.atlas import atlas
from alphagenome.data import genome

p = argparse.ArgumentParser()
p.add_argument("dirs", nargs="+")
p.add_argument("--workers", type=int, default=8)
p.add_argument("--attempts", type=int, default=6)
a = p.parse_args()
client = atlas.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)
SCORERS = ["AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE"]


def fetch(v, features):
    err = ""
    for k in range(a.attempts):
        try:
            s = client.query_variant(genome.Variant.from_str(v), requested_scorers=SCORERS)
            avi, fi = s["AVI_SCORE"], s["AVI_SCORE_FEATURE_IMPORTANCE"]
            returned = str(avi.obs["variant"].iloc[0]) if "variant" in avi.obs else v
            if returned.replace("chr", "") != v.replace("chr", ""):
                return v, None, f"returned variant {returned}"
            d = dict(zip(fi.var["name"], (float(x) for x in fi.X[0])))
            return v, [float(avi.X[0, 0]), float(avi.layers["quantiles"][0, 0])] + [d.get(n, "") for n in features], ""
        except Exception as e:
            err = f"{type(e).__name__}: {str(e)[:200]}"
            time.sleep(min(60, 5 * 2 ** k) if "RESOURCE_EXHAUSTED" in err else 2)
    return v, None, err


with ThreadPoolExecutor(max_workers=a.workers) as ex:
    for d in a.dirs:
        n_failed = n_fixed = 0
        for f in sorted(glob.glob(os.path.join(d, "chunk_*.tsv.gz"))):
            with gzip.open(f, "rt", newline="") as fh:
                rows = list(csv.reader(fh, delimiter="\t"))
            head, body = rows[0], rows[1:]
            features = [h[3:] for h in head if h.startswith("fi_")]
            bad = [i for i, r in enumerate(body) if r[-1] and not r[-1].startswith("returned")]
            if not bad:
                continue
            n_failed += len(bad)
            res = {v: (vals, err) for v, vals, err in ex.map(lambda v: fetch(v, features), [body[i][0] for i in bad])}
            for i in bad:
                vals, err = res[body[i][0]]
                if vals:
                    body[i] = [body[i][0]] + vals + [""]; n_fixed += 1
                else:
                    body[i][-1] = err
            tmp = f + ".part"
            with gzip.open(tmp, "wt", newline="") as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(head); w.writerows(body)
            os.replace(tmp, f)
        with open(os.path.join(d, "repair.tsv"), "a") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M')}\tfailed {n_failed}\trecovered {n_fixed}\n")
        print(d, "failed", n_failed, "recovered", n_fixed)
