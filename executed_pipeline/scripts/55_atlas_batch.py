"""Paper 2: resumable batch retrieval of AVI scores and the 18 SHAP attributions from the AlphaGenome Atlas API.

Splits the input (TSV with a column `variant`, GRCh38 chr:pos:ref>alt) into chunks of --chunk variants and
writes one gzipped TSV per chunk to OUTDIR (chunk_000001.tsv.gz, ...). Chunks already written are skipped,
so an interrupted run resumes where it stopped. Columns: variant, avi_raw, avi_cdf, fi_<18 features>, error.
Per-variant failures are retried (--retries passes) and otherwise recorded in `error`, never fatal.
Progress is appended to OUTDIR/progress.tsv. The API key is read only from ALPHAGENOME_API_KEY.
Outputs are under AlphaGenome non-commercial terms and are not versioned (see .gitignore).

Usage: uv run --python 3.12 --with alphagenome python -I 55_atlas_batch.py IN.tsv OUTDIR [--workers 16]
       [--chunk 10000] [--retries 2]
"""
import argparse, csv, gzip, os, time
from concurrent.futures import ThreadPoolExecutor

from alphagenome.atlas import atlas
from alphagenome.data import genome

FEATURES = ["MAX_ABS_ATAC", "MAX_ABS_CAGE", "MAX_ABS_CHIP_HISTONE", "MAX_ABS_CHIP_TF", "MAX_ABS_CONTACT_MAPS",
            "MAX_ABS_DNASE", "MAX_ABS_POLYADENYLATION", "MAX_ABS_PROCAP", "MAX_ABS_RNA_SEQ", "MERGED_SPLICING",
            "ALPHAMISSENSE", "CACTUS_241_WAY", "PHASTCONS_470_WAY", "PROTEIN_TERMINATION", "START_LOST",
            "STOP_LOST", "IS_INSERTION", "IS_DELETION"]
SCORERS = ["AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE"]

p = argparse.ArgumentParser()
p.add_argument("input"); p.add_argument("outdir")
p.add_argument("--workers", type=int, default=16)
p.add_argument("--chunk", type=int, default=10000)
p.add_argument("--retries", type=int, default=2)
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)

with open(a.input, newline="") as f:
    variants = [r["variant"] for r in csv.DictReader(f, delimiter="\t")]
client = atlas.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)


def fetch(v):
    try:
        s = client.query_variant(genome.Variant.from_str(v), requested_scorers=SCORERS)
        avi, fi = s["AVI_SCORE"], s["AVI_SCORE_FEATURE_IMPORTANCE"]
        returned = str(avi.obs["variant"].iloc[0]) if "variant" in avi.obs else v
        if returned.replace("chr", "") != v.replace("chr", ""):
            return v, None, f"returned variant {returned}"
        d = dict(zip(fi.var["name"], (float(x) for x in fi.X[0])))
        return v, [float(avi.X[0, 0]), float(avi.layers["quantiles"][0, 0])] + [d.get(n, "") for n in FEATURES], ""
    except Exception as e:  # recorded per variant, never fatal
        return v, None, f"{type(e).__name__}: {str(e)[:200]}"


n_chunks = (len(variants) + a.chunk - 1) // a.chunk
prog = os.path.join(a.outdir, "progress.tsv")
if not os.path.exists(prog):
    with open(prog, "w") as f:
        f.write("chunk\tvariants\terrors\tseconds\tvariants_per_second\tfinished\n")
with ThreadPoolExecutor(max_workers=a.workers) as ex:
    for i in range(n_chunks):
        out = os.path.join(a.outdir, f"chunk_{i + 1:06d}.tsv.gz")
        if os.path.exists(out):
            continue
        batch = variants[i * a.chunk:(i + 1) * a.chunk]
        t0 = time.time()
        res = {v: (vals, err) for v, vals, err in ex.map(fetch, batch)}
        for _ in range(a.retries):
            failed = [v for v, (vals, err) in res.items() if vals is None and not err.startswith("returned")]
            if not failed:
                break
            time.sleep(5)
            res.update({v: (vals, err) for v, vals, err in ex.map(fetch, failed)})
        tmp = out + ".part"
        with gzip.open(tmp, "wt", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["variant", "avi_raw", "avi_cdf"] + [f"fi_{n}" for n in FEATURES] + ["error"])
            for v in batch:
                vals, err = res[v]
                w.writerow([v] + (vals if vals else [""] * (2 + len(FEATURES))) + [err])
        os.replace(tmp, out)
        dt = time.time() - t0
        n_err = sum(1 for vals, err in res.values() if vals is None)
        with open(prog, "a") as f:
            f.write(f"{i + 1}/{n_chunks}\t{len(batch)}\t{n_err}\t{dt:.1f}\t{len(batch) / dt:.1f}\t{time.strftime('%Y-%m-%d %H:%M')}\n")
print("done", n_chunks, "chunks")
