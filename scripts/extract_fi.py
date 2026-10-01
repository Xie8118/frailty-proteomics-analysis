# ============================================================================
# Extract FI summary statistics around candidate genes
# ============================================================================
# Purpose: Extract FI summary statistics around candidate genes.
# Inputs: Candidate gene boundaries and compressed FI GWAS.
# Outputs: fi_regions.tsv and extraction metadata.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import csv, gzip, json, hashlib, time

root = Path(".")
w = root / "work"
regions = {}
for r in csv.DictReader(
    (w / "inventory_v2/fenland_candidates.tsv").open(), delimiter="\t"
):
    regions.setdefault(r["Chr"], []).append(
        (int(r["start"]) - 3000000, int(r["end"]) + 3000000)
    )
src = Path("data/fi/fi_summary_GRCh37.txt.gz")
out = w / "fi_regions.tsv"
n = 0
k = 0
t = time.time()
with gzip.open(src, "rt") as f, out.open("x") as g:
    h = f.readline().split()
    ix = {x: i for i, x in enumerate(h)}
    writer = csv.writer(g, delimiter="\t")
    writer.writerow(h)
    for line in f:
        n += 1
        a = line.split()
        if len(a) != len(h):
            raise ValueError(n)
        c = a[ix["CHR"]]
        p = int(a[ix["POS"]])
        if c in regions and any(lo <= p <= hi for lo, hi in regions[c]):
            writer.writerow(a)
            k += 1
(w / "fi_regions.json").write_text(
    json.dumps(
        dict(
            rows=n,
            selected=k,
            seconds=time.time() - t,
            sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
        ),
        indent=2,
    )
)
print("FI_EXTRACT_COMPLETE", n, k)
