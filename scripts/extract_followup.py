# ============================================================================
# Extract follow-up pQTL regions
# ============================================================================
# Purpose: Extract follow-up pQTL regions.
# Inputs: followup_tasks.json and task-specific compressed pQTL input.
# Outputs: Follow-up regional TSV and extraction metadata.
# Arguments: one-based task ID and mode (pilot or full).
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import csv, json, sys, gzip, hashlib, subprocess, time, os

root = Path(".")
w = root / "work"
tid = int(sys.argv[1])
mode = sys.argv[2]
t = json.loads((w / "followup_tasks.json").read_text())[tid - 1]
out = w / ("followup_extract_pilot" if mode == "pilot" else "followup_extracted")
out.mkdir(exist_ok=True)
target = out / f"task_{tid}.tsv"
assert not target.exists()
source = Path(t["input_file"])
start = time.time()
n = 0
rows = []
bad = 0
with gzip.open(source, "rt") as f:
    names = f.readline().split()
    ix = {k: i for i, k in enumerate(names)}
    platform = t["platform"]
    ck = (
        "chr"
        if platform == "Fenland"
        else ("Chrom" if platform == "deCODE" else "CHROM")
    )
    pk = (
        "pos" if platform == "Fenland" else ("Pos" if platform == "deCODE" else "POS19")
    )
    for line in f:
        n += 1
        a = line.split()
        if len(a) != len(names):
            bad += 1
            continue
        if a[ix[ck]].replace("chr", "") != t["chr"]:
            continue
        p = int(a[ix[pk]])
        if not (t["start"] <= p <= t["end"]):
            continue

        def g(k):
            return a[ix[k]]

        if platform == "Fenland":
            r = [
                g("rsid"),
                t["chr"],
                p,
                g("Allele1").upper(),
                g("Allele2").upper(),
                g("Effect"),
                g("StdErr"),
                g("Pvalue"),
                g("Freq1"),
                g("TotalSampleSize"),
            ]
        elif platform == "deCODE":
            r = [
                g("rsids"),
                t["chr"],
                p,
                g("effectAllele").upper(),
                g("otherAllele").upper(),
                g("Beta"),
                g("SE"),
                g("Pval"),
                g("ImpMAF"),
                g("N"),
            ]
        else:
            r = [
                g("rsid"),
                t["chr"],
                p,
                g("ALLELE1").upper(),
                g("ALLELE0").upper(),
                g("BETA"),
                g("SE"),
                g("P"),
                g("A1FREQ"),
                g("N"),
            ]
        rows.append(r)
assert bad == 0, (bad, n)
if t["platform"] == "deCODE":
    bed = out / f"task_{tid}.hg38.bed"
    mapped = out / f"task_{tid}.hg19.bed"
    unmapped = out / f"task_{tid}.unmapped.bed"
    with bed.open("x") as f:
        for i, r in enumerate(rows):
            f.write(f"chr{r[1]}\t{int(r[2])-1}\t{r[2]}\t{i}\n")
    cmd = [
        "bin/liftOver",
        str(bed),
        "data/reference/UCSC/liftOver/hg38ToHg19/hg38ToHg19.over.chain.gz",
        str(mapped),
        str(unmapped),
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    m = {}
    for line in mapped.open():
        a = line.split()
        m[int(a[3])] = (a[0].replace("chr", ""), int(a[1]) + 1)
    retained = []
    for i, r in enumerate(rows):
        if i in m:
            r[1], r[2] = m[i]
            retained.append(r)
    rows = retained
with target.open("x") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["SNP", "CHR", "POS", "EA", "OA", "BETA", "SE", "P", "FREQ", "N"])
    writer.writerows(rows)
h = hashlib.sha256()
with source.open("rb") as f:
    for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
        h.update(b)
meta = dict(
    task=t,
    mode=mode,
    rows_scanned=n,
    rows_selected=len(rows),
    bad_rows=bad,
    input_sha256=h.hexdigest(),
    output_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
    seconds=time.time() - start,
    frequency_semantics=(
        "minor_allele_frequency"
        if t["platform"] == "deCODE"
        else "effect_allele_frequency"
    ),
)
(out / f"task_{tid}.json").write_text(json.dumps(meta, indent=2))
(out / f"task_{tid}.done").write_text("EXTRACTION_COMPLETE\n")
print(json.dumps(meta))
