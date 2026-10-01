# ============================================================================
# Estimate cis-instrument Mendelian randomization for one assay
# ============================================================================
# Purpose: Estimate cis-instrument Mendelian randomization for one assay.
# Inputs: full_mr_v4_tasks.json, pQTL input, FI SQLite index and LD reference.
# Outputs: Instrument table, effect estimates, exclusions and provenance.
# Arguments: one-based task ID and mode (pilot or full).
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import csv, gzip, json, math, sys, sqlite3, subprocess, hashlib, time

w = Path("work")
tid = int(sys.argv[1])
mode = sys.argv[2]
t = json.loads((w / "full_mr_v4_tasks.json").read_text())[tid - 1]
out = w / ("full_mr_pilot" if mode == "pilot" else "full_mr") / f"task_{tid}"
out.mkdir(parents=True, exist_ok=False)
src = Path(t["input_file"])
start = time.time()
rows = []
n = 0
invalid_numeric = 0

# Read source statistics and apply instrument significance and strength filters.
with gzip.open(src, "rt") as f:
    h = f.readline().split()
    ix = {a: i for i, a in enumerate(h)}
    fen = t["platform"] == "Fenland"
    for line in f:
        n += 1
        a = line.split()
        pk = "Pvalue" if fen else "Pval"
        try:
            p = float(a[ix[pk]])
        except ValueError:
            continue
        if not (0 <= p < 5e-8):
            continue
        ck = "chr" if fen else "Chrom"
        posk = "pos" if fen else "Pos"
        chrom = a[ix[ck]].replace("chr", "")
        pos = int(a[ix[posk]])
        if (
            chrom != t["chr"]
            or not t["start_native"] - 2000000 <= pos <= t["end_native"] + 2000000
        ):
            continue

        def v(k):
            return a[ix[k]]

        try:
            r = dict(
                SNP=v("rsid" if fen else "rsids"),
                CHR=chrom,
                POS=pos,
                EA=v("Allele1" if fen else "effectAllele").upper(),
                OA=v("Allele2" if fen else "otherAllele").upper(),
                BETA=float(v("Effect" if fen else "Beta")),
                SE=float(v("StdErr" if fen else "SE")),
                P=p,
                FREQ=float(v("Freq1" if fen else "ImpMAF")),
                N=float(v("TotalSampleSize" if fen else "N")),
            )
        except (ValueError, IndexError):
            invalid_numeric += 1
            continue
        if not all(math.isfinite(r[k]) for k in ["BETA", "SE", "P", "FREQ", "N"]):
            continue
        if (
            not r["SNP"].startswith("rs")
            or r["SE"] <= 0
            or not 0 < r["FREQ"] < 1
            or r["N"] <= 2
        ):
            continue
        if (r["BETA"] / r["SE"]) ** 2 <= 10:
            continue
        rows.append(r)

# Convert deCODE coordinates to GRCh37 before matching the FI outcome.
if not fen and rows:
    bed = out / "variants.hg38.bed"
    mapped = out / "variants.hg19.bed"
    unmapped = out / "unmapped.bed"
    with bed.open("w") as f:
        for i, r in enumerate(rows):
            f.write(f"chr{r['CHR']}\t{r['POS']-1}\t{r['POS']}\t{i}\n")
    subprocess.run(
        [
            "bin/liftOver",
            str(bed),
            "data/reference/UCSC/liftOver/hg38ToHg19/hg38ToHg19.over.chain.gz",
            str(mapped),
            str(unmapped),
        ],
        check=True,
        stdout=(out / "lift.stdout").open("w"),
        stderr=(out / "lift.stderr").open("w"),
    )
    m = {}
    for line in mapped.open():
        a = line.split()
        m[int(a[3])] = (a[0].replace("chr", ""), int(a[1]) + 1)
    rows = [dict(r, CHR=m[i][0], POS=m[i][1]) for i, r in enumerate(rows) if i in m]
rows = [
    r
    for r in rows
    if r["CHR"] == t["chr"]
    and t["start37"] - 1000000 <= r["POS"] <= t["end37"] + 1000000
    and not (r["CHR"] == "6" and 28477797 <= r["POS"] <= 33448354)
]
# Duplicate variant IDs are excluded rather than silently selecting an alternate allele.
counts = {}
for r in rows:
    counts[r["SNP"]] = counts.get(r["SNP"], 0) + 1
rows = [r for r in rows if counts[r["SNP"]] == 1]
fields = ["SNP", "CHR", "POS", "EA", "OA", "BETA", "SE", "P", "FREQ", "N"]
with (out / "preclump.tsv").open("w") as f:
    d = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    d.writeheader()
    d.writerows(rows)

# LD-clump eligible instruments using the fixed reference and thresholds.
leads = set()
if rows:
    cmd = [
        "bin/plink",
        "--bfile",
        "data/reference/EUR_1000Genome_phase3_all",
        "--clump",
        str(out / "preclump.tsv"),
        "--clump-snp-field",
        "SNP",
        "--clump-field",
        "P",
        "--clump-p1",
        "5e-8",
        "--clump-p2",
        "1",
        "--clump-r2",
        "0.01",
        "--clump-kb",
        "10000",
        "--threads",
        "1",
        "--memory",
        "4000",
        "--out",
        str(out / "clump"),
    ]
    subprocess.run(
        cmd,
        check=True,
        stdout=(out / "plink.stdout").open("w"),
        stderr=(out / "plink.stderr").open("w"),
    )
    if (out / "clump.clumped").exists():
        for line in (out / "clump.clumped").read_text().splitlines()[1:]:
            a = line.split()
            if a:
                leads.add(a[2])

# Match outcome alleles; unresolved palindromic variants are excluded.
con = sqlite3.connect(f"file:{w}/fi_lookup.sqlite?mode=ro", uri=True)
kept = []
excluded = []
complement = str.maketrans("ACGT", "TGCA")
for r in rows:
    if r["SNP"] not in leads:
        continue
    ys = con.execute(
        "SELECT chr,pos,ea,oa,beta,se,p,n,freq FROM fi WHERE snp=?", (r["SNP"],)
    ).fetchall()
    matches = []
    for ch, pos, ea, oa, b, se, p, ns, fr in ys:
        ea = ea.upper()
        oa = oa.upper()
        a = r["EA"]
        o = r["OA"]
        same = (a == ea and o == oa) or (
            len(a) == len(o) == 1
            and a.translate(complement) == ea
            and o.translate(complement) == oa
        )
        flip = (a == oa and o == ea) or (
            len(a) == len(o) == 1
            and a.translate(complement) == oa
            and o.translate(complement) == ea
        )
        if (
            str(ch) != r["CHR"]
            or pos != r["POS"]
            or not (same or flip)
            or not set(a + o) <= set("ACGT")
        ):
            continue
        if (
            not all(math.isfinite(x) for x in [b, se, p, ns, fr])
            or se <= 0
            or ns <= 2
            or not 0 < fr < 1
        ):
            continue
        sign = 1 if same else -1
        pal = {a, o} in [{"A", "T"}, {"C", "G"}]
        if pal:
            if (
                not fen
                or min(r["FREQ"], 1 - r["FREQ"]) > 0.42
                or min(fr, 1 - fr) > 0.42
            ):
                continue
            ds = abs(r["FREQ"] - fr)
            df = abs(r["FREQ"] - (1 - fr))
            if min(ds, df) > 0.1:
                continue
            sign = 1 if ds < df else -1
        matches.append(
            dict(
                r,
                BETA_Y=b * sign,
                SE_Y=se,
                P_Y=p,
                N_Y=ns,
                EAF_Y=fr,
                F=(r["BETA"] / r["SE"]) ** 2,
                R2_X=r["BETA"] ** 2 / (r["BETA"] ** 2 + r["N"] * r["SE"] ** 2),
                R2_Y=b * b / (b * b + ns * se * se),
            )
        )
    if len(matches) == 1:
        kept.append(matches[0])
    else:
        excluded.append(
            dict(
                SNP=r["SNP"],
                reason="missing_ambiguous_or_unresolved_alleles",
                n_outcome_records=len(ys),
            )
        )
con.close()
with (out / "instruments.tsv").open("w") as f:
    d = csv.DictWriter(
        f,
        fieldnames=fields
        + ["BETA_Y", "SE_Y", "P_Y", "N_Y", "EAF_Y", "F", "R2_X", "R2_Y"],
        delimiter="\t",
    )
    d.writeheader()
    d.writerows(kept)
r = dict(
    task_id=tid,
    platform=t["platform"],
    gene=t["gene"],
    assay=t["assay"],
    mode=mode,
    preclump=len(rows),
    clumped=len(leads),
    nsnp=len(kept),
    status="not_estimable_after_filters",
)
if kept:
    denom = sum(a["BETA"] ** 2 / a["SE_Y"] ** 2 for a in kept)
    beta = sum(a["BETA"] * a["BETA_Y"] / a["SE_Y"] ** 2 for a in kept) / denom
    se = math.sqrt(1 / denom)
    q = sum((a["BETA_Y"] - beta * a["BETA"]) ** 2 / a["SE_Y"] ** 2 for a in kept)
    ser = se * math.sqrt(max(1, q / (len(kept) - 1))) if len(kept) > 1 else se
    r.update(
        status="estimable",
        beta=beta,
        se=se,
        p=math.erfc(abs(beta / se) / math.sqrt(2)),
        lower=beta - 1.96 * se,
        upper=beta + 1.96 * se,
        random_se=ser,
        random_p=math.erfc(abs(beta / ser) / math.sqrt(2)),
        Q=q,
        Fmin=min(a["F"] for a in kept),
        N_X_min=min(a["N"] for a in kept),
        N_X_max=max(a["N"] for a in kept),
        PVE_X_approx=sum(a["R2_X"] for a in kept),
        PVE_Y_approx=sum(a["R2_Y"] for a in kept),
    )
    r["steiger_direction_approx"] = r["PVE_X_approx"] > r["PVE_Y_approx"]
h = hashlib.sha256()
with src.open("rb") as f:
    for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
        h.update(block)
r.update(
    source_sha256=h.hexdigest(),
    seconds=time.time() - start,
    rows_scanned=n,
    invalid_numeric=invalid_numeric,
    original_source=str(src),
    excluded_leads=excluded,
)
(out / "result.json").write_text(json.dumps(r, indent=2, allow_nan=False))
(out / "done").write_text("MR_TASK_COMPLETE\n")
print(
    json.dumps(
        {k: r[k] for k in ["task_id", "gene", "platform", "nsnp", "status", "seconds"]}
    )
)
