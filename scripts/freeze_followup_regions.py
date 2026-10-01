# ============================================================================
# Define fixed follow-up locus anchors
# ============================================================================
# Purpose: Define fixed follow-up locus anchors.
# Inputs: followup_tasks.json, extracted statistics and FI lookup database.
# Outputs: Frozen regions, gene boundaries and follow-up FI regions.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import csv, json, math, sqlite3, datetime, hashlib

w = Path("work")
tasks = json.loads((w / "followup_tasks.json").read_text())
genes = {}
regions = []

# Process each declared task; preserve its original identifiers.
for t in tasks:
    assert Path(t["output_file"]).with_suffix(".done").exists()
    if t["platform"] == "Fenland" and (
        t["gene"] not in genes
        or (int(t["assay"].split("_")[0]), t["assay"])
        < (int(genes[t["gene"]]["assay"].split("_")[0]), genes[t["gene"]]["assay"])
    ):
        genes[t["gene"]] = t
for gene, t in sorted(genes.items()):
    rows = []
    for a in csv.DictReader(Path(t["output_file"]).open(), delimiter="\t"):
        try:
            ch = int(a["CHR"])
            pos = int(a["POS"])
            p = float(a["P"])
            b = float(a["BETA"])
            se = float(a["SE"])
        except ValueError:
            continue
        if (
            ch != int(t["chr"])
            or not t["start37"] - 1000000 <= pos <= t["end37"] + 1000000
            or not math.isfinite(p)
            or not math.isfinite(b)
            or not se > 0
            or not a["SNP"].startswith("rs")
        ):
            continue
        rows.append((p, -abs(b / se), pos, a["SNP"]))
    assert rows
    pv, z, anchor, snp = min(rows)
    regions.append(
        dict(
            gene=gene,
            anchor_gene=gene,
            source="Fenland",
            assay=t["assay"],
            SNP=snp,
            CHR=t["chr"],
            anchor=anchor,
            p=pv,
            main_start=anchor - 500000,
            main_end=anchor + 500000,
            sensitivity_start=anchor - 1000000,
            sensitivity_end=anchor + 1000000,
            build="GRCh37",
            frozen_at=datetime.datetime.now().astimezone().isoformat(),
        )
    )
with (w / "followup_frozen_regions.tsv").open("x") as f:
    z = csv.DictWriter(f, fieldnames=list(regions[0]), delimiter="\t")
    z.writeheader()
    z.writerows(regions)
con = sqlite3.connect(f"file:{w}/fi_lookup.sqlite?mode=ro", uri=True)
clauses = []
params = []
for r in regions:
    clauses.append("(chr=? AND pos BETWEEN ? AND ?)")
    params += [int(r["CHR"]), r["sensitivity_start"], r["sensitivity_end"]]
rows = con.execute(
    "SELECT snp,chr,pos,ea,oa,beta,se,p,n,freq FROM fi WHERE " + " OR ".join(clauses),
    params,
).fetchall()
con.close()
with (w / "followup_fi_regions.tsv").open("x") as f:
    z = csv.writer(f, delimiter="\t")
    z.writerow(
        ["SNP", "CHR", "POS", "INC_ALLELE", "DEC_ALLELE", "BETA", "SE", "P", "N", "MAF"]
    )
    z.writerows(rows)
(w / "followup_region_freeze.json").write_text(
    json.dumps(
        dict(
            status="frozen_before_new_coloc",
            regions=regions,
            FI_rows=len(rows),
            rules="analysis_rules_v1",
            manifest_sha256=hashlib.sha256(
                (w / "followup_tasks.json").read_bytes()
            ).hexdigest(),
        ),
        indent=2,
    )
)
print(json.dumps({"regions": len(regions), "FI_rows": len(rows)}))
