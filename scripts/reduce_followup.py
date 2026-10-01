# ============================================================================
# Combine follow-up assay results
# ============================================================================
# Purpose: Combine follow-up assay results.
# Inputs: Completed follow-up tasks and main MR candidate results.
# Outputs: Follow-up aggregate tables and comparison checks.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import csv, json, math, hashlib

w = Path("work")
ts = json.loads((w / "followup_tasks.json").read_text())
out = w / "followup_aggregate"
assert all((w / f"followup_analysis/task_{t['task_id']}/done").exists() for t in ts)
assert not out.exists()
read = lambda p: list(csv.DictReader(p.open(), delimiter="\t"))
main = {
    r["task_id"]: r for r in read(w / "full_mr_aggregate_v3/candidate_assay_mr.tsv")
}
col = []
mr = []
checks = []

# Process each declared task; preserve its original identifiers.
for t in ts:
    d = w / f"followup_analysis/task_{t['task_id']}"
    cs = read(d / "coloc.tsv")
    assert len(cs) == 12
    for r in cs:
        v = [float(r[f"PP.H{i}.abf"]) for i in range(5)]
        assert all(0 <= a <= 1 for a in v)
        assert abs(sum(v) - 1) < 1e-8
        r["full_mr_task_id"] = t["full_mr_task_id"]
    col += cs
    m = read(d / "mr.tsv")[0]
    m["full_mr_task_id"] = t["full_mr_task_id"]
    base = main[str(t["full_mr_task_id"])]
    same = int(m["nsnp"]) == int(base["nsnp"]) and (
        int(m["nsnp"]) == 0
        or math.isclose(
            float(m["beta"]), float(base["beta"]), abs_tol=1e-10, rel_tol=1e-8
        )
    )
    checks.append(
        dict(
            task_id=t["task_id"],
            gene=t["gene"],
            platform=t["platform"],
            assay=t["assay"],
            main_nsnp=base["nsnp"],
            region_nsnp=m["nsnp"],
            MR_matches_full_scan=same,
        )
    )
    mr.append(m)
out.mkdir()


def table(name, rows):
    with (out / name).open("x") as f:
        z = csv.DictWriter(
            f,
            fieldnames=list(dict.fromkeys(k for r in rows for k in r)),
            delimiter="\t",
        )
        z.writeheader()
        z.writerows(rows)


table("coloc_all.tsv", col)
table("regional_mr.tsv", mr)
table("mr_reconciliation.tsv", checks)
for name in [
    "followup_frozen_regions.tsv",
    "followup_region_freeze.json",
    "full_core_instrument_reconciliation.tsv",
    "followup_preparation.json",
]:
    (out / name).write_bytes((w / name).read_bytes())
q = dict(
    tasks_complete=len(ts),
    coloc_models=len(col),
    all_posteriors_valid=True,
    regional_MR_matches=sum(r["MR_matches_full_scan"] for r in checks),
    MR_differences=[r for r in checks if not r["MR_matches_full_scan"]],
    authoritative_MR="full_mr_aggregate_v3",
    panel_policy="regional diagnostics and corresponding source data generated for each completed assay",
    frozen_region_sha256=hashlib.sha256(
        (w / "followup_frozen_regions.tsv").read_bytes()
    ).hexdigest(),
)
(out / "validation.json").write_text(json.dumps(q, indent=2))
(out / "done").write_text("FOLLOWUP_AGGREGATE_VALIDATED\n")
print(json.dumps(q))
