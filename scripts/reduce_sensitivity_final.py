# ============================================================================
# Combine completed MR sensitivity tasks
# ============================================================================
# Purpose: Combine completed MR sensitivity tasks.
# Inputs: 46 candidate tasks and their completed sensitivity results.
# Outputs: Combined sensitivity, leave-one-out and directionality tables.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import csv, json, math

w = Path("work")
tasks = list(
    csv.DictReader((w / "final_candidate_inputs/tasks.tsv").open(), delimiter="\t")
)
out = w / "sensitivity_final_aggregate"
assert len(tasks) == 46
missing = [
    t["task_id"]
    for t in tasks
    if not (w / f"sensitivity_final/task_{t['task_id']}/done").exists()
]
assert not missing, missing
assert not out.exists()
rows = []
loo = []
ste = []
checks = []

# Process each declared task; preserve its original identifiers.
for t in tasks:
    d = w / f"sensitivity_final/task_{t['task_id']}"
    r = list(csv.DictReader((d / "sensitivity.tsv").open(), delimiter="\t"))[0]
    baseline = t
    assert int(r["nsnp"]) == int(baseline["nsnp"])
    if int(r["nsnp"]):
        for k in ["beta"]:
            assert math.isclose(
                float(r[k]), float(baseline[k]), rel_tol=1e-9, abs_tol=1e-12
            )
    if r["egger_status"] == "computed":
        assert all(
            math.isfinite(float(r[k]))
            for k in ["egger_beta", "egger_se", "egger_intercept"]
        )
    for name, target in [("leave_one_out.tsv", loo), ("steiger_variants.tsv", ste)]:
        if (d / name).exists():
            for a in csv.DictReader((d / name).open(), delimiter="\t"):
                target.append(
                    dict(
                        {k: r[k] for k in ["task_id", "platform", "gene", "assay"]}, **a
                    )
                )
    rows.append(r)
    checks.append(json.loads((d / "validation.json").read_text()))
out.mkdir()


def table(name, z):
    with (out / name).open("x") as f:
        wr = csv.DictWriter(
            f, fieldnames=list(dict.fromkeys(k for r in z for k in r)), delimiter="\t"
        )
        wr.writeheader()
        wr.writerows(z)


table("candidate_sensitivity.tsv", rows)
table("leave_one_out.tsv", loo)
table("steiger_variants.tsv", ste)
qc = {
    "tasks_complete": len(rows),
    "independent_egger_checks": sum(
        x.get("independent_egger_check") is True for x in checks
    ),
    "all_primary_estimates_match": True,
    "egger_estimable": sum(x["egger_status"] == "computed" for x in rows),
    "leave_one_out_assays": sum(x["loo_status"] == "computed" for x in rows),
    "Steiger_inference": "Approximate direction only, no formal p with unknown overlap",
    "status": "validated",
}
(out / "validation.json").write_text(json.dumps(qc, indent=2))
(out / "done").write_text("SENSITIVITY_AGGREGATE_COMPLETE\n")
print(json.dumps(qc))
