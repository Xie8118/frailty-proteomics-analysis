# ============================================================================
# Assemble harmonized candidate instruments for sensitivity analyses
# ============================================================================
# Purpose: Assemble harmonized candidate instruments for sensitivity analyses.
# Inputs: Main-scan and regional MR candidate/instrument tables.
# Outputs: final_candidate_inputs task and instrument tables.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import csv, json, math, hashlib

w = Path("work")
out = w / "final_candidate_inputs"
out.mkdir(exist_ok=False)
read = lambda p: list(csv.DictReader(p.open(), delimiter="\t"))
main = read(w / "full_mr_aggregate_v3/candidate_assay_mr.tsv")
core = read(w / "aggregate_v3/mr_candidate_assays.tsv")
mi = read(w / "full_mr_aggregate_v3/all_instruments.tsv")
ci = read(w / "aggregate_v3/mr_all_instruments.tsv")
key = lambda a: (a["platform"], a["gene"], a["assay"])
keys = {key(a) for a in main}
tasks = []
ivs = []
for a in main + [a for a in core if key(a) not in keys]:
    ismain = key(a) in keys
    tid = len(tasks) + 1
    src = mi if ismain else ci
    rows = [
        r
        for r in src
        if r["task_id"] == a["task_id"]
        and (ismain or r.get("mr_keep") in ["TRUE", "True", "1"])
    ]
    z = []
    for r in rows:
        if ismain:
            q = dict(
                SNP=r["SNP"],
                **{
                    "BETA.x": r["BETA"],
                    "SE.x": r["SE"],
                    "BETA.y.aligned": r["BETA_Y"],
                    "SE.y": r["SE_Y"],
                    "N.x": r["N"],
                    "N.y": r["N_Y"],
                },
            )
        else:
            q = {
                k: r[k]
                for k in [
                    "SNP",
                    "BETA.x",
                    "SE.x",
                    "BETA.y.aligned",
                    "SE.y",
                    "N.x",
                    "N.y",
                ]
            }
        q.update(
            task_id=tid,
            platform=a["platform"],
            gene=a["gene"],
            assay=a["assay"],
            mr_keep="TRUE",
        )
        z.append(q)
    if z:
        den = sum(float(r["BETA.x"]) ** 2 / float(r["SE.y"]) ** 2 for r in z)
        b = (
            sum(
                float(r["BETA.x"]) * float(r["BETA.y.aligned"]) / float(r["SE.y"]) ** 2
                for r in z
            )
            / den
        )
        assert math.isclose(b, float(a["beta"]), rel_tol=1e-8, abs_tol=1e-10)
    assert len(z) == int(a["nsnp"])
    ivs += z
    task = dict(
        task_id=tid,
        platform=a["platform"],
        gene=a["gene"],
        assay=a["assay"],
        source="full_mr_aggregate_v3" if ismain else "aggregate_v3",
        source_task_id=a["task_id"],
        nsnp=len(z),
        beta=a.get("beta", ""),
        p=a.get("p", ""),
        q_platform_assays=a.get("q_platform_assays", ""),
        fdr_family_n=a.get("fdr_family_n", ""),
        output_file=str(w / f"sensitivity_final/task_{tid}/sensitivity.tsv"),
    )
    tasks.append(task)
fields = [
    "SNP",
    "BETA.x",
    "SE.x",
    "BETA.y.aligned",
    "SE.y",
    "N.x",
    "N.y",
    "task_id",
    "platform",
    "gene",
    "assay",
    "mr_keep",
]


def write(name, rows, columns):
    with (out / name).open("x") as f:
        d = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        d.writeheader()
        d.writerows(rows)


write("tasks.tsv", tasks, list(tasks[0]))
write("instruments.tsv", ivs, fields)
(out / "validation.json").write_text(
    json.dumps(
        dict(
            tasks=len(tasks),
            instruments=len(ivs),
            full_main_assays=len(main),
            additional_core_assays=len(tasks) - len(main),
            all_beta_nsnp_matched=True,
            source_sha256=hashlib.sha256(
                (w / "full_mr_aggregate_v3/all_instruments.tsv").read_bytes()
            ).hexdigest(),
        ),
        indent=2,
    )
)
print((out / "validation.json").read_text())
