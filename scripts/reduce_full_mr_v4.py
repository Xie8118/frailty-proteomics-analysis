# ============================================================================
# Aggregate the full MR scan and adjust for multiple testing
# ============================================================================
# Purpose: Aggregate the full MR scan and adjust for multiple testing.
# Inputs: Completed full MR tasks, previous MR table and PAV coverage list.
# Outputs: Assay/gene summaries, selection flow and follow-up inputs.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import json, csv, math, collections, hashlib

w = Path("work")
tasks = json.loads((w / "full_mr_v4_tasks.json").read_text())
out = w / "full_mr_aggregate_v3"
out.mkdir(exist_ok=False)
missing = [
    t["task_id"]
    for t in tasks
    if not (w / f"full_mr/task_{t['task_id']}/done").exists()
]
if missing:
    (out / "incomplete.json").write_text(json.dumps(dict(missing_tasks=missing)))
    raise RuntimeError(f"{len(missing)} tasks incomplete")
results = []
inst = []
excluded = []
seen = {}
duplicate_checks = []
# Restore gene identity from original UniProt and existing org.Hs.eg.db3.18.0, outcome-independent.
identity = {354: "SELENOM", 668: "HDGFL3"}
identity_audit = []

# Process each declared task; preserve its original identifiers.
for t in tasks:
    p = w / f"full_mr/task_{t['task_id']}"
    key = (t["platform"], identity.get(t["task_id"], t["gene"]), t["assay"])
    if key in seen:
        prior = seen[key]
        assert (p / "instruments.tsv").read_bytes() == (
            prior / "instruments.tsv"
        ).read_bytes(), (
            "Duplicate annotation produced different instrument sets; requires review"
        )
        duplicate_checks.append(
            dict(
                platform=t["platform"],
                gene=t["gene"],
                assay=t["assay"],
                duplicate_task=t["task_id"],
                retained_task=int(prior.name.split("_")[1]),
                identical_instruments=True,
            )
        )
        continue
    seen[key] = p
    r = json.loads((p / "result.json").read_text())
    assert all(r[k] == t[k] for k in ["task_id", "platform", "gene", "assay"])
    assert len(r["source_sha256"]) == 64
    if t["task_id"] in identity:
        assert t["gene"] == ""
        r["gene"] = identity[t["task_id"]]
        identity_audit.append(
            dict(
                task_id=t["task_id"],
                assay=t["assay"],
                old_gene=t["gene"],
                gene=r["gene"],
            )
        )
    assert r["gene"].strip()
    z = list(csv.DictReader((p / "instruments.tsv").open(), delimiter="\t"))
    assert len(z) == r["nsnp"]
    assert len({a["SNP"] for a in z}) == len(z)
    if z:
        den = sum(float(a["BETA"]) ** 2 / float(a["SE_Y"]) ** 2 for a in z)
        b = (
            sum(
                float(a["BETA"]) * float(a["BETA_Y"]) / float(a["SE_Y"]) ** 2 for a in z
            )
            / den
        )
        se = math.sqrt(1 / den)
        assert math.isclose(b, r["beta"], abs_tol=1e-12) and math.isclose(
            se, r["se"], rel_tol=1e-10
        )
        for a in z:
            assert (
                float(a["N"]) > 2
                and float(a["N_Y"]) > 2
                and float(a["F"]) > 10
                and float(a["SE"]) > 0
            )
            a.update({k: r[k] for k in ["task_id", "platform", "gene", "assay"]})
            inst.append(a)
    for a in r.pop("excluded_leads"):
        a.update({k: r[k] for k in ["task_id", "platform", "gene", "assay"]})
        excluded.append(a)
    results.append(r)


# Benjamini-Hochberg adjustment; preserve the platform-specific testing families.
def bh(vals):
    n = len(vals)
    order = sorted(range(n), key=lambda i: vals[i])
    out = [1.0] * n
    q = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        q = min(q, vals[i] * n / rank)
        out[i] = q
    return out


for platform in ["Fenland", "deCODE"]:
    rr = [
        r for r in results if r["platform"] == platform and r["status"] == "estimable"
    ]
    for r, q in zip(rr, bh([r["p"] for r in rr])):
        r["q_platform_assays"] = q
        r["fdr_family_n"] = len(rr)
    # All mapped assays, including non-estimable represented by p=1, provide a conservative family sensitivity.
    allr = [r for r in results if r["platform"] == platform]
    for r, q in zip(allr, bh([r.get("p", 1.0) for r in allr])):
        r["q_all_mapped_assays_sensitivity"] = q

# Summarize multiple assays per gene before gene-level adjustment.
bygene = collections.defaultdict(list)
for r in results:
    bygene[(r["platform"], r["gene"])].append(r)
genes = []
for (platform, gene), rs in sorted(bygene.items()):
    est = [r for r in rs if r["status"] == "estimable"]
    p = min(1.0, min([r["p"] for r in est], default=1.0) * len(rs))
    genes.append(
        dict(
            platform=platform,
            gene=gene,
            total_assays=len(rs),
            estimable_assays=len(est),
            assays=";".join(r["assay"] for r in rs),
            within_gene_bonferroni_p=p,
            min_assay_p=min([r["p"] for r in est], default=1.0),
            status="estimable" if est else "not_estimable",
        )
    )
for platform in ["Fenland", "deCODE"]:
    rows = [
        g for g in genes if g["platform"] == platform and g["status"] == "estimable"
    ]
    for r, q in zip(rows, bh([r["within_gene_bonferroni_p"] for r in rows])):
        r["q_platform_genes"] = q
        r["fdr_gene_family_n"] = len(rows)


def table(name, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    with (out / name).open("x") as f:
        d = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        d.writeheader()
        d.writerows(rows)


table("all_assay_mr.tsv", results)
table("gene_multiplicity_summary.tsv", genes)
table("all_instruments.tsv", inst)
table(
    "excluded_leads.tsv",
    excluded,
    ["SNP", "reason", "n_outcome_records", "task_id", "platform", "gene", "assay"],
)

# Compare target coverage with the previous MR table.
old = list(csv.DictReader((Path("data/previous_mr.tsv")).open(), delimiter="\t"))
wanted = {p: set() for p in ["Fenland", "deCODE"]}
for r in old:
    wanted["Fenland" if r["Group"].lower() == "fenland" else "deCODE"].add(
        r["Exp"].lstrip("_")
    )
# The old table contains one blank label; two original probes recover distinct symbols.
raw_label_counts = {p: len(wanted[p]) for p in wanted}
old_blank = {p: ("" in wanted[p]) for p in wanted}
wanted["Fenland"].discard("")
wanted["Fenland"].update(identity.values())
flow = []
for p in wanted:
    rr = [r for r in results if r["platform"] == p]
    gs = [g for g in genes if g["platform"] == p]
    flow.append(
        dict(
            platform=p,
            original_table_labels=raw_label_counts[p],
            resolved_original_targets=len(wanted[p]),
            mapped_genes=len(gs),
            unmapped_genes=";".join(sorted(wanted[p] - {g["gene"] for g in gs})),
            raw_assays=len(rr),
            cis_significant_assays=sum(r["preclump"] > 0 for r in rr),
            ld_clumped_assays=sum(r["clumped"] > 0 for r in rr),
            estimable_assays=sum(r["nsnp"] > 0 for r in rr),
            estimable_genes=sum(g["estimable_assays"] > 0 for g in gs),
            significant_assays=sum(r.get("q_platform_assays", 1) < 0.05 for r in rr),
            significant_genes=sum(g.get("q_platform_genes", 1) < 0.05 for g in gs),
            retained_instruments=sum(r["nsnp"] for r in rr),
        )
    )
table("selection_flow.tsv", flow)
table(
    "gene_identity_corrections.tsv",
    identity_audit,
    ["task_id", "assay", "old_gene", "gene"],
)

# Prepare follow-up tasks from significant and prespecified core candidates.
coregenes = set(
    r["gene"]
    for r in csv.DictReader((w / "extraction_tasks.tsv").open(), delimiter="\t")
)
selected = {
    r["gene"] for r in results if r.get("q_platform_assays", 1) < 0.05
} | coregenes
follow = []
follow_seen = set()

# Process each declared task; preserve its original identifiers.
for t in tasks:
    gene = identity.get(t["task_id"], t["gene"])
    key = (t["platform"], gene, t["assay"])
    if gene in selected and key not in follow_seen:
        follow.append(dict(t, gene=gene))
        follow_seen.add(key)
(out / "candidate_followup_tasks.json").write_text(json.dumps(follow, indent=2))
table("candidate_assay_mr.tsv", [r for r in results if r["gene"] in selected])
(out / "candidate_genes.txt").write_text("\n".join(sorted(selected)) + "\n")
coverage = set((w / "pav_ld/all_instruments.txt").read_text().splitlines())
extra = sorted({r["SNP"] for r in inst} - coverage)
(out / "additional_pav_instruments.txt").write_text(
    "\n".join(extra) + ("\n" if extra else "")
)
for p in wanted:
    (out / f"{p}_background_genes.txt").write_text(
        "\n".join(
            sorted(
                g["gene"]
                for g in genes
                if g["platform"] == p and g["status"] == "estimable"
            )
        )
        + "\n"
    )
    (out / f"{p}_foreground_genes.txt").write_text(
        "\n".join(
            sorted(
                g["gene"]
                for g in genes
                if g["platform"] == p and g.get("q_platform_genes", 1) < 0.05
            )
        )
        + "\n"
    )
table(
    "duplicate_annotation_checks.tsv",
    duplicate_checks,
    [
        "platform",
        "gene",
        "assay",
        "duplicate_task",
        "retained_task",
        "identical_instruments",
    ],
)
q = dict(
    unique_assays=len(results),
    duplicate_annotation_rows=len(duplicate_checks),
    status="validated_main_MR_aggregate",
    tasks=len(tasks),
    all_tasks_complete=True,
    independent_estimator_check=True,
    probes_separate=True,
    gene_summary="within gene Bonferroni across all mapped probes, then platform BH among estimable genes",
    additional_PAV_instruments=len(extra),
    followup_assays=len(follow),
    gene_identity_restored=identity_audit,
    old_blank_labels=old_blank,
    flow=flow,
    task_manifest_sha256=hashlib.sha256(
        (w / "full_mr_v4_tasks.json").read_bytes()
    ).hexdigest(),
)
(out / "validation.json").write_text(json.dumps(q, indent=2))
(out / "done").write_text("FULL_MR_AGGREGATION_COMPLETE\n")
print(json.dumps(q))
