# ============================================================================
# Combine allele-matched PAV annotations and LD proxies
# ============================================================================
# Purpose: Combine allele-matched PAV annotations and LD proxies.
# Inputs: VEP batches, LD reference, LD edges and candidate instruments.
# Outputs: Variant annotations, proxy flags and PAV sensitivity results.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import json, csv, hashlib, collections, math

w = Path("work")
tasks = json.loads((w / "pav_ld_extra/combined_vep_tasks.json").read_text())
out = w / "pav_aggregate_final"
out.mkdir(exist_ok=False)
missing = [
    t["task_id"]
    for t in tasks
    if not Path(t["output_file"]).with_suffix(".done").exists()
]
if missing:
    (out / "incomplete.json").write_text(json.dumps({"missing_batches": missing}))
    raise RuntimeError(f"Missing {len(missing)} batches; no final aggregation")
ids = set(i for t in tasks for i in t["ids"])
bim = {}
dups = set()
for line in open("data/reference/EUR_1000Genome_phase3_all.bim"):
    a = line.split()
    s = a[1]
    if s not in ids:
        continue
    if s in bim:
        dups.add(s)
    bim[s] = (a[0], int(a[3]), a[4], a[5])
for s in dups:
    bim.pop(s)

# Define coding and essential splice consequences used for PAV exclusion.
pavterms = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
}
comp = str.maketrans("ACGT", "TGCA")
ann = {}
details = []
nrecords = 0
missing_response = set()

# Process each declared task; preserve its original identifiers.
for t in tasks:
    p = Path(t["output_file"])
    raw = p.read_bytes()
    meta = json.loads(p.with_suffix(".meta.json").read_text())
    assert hashlib.sha256(raw).hexdigest() == meta["sha256"]
    recs = json.loads(raw)
    returned = set()
    for r in recs:
        s = r.get("input", r.get("id"))
        returned.add(s)
        nrecords += 1
        if s not in ids:
            raise RuntimeError("Unexpected VEP identifier " + s)
        z = ann.setdefault(
            s,
            dict(
                matched=0,
                pav_genes=set(),
                consequences=set(),
                transcripts=set(),
                splice_region_genes=set(),
                records=0,
            ),
        )
        z["records"] += 1
        if s not in bim:
            continue
        ch, pos, a1, a2 = bim[s]
        alleles = r.get("allele_string", "").split("/")
        if (
            r.get("assembly_name") != "GRCh37"
            or str(r.get("seq_region_name")) != ch
            or r.get("start") != pos
        ):
            continue
        # SNP and explicit sequence allele pairs only; indel representation mismatch remains unknown.
        pair = {a1, a2}
        direct = pair <= set(alleles) and alleles[0] in pair
        cpair = {a.translate(comp) for a in pair}
        reverse = (
            all(len(a) == 1 and a in "ACGT" for a in [a1, a2])
            and cpair <= set(alleles)
            and alleles[0] in cpair
        )
        if len(alleles) < 2 or not (direct or reverse):
            continue
        observed_alt = (pair if direct else cpair) - {alleles[0]}
        if len(observed_alt) != 1:
            continue
        z["matched"] += 1
        for tr in r.get("transcript_consequences", []):
            if tr.get("variant_allele") not in observed_alt:
                continue
            terms = set(tr.get("consequence_terms", []))
            gene = tr.get("gene_symbol", tr.get("gene_id", "unknown"))
            tx = tr.get("transcript_id", "")
            z["consequences"] |= terms
            hit = bool(terms & pavterms)
            if hit:
                z["pav_genes"].add(gene)
                z["transcripts"].add(tx)
            if "splice_region_variant" in terms:
                z["splice_region_genes"].add(gene)
            if hit or "splice_region_variant" in terms:
                details.append(
                    dict(
                        SNP=s,
                        CHR=ch,
                        POS=pos,
                        alleles="/".join(alleles),
                        gene=gene,
                        transcript=tx,
                        consequences=";".join(sorted(terms)),
                        PAV=hit,
                        canonical=tr.get("canonical", 0),
                    )
                )
    missing_response |= set(t["ids"]) - returned

# Keep missing or allele-mismatched annotation distinct from a negative result.
variant = []
for s in sorted(ids):
    z = ann.get(s, {})
    known = z.get("matched", 0) > 0
    variant.append(
        dict(
            SNP=s,
            annotation_status=(
                "matched"
                if known
                else (
                    "absent_from_LD_reference"
                    if s not in bim
                    else "no_matching_allele_annotation"
                )
            ),
            PAV=bool(z.get("pav_genes")),
            pav_genes=";".join(sorted(z.get("pav_genes", []))),
            consequences=";".join(sorted(z.get("consequences", []))),
            transcripts=";".join(sorted(z.get("transcripts", []))),
            splice_region_genes=";".join(sorted(z.get("splice_region_genes", []))),
        )
    )


def table(name, rows, fields=None):
    with (out / name).open("x") as f:
        d = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter="\t")
        d.writeheader()
        d.writerows(rows)


table("variant_annotations.tsv", variant)
table(
    "pav_transcript_details.tsv",
    details,
    [
        "SNP",
        "CHR",
        "POS",
        "alleles",
        "gene",
        "transcript",
        "consequences",
        "PAV",
        "canonical",
    ],
)
lookup = {x["SNP"]: x for x in variant}
instruments = (w / "pav_ld_extra/all_instruments.txt").read_text().splitlines()
neighbours = {s: {s: 1.0} for s in instruments}
for ld_file in ["pav_ld/all.ld", "pav_ld_extra/all.ld"]:
    with (w / ld_file).open() as f:
        h = f.readline().split()
        for line in f:
            a = line.split()
            s = a[2]
            p = a[5]
            r2 = float(a[6])
            if s in neighbours and r2 >= 0.6:
                neighbours[s][p] = max(r2, neighbours[s].get(p, 0))

# Evaluate the main LD threshold (0.8) and sensitivity threshold (0.6).
flags = []
proxy_rows = []
for s in instruments:
    ns = neighbours[s]
    for threshold in [0.8, 0.6]:
        relevant = {p: r for p, r in ns.items() if r >= threshold}
        hits = [p for p in relevant if lookup.get(p, {}).get("PAV", False)]
        unknown = [
            p
            for p in relevant
            if lookup.get(p, {}).get("annotation_status") != "matched"
        ]
        target_genes = set(
            g for p in hits for g in lookup[p]["pav_genes"].split(";") if g
        )
        flags.append(
            dict(
                SNP=s,
                ld_threshold=threshold,
                ld_reference_covered=s in bim,
                neighbours=len(relevant),
                annotation_unknown=len(unknown),
                PAV_flag=bool(hits),
                direct_PAV=s in hits,
                PAV_proxies=";".join(sorted(hits)),
                PAV_genes=";".join(sorted(target_genes)),
                evidence_status=(
                    "PAV_linked"
                    if hits
                    else (
                        "annotation_incomplete"
                        if unknown or s not in bim
                        else "no_PAV_in_assessed_region"
                    )
                ),
            )
        )
        for p in hits:
            proxy_rows.append(
                dict(
                    SNP=s,
                    proxy=p,
                    r2=relevant[p],
                    threshold=threshold,
                    genes=lookup[p]["pav_genes"],
                    consequences=lookup[p]["consequences"],
                )
            )
table("instrument_pav_flags.tsv", flags)
table(
    "pav_proxy_evidence.tsv",
    proxy_rows,
    ["SNP", "proxy", "r2", "threshold", "genes", "consequences"],
)
core = list(
    csv.DictReader(
        (w / "final_candidate_inputs/instruments.tsv").open(), delimiter="\t"
    )
)
assays = list(
    csv.DictReader((w / "final_candidate_inputs/tasks.tsv").open(), delimiter="\t")
)
fl = {(a["SNP"], a["ld_threshold"]): a for a in flags}
results = []
for assay in assays:
    rows = [r for r in core if r["task_id"] == assay["task_id"]]
    usable = [r for r in rows if r["mr_keep"] in ["TRUE", "True", "1"]]
    for threshold in [0.8, 0.6]:
        kept = [
            r
            for r in usable
            if not fl.get((r["SNP"], threshold), {}).get("PAV_flag", False)
        ]
        unknown = sum(
            fl.get((r["SNP"], threshold), {}).get(
                "evidence_status", "annotation_incomplete"
            )
            == "annotation_incomplete"
            for r in kept
        )
        r = {k: assay[k] for k in ["task_id", "platform", "gene", "assay"]}
        r.update(
            ld_threshold=threshold,
            original_nsnp=len(usable),
            remaining_nsnp=len(kept),
            excluded_PAV=len(usable) - len(kept),
            remaining_unknown=unknown,
            status="estimable" if kept else "not_estimable",
        )
        if kept:
            den = sum(float(z["BETA.x"]) ** 2 / float(z["SE.y"]) ** 2 for z in kept)
            beta = (
                sum(
                    float(z["BETA.x"])
                    * float(z["BETA.y.aligned"])
                    / float(z["SE.y"]) ** 2
                    for z in kept
                )
                / den
            )
            se = math.sqrt(1 / den)
            r.update(
                beta=beta,
                se=se,
                p=math.erfc(abs(beta / se) / math.sqrt(2)),
                lower=beta - 1.96 * se,
                upper=beta + 1.96 * se,
            )
        results.append(r)
fields = list(dict.fromkeys(k for r in results for k in r))
table("candidate_pav_sensitivity.tsv", results, fields)
q = dict(
    requested_batches=len(tasks),
    records=nrecords,
    requested_variants=len(ids),
    missing_response_ids=sorted(missing_response),
    matched_alleles=sum(x["annotation_status"] == "matched" for x in variant),
    PAV_variants=sum(x["PAV"] for x in variant),
    instruments=len(instruments),
    core_assays=len(assays),
    LD_reference="1000 Genomes phase3 EUR 503 participants",
    LD_window_kb=1000,
    thresholds=[0.8, 0.6],
    PAV_terms=sorted(pavterms),
    splice_region_policy="reported separately; essential donor/acceptor included in exclusion",
    unknown_policy="retained but explicitly flagged; not counted as annotation negative",
    source="Ensembl GRCh37 REST15.12; coordinate/ref/observed ALT matched to LD reference, unrelated alternate alleles ignored",
)
(out / "validation.json").write_text(json.dumps(q, indent=2))
(out / "done").write_text("PAV_CORE_AGGREGATION_COMPLETE\n")
print(json.dumps({k: v for k, v in q.items() if k != "missing_response_ids"}))
