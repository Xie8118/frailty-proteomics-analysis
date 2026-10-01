# ============================================================================
# Extract STRING interactions among foreground proteins
# ============================================================================
# Purpose: Extract STRING interactions among foreground proteins.
# Inputs: Enrichment foreground genes and STRING v12.0 reference files.
# Outputs: Platform-specific interaction tables and summary metadata.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import gzip, csv, json, collections, hashlib

w = Path("work")
out = w / "local_enrichment_v2"
assert (out / "done").exists()
pub = Path("data/reference/STRING/v12.0")
info = pub / "9606.protein.info.v12.0.txt.gz"
links = pub / "9606.protein.links.detailed.v12.0.txt.gz"
bygene = collections.defaultdict(set)
with gzip.open(info, "rt") as f:
    header = f.readline().strip().split("\t")
    for line in f:
        a = line.rstrip().split("\t")
        bygene[a[1]].add(a[0])
fg = {
    p: set((out / f"{p}_PPI_foreground.txt").read_text().split())
    for p in ["Fenland", "deCODE"]
}
mapping = {
    p: {g: next(iter(bygene[g])) for g in fg[p] if len(bygene.get(g, [])) == 1}
    for p in fg
}
ids = {p: set(mapping[p].values()) for p in fg}
edges = {p: [] for p in fg}
with gzip.open(links, "rt") as f:
    h = f.readline().split()
    score = h.index("combined_score")
    for line in f:
        a = line.split()
        if int(a[score]) < 700 or a[0] >= a[1]:
            continue
        for p in fg:
            if a[0] in ids[p] and a[1] in ids[p]:
                edges[p].append(dict(zip(h, a)))
qc = []
for p in fg:
    inv = {v: k for k, v in mapping[p].items()}
    with (out / f"{p}_PPI_nodes.tsv").open("x") as f:
        d = csv.writer(f, delimiter="\t")
        d.writerow(["gene", "string_id"])
        d.writerows(sorted(mapping[p].items()))
    with (out / f"{p}_PPI_edges.tsv").open("x") as f:
        d = csv.DictWriter(f, fieldnames=h + ["gene_A", "gene_B"], delimiter="\t")
        d.writeheader()
        d.writerows(
            dict(e, gene_A=inv[e["protein1"]], gene_B=inv[e["protein2"]])
            for e in edges[p]
        )
    qc.append(
        dict(
            platform=p,
            foreground=len(fg[p]),
            mapped_nodes=len(ids[p]),
            unmapped_or_ambiguous=sorted(fg[p] - set(mapping[p])),
            unique_edges=len(edges[p]),
        )
    )
(out / "PPI_validation.json").write_text(
    json.dumps(
        dict(
            STRING_version="12.0",
            species=9606,
            score_threshold=700,
            added_nodes=0,
            network="functional",
            channel_evidence_retained=True,
            no_external_candidate_queries=True,
            platforms=qc,
            reference_manifest_files=[str(x) + ".source.json" for x in [info, links]],
        ),
        indent=2,
    )
)
(out / "ppi_done").write_text("LOCAL_PPI_COMPLETE\n")
print(json.dumps(qc))
