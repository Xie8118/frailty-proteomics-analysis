# ============================================================================
# Fetch variant consequences from Ensembl GRCh37 REST
# ============================================================================
# Purpose: Fetch variant consequences from Ensembl GRCh37 REST.
# Inputs: pav_ld/vep_tasks.json; one-based task ID.
# Outputs: Raw annotation JSON, retrieval metadata and completion marker.
# Arguments: one-based task ID.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import json, sys, urllib.request, urllib.error, time, hashlib, datetime

w = Path("work")
tid = int(sys.argv[1])
t = json.loads((w / "pav_ld/vep_tasks.json").read_text())[tid - 1]
p = Path(t["output_file"])
if p.exists():
    raise RuntimeError("Output already exists; validate before reuse")
u = "https://grch37.rest.ensembl.org/vep/human/id?canonical=1"
payload = json.dumps({"ids": t["ids"]}).encode()
r = None

# Retry only transient service/network failures; propagate permanent errors.
for attempt in range(3):
    try:
        req = urllib.request.Request(
            u,
            data=payload,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=90) as res:
            raw = res.read()
        r = json.loads(raw)
        break
    except urllib.error.HTTPError as e:
        if e.code not in [429, 500, 502, 503, 504]:
            raise
        if attempt == 2:
            raise
        time.sleep(min(30, int(e.headers.get("Retry-After", "5"))))
    except (TimeoutError, urllib.error.URLError):
        if attempt == 2:
            raise
        time.sleep(5)
assert isinstance(r, list)
with p.open("xb") as f:
    f.write(raw)
p.with_suffix(".meta.json").write_text(
    json.dumps(
        dict(
            task_id=tid,
            url=u,
            requested=t["ids"],
            records=len(r),
            retrieved_at=datetime.datetime.now().isoformat(),
            sha256=hashlib.sha256(raw).hexdigest(),
        ),
        indent=2,
    )
)
p.with_suffix(".done").write_text("ANNOTATION_FETCHED\n")
print(tid, len(r))
