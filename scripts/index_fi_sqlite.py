# ============================================================================
# Build an indexed FI GWAS lookup database
# ============================================================================
# Purpose: Build an indexed FI GWAS lookup database.
# Inputs: data/fi/fi_summary_GRCh37.txt.gz.
# Outputs: work/fi_lookup.sqlite and completion metadata.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

from pathlib import Path
import sqlite3, gzip, csv, json, hashlib, time

w = Path("work")
p = w / "fi_lookup.sqlite"
assert not p.exists()
src = Path("data/fi/fi_summary_GRCh37.txt.gz")
start = time.time()
con = sqlite3.connect(p)
con.execute("PRAGMA journal_mode=OFF")
con.execute("PRAGMA synchronous=OFF")
con.execute("PRAGMA temp_store=MEMORY")
con.execute(
    "CREATE TABLE fi(snp TEXT,chr INTEGER,pos INTEGER,ea TEXT,oa TEXT,beta REAL,se REAL,p REAL,n REAL,freq REAL)"
)
batch = []
n = 0
with gzip.open(src, "rt") as f:
    h = f.readline().split()
    ix = {a: i for i, a in enumerate(h)}
    for line in f:
        a = line.split()
        r = [
            a[ix[k]]
            for k in [
                "SNP",
                "CHR",
                "POS",
                "INC_ALLELE",
                "DEC_ALLELE",
                "BETA",
                "SE",
                "P",
                "N",
                "MAF",
            ]
        ]
        batch.append(r)
        n += 1
        if len(batch) == 20000:
            con.executemany("INSERT INTO fi VALUES (?,?,?,?,?,?,?,?,?,?)", batch)
            batch = []
if batch:
    con.executemany("INSERT INTO fi VALUES (?,?,?,?,?,?,?,?,?,?)", batch)
con.execute("CREATE INDEX fi_snp ON fi(snp)")
con.commit()
assert con.execute("SELECT COUNT(*) FROM fi").fetchone()[0] == n
con.close()
(w / "fi_lookup.done").write_text(
    json.dumps(
        dict(
            rows=n,
            seconds=time.time() - start,
            source=str(src),
            source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
        ),
        indent=2,
    )
)
print("FI_INDEX_COMPLETE", n)
