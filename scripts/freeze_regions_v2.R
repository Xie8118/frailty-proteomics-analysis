# ============================================================================
# Define initial locus anchors from Fenland discovery signals
# ============================================================================
# Purpose: Define initial locus anchors from Fenland discovery signals.
# Inputs: extraction_tasks.tsv, candidate boundaries and extracted pQTLs.
# Outputs: frozen_regions.tsv and anchor provenance.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

library(data.table)
library(jsonlite)
w <- 'work'
ts <- fread(file.path(
  w,
  'extraction_tasks.tsv'
))
an <- fread(file.path(
  w,
  'inventory_v2/fenland_candidates.tsv'
))
res <- list()
for (g in unique(ts$gene)) {
  t <- ts[gene == g & platform == 'Fenland']
  if (!nrow(t)) next
  t[, assay_number := as.integer(sub('_.*', '', assay))]
  setorder(t, assay_number, assay)
  t <- t[1]
  x <- fread(t$output_file)
  x[
    ,
    P := as.numeric(P)
  ]
  a <- an[gene == g][1]
  x <- x[CHR == as.integer(a$Chr) & POS >= a$start - 1e6 & POS <= a$end + 1e6 & is.finite(P) & is.finite(BETA) & SE > 0 & grepl(
    '^rs[0-9]+$',
    SNP
  )]
  x[, rankz := abs(BETA / SE)]
  setorder(x, P, -rankz, POS, SNP)
  stopifnot(nrow(x) > 0)
  x <- x[1]
  res[[g]] <- data.table(
    gene = g,
    anchor_gene = g,
    source = 'Fenland',
    assay = t$assay,
    SNP = x$SNP,
    CHR = x$CHR,
    anchor = x$POS,
    p = x$P
  )
}
z <- rbindlist(res)
a <- z[gene == 'NMT1']
for (g in c(
  'HEXIM1',
  'HEXIM2'
)) z[
  gene == g,
  c(
    'anchor_gene',
    'assay',
    'SNP',
    'CHR',
    'anchor',
    'p'
  ) := a[
    ,
    .(
      anchor_gene,
      assay,
      SNP,
      CHR,
      anchor,
      p
    )
  ]
]
z[
  ,
  `:=`(
    main_start = anchor - 5e5,
    main_end = anchor + 5e5,
    sensitivity_start = anchor - 1e6,
    sensitivity_end = anchor + 1e6,
    build = 'GRCh37',
    frozen_at = as.character(Sys.time())
  )
]
fwrite(z,
  file.path(
    w,
    'frozen_regions.tsv'
  ),
  sep = '\t'
)
write_json(
  list(
    status = 'frozen_before_updated_coloc',
    discovery_only = TRUE,
    n = nrow(z)
  ),
  file.path(
    w,
    'region_freeze.json'
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)
cat('REGIONS_FROZEN\n')
