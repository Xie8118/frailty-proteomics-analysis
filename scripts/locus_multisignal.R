# ============================================================================
# Fine-map the shared NMT1/HEXIM region
# ============================================================================
# Purpose: Fine-map the shared NMT1/HEXIM region.
# Inputs: Core pQTL/FI regions, fixed locus anchor and PLINK LD reference.
# Outputs: Trait-level SuSiE results and eligible paired colocalization results.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(susieR)
  library(jsonlite)
})
setDTthreads(1)
set.seed(20260916)
w <- 'work'
out <- file.path(
  w,
  'multisignal'
)
dir.create(out,
  showWarnings = FALSE
)
ts <- fread(file.path(
  w,
  'extraction_tasks.tsv'
))
reg <- fread(file.path(
  w,
  'frozen_regions.tsv'
))[gene == 'NMT1']
ts <- ts[gene %in% c(
  'NMT1',
  'HEXIM1',
  'HEXIM2'
)]

# Harmonize numeric fields and remove ambiguous duplicate variant IDs.
clean <- function(z) {
  for (nm in c(
    'BETA',
    'SE',
    'P',
    'FREQ',
    'N',
    'POS',
    'CHR'
  )) set(z,
    j = nm,
    value = as.numeric(z[[nm]])
  )
  z <- z[CHR == 17 & abs(POS - reg$anchor) <= 5e5 & grepl(
    '^rs[0-9]+$',
    SNP
  ) & is.finite(BETA) & SE > 0 & FREQ > 0 & FREQ < 1]
  z <- unique(z)
  z <- z[!SNP %in% z[,
    .(n = .N),
    by = SNP
  ][
    n > 1,
    SNP
  ]]
  z[grepl(
    '^[ACGT]$',
    EA
  ) & grepl(
    '^[ACGT]$',
    OA
  ) & !((EA == 'A' & OA == 'T') | (EA == 'T' & OA == 'A') | (EA == 'C' & OA == 'G') | (EA == 'G' & OA == 'C'))]
}
y <- fread(file.path(
  w,
  'fi_regions.tsv'
))
setnames(
  y,
  c(
    'INC_ALLELE',
    'DEC_ALLELE',
    'MAF'
  ),
  c(
    'EA',
    'OA',
    'FREQ'
  )
)
all <- list(FI = clean(y))
for (i in seq_len(nrow(ts))) {
  t <- ts[i]
  all[[paste(t$platform, t$gene, t$assay, sep = '_')]] <- clean(fread(t$output_file))
}

# Use a common set of variants across the locus traits.
ids <- Reduce(
  intersect,
  lapply(
    all,
    function(z) z$SNP
  )
)
stopifnot(length(ids) > 100)
writeLines(
  ids,
  file.path(
    out,
    'common_rsids.txt'
  )
)
bin <- 'bin/plink'
ref <- 'data/reference/EUR_1000Genome_phase3_all'
pref <- file.path(out, 'locus')
cmd <- c(
  '--bfile',
  ref,
  '--extract',
  file.path(
    out,
    'common_rsids.txt'
  ),
  '--keep-allele-order',
  '--make-bed',
  '--threads',
  '1',
  '--memory',
  '6000',
  '--out',
  pref
)
stopifnot(system2(bin,
  cmd,
  stdout = file.path(
    out,
    'extract.log'
  ),
  stderr = file.path(
    out,
    'extract.err'
  )
) == 0)
stopifnot(system2(bin,
  c(
    '--bfile',
    pref,
    '--keep-allele-order',
    '--r',
    'square',
    'gz',
    '--threads',
    '1',
    '--memory',
    '6000',
    '--out',
    pref
  ),
  stdout = file.path(
    out,
    'ld.log'
  ),
  stderr = file.path(
    out,
    'ld.err'
  )
) == 0)
b <- fread(
  paste0(
    pref,
    '.bim'
  ),
  col.names = c(
    'CHR',
    'SNP',
    'CM',
    'POS',
    'A1',
    'A2'
  )
)
R <- as.matrix(fread(cmd = paste(
  'gzip -dc',
  paste0(
    pref,
    '.ld.gz'
  )
)))
stopifnot(
  nrow(R) == nrow(b),
  max(abs(R - t(R))) < 1e-6,
  all(is.finite(R)),
  all(abs(diag(R) - 1) < 1e-5)
)
dimnames(R) <- list(
  b$SNP,
  b$SNP
)
fits <- list()
qc <- list()
comp <- function(a) chartr('ACGT', 'TGCA', a)
for (nm in names(all)) {
  z <- all[[nm]][match(
    b$SNP,
    SNP
  )]
  stopifnot(all(z$SNP == b$SNP))
  same <- (z$EA == b$A1 & z$OA == b$A2) | (comp(z$EA) == b$A1 & comp(z$OA) == b$A2)
  flip <- (z$EA == b$A2 & z$OA == b$A1) | (comp(z$EA) == b$A2 & comp(z$OA) == b$A1)
  if (!all(same | flip)) stop('allele mismatch ', nm)
  zscore <- z$BETA / z$SE * ifelse(same, 1, -1)
  n <- median(z$N)
  s <- estimate_s_rss(z = zscore, R = R, n = n)
  qc[[nm]] <- data.table(trait = nm, snps = nrow(b), N = n, regularization_s = s, qualified = s <= 0.2)
  if (s > 0.2) next
  fit <- susie_rss(
    z = zscore,
    R = R,
    n = n,
    L = 5,
    coverage = .95,
    estimate_residual_variance = FALSE,
    max_iter = 1000,
    check_prior = TRUE
  )
  # coloc.susie expects named SNPs in log-BF columns.
  colnames(fit$lbf_variable) <- b$SNP
  names(fit$pip) <- b$SNP
  fit$sets$cs_index <- fit$sets$cs_index
  fits[[nm]] <- fit
  fwrite(
    data.table(
      SNP = b$SNP,
      POS = b$POS,
      pip = fit$pip
    ),
    file.path(
      out,
      paste0(
        nm,
        '_pip.tsv'
      )
    ),
    sep = '\t'
  )
  saveRDS(
    fit,
    file.path(
      out,
      paste0(
        nm,
        '.rds'
      )
    )
  )
  qc[[nm]][, `:=`(converged = fit$converged, credible_sets = length(fit$sets$cs))]
}
fwrite(rbindlist(qc, fill = TRUE), file.path(out, 'ld_model_qc.tsv'), sep = '\t')
rs <- list()
if ('FI' %in% names(fits)) for (nm in setdiff(names(fits), 'FI')) {
  if (!length(fits[[nm]]$sets$cs) || !length(fits$FI$sets$cs)) next
  rr <- coloc.susie(fits[[nm]], fits$FI, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  if (!is.null(rr$summary)) {
    z <- as.data.table(rr$summary)
    z[, trait := nm]
    rs[[nm]] <- z
  }
}
if (length(rs)) fwrite(rbindlist(rs, fill = TRUE), file.path(out, 'coloc_susie.tsv'), sep = '\t')
write_json(
  list(
    n_traits = length(all),
    n_qualified = length(fits),
    n_snps = nrow(b),
    L = 5,
    coverage = .95,
    regularization_limit = .2,
    scope = 'shared_locus_resolution_support_not_causal_protein_proof'
  ),
  file.path(
    out,
    'validation.json'
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)
cat('MULTISIGNAL_COMPLETE\n')
