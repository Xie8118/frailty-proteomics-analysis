# ============================================================================
# Core-assay MR and quantitative-trait colocalization
# ============================================================================
# Purpose: Core-assay MR and quantitative-trait colocalization.
# Inputs: extraction_tasks.tsv, extracted pQTLs, FI regions and frozen regions.
# Outputs: Per-task MR, colocalization, instrument tables and diagnostic plots.
# Arguments: one-based task ID and mode (pilot or full).
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(jsonlite)
  library(ggplot2)
})
setDTthreads(1)
args <- commandArgs(TRUE)
tid <- as.integer(args[1])
mode <- args[2]
w <- 'work'
ts <- fread(file.path(
  w,
  'extraction_tasks.tsv'
))
task <- ts[task_id == tid]
stopifnot(nrow(task) == 1)
out <- file.path(
  w,
  if (mode == 'pilot') 'analysis_pilot_v3' else 'analysis_v3',
  paste0(
    'task_',
    tid
  )
)
stopifnot(!dir.exists(out))
dir.create(out,
  recursive = TRUE
)
x <- fread(task$output_file)
y <- fread(file.path(
  w,
  'fi_regions.tsv'
))
reg <- fread(file.path(
  w,
  'frozen_regions.tsv'
))[gene == task$gene]
an <- fread(file.path(
  w,
  'inventory_v2/fenland_candidates.tsv'
))[gene == task$gene][1]

# Harmonize numeric fields and remove ambiguous duplicate variant IDs.
clean <- function(z) {
  for (nm in c('BETA', 'SE', 'P', 'FREQ', 'N', 'POS', 'CHR')) set(z, j = nm, value = as.numeric(z[[nm]]))
  z <- z[grepl(
    '^rs[0-9]+$',
    SNP
  ) & is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(FREQ) & FREQ > 0 & FREQ < 1 & is.finite(N) & N > 2]
  z <- unique(z)
  dup <- z[, .(n = .N), by = SNP][n > 1, SNP]
  z[!SNP %in% dup]
}
rawnx <- nrow(x)
x <- clean(x)
setnames(y, c('INC_ALLELE', 'DEC_ALLELE', 'MAF'), c('EA', 'OA', 'FREQ'))
y <- clean(y)
x[, MAF := pmin(FREQ, 1 - FREQ)]
y[, MAF := pmin(FREQ, 1 - FREQ)]
# Both single-signal ABFs and robust effect comparisons require allele-compatible variants.
j <- merge(x,
  y,
  by = 'SNP',
  suffixes = c(
    '.x',
    '.y'
  )
)
j[
  ,
  `:=`(
    EA.x = toupper(EA.x),
    OA.x = toupper(OA.x),
    EA.y = toupper(EA.y),
    OA.y = toupper(OA.y)
  )
]
comp <- function(a) chartr('ACGT', 'TGCA', a)
j[, same := (EA.x == EA.y & OA.x == OA.y)]
j[, swap := (EA.x == OA.y & OA.x == EA.y)]
j[, csame := (nchar(EA.x) == 1 & nchar(OA.x) == 1 & comp(EA.x) == EA.y & comp(OA.x) == OA.y)]
j[, cswap := (nchar(EA.x) == 1 & nchar(OA.x) == 1 & comp(EA.x) == OA.y & comp(OA.x) == EA.y)]
j[
  ,
  compatible := (same | swap | csame | cswap) & CHR.x == CHR.y & POS.x == POS.y & grepl(
    '^[ACGT]+$',
    EA.x
  ) & grepl(
    '^[ACGT]+$',
    OA.x
  )
]
fwrite(
  j[
    compatible == FALSE,
    .(
      SNP,
      CHR.x,
      POS.x,
      CHR.y,
      POS.y,
      EA.x,
      OA.x,
      EA.y,
      OA.y
    )
  ],
  file.path(
    out,
    'mismatch.tsv'
  ),
  sep = '\t'
)
j <- j[compatible == TRUE]
j[
  ,
  palindrome := (EA.x == 'A' & OA.x == 'T') | (EA.x == 'T' & OA.x == 'A') | (EA.x == 'C' & OA.x == 'G') | (EA.x == 'G' & OA.x == 'C')
]
j[, alignment := fifelse(same | csame, 1, -1)]
j[, mr_keep := TRUE]
# deCODE ImpMAF is not EAF; frequency-based orientation is unavailable for palindromes.
if (task$platform == 'deCODE') j[palindrome == TRUE, mr_keep := FALSE] else {
  j[
    palindrome == TRUE,
    mr_keep := MAF.x <= 0.42 & MAF.y <= 0.42 & (abs(FREQ.x - FREQ.y) <= 0.1 | abs(FREQ.x - (1 - FREQ.y)) <= 0.1)
  ]
  j[palindrome & mr_keep, alignment := fifelse(abs(FREQ.x - FREQ.y) < abs(FREQ.x - (1 - FREQ.y)), 1, -1)]
}
j[, BETA.y.aligned := BETA.y * alignment]

# Estimate quantitative-trait SD from effect variance, frequency and sample size.
est_sd <- function(v, f, n) {
  u <- 1 / v
  a <- 2 * n * f * (1 - f)
  sqrt(sum(u * a) / sum(u * u))
}
ans <- list()
qq <- 0

# Evaluate prespecified region widths, trait-scale assumptions and priors.
for (width in c(500000, 1000000)) {
  z <- j[CHR.y == reg$CHR & abs(POS.y - reg$anchor) <= width]
  if (nrow(z) < 100) stop(
    'Inadequate shared variant coverage: ',
    nrow(z)
  )
  sx <- est_sd(
    z$SE.x^2,
    z$MAF.x,
    z$N.x
  )
  sy <- est_sd(
    z$SE.y^2,
    z$MAF.y,
    z$N.y
  )
  stopifnot(
    is.finite(sx),
    is.finite(sy),
    sx > 0,
    sy > 0
  )
  for (sdmode in c('estimated', 'unit_variance_sensitivity')) for (prior in c(1e-5, 1e-6, 1e-4)) {
    dx <- list(
      beta = z$BETA.x,
      varbeta = z$SE.x^2,
      type = 'quant',
      sdY = if (sdmode == 'estimated') sx else 1,
      snp = z$SNP,
      position = z$POS.y
    )
    dy <- list(
      beta = z$BETA.y,
      varbeta = z$SE.y^2,
      type = 'quant',
      sdY = if (sdmode == 'estimated') sy else 1,
      snp = z$SNP,
      position = z$POS.y
    )
    r <- suppressWarnings(coloc.abf(dx,
      dy,
      p1 = 1e-4,
      p2 = 1e-4,
      p12 = prior
    ))
    q <- as.list(r$summary)
    q <- c(
      list(
        task_id = tid,
        platform = task$platform,
        gene = task$gene,
        assay = task$assay,
        width = width,
        p12 = prior,
        sdmode = sdmode,
        sdY_protein = sx,
        sdY_FI = sy,
        N_protein_min = min(z$N.x),
        N_protein_max = max(z$N.x),
        N_FI_min = min(z$N.y),
        N_FI_max = max(z$N.y)
      ),
      q
    )
    qq <- qq + 1
    ans[[qq]] <- as.data.table(q)
  }
}
fwrite(rbindlist(ans, fill = TRUE), file.path(out, 'coloc.tsv'), sep = '\t')
# Rebuild cis instruments independently of FI outcome significance.
iv <- x[CHR == an$Chr & POS >= an$start - 1e6 & POS <= an$end + 1e6 & P < 5e-8 & (BETA / SE)^2 > 10 & !(CHR == 6 & POS >= 28477797 & POS <= 33448354)]
fwrite(iv, file.path(out, 'cis_preclump.tsv'), sep = '\t')
snps <- character()
if (nrow(iv)) {
  pref <- file.path(
    out,
    'clump'
  )
  cmd <- c(
    '--bfile',
    'data/reference/EUR_1000Genome_phase3_all',
    '--clump',
    file.path(
      out,
      'cis_preclump.tsv'
    ),
    '--clump-snp-field',
    'SNP',
    '--clump-field',
    'P',
    '--clump-p1',
    '5e-8',
    '--clump-p2',
    '1',
    '--clump-r2',
    '0.01',
    '--clump-kb',
    '10000',
    '--threads',
    '1',
    '--memory',
    '4000',
    '--out',
    pref
  )
  status <- system2('bin/plink',
    cmd,
    stdout = file.path(
      out,
      'plink.stdout'
    ),
    stderr = file.path(
      out,
      'plink.stderr'
    )
  )
  stopifnot(status == 0)
  if (file.exists(paste0(pref, '.clumped'))) snps <- fread(paste0(pref, '.clumped'))$SNP
}
z <- j[SNP %in% snps & mr_keep]
z[, `:=`(F = (BETA.x / SE.x)^2, PVE = BETA.x^2 / (BETA.x^2 + N.x * SE.x^2))]
fwrite(z, file.path(out, 'mr_instruments.tsv'), sep = '\t')
if (nrow(z)) {
  b <- sum(z$BETA.x * z$BETA.y.aligned / z$SE.y^2) / sum(z$BETA.x^2 / z$SE.y^2)
  se <- sqrt(1 / sum(z$BETA.x^2 / z$SE.y^2))
  p <- 2 * pnorm(-abs(b / se))
  Q <- sum((z$BETA.y.aligned - b * z$BETA.x)^2 / z$SE.y^2)
  ser <- se * sqrt(if (nrow(z) > 1) max(
    1,
    Q / (nrow(z) - 1)
  ) else 1)
  independent_b <- unname(coef(lm(z$BETA.y.aligned ~ z$BETA.x - 1,
    weights = 1 / z$SE.y^2
  ))[1])
  stopifnot(abs(b - independent_b) < 1e-10)
  mr <- data.table(
    task_id = tid,
    platform = task$platform,
    gene = task$gene,
    assay = task$assay,
    status = 'estimable',
    method = if (nrow(z) == 1) 'Wald' else 'IVW_fixed',
    nsnp = nrow(z),
    beta = b,
    se = se,
    lower = b - 1.96 * se,
    upper = b + 1.96 * se,
    p = p,
    random_se = ser,
    random_p = 2 * pnorm(-abs(b / ser)),
    Q = Q,
    Qp = if (nrow(z) > 1) pchisq(Q,
      nrow(z) - 1,
      lower.tail = FALSE
    ) else NA_real_,
    Fmin = min(z$F),
    PVE_sum_approx = sum(z$PVE)
  )
} else mr <- data.table(
  task_id = tid,
  platform = task$platform,
  gene = task$gene,
  assay = task$assay,
  status = 'not_estimable_after_filters',
  nsnp = 0
)
fwrite(mr, file.path(out, 'mr.tsv'), sep = '\t')
source <- j[
  CHR.y == reg$CHR & abs(POS.y - reg$anchor) <= 1000000,
  .(SNP,
    POS = POS.y,
    p_protein = pmax(
      P.x,
      1e-300
    ),
    p_FI = pmax(
      P.y,
      1e-300
    ),
    beta_protein = BETA.x,
    beta_FI = BETA.y.aligned
  )
]
fwrite(source,
  file.path(
    out,
    'regional_source.tsv'
  ),
  sep = '\t'
)
plotdata <- melt(
  source[
    ,
    .(
      POS,
      p_protein,
      p_FI
    )
  ],
  id.vars = 'POS',
  variable.name = 'trait',
  value.name = 'p'
)
g <- ggplot(plotdata, aes(POS / 1e6, -log10(p), color = trait)) +
  geom_point(size = .5, alpha = .6) +
  facet_wrap(~trait, ncol = 1, scales = 'free_y') +
  theme_bw(base_size = 10) +
  labs(
    title = paste(
      task$gene,
      task$platform,
      task$assay
    ),
    x = 'GRCh37 position (Mb)',
    y = expression(-log[10](P))
  ) +
  theme(legend.position = 'none')
ggsave(file.path(out, 'regional_candidate.png'), g, width = 7, height = 4.5, dpi = 180)
write_json(
  list(
    task_id = tid,
    mode = mode,
    raw_region_rows = rawnx,
    valid_exposure_rows = nrow(x),
    harmonized_coloc_rows = nrow(j),
    preclump = nrow(iv),
    clumped = length(snps),
    mr_instruments = nrow(z),
    rules = 'analysis_rules_v1',
    status = 'computed_validated',
    frequency_caveat = 'deCODE ImpMAF not used as EAF; unresolved MR palindromes excluded'
  ),
  file.path(
    out,
    'validation.json'
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)
writeLines('ANALYSIS_COMPLETE', file.path(out, 'done'))
cat('ASSAY_COMPLETE', tid, '\n')
