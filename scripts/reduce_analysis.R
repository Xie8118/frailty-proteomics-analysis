# ============================================================================
# Combine and validate core-assay results
# ============================================================================
# Purpose: Combine and validate core-assay results.
# Inputs: Completed analysis_v3 task directories.
# Outputs: Combined MR, colocalization and instrument tables.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
})
w <- 'work'
ts <- fread(file.path(
  w,
  'extraction_tasks.tsv'
))
out <- file.path(
  w,
  'aggregate_v3'
)
dir.create(out,
  showWarnings = FALSE
)
cs <- list()
ms <- list()
ivs <- list()
val <- list()
for (i in ts$task_id) {
  d <- file.path(w, 'analysis_v3', paste0('task_', i))
  stopifnot(file.exists(file.path(d, 'done')))
  c <- fread(file.path(
    d,
    'coloc.tsv'
  ))
  stopifnot(
    nrow(c) == 12,
    all(abs(rowSums(c[,
      grep(
        '^PP.H',
        names(c)
      ),
      with = FALSE
    ]) - 1) < 1e-8)
  )
  cs[[i]] <- c
  m <- fread(file.path(d, 'mr.tsv'))
  stopifnot(nrow(m) == 1)
  ms[[i]] <- m
  z <- fread(file.path(d, 'mr_instruments.tsv'))
  if (nrow(z)) {
    z[
      ,
      `:=`(
        task_id = i,
        platform = ts[
          task_id == i,
          platform
        ],
        gene = ts[
          task_id == i,
          gene
        ],
        assay = ts[
          task_id == i,
          assay
        ]
      )
    ]
    ivs[[length(ivs) + 1]] <- z
  }
  val[[i]] <- fromJSON(file.path(d, 'validation.json'))
}
c <- rbindlist(cs)
m <- rbindlist(ms, fill = TRUE)
z <- rbindlist(ivs, fill = TRUE)
m[,
  q_candidate_family := p.adjust(
    p,
    'BH'
  ),
  by = platform
]
m[
  ,
  multiplicity_scope := 'preselected_candidate_assays_not_proteome_discovery'
]
m[
  platform == 'Olink',
  multiplicity_scope := 'prespecified_Olink_validation_family'
]
fwrite(c,
  file.path(
    out,
    'coloc_all_sensitivities.tsv'
  ),
  sep = '\t'
)
fwrite(c[width == 5e5 & p12 == 1e-5 & sdmode == 'estimated'],
  file.path(
    out,
    'coloc_primary.tsv'
  ),
  sep = '\t'
)
fwrite(m,
  file.path(
    out,
    'mr_candidate_assays.tsv'
  ),
  sep = '\t'
)
fwrite(z,
  file.path(
    out,
    'mr_all_instruments.tsv'
  ),
  sep = '\t'
)
fwrite(
  unique(z[
    ,
    .(SNP)
  ]),
  file.path(
    out,
    'instrument_rsids.txt'
  ),
  col.names = FALSE
)
write_json(val,
  file.path(
    out,
    'task_validations.json'
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)
core <- m[gene %in% c(
  'NMT1',
  'HDGF',
  'HEXIM1',
  'HEXIM2'
) & is.finite(beta)]
core[
  ,
  label := paste(gene,
    assay,
    platform,
    sep = ' / '
  )
]
core[
  ,
  label := factor(label,
    levels = rev(unique(label))
  )
]
fwrite(core,
  file.path(
    out,
    'forest_source.tsv'
  ),
  sep = '\t'
)
p <- ggplot(core, aes(beta, label, color = platform)) +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
  geom_segment(aes(x = lower, xend = upper, yend = label), linewidth = .6) +
  geom_point(size = 2) +
  theme_bw(base_size = 10) +
  labs(x = 'MR beta (continuous FI)', y = NULL) +
  theme(legend.position = 'bottom')
ggsave(
  file.path(
    out,
    'core_forest_candidate.png'
  ),
  p,
  width = 8,
  height = 6,
  dpi = 180
)
cp <- c[gene %in% c(
  'NMT1',
  'HDGF',
  'HEXIM1',
  'HEXIM2'
) & width == 5e5 & sdmode == 'estimated']
cp[
  ,
  label := paste(gene,
    assay,
    platform,
    sep = ' / '
  )
]
fwrite(cp,
  file.path(
    out,
    'coloc_sensitivity_source.tsv'
  ),
  sep = '\t'
)
p <- ggplot(cp, aes(p12, PP.H4.abf, color = label)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  ylim(0, 1) +
  theme_bw(base_size = 10) +
  labs(
    x = 'Shared prior p12',
    y = 'Posterior probability H4',
    color = 'Assay / platform'
  )
ggsave(
  file.path(
    out,
    'coloc_sensitivity_candidate.png'
  ),
  p,
  width = 9,
  height = 5,
  dpi = 180
)
write_json(
  list(
    assay_tasks = nrow(ts),
    completed = length(cs),
    coloc_rows = nrow(c),
    mr_estimable = sum(is.finite(m$beta)),
    posterior_sums_checked = TRUE,
    all_completion_markers = TRUE,
    scope = '15_prior_candidates_and_available_core_Olink_assays_not_full_discovery'
  ),
  file.path(
    out,
    'validation.json'
  ),
  auto_unbox = TRUE,
  pretty = TRUE
)
cat('REDUCER_COMPLETE\n')
