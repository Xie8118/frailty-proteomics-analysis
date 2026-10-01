# ============================================================================
# Test enrichment using platform-specific tested-gene backgrounds
# ============================================================================
# Purpose: Test enrichment using platform-specific tested-gene backgrounds.
# Inputs: MR foreground/background genes, GO and Reactome annotations.
# Outputs: Mapping tables, enrichment statistics and diagnostic plot.
# Arguments: none; inputs are read from the paths below.
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(jsonlite)
  library(ggplot2)
})
setDTthreads(1)
w <- 'work'
out <- file.path(
  w,
  'local_enrichment_v2'
)
stopifnot(
  file.exists(file.path(
    w,
    'full_mr_aggregate_v3/done'
  )),
  !dir.exists(out)
)
dir.create(out)
pub <- 'data/reference'
gmt <- file.path(
  pub,
  'Reactome/ReactomePathways.gmt'
)
paths <- lapply(
  strsplit(readLines(gmt),
    '\t',
    fixed = TRUE
  ),
  function(a) list(
    name = a[1],
    id = a[2],
    genes = unique(a[-(1:2)])
  )
)
paths <- Filter(
  function(z) grepl(
    '^R-HSA-',
    z$id
  ),
  paths
)
# GO names are display annotations only; term-to-gene mappings use the frozen local org.Hs.eg.db package.
ontology_labels <- file.path(
  pub,
  'GO/go-edit.obo'
)
obo <- if (file.exists(ontology_labels)) readLines(ontology_labels) else character()
namesGO <- list()
cur <- NULL
for (s in obo) {
  if (startsWith(
    s,
    'id: GO:'
  )) cur <- substring(
    s,
    5
  )
  if (startsWith(
    s,
    'name: '
  ) && !is.null(cur)) namesGO[[cur]] <- substring(
    s,
    7
  )
  if (s == '') cur <- NULL
}
metadata <- as.data.table(AnnotationDbi::metadata(org.Hs.eg.db))
fwrite(metadata,
  file.path(
    out,
    'GO_annotation_metadata.tsv'
  ),
  sep = '\t'
)
allres <- list()
qc <- list()

# Use each platform's tested genes as its enrichment background.
for (plat in c('Fenland', 'deCODE')) {
  bg <- unique(scan(
    file.path(
      w,
      'full_mr_aggregate_v3',
      paste0(
        plat,
        '_background_genes.txt'
      )
    ),
    what = '',
    quiet = TRUE
  ))
  fg <- unique(scan(
    file.path(
      w,
      'full_mr_aggregate_v3',
      paste0(
        plat,
        '_foreground_genes.txt'
      )
    ),
    what = '',
    quiet = TRUE
  ))
  stopifnot(all(fg %in% bg))
  valid <- intersect(
    bg,
    AnnotationDbi::keys(org.Hs.eg.db,
      keytype = 'SYMBOL'
    )
  )
  m <- as.data.table(AnnotationDbi::select(org.Hs.eg.db,
    keys = valid,
    keytype = 'SYMBOL',
    columns = 'ENTREZID'
  ))
  m[
    ,
    query := SYMBOL
  ]
  missing <- setdiff(
    bg,
    valid
  )
  aliases <- intersect(
    missing,
    AnnotationDbi::keys(org.Hs.eg.db,
      keytype = 'ALIAS'
    )
  )
  if (length(aliases)) {
    a <- as.data.table(AnnotationDbi::select(org.Hs.eg.db,
      keys = aliases,
      keytype = 'ALIAS',
      columns = c(
        'SYMBOL',
        'ENTREZID'
      )
    ))
    a[
      ,
      query := ALIAS
    ]
    m <- rbind(
      m,
      a[
        ,
        .(
          SYMBOL,
          ENTREZID,
          query
        )
      ]
    )
  }
  m <- unique(m[!is.na(ENTREZID)])
  amb <- m[,
    .(n = uniqueN(ENTREZID)),
    by = query
  ][
    n > 1,
    query
  ]
  m <- m[!query %in% amb]
  m <- unique(m,
    by = 'query'
  )
  map <- merge(data.table(query = bg),
    m,
    by = 'query',
    all.x = TRUE
  )
  map[
    ,
    analysis_id := fifelse(
      is.na(ENTREZID),
      paste0(
        'unmapped:',
        query
      ),
      ENTREZID
    )
  ]
  map[
    ,
    canonical_symbol := fifelse(
      is.na(SYMBOL),
      query,
      SYMBOL
    )
  ]
  map[
    ,
    foreground := query %in% fg
  ]
  fwrite(map,
    file.path(
      out,
      paste0(
        plat,
        '_gene_mapping.tsv'
      )
    ),
    sep = '\t'
  )
  B <- unique(map$analysis_id)
  F <- unique(map[
    foreground == TRUE,
    analysis_id
  ])
  BS <- unique(map$canonical_symbol)
  FS <- unique(map[
    foreground == TRUE,
    canonical_symbol
  ])
  validids <- unique(na.omit(map$ENTREZID))
  go <- as.data.table(AnnotationDbi::select(org.Hs.eg.db,
    keys = validids,
    keytype = 'ENTREZID',
    columns = c(
      'GOALL',
      'ONTOLOGYALL'
    )
  ))
  go <- unique(go[
    ONTOLOGYALL == 'BP' & !is.na(GOALL),
    .(
      id = GOALL,
      gene = ENTREZID
    )
  ])
  goSet <- split(
    go$gene,
    go$id
  )

  # Hypergeometric enrichment retains all background-represented terms.
  test <- function(id, label, genes, back, front, source) {
    genes <- intersect(
      unique(genes),
      back
    )
    K <- length(genes)
    if (!K) return(NULL)
    hits <- intersect(
      front,
      genes
    )
    data.table(
      platform = plat,
      source = source,
      id = id,
      name = label,
      background_n = length(back),
      foreground_n = length(front),
      term_in_background = K,
      overlap = length(hits),
      genes = paste(sort(hits),
        collapse = ';'
      ),
      p = phyper(length(hits) - 1,
        K,
        length(back) - K,
        length(front),
        lower.tail = FALSE
      )
    )
  }
  g <- rbindlist(lapply(
    names(goSet),
    function(id) test(
      id,
      if (is.null(namesGO[[id]])) id else namesGO[[id]],
      goSet[[id]],
      B,
      F,
      'GO_BP'
    )
  ))
  r <- rbindlist(lapply(
    paths,
    function(z) test(
      z$id,
      z$name,
      z$genes,
      BS,
      FS,
      'Reactome'
    )
  ))
  res <- rbind(g,
    r,
    fill = TRUE
  )
  res[,
    q := p.adjust(
      p,
      'BH'
    ),
    by = source
  ]
  stopifnot(
    all(is.finite(res$p)),
    all(res$p >= 0 & res$p <= 1),
    all(res$q >= res$p - 1e-12)
  )
  fwrite(res,
    file.path(
      out,
      paste0(
        plat,
        '_enrichment.tsv'
      )
    ),
    sep = '\t'
  )
  allres[[plat]] <- res
  qc[[plat]] <- list(
    platform = plat,
    background_raw = length(bg),
    background_unique_ids = length(B),
    foreground_raw = length(fg),
    foreground_unique_ids = length(F),
    GO_background_annotated = uniqueN(go$gene),
    unmapped = map[
      is.na(ENTREZID),
      query
    ],
    ambiguous_alias = amb,
    GO_terms = nrow(g),
    Reactome_terms = nrow(r),
    significant = res[q < .05,
      .N,
      by = source
    ]
  )
  writeLines(
    sort(FS),
    file.path(
      out,
      paste0(
        plat,
        '_PPI_foreground.txt'
      )
    )
  )
}

# Export descriptive panels; significance remains defined by adjusted P values.
res <- rbindlist(allres)
top <- res[overlap > 0][order(p)][,
  head(
    .SD,
    8
  ),
  by = .(
    platform,
    source
  )
]
fwrite(top,
  file.path(
    out,
    'enrichment_panel_source.tsv'
  ),
  sep = '\t'
)
if (nrow(top)) {
  top[
    ,
    label := paste0(
      substr(
        name,
        1,
        65
      ),
      ' (',
      overlap,
      ')'
    )
  ]
  top[
    ,
    label := factor(label,
      levels = rev(unique(label))
    )
  ]
  g <- ggplot(
    top,
    aes(
      -log10(pmax(
        q,
        1e-300
      )),
      label,
      colour = q < .05
    )
  ) +
    geom_point(size = 2.3) +
    facet_wrap(~ platform + source, ncol = 2, scales = 'free_y') +
    scale_colour_manual(values = c('FALSE' = '#888888', 'TRUE' = '#0072B2'), name = 'FDR < 0.05') +
    theme_bw(base_size = 9) +
    labs(
      x = expression(-log[10](FDR)),
      y = NULL,
      title = 'Enrichment of primary cis MR candidates',
      caption = 'Platform-specific tested-gene background; BH within annotation family. Grey terms do not meet FDR <0.05.'
    )
  ggsave(file.path(out, 'enrichment_candidate.png'), g, width = 13, height = 11, dpi = 220)
} else writeLines(
  'No foreground overlaps with tested annotation terms.',
  file.path(
    out,
    'no_panel_rationale.txt'
  )
)
write_json(
  list(
    platforms = qc,
    GO_package = as.character(packageVersion('org.Hs.eg.db')),
    GO_display_labels = if (length(namesGO)) 'official_ontology' else 'stable_GO_identifiers_term_name_source_unavailable',
    Reactome_version = 97,
    no_external_candidate_queries = TRUE,
    methods = 'hypergeometric; all background-represented terms including zero overlaps; BH by platform and source; unmapped background entries counted but unannotated'
  ),
  file.path(
    out,
    'validation.json'
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)
writeLines(
  'LOCAL_ENRICHMENT_COMPLETE',
  file.path(
    out,
    'done'
  )
)
