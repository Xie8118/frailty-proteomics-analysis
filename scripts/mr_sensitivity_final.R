# ============================================================================
# Assess MR heterogeneity, directionality and influence
# ============================================================================
# Purpose: Assess MR heterogeneity, directionality and influence.
# Inputs: final_candidate_inputs/tasks.tsv and instruments.tsv.
# Outputs: Per-task sensitivity statistics, leave-one-out and directionality tables.
# Arguments: one-based task ID and mode (pilot or full).
# Run from the repository root; see README.md for prerequisites.
# This script requires prepared study inputs; these are not bundled.
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
  library(jsonlite)
})
setDTthreads(1)
args <- commandArgs(TRUE)
tid <- as.integer(args[1])
mode <- args[2]
w <- 'work'
task <- fread(file.path(w, 'final_candidate_inputs/tasks.tsv'))[task_id == tid]
stopifnot(nrow(task) == 1)
out <- file.path(
  w,
  if (mode == 'pilot') 'sensitivity_final_pilot' else 'sensitivity_final',
  paste0(
    'task_',
    tid
  )
)
stopifnot(!dir.exists(out))
dir.create(out,
  recursive = TRUE
)
z <- fread(file.path(w, 'final_candidate_inputs/instruments.tsv'))[task_id == tid]
n <- nrow(z)
res <- data.table(
  task_id = tid,
  platform = task$platform,
  gene = task$gene,
  assay = task$assay,
  nsnp = n
)
res[
  ,
  `:=`(
    egger_status = if (n >= 3) 'pending' else 'fewer_than_3_instruments',
    loo_status = if (n >= 2) 'pending' else 'fewer_than_2_instruments'
  )
]

# Compute fixed-effect IVW and heterogeneity from harmonized instruments.
if (n) {
  bx <- z$BETA.x
  by <- z$BETA.y.aligned
  sx <- z$SE.x
  sy <- z$SE.y
  stopifnot(all(is.finite(bx)), all(sy > 0), all(sx > 0), all(z$N.x > 3), all(z$N.y > 3))
  den <- sum(bx^2 / sy^2)
  num <- sum(bx * by / sy^2)
  b <- num / den
  se <- sqrt(1 / den)
  Q <- sum((by - b * bx)^2 / sy^2)
  res[
    ,
    `:=`(
      beta = b,
      se = se,
      Q = Q,
      Qp = if (n > 1) pchisq(Q,
        n - 1,
        lower.tail = FALSE
      ) else NA_real_,
      Qdf = n - 1
    )
  ]
  # Correlation from t statistics; N-2 approximation treats covariate-adjusted GWAS as marginal summary data.
  rx2 <- bx^2 / (bx^2 + (z$N.x - 2) * sx^2)
  ry2 <- by^2 / (by^2 + (z$N.y - 2) * sy^2)
  steiger <- data.table(
    SNP = z$SNP,
    R2_X_approx = rx2,
    R2_Y_approx = ry2,
    N_X = z$N.x,
    N_Y = z$N.y,
    direction_X_to_Y = rx2 > ry2
  )
  fwrite(steiger, file.path(out, 'steiger_variants.tsv'), sep = '\t')
  res[
    ,
    `:=`(
      R2_X_sum_approx = sum(rx2),
      R2_Y_sum_approx = sum(ry2),
      steiger_direction_approx = sum(rx2) > sum(ry2),
      steiger_p_status = 'not_reported_unknown_overlap_and_summary_R2_approximation'
    )
  ]

  # Leave-one-out estimates require at least two instruments.
  if (n >= 2) {
    # Closed-form deletion diagnostics, not stochastic resampling.
    db <- den - bx^2 / sy^2
    nb <- num - bx * by / sy^2
    stopifnot(all(db > 0))
    loo <- data.table(
      removed_SNP = z$SNP,
      beta = nb / db,
      se = sqrt(1 / db)
    )
    loo[
      ,
      `:=`(
        lower = beta - 1.96 * se,
        upper = beta + 1.96 * se,
        p = 2 * pnorm(-abs(beta / se))
      )
    ]
    fwrite(loo,
      file.path(
        out,
        'leave_one_out.tsv'
      ),
      sep = '\t'
    )
    check <- sum(bx[-1] * by[-1] / sy[-1]^2) / sum(bx[-1]^2 / sy[-1]^2)
    stopifnot(abs(check - loo$beta[1]) < 1e-10)
    res[
      ,
      `:=`(
        loo_status = 'computed',
        loo_beta_min = min(loo$beta),
        loo_beta_max = max(loo$beta),
        loo_any_direction_change = any(sign(loo$beta) != sign(b))
      )
    ]
  }

  # MR-Egger requires enough instruments and variation in exposure effects.
  if (n >= 3 && var(abs(bx)) > 1e-12) {
    obj <- mr_input(bx = bx, bxse = sx, by = by, byse = sy, snps = z$SNP)
    eg <- mr_egger(obj,
      distribution = 'normal'
    )
    xx <- abs(bx)
    yy <- by * sign(bx)
    fit <- lm(yy ~ xx,
      weights = 1 / sy^2
    )
    stopifnot(
      abs(unname(coef(fit)[2]) - eg@Estimate) < 1e-9,
      abs(unname(coef(fit)[1]) - eg@Intercept) < 1e-9
    )
    res[
      ,
      `:=`(
        egger_status = 'computed',
        egger_beta = eg@Estimate,
        egger_se = eg@StdError.Est,
        egger_p = eg@Pvalue.Est,
        egger_lower = eg@CILower.Est,
        egger_upper = eg@CIUpper.Est,
        egger_intercept = eg@Intercept,
        egger_intercept_se = eg@StdError.Int,
        egger_intercept_p = eg@Pvalue.Int,
        egger_I2GX = eg@I.sq
      )
    ]
  } else if (n >= 3) res[, egger_status := 'insufficient_exposure_effect_variation']
}
fwrite(res, file.path(out, 'sensitivity.tsv'), sep = '\t')
write_json(
  list(
    task_id = tid,
    nsnp = n,
    package = as.character(packageVersion('MendelianRandomization')),
    independent_egger_check = if (n >= 3 && var(abs(z$BETA.x)) > 1e-12) TRUE else NA,
    deterministic = TRUE,
    seed_policy = 'no_simulation_or_bootstrap',
    role = 'sensitivity_not_primary_selection'
  ),
  file.path(
    out,
    'validation.json'
  ),
  auto_unbox = TRUE,
  pretty = TRUE
)
writeLines('SENSITIVITY_COMPLETE', file.path(out, 'done'))
