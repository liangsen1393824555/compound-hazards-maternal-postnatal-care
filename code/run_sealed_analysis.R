args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: run_sealed_analysis.R <restricted-dir> <sealed-dir>")

Sys.setlocale("LC_ALL", "C")
Sys.setenv(TZ = "Asia/Shanghai")
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(jsonlite))

input_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
sealed_dir <- args[[2]]
dir.create(sealed_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260911)
ck_levels <- c("00", "10", "01", "11")

as_num_list <- function(x) {
  y <- as.list(as.numeric(x)); names(y) <- names(x); y
}

quantiles <- function(x) {
  values <- stats::quantile(as.numeric(x), probs = c(0, .01, .25, .5, .75, .99, 1), na.rm = TRUE, names = FALSE)
  names(values) <- c("min", "p01", "p25", "p50", "p75", "p99", "max")
  as_num_list(values)
}

prepare_data <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  for (v in c("CK_A25", "CK_A50", "CK_A100")) d[[v]] <- sprintf("%02d", as.integer(d[[v]]))
  d$age_group <- factor(d$age_group, levels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"))
  d$pord_group <- factor(d$pord_group, levels = c("1", "2-4", "5+"))
  d$outcome_year <- factor(d$outcome_year)
  d$calendar_month <- factor(sprintf("%02d", as.integer(d$calendar_month)), levels = sprintf("%02d", 1:12))
  d$outcome_month <- factor(d$outcome_month)
  d$v024 <- factor(d$v024)
  d$v025 <- factor(d$v025)
  d$v022_fe <- factor(d$v022)
  d$v106 <- factor(d$v106)
  d$v190 <- factor(d$v190)
  d
}

design_full <- function(d, weight_name = "weight", id_name = "v021", strata_name = "v022") {
  svydesign(
    ids = as.formula(paste0("~", id_name)),
    strata = as.formula(paste0("~", strata_name)),
    weights = as.formula(paste0("~", weight_name)),
    data = d, nest = TRUE
  )
}

support_gate <- function(d, ck_name, extra = rep(TRUE, nrow(d))) {
  use <- !is.na(d$Y) & extra
  groups <- lapply(ck_levels, function(level) {
    z <- d[use & d[[ck_name]] == level, , drop = FALSE]
    counts <- table(factor(z$Y, levels = c(0, 1)))
    list(
      ck = level,
      outcome_known_n = nrow(z),
      independent_clusters = length(unique(z$v021)),
      smaller_outcome_category_n = min(as.integer(counts)),
      pass = nrow(z) >= 100L && length(unique(z$v021)) >= 50L && min(as.integer(counts)) >= 30L
    )
  })
  names(groups) <- ck_levels
  list(pass = all(vapply(groups, function(x) x$pass, logical(1))), cells = groups)
}

overlap_gate <- function(d, ck_name, dimensions, use) {
  total <- sum(d$weight[use])
  by_dimension <- lapply(dimensions, function(v) {
    levels_v <- unique(d[[v]][use & !is.na(d[[v]])])
    supported <- levels_v[vapply(levels_v, function(a) all(ck_levels %in% unique(d[[ck_name]][use & d[[v]] == a])), logical(1))]
    coverage <- sum(d$weight[use & d[[v]] %in% supported]) / total
    list(levels = length(levels_v), levels_with_all_four_ck = length(supported), weighted_coverage = unname(coverage), pass = coverage >= .80)
  })
  names(by_dimension) <- dimensions
  list(pass = all(vapply(by_dimension, function(x) x$pass, logical(1))), threshold = .80, dimensions = by_dimension)
}

model_formula <- function(strict = FALSE, wealth = FALSE) {
  rhs <- if (strict) {
    c("CK", "v022_fe", "outcome_month", "age_group", "pord_group", "v106")
  } else {
    c("CK", "outcome_year", "calendar_month", "v024", "v025", "age_group", "pord_group", "v106")
  }
  if (wealth) rhs <- c(rhs, "v190")
  as.formula(paste("Y ~", paste(rhs, collapse = " + ")))
}

standardize_levels <- function(fit, target, variable, levels_value, link = "logit") {
  beta <- coef(fit)
  variance <- vcov(fit)
  w <- target$weight
  risks <- numeric(length(levels_value))
  jac <- matrix(NA_real_, nrow = length(levels_value), ncol = length(beta), dimnames = list(levels_value, names(beta)))
  for (i in seq_along(levels_value)) {
    newdata <- target
    template <- target[[variable]]
    if (is.factor(template)) newdata[[variable]] <- factor(levels_value[[i]], levels = levels(template)) else newdata[[variable]] <- levels_value[[i]]
    mm <- model.matrix(delete.response(terms(fit)), newdata, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
    mm <- mm[, names(beta), drop = FALSE]
    eta <- drop(mm %*% beta)
    prediction <- if (link == "logit") plogis(eta) else eta
    risks[[i]] <- sum(w * prediction) / sum(w)
    slope <- if (link == "logit") prediction * (1 - prediction) else rep(1, length(prediction))
    jac[i, ] <- colSums(mm * (w * slope)) / sum(w)
  }
  names(risks) <- levels_value
  risk_cov <- jac %*% variance %*% t(jac)
  list(risks = risks, jacobian = jac, covariance = risk_cov)
}

estimate_linear <- function(std, weights, df, ci_level = .95) {
  est <- sum(weights * std$risks[names(weights)])
  grad <- drop(weights %*% std$jacobian[names(weights), , drop = FALSE])
  se <- sqrt(drop(t(grad) %*% vcov_current %*% grad))
  crit <- qt((1 + ci_level) / 2, df = df)
  p <- 2 * pt(-abs(est / se), df = df)
  list(estimate = unname(est), standard_error = unname(se), ci_level = ci_level, conf_low = unname(est - crit * se), conf_high = unname(est + crit * se), p_value = unname(p))
}

estimate_logratio <- function(std, numerator, denominator, df) {
  rn <- std$risks[[numerator]]; rd <- std$risks[[denominator]]
  est <- rn / rd
  grad_log <- std$jacobian[numerator, ] / rn - std$jacobian[denominator, ] / rd
  se_log <- sqrt(drop(t(grad_log) %*% vcov_current %*% grad_log))
  crit <- qt(.975, df = df)
  list(estimate = unname(est), standard_error_log = unname(se_log), conf_low = unname(exp(log(est) - crit * se_log)), conf_high = unname(exp(log(est) + crit * se_log)))
}

estimate_multiplicative_interaction <- function(std, df) {
  r <- std$risks
  est <- r[["11"]] * r[["00"]] / (r[["10"]] * r[["01"]])
  grad_log <- std$jacobian["11", ] / r[["11"]] + std$jacobian["00", ] / r[["00"]] - std$jacobian["10", ] / r[["10"]] - std$jacobian["01", ] / r[["01"]]
  se_log <- sqrt(drop(t(grad_log) %*% vcov_current %*% grad_log))
  crit <- qt(.975, df = df)
  list(estimate = unname(est), standard_error_log = unname(se_log), conf_low = unname(exp(log(est) - crit * se_log)), conf_high = unname(exp(log(est) + crit * se_log)))
}

diagnostics <- function(fit, target, design, warning_text, std) {
  mm <- model.matrix(fit)
  infl <- attr(fit, "influence")
  psu_influence <- NULL
  if (!is.null(infl)) {
    psu_codes <- target$v021
    aggregated <- rowsum(infl, group = psu_codes, reorder = FALSE)
    psu_influence <- quantiles(sqrt(rowSums(aggregated^2)))
  }
  residual <- residuals(fit, type = "pearson")
  list(
    converged = isTRUE(fit$converged), boundary = isTRUE(fit$boundary), iterations = fit$iter,
    model_columns = ncol(mm), model_rank = fit$rank, full_rank = fit$rank == ncol(mm),
    covariance_rank = qr(vcov(fit))$rank, coefficient_count = length(coef(fit)),
    warnings = warning_text, coefficient_abs_max = max(abs(coef(fit))), standard_error_max = max(sqrt(diag(vcov(fit)))),
    standardized_prediction_min = min(std$risks), standardized_prediction_max = max(std$risks),
    pearson_residual_quantiles = quantiles(residual), anonymous_psu_influence_norm_quantiles = psu_influence,
    design_degrees_of_freedom = degf(design)
  )
}

run_ck_model <- function(d, ck_name, tag, strict = FALSE, wealth = FALSE, y_override = NULL, subset_extra = rep(TRUE, nrow(d))) {
  if (!is.null(y_override)) d$Y <- y_override
  d$CK <- factor(d[[ck_name]], levels = ck_levels)
  required <- c("Y", "CK", "weight", "v021", "v022", "age_group", "pord_group", "v106")
  if (strict) required <- c(required, "v022_fe", "outcome_month") else required <- c(required, "outcome_year", "calendar_month", "v024", "v025")
  if (wealth) required <- c(required, "v190")
  use <- complete.cases(d[, unique(required), drop = FALSE]) & subset_extra
  full <- design_full(d)
  des <- subset(full, use)
  target <- d[use, , drop = FALSE]
  form <- model_formula(strict, wealth)
  warning_text <- character()
  fit <- withCallingHandlers(
    tryCatch(svyglm(form, design = des, family = quasibinomial(), influence = TRUE), error = function(e) e),
    warning = function(w) { warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  logistic_ok <- !inherits(fit, "error") && isTRUE(fit$converged) && fit$rank == length(coef(fit)) && all(is.finite(coef(fit))) && all(is.finite(vcov(fit))) && !any(grepl("0 or 1|did not converge|singular", warning_text, ignore.case = TRUE))
  model_type <- "survey_weighted_logistic"
  if (!logistic_ok) {
    warning_text <- c(warning_text, if (inherits(fit, "error")) conditionMessage(fit) else "logistic_failure_rule_triggered")
    fit <- withCallingHandlers(
      tryCatch(svyglm(form, design = des, family = gaussian(), influence = TRUE), error = function(e) e),
      warning = function(w) { warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning") }
    )
    if (inherits(fit, "error") || !isTRUE(fit$converged) || fit$rank != length(coef(fit)) || degf(des) <= 0) stop(paste("LPM fallback failed for", tag))
    model_type <- "survey_weighted_lpm_fallback"
  }
  std <- standardize_levels(fit, target, "CK", ck_levels, if (model_type == "survey_weighted_logistic") "logit" else "identity")
  assign("vcov_current", vcov(fit), envir = .GlobalEnv)
  df <- degf(des)
  risk_items <- lapply(ck_levels, function(level) {
    se <- sqrt(std$covariance[level, level]); crit <- qt(.975, df)
    list(estimate = unname(std$risks[[level]]), standard_error = unname(se), conf_low = unname(std$risks[[level]] - crit * se), conf_high = unname(std$risks[[level]] + crit * se))
  }); names(risk_items) <- ck_levels
  rd_weights <- list(
    RD_10_00 = c("00" = -1, "10" = 1, "01" = 0, "11" = 0),
    RD_01_00 = c("00" = -1, "10" = 0, "01" = 1, "11" = 0),
    RD_11_00 = c("00" = -1, "10" = 0, "01" = 0, "11" = 1),
    IC_add = c("00" = 1, "10" = -1, "01" = -1, "11" = 1)
  )
  contrasts <- lapply(rd_weights, function(w) estimate_linear(std, w, df, .95))
  primary <- contrasts[c("RD_11_00", "IC_add")]
  p_raw <- vapply(primary, function(x) x$p_value, numeric(1))
  p_holm <- p.adjust(p_raw, method = "holm")
  for (i in seq_along(primary)) {
    primary[[i]]$p_holm <- unname(p_holm[[i]])
    simultaneous <- estimate_linear(std, rd_weights[[names(primary)[i]]], df, .975)
    primary[[i]]$simultaneous_97_5_low <- simultaneous$conf_low
    primary[[i]]$simultaneous_97_5_high <- simultaneous$conf_high
  }
  contrasts[c("RD_11_00", "IC_add")] <- primary
  ratios <- NULL
  if (model_type == "survey_weighted_logistic") {
    ratios <- list(
      RR_10_00 = estimate_logratio(std, "10", "00", df),
      RR_01_00 = estimate_logratio(std, "01", "00", df),
      RR_11_00 = estimate_logratio(std, "11", "00", df),
      ratio_scale_interaction = estimate_multiplicative_interaction(std, df)
    )
  }
  coef_table <- data.frame(tag = tag, term = names(coef(fit)), estimate = unname(coef(fit)), standard_error = sqrt(diag(vcov(fit))), stringsAsFactors = FALSE)
  cov_table <- as.data.frame(as.table(vcov(fit)), stringsAsFactors = FALSE)
  names(cov_table) <- c("row_term", "column_term", "covariance")
  cov_table$tag <- tag
  cov_table <- cov_table[, c("tag", "row_term", "column_term", "covariance")]
  list(
    result = list(tag = tag, model_type = model_type, n = nrow(target), risks = risk_items, contrasts = contrasts, ratios = ratios,
                  diagnostics = diagnostics(fit, target, des, unique(warning_text), std)),
    coefficients = coef_table,
    covariance = cov_table,
    estimates = data.frame(tag = tag, estimand = c(paste0("R_", ck_levels), names(contrasts)), estimate = c(std$risks, vapply(contrasts, function(x) x$estimate, numeric(1))), stringsAsFactors = FALSE)
  )
}

run_direct_cells <- function(d, country_tag) {
  d$CK <- factor(d$CK_A50, levels = ck_levels)
  use <- !is.na(d$Y) & !is.na(d$CK)
  des <- subset(design_full(d), use)
  target <- d[use, , drop = FALSE]
  warning_text <- character()
  fit <- withCallingHandlers(svyglm(Y ~ 0 + CK, design = des, family = quasibinomial()), warning = function(w) { warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning") })
  beta <- coef(fit); vv <- vcov(fit); risks <- plogis(beta); deriv <- risks * (1 - risks); rcov <- diag(deriv) %*% vv %*% diag(deriv)
  names(risks) <- sub("^CK", "", names(beta)); rownames(rcov) <- colnames(rcov) <- names(risks)
  df <- degf(des); crit <- qt(.975, df)
  cells <- lapply(ck_levels, function(level) {se <- sqrt(rcov[level, level]); list(estimate = unname(risks[[level]]), standard_error = unname(se), conf_low = unname(risks[[level]] - crit * se), conf_high = unname(risks[[level]] + crit * se))}); names(cells) <- ck_levels
  w <- c("00" = -1, "10" = 0, "01" = 0, "11" = 1); est <- sum(w * risks); se <- sqrt(drop(t(w) %*% rcov %*% w))
  list(country = country_tag, method = "direct_survey_weighted_cell_means_via_saturated_cells", n = nrow(target), cells = cells,
       rd_11_00 = list(estimate = unname(est), standard_error = unname(se), conf_low = unname(est - crit * se), conf_high = unname(est + crit * se), p_value_reported = FALSE),
       design_degrees_of_freedom = df, warnings = unique(warning_text))
}

run_bf_reduced <- function(d) {
  d$C <- factor(d$C_A, levels = c(0, 1)); d$K <- factor(d$K_A50, levels = c(0, 1))
  required <- c("Y", "C", "K", "weight", "v021", "v022", "outcome_year", "calendar_month", "v024", "v025", "age_group", "pord_group", "v106")
  use <- complete.cases(d[, required])
  des <- subset(design_full(d), use); target <- d[use, , drop = FALSE]
  form <- Y ~ C + K + outcome_year + calendar_month + v024 + v025 + age_group + pord_group + v106
  warning_text <- character()
  fit <- withCallingHandlers(svyglm(form, design = des, family = quasibinomial()), warning = function(w) { warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (!isTRUE(fit$converged) || fit$rank != length(coef(fit))) stop("Burkina Faso reduced C+K model failed")
  assign("vcov_current", vcov(fit), envir = .GlobalEnv)
  df <- degf(des)
  summarize_binary <- function(variable) {
    std <- standardize_levels(fit, target, variable, c("0", "1"), "logit")
    names(std$risks) <- rownames(std$jacobian) <- c("0", "1")
    rd <- estimate_linear(std, c("0" = -1, "1" = 1), df, .95)
    rr <- estimate_logratio(std, "1", "0", df)
    list(risk_0 = unname(std$risks[["0"]]), risk_1 = unname(std$risks[["1"]]), risk_difference = rd, risk_ratio = rr)
  }
  list(model = "survey_weighted_logistic_C_plus_K_no_interaction", n = nrow(target), C = summarize_binary("C"), K = summarize_binary("K"),
       diagnostics = list(converged = isTRUE(fit$converged), model_rank = fit$rank, model_columns = ncol(model.matrix(fit)), covariance_rank = qr(vcov(fit))$rank, design_degrees_of_freedom = df, warnings = unique(warning_text)),
       coefficients = data.frame(tag = "burkina_reduced_C_plus_K", term = names(coef(fit)), estimate = unname(coef(fit)), standard_error = sqrt(diag(vcov(fit))), stringsAsFactors = FALSE),
       covariance = transform(as.data.frame(as.table(vcov(fit)), stringsAsFactors = FALSE), tag = "burkina_reduced_C_plus_K"))
}

run_pooled_direct <- function(ng, bf) {
  d <- rbind(ng, bf)
  d$CK <- factor(d$CK_A50, levels = ck_levels); d$country <- factor(d$country, levels = c("nigeria", "burkina_faso"))
  use <- !is.na(d$Y) & !is.na(d$CK)
  totals <- tapply(d$weight[use], d$country[use], sum)
  d$pooled_weight <- .5 * d$weight / totals[as.character(d$country)]
  d$pooled_psu <- interaction(d$country, d$v021, drop = TRUE)
  d$pooled_strata <- interaction(d$country, d$v022, drop = TRUE)
  des <- subset(design_full(d, "pooled_weight", "pooled_psu", "pooled_strata"), use)
  fit <- svyglm(Y ~ 0 + country:CK, design = des, family = quasibinomial())
  beta <- coef(fit); vv <- vcov(fit); cell_risk <- plogis(beta); deriv <- cell_risk * (1-cell_risk); rcov <- diag(deriv) %*% vv %*% diag(deriv)
  risk_names <- names(beta)
  avg <- numeric(4); j <- matrix(0, 4, length(beta), dimnames = list(ck_levels, risk_names))
  for (i in seq_along(ck_levels)) {
    hits <- grep(paste0("CK", ck_levels[[i]], "$"), risk_names)
    avg[[i]] <- mean(cell_risk[hits]); j[i, hits] <- .5
  }
  names(avg) <- ck_levels; avg_cov <- j %*% rcov %*% t(j); df <- degf(des); crit <- qt(.975, df)
  cells <- lapply(ck_levels, function(level) {se <- sqrt(avg_cov[level, level]); list(estimate = unname(avg[[level]]), standard_error = unname(se), conf_low = unname(avg[[level]] - crit*se), conf_high = unname(avg[[level]] + crit*se))}); names(cells) <- ck_levels
  w <- c("00"=-1,"10"=0,"01"=0,"11"=1); rd <- sum(w*avg); se <- sqrt(drop(t(w)%*%avg_cov%*%w))
  list(method = "equal_country_weight_saturated_direct_CK", hypothesis_tests = FALSE, country_total_weights = as_num_list(tapply(d$pooled_weight[use], d$country[use], sum)),
       cells = cells, rd_11_00 = list(estimate=unname(rd), standard_error=unname(se), conf_low=unname(rd-crit*se), conf_high=unname(rd+crit*se)),
       design_degrees_of_freedom = df, model_rank = fit$rank, model_columns = ncol(model.matrix(fit)))
}

ng <- prepare_data(file.path(input_dir, "restricted_nigeria.csv"))
bf <- prepare_data(file.path(input_dir, "restricted_burkina_faso.csv"))
sample_flow <- list(
  nigeria = list(eligible_n = nrow(ng), outcome_known_n = sum(!is.na(ng$Y)), outcome_unknown_n = sum(is.na(ng$Y)), complete_case_n = sum(complete.cases(ng[, c("Y","CK_A50","weight","v021","v022","v024","v025","age_group","pord_group","v106")]))),
  burkina_faso = list(eligible_n = nrow(bf), outcome_known_n = sum(!is.na(bf$Y)), outcome_unknown_n = sum(is.na(bf$Y)), complete_case_n = sum(complete.cases(bf[, c("Y","CK_A50","weight","v021","v022","v024","v025","age_group","pord_group","v106")])) )
)

gates <- list(
  nigeria_primary = support_gate(ng, "CK_A50"),
  burkina_primary = support_gate(bf, "CK_A50"),
  nigeria_radius_25 = support_gate(ng, "CK_A25"),
  nigeria_radius_100 = support_gate(ng, "CK_A100"),
  nigeria_live_birth = support_gate(ng, "CK_A50", ng$m80 == 1),
  nigeria_strict_FE = support_gate(ng, "CK_A50"),
  nigeria_v190 = support_gate(ng, "CK_A50")
)
if (!gates$nigeria_primary$pass) stop("Primary Nigeria support gate failed")

fits <- list(); coef_tables <- list(); covariance_tables <- list(); estimate_tables <- list()
add_fit <- function(name, value) {
  fits[[name]] <<- value$result
  coef_tables[[name]] <<- value$coefficients
  covariance_tables[[name]] <<- value$covariance
  estimate_tables[[name]] <<- value$estimates
}
add_fit("nigeria_primary", run_ck_model(ng, "CK_A50", "nigeria_primary"))

unknown_one <- ng$Y; unknown_one[is.na(unknown_one)] <- 1
unknown_zero <- ng$Y; unknown_zero[is.na(unknown_zero)] <- 0
add_fit("nigeria_unknown_as_1", run_ck_model(ng, "CK_A50", "nigeria_unknown_as_1", y_override = unknown_one))
add_fit("nigeria_unknown_as_0", run_ck_model(ng, "CK_A50", "nigeria_unknown_as_0", y_override = unknown_zero))

if (gates$nigeria_radius_25$pass) add_fit("nigeria_radius_25", run_ck_model(ng, "CK_A25", "nigeria_radius_25"))
if (gates$nigeria_radius_100$pass) add_fit("nigeria_radius_100", run_ck_model(ng, "CK_A100", "nigeria_radius_100"))
if (gates$nigeria_live_birth$pass) add_fit("nigeria_live_birth", run_ck_model(ng, "CK_A50", "nigeria_live_birth", subset_extra = ng$m80 == 1))
if (gates$nigeria_strict_FE$pass) add_fit("nigeria_strict_FE", run_ck_model(ng, "CK_A50", "nigeria_strict_FE", strict = TRUE))
if (gates$nigeria_v190$pass) add_fit("nigeria_v190", run_ck_model(ng, "CK_A50", "nigeria_v190", wealth = TRUE))

bf_direct <- run_direct_cells(transform(bf, CK = factor(CK_A50, levels = ck_levels)), "burkina_faso")
bf_reduced <- run_bf_reduced(bf)
coef_tables$burkina_reduced_C_plus_K <- bf_reduced$coefficients
names(bf_reduced$covariance)[1:3] <- c("row_term", "column_term", "covariance")
bf_reduced$covariance <- bf_reduced$covariance[, c("tag", "row_term", "column_term", "covariance")]
covariance_tables$burkina_reduced_C_plus_K <- bf_reduced$covariance
bf_reduced$coefficients <- NULL
bf_reduced$covariance <- NULL
pooled <- run_pooled_direct(ng, bf)

design_diagnostics <- lapply(list(nigeria = ng, burkina_faso = bf), function(d) {
  use <- complete.cases(d[, c("Y","CK_A50","weight","v021","v022","v024","v025","age_group","pord_group","v106")])
  z <- d[use,]
  psu_by_strata <- tapply(z$v021, z$v022, function(x) length(unique(x)))
  list(n = nrow(z), psus = length(unique(z$v021)), strata = length(unique(z$v022)), singleton_strata = sum(psu_by_strata == 1), design_degrees_of_freedom = length(unique(z$v021))-length(unique(z$v022)),
       weight_quantiles = quantiles(z$weight), effective_sample_size = unname(sum(z$weight)^2/sum(z$weight^2)),
       overlap_primary = overlap_gate(d, "CK_A50", c("v024","outcome_year","calendar_month"), use))
})

derivation_manifest <- fromJSON(file.path(input_dir, "derivation_manifest.json"), simplifyVector = FALSE)
payload <- list(
  state = "MODEL_OUTPUT_SEALED", seed = 20260911,
  survey_options = list(lonely_psu = getOption("survey.lonely.psu"), adjust_domain_lonely = isTRUE(getOption("survey.adjust.domain.lonely"))),
  sample_flow = sample_flow, support_gates = gates, design_diagnostics = design_diagnostics,
  nigeria_models = fits, burkina_direct = bf_direct, burkina_reduced_C_plus_K = bf_reduced,
  pooled_equal_country_direct = pooled,
  scheme_b_sensitivity = derivation_manifest$scheme_b_sensitivity,
  displacement_sensitivity = derivation_manifest$displacement_sensitivity,
  restrictions = list(burkina_adjusted_CxK_fitted = FALSE, pooled_adjusted_common_slope_fitted = FALSE, ordinary_robust_SE_used = FALSE, stepwise_or_significance_selection_used = FALSE, high_influence_PSU_deleted = FALSE)
)
write_json(payload, file.path(sealed_dir, "sealed_results.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16, na = "null")
write.csv(do.call(rbind, coef_tables), file.path(sealed_dir, "sealed_coefficients.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(do.call(rbind, covariance_tables), file.path(sealed_dir, "sealed_covariances.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(do.call(rbind, estimate_tables), file.path(sealed_dir, "sealed_estimates.csv"), row.names = FALSE, quote = TRUE, na = "")
writeLines(capture.output(sessionInfo()), file.path(sealed_dir, "R_sessionInfo.txt"), useBytes = TRUE)
