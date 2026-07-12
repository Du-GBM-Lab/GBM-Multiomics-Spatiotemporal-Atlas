# Shared helpers for R4 STOP5.
# Cohorts are scored independently; never merge cohorts before ssGSEA.

suppressPackageStartupMessages({
  library(GSVA)
  library(HGNChelper)
  library(survival)
  library(ggplot2)
})

PROGS <- c("NPC_P", "OPC_M", "MES_V", "MES_I")

to_current <- function(x) {
  res <- suppressWarnings(suppressMessages(
    HGNChelper::checkGeneSymbols(as.character(x), species = "human")
  ))
  s <- as.character(res$Suggested.Symbol)
  s <- sub(" ///.*$", "", s)
  s[is.na(s) | s == ""] <- NA_character_
  s
}

collapse_duplicate_rows <- function(expr) {
  expr <- as.matrix(expr)
  mode(expr) <- "numeric"
  if (!any(duplicated(rownames(expr)))) return(expr)
  sums <- rowsum(expr, group = rownames(expr), reorder = FALSE)
  counts <- as.vector(table(factor(rownames(expr), levels = rownames(sums))))
  sweep(sums, 1, counts, "/")
}

score_cohort <- function(expr, gene_sets, tag = "cohort", out_dir = NULL) {
  expr <- collapse_duplicate_rows(expr)
  expr_cur <- to_current(rownames(expr))
  keep <- !is.na(expr_cur) & expr_cur != ""
  expr <- expr[keep, , drop = FALSE]
  rownames(expr) <- expr_cur[keep]
  expr <- collapse_duplicate_rows(expr)
  matched <- lapply(gene_sets, function(g) {
    g_cur <- unique(to_current(g))
    g_cur <- g_cur[!is.na(g_cur)]
    intersect(g_cur, rownames(expr))
  })
  input_n <- vapply(gene_sets, function(g) length(unique(na.omit(to_current(g)))), integer(1))

  cov <- data.frame(
    gene_set = names(gene_sets),
    n_input = input_n,
    n_matched = lengths(matched),
    coverage = round(lengths(matched) / input_n, 3),
    stringsAsFactors = FALSE
  )
  cat("[", tag, "] coverage: ",
      paste(sprintf("%s=%.0f%%", cov$gene_set, 100 * cov$coverage), collapse = " "),
      "\n", sep = "")
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(cov, file.path(out_dir, paste0("STOP5_", tag, "_score_coverage.csv")), row.names = FALSE)
  }

  if (exists("ssgseaParam", where = asNamespace("GSVA"), mode = "function")) {
    es <- GSVA::gsva(GSVA::ssgseaParam(expr, matched, normalize = TRUE), verbose = FALSE)
  } else {
    es <- GSVA::gsva(expr, matched, method = "ssgsea", verbose = FALSE)
  }
  t(as.matrix(es))
}

km_curve_data <- function(fit) {
  s <- summary(fit)
  strata <- if (is.null(s$strata)) rep("all", length(s$time)) else as.character(s$strata)
  data.frame(time = s$time, surv = s$surv, strata = strata)
}

save_km_plot <- function(d, group_col, tag, out_dir, title) {
  f <- as.formula(paste("Surv(OS.time, OS.event) ~", group_col))
  fit <- survival::survfit(f, data = d)
  lr <- survival::survdiff(f, data = d)
  pval <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  km <- km_curve_data(fit)
  p <- ggplot(km, aes(time / 365, surv, color = strata)) +
    geom_step(linewidth = 0.8) +
    labs(x = "Time (years)", y = "Overall survival", color = NULL,
         title = title, subtitle = paste0("log-rank p = ", signif(pval, 3))) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 11))
  ggsave(file.path(out_dir, paste0("STOP5_KM_", tag, ".pdf")), p, width = 6, height = 4.3)
  data.frame(tag = tag, chisq = lr$chisq, df = length(lr$n) - 1, p = pval)
}

mes_landscape <- function(es, surv, tag, out_dir) {
  es <- as.matrix(es)
  if (!all(PROGS %in% colnames(es))) stop("mes_landscape expects sample x 4 score matrix.")
  surv <- surv[surv$id %in% rownames(es), , drop = FALSE]
  d <- data.frame(es[surv$id, PROGS, drop = FALSE],
                  OS.time = as.numeric(surv$OS.time),
                  OS.event = as.numeric(surv$OS.event))
  d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS.event) & d$OS.event %in% c(0, 1), ]

  d$z <- as.numeric(scale(d$MES_V))
  m <- survival::coxph(Surv(OS.time, OS.event) ~ z, data = d)
  sm <- summary(m)
  zph <- tryCatch(survival::cox.zph(m)$table["z", "p"], error = function(e) NA_real_)

  zall <- scale(as.matrix(d[, PROGS]))
  d$dominant_program <- factor(PROGS[max.col(zall, ties.method = "first")], levels = PROGS)
  dom_tab <- table(d$dominant_program, useNA = "ifany")
  write.csv(data.frame(program = names(dom_tab), n = as.integer(dom_tab)),
            file.path(out_dir, paste0("STOP5_", tag, "_argmax_counts.csv")), row.names = FALSE)
  dom_median <- aggregate(OS.time ~ dominant_program, data = d[d$OS.event == 1, , drop = FALSE], median)
  write.csv(dom_median, file.path(out_dir, paste0("STOP5_", tag, "_argmax_event_median_OS.csv")),
            row.names = FALSE)
  km_argmax <- save_km_plot(d, "dominant_program", paste0("argmax_", tag), out_dir,
                            paste0("Dominant program: ", tag))

  d$MES_V_median_group <- factor(ifelse(d$MES_V > median(d$MES_V, na.rm = TRUE),
                                        "MES_V-high", "MES_V-low"),
                                 levels = c("MES_V-low", "MES_V-high"))
  km_median <- save_km_plot(d, "MES_V_median_group", paste0("MESV_median_", tag), out_dir,
                            paste0("MES_V median split: ", tag))

  cat(sprintf("[%s] n=%d events=%d | MES_V perSD HR=%.3f (%.3f-%.3f) p=%.2e | argmax p=%.2e | medianKM p=%.2e\n",
              tag, nrow(d), sum(d$OS.event == 1),
              sm$coef["z", "exp(coef)"],
              sm$conf.int["z", "lower .95"],
              sm$conf.int["z", "upper .95"],
              sm$coef["z", "Pr(>|z|)"],
              km_argmax$p, km_median$p))

  data.frame(
    cohort = tag,
    n = nrow(d),
    events = sum(d$OS.event == 1),
    MESV_HR = sm$coef["z", "exp(coef)"],
    lo = sm$conf.int["z", "lower .95"],
    hi = sm$conf.int["z", "upper .95"],
    MESV_p = sm$coef["z", "Pr(>|z|)"],
    MESV_zph_p = zph,
    argmax_p = km_argmax$p,
    medianKM_p = km_median$p,
    dir_MESV = ifelse(sm$coef["z", "exp(coef)"] > 1, "adverse", "protective"),
    stringsAsFactors = FALSE
  )
}
