# =============================================================================
# drilldown_html_utils.R
# Shared helpers + entry points for step 08 HTML drilldown outputs.
#
# Exports two top-level functions called by 08_aggregate_results.R:
#   generate_positive_drilldown(contrast_dir, contrast_id, strata, cfg,
#                                pipeline_root, summary_dir)
#   generate_negative_drilldown(contrast_dir, contrast_id, strata, cfg,
#                                pipeline_root, summary_dir)
# =============================================================================

suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(fgsea)
  library(reactable)
  library(htmltools)
  library(htmlwidgets)
})

# =============================================================================
# Shared: load per-stratum DE results (RUV preferred, base fallback)
# =============================================================================

load_stratum_de <- function(contrast_dir, stratum) {
  finals <- file.path(contrast_dir, stratum, "finals")
  ruv_f  <- file.path(finals, "results_ruv.tsv")
  base_f <- file.path(finals, "results_base.tsv")
  ruv_usable <- file.exists(ruv_f) &&
    !identical(trimws(readLines(ruv_f, n = 1, warn = FALSE)), "skipped=TRUE")
  f <- if (ruv_usable) ruv_f else if (file.exists(base_f)) base_f else NULL
  if (is.null(f)) return(NULL)
  res <- fread(f, quote = '"')
  if (!all(c("gene", "stat") %in% names(res))) return(NULL)
  res <- res[!is.na(stat) & !is.na(gene)]
  res[, absstat := abs(stat)]
  res <- res[order(-absstat)]
  res <- res[!duplicated(gene)]
  list(
    ranks = setNames(res$stat, res$gene),
    de    = res[, .(gene, log2FoldChange, padj)]
  )
}

# =============================================================================
# Shared: NES color scale (diverging red/white/blue, clamped at |NES| = 3)
# =============================================================================

nes_color <- function(value) {
  if (is.na(value)) return("#f2f2f2")
  v <- max(min(value / 3, 1), -1)
  if (v >= 0) {
    r <- 255; g <- round(255 * (1 - v)); b <- round(255 * (1 - v))
  } else {
    r <- round(255 * (1 + v)); g <- round(255 * (1 + v)); b <- 255
  }
  sprintf("rgb(%d,%d,%d)", r, g, b)
}

nes_cell <- function(value, padj) {
  if (is.na(value)) return(div(style = "color:#ccc;font-size:11px;", "\u2014"))
  sig  <- !is.na(padj) && padj < 0.05
  star <- if (!is.na(padj) && padj < 0.001) "***"
          else if (!is.na(padj) && padj < 0.01) "**"
          else if (!is.na(padj) && padj < 0.05) "*"
          else ""
  cls <- if (sig) "nes-cell-sig" else "nes-cell-nonsig"
  div(class = cls, onclick = "event.stopPropagation()",
      style = paste0("background:", nes_color(value),
                     ";padding:4px 6px;text-align:center;",
                     "border-radius:3px;font-variant-numeric:tabular-nums;",
                     "color:#222;font-weight:400;"),
      paste0(sprintf("%.2f", value), star))
}

make_nes_cols <- function(tbl, strata) {
  out <- list()
  for (s in strata) {
    padj_col <- paste0(s, "__padj")
    local({
      s_l <- s; p_l <- padj_col
      out[[s_l]] <<- colDef(
        name = s_l, minWidth = 70, align = "center",
        cell = function(value, index) nes_cell(value, tbl[[p_l]][index])
      )
    })
  }
  out
}

# =============================================================================
# Shared: L3 gene-level reactable (leading-edge genes × strata)
# =============================================================================

build_gene_reactable <- function(le_tbl) {
  if (is.null(le_tbl) || nrow(le_tbl) == 0)
    return(div(style = "padding:10px;color:#999;", "No leading-edge genes."))

  lfc_color <- function(value) {
    v <- max(min(value / 2, 1), -1)
    if (v >= 0) sprintf("rgb(255,%d,%d)", round(255*(1-v)), round(255*(1-v)))
    else        sprintf("rgb(%d,%d,255)", round(255*(1+v)), round(255*(1+v)))
  }
  lfc_cell <- function(value, padj) {
    if (is.na(value)) return(div(style = "color:#ccc;font-size:11px;", "\u2014"))
    sig  <- !is.na(padj) && padj < 0.05
    star <- if (!is.na(padj) && padj < 0.001) "***"
            else if (!is.na(padj) && padj < 0.01) "**"
            else if (!is.na(padj) && padj < 0.05) "*"
            else ""
    cls <- if (sig) "lfc-cell-sig" else "lfc-cell-nonsig"
    div(class = cls, onclick = "event.stopPropagation()",
        style = paste0("background:", lfc_color(value), ";padding:3px 5px;",
                       "text-align:center;border-radius:3px;",
                       "font-variant-numeric:tabular-nums;color:#222;font-weight:400;"),
        paste0(sprintf("%.2f", value), star))
  }

  lfc_wide  <- dcast(le_tbl, gene ~ stratum, value.var = "log2FoldChange")
  padj_wide <- dcast(le_tbl, gene ~ stratum, value.var = "padj")
  num_strata <- setdiff(names(lfc_wide), "gene")
  setnames(padj_wide, num_strata, paste0(num_strata, "__padj"))
  g_tbl <- merge(lfc_wide, padj_wide, by = "gene")
  setDT(g_tbl)

  g_tbl[, max_abs_lfc := apply(.SD, 1, function(v) max(abs(v), na.rm = TRUE)),
        .SDcols = num_strata]
  padj_cols_g <- paste0(num_strata, "__padj")
  g_tbl[, sig_any := apply(.SD, 1, function(x) any(!is.na(x) & x < 0.05)),
        .SDcols = padj_cols_g]
  setorder(g_tbl, -sig_any, -max_abs_lfc)
  n_sig  <- sum(g_tbl$sig_any)
  n_keep <- min(nrow(g_tbl), max(n_sig, 10L))
  g_tbl[, row_class := "nonsig-l3"]
  g_tbl[seq_len(n_keep), row_class := ""]
  g_tbl[, c("max_abs_lfc", "sig_any") := NULL]

  cols <- list(
    gene      = colDef(name = "gene", width = 110,
                       style = list(fontFamily = "monospace", fontWeight = "500")),
    row_class = colDef(show = FALSE)
  )
  for (s in num_strata) {
    padj_col <- paste0(s, "__padj")
    local({
      s_l <- s; p_l <- padj_col
      cols[[s_l]] <<- colDef(
        name = s_l, minWidth = 80, align = "center",
        cell = function(value, index) lfc_cell(value, g_tbl[[p_l]][index])
      )
    })
  }
  for (s in num_strata) cols[[paste0(s, "__padj")]] <- colDef(show = FALSE)

  reactable(g_tbl, columns = cols,
            rowClass = function(index) g_tbl$row_class[index],
            defaultPageSize = max(20L, nrow(g_tbl)),
            searchable = TRUE, compact = TRUE, bordered = TRUE,
            style = list(fontSize = "12px"), fullWidth = TRUE)
}

# =============================================================================
# Shared: build leading-edge lookup table for a set of gene sets + strata
# =============================================================================

build_le_lookup <- function(set_ids, fg_dt, de_by_stratum) {
  setNames(lapply(set_ids, function(pid) {
    sub <- fg_dt[pathway == pid & !is.na(leadingEdge) & leadingEdge != ""]
    if (nrow(sub) == 0) return(data.table())
    rows <- lapply(seq_len(nrow(sub)), function(i) {
      genes <- unlist(strsplit(sub$leadingEdge[i], "\\|"))
      s     <- sub$stratum[i]
      de    <- de_by_stratum[[s]]
      if (is.null(de)) return(NULL)
      le <- de[gene %in% genes]
      if (nrow(le) == 0) return(NULL)
      le[, stratum := s]
      le[, .(gene, stratum, log2FoldChange, padj)]
    })
    rbindlist(rows, fill = TRUE)
  }), set_ids)
}

# =============================================================================
# Shared: HTML page shell with toggle bar
# =============================================================================

.page_css <- "
  html, body { margin: 0; padding: 0; }
  body { font-family: -apple-system, system-ui, sans-serif;
         padding: 20px; box-sizing: border-box; }
  h1 { font-size: 20px; margin-bottom: 4px; }
  .rt-table { font-variant-numeric: tabular-nums; }
  .rt-th, .rt-td { padding: 6px 8px !important; }
  .rt-th { white-space: nowrap !important; }
  .reactable { width: 100% !important; }
  body.hide-nonsig-l1 .rt-tr.nonsig-l1 { display: none; }
  body.hide-nonsig-l2 .rt-tr.nonsig-l2 { display: none; }
  body.hide-nonsig-l3 .rt-tr.nonsig-l3 { display: none; }
  body.highlight-sig .nes-cell-sig {
    color: white !important; font-weight: 700 !important;
    box-shadow: inset 0 0 0 1.5px rgba(0,0,0,0.22) !important; }
  body.highlight-sig .nes-cell-nonsig {
    filter: saturate(0.12) !important; color: #bbb !important; }
  body.highlight-sig .lfc-cell-sig {
    color: white !important; font-weight: 700 !important;
    box-shadow: inset 0 0 0 1.5px rgba(0,0,0,0.22) !important; }
  body.highlight-sig .lfc-cell-nonsig {
    filter: saturate(0.12) !important; color: #bbb !important; }
  .sig-toggles { display: flex; gap: 22px; align-items: center;
    margin: 10px 0 14px; padding: 8px 12px; background: #f0f4f8;
    border-radius: 6px; font-size: 13px; flex-wrap: wrap; }
  .sig-toggles label { display: flex; align-items: center; gap: 6px;
    cursor: pointer; user-select: none; }
  .sig-toggles input[type=checkbox] { width: 15px; height: 15px; cursor: pointer; }
"

assemble_page <- function(contrast_id, l1_tbl, legend, how_to, title_str,
                          extra_css = "", default_highlight = FALSE) {
  hl_script <- if (default_highlight)
    tags$script(HTML("document.body.classList.add('highlight-sig');")) else NULL
  hl_checked <- if (default_highlight)
    HTML('<input type="checkbox" checked
           onclick="document.body.classList.toggle(\'highlight-sig\')"
           style="width:15px;height:15px;cursor:pointer;">')
  else
    tags$input(type = "checkbox",
               onclick = "document.body.classList.toggle('highlight-sig')")

  tags$html(
    tags$head(
      tags$title(title_str),
      tags$style(HTML(paste0(.page_css, extra_css)))
    ),
    tags$body(
      hl_script,
      h1(title_str),
      legend, how_to,
      tags$div(class = "sig-toggles",
        tags$b(style = "color:#333;", "Filter:"),
        tags$label(
          tags$input(type = "checkbox",
                     onclick = "document.body.classList.toggle('hide-nonsig-l1')"),
          "L1: only significant in \u22651 cell type"),
        tags$label(
          tags$input(type = "checkbox",
                     onclick = "document.body.classList.toggle('hide-nonsig-l2')"),
          "L2: only significant in \u22651 cell type"),
        tags$label(
          tags$input(type = "checkbox",
                     onclick = "document.body.classList.toggle('hide-nonsig-l3')"),
          "Genes: only significant (padj<0.05)"),
        tags$span(style = "border-left:1px solid #ccc;margin:0 4px;height:16px;display:inline-block;"),
        tags$b(style = "color:#333;", "Display:"),
        tags$label(hl_checked, "Highlight significant (padj<0.05)")
      ),
      l1_tbl
    )
  )
}

# =============================================================================
# POSITIVE: pathway drilldown (mega-sets → member pathways → leading-edge genes)
# =============================================================================

source_of <- function(x) {
  ifelse(grepl("^REACTOME_", x), "Reactome",
  ifelse(grepl("^HALLMARK_", x), "Hallmark", "KEGG"))
}
shorten <- function(x) {
  s <- gsub("^(REACTOME|KEGG_MEDICUS_REFERENCE|KEGG_MEDICUS|KEGG|HALLMARK)_", "", x)
  gsub("_", " ", s)
}

generate_positive_drilldown <- function(contrast_dir, contrast_id, strata,
                                        cfg, pipeline_root, summary_dir) {
  pd <- cfg$pathway_drilldown %||% list()
  mega_f   <- file.path(pipeline_root, pd$mega_sets   %||% "resources/pathway_drilldown/mega_sets.tsv")
  hall_f   <- file.path(pipeline_root, pd$hallmark_gmt %||% "resources/gsea_files/hallmark.gmt")
  rk_f     <- file.path(pipeline_root, cfg$fgsea$gene_sets_file %||% "resources/gsea_files/reactome_kegg.gmt.txt")
  fgsea_f  <- file.path(summary_dir, sprintf("%s_fgsea_final.tsv", contrast_id))
  out_html <- file.path(summary_dir, "positive_benchmarking.html")

  for (f in c(mega_f, hall_f, rk_f, fgsea_f)) {
    if (!file.exists(f)) {
      cat(sprintf("  [WARN] positive drilldown skipped: missing %s\n", f))
      return(invisible(NULL))
    }
  }

  cat("  Generating positive pathway drilldown HTML...\n")

  suppressWarnings({
    rk_sets <- gmtPathways(rk_f)
    h_sets  <- gmtPathways(hall_f)
  })
  all_sets <- c(rk_sets, h_sets)

  # Load per-stratum DE
  cache <- lapply(strata, function(s) load_stratum_de(contrast_dir, s))
  names(cache) <- strata
  cache <- cache[!sapply(cache, is.null)]
  strata_ok <- names(cache)
  de_by_stratum <- lapply(cache, `[[`, "de")

  # Pipeline fGSEA (Reactome/KEGG) from _fgsea_final
  fg <- fread(fgsea_f)
  fg <- fg[stratum %in% strata_ok]

  # Hallmark fGSEA (recompute — not in pipeline GMT)
  hall_fg <- rbindlist(lapply(strata_ok, function(s) {
    set.seed(1)
    r <- fgsea(pathways = h_sets, stats = cache[[s]]$ranks,
               minSize = 10, maxSize = 500)
    r$stratum <- s
    r
  }), fill = TRUE)
  hall_fg[, leadingEdge := sapply(leadingEdge, paste, collapse = "|")]

  individual_fg <- rbind(
    fg[, .(pathway, stratum, NES, padj, size, leadingEdge)],
    hall_fg[, .(pathway, stratum, NES, padj, size, leadingEdge)]
  )

  # Mega-set fGSEA
  mega <- fread(mega_f)
  mega_genes <- lapply(mega$member_pathways, function(s) {
    hits <- intersect(unlist(strsplit(s, "\\|")), names(all_sets))
    unique(unlist(all_sets[hits]))
  })
  names(mega_genes) <- mega$mega_set_id
  mega$num_genes  <- sapply(mega_genes, length)
  mega$n_resolved <- sapply(mega$member_pathways, function(s)
    length(intersect(unlist(strsplit(s, "\\|")), names(all_sets))))
  mega$sole_member <- sapply(mega$member_pathways, function(s) {
    hits <- intersect(unlist(strsplit(s, "\\|")), names(all_sets))
    if (length(hits) == 1) hits else NA_character_
  })

  multi_ids  <- mega$mega_set_id[mega$n_resolved >  1]
  single_ids <- mega$mega_set_id[mega$n_resolved == 1]

  mega_fg_multi <- rbindlist(lapply(strata_ok, function(s) {
    set.seed(1)
    r <- fgsea(pathways = mega_genes[multi_ids], stats = cache[[s]]$ranks,
               minSize = 5, maxSize = 5000)
    r$stratum <- s; r
  }), fill = TRUE)
  mega_fg_multi[, leadingEdge := sapply(leadingEdge, paste, collapse = "|")]

  single_lookup <- mega[mega_set_id %in% single_ids, .(mega_set_id, sole_member)]
  mega_fg_single <- merge(single_lookup,
    individual_fg[, .(pathway, stratum, NES, padj, size, leadingEdge)],
    by.x = "sole_member", by.y = "pathway")
  mega_fg_single[, pathway := mega_set_id]
  mega_fg_single[, c("mega_set_id", "sole_member") := NULL]
  mega_fg <- rbind(mega_fg_multi, mega_fg_single, fill = TRUE)

  # L1 table
  mega_nes  <- dcast(mega_fg, pathway ~ stratum, value.var = "NES")
  mega_padj <- dcast(mega_fg, pathway ~ stratum, value.var = "padj")
  setnames(mega_padj, setdiff(names(mega_padj), "pathway"),
           paste0(setdiff(names(mega_padj), "pathway"), "__padj"))
  mega_table <- merge(mega[, .(mega_set_id, category, subcategory, description, num_genes)],
                      mega_nes,  by.x = "mega_set_id", by.y = "pathway", all.x = TRUE)
  mega_table <- merge(mega_table, mega_padj,
                      by.x = "mega_set_id", by.y = "pathway", all.x = TRUE)
  setorder(mega_table, category, subcategory, mega_set_id)

  strata_cols <- intersect(strata_ok, names(mega_table))
  padj_cols_m <- paste0(strata_cols, "__padj")
  padj_cols_m <- padj_cols_m[padj_cols_m %in% names(mega_table)]
  mega_table[, sig_any := apply(.SD, 1, function(x) any(!is.na(x) & x < 0.05)),
             .SDcols = padj_cols_m]

  # L2 member table
  member_long <- mega[, .(member_pathway = unlist(strsplit(member_pathways, "\\|"))),
                      by = .(mega_set_id)]
  member_long <- member_long[member_pathway %in% names(all_sets)]
  member_long[, source     := source_of(member_pathway)]
  member_long[, short_name := shorten(member_pathway)]
  member_long[, n_genes    := sapply(member_pathway, function(p) length(all_sets[[p]]))]

  member_nes  <- dcast(individual_fg, pathway ~ stratum, value.var = "NES")
  member_padj <- dcast(individual_fg, pathway ~ stratum, value.var = "padj")
  setnames(member_padj, setdiff(names(member_padj), "pathway"),
           paste0(setdiff(names(member_padj), "pathway"), "__padj"))
  member_table <- merge(member_long, member_nes,
                        by.x = "member_pathway", by.y = "pathway", all.x = TRUE)
  member_table <- merge(member_table, member_padj,
                        by.x = "member_pathway", by.y = "pathway", all.x = TRUE)

  # LE lookup for member pathways
  le_lookup <- build_le_lookup(unique(member_long$member_pathway),
                               individual_fg, de_by_stratum)

  # L2 builder
  build_member_reactable <- function(mega_id) {
    sub <- member_table[mega_set_id == mega_id]
    if (nrow(sub) == 0)
      return(div(style = "padding:10px;color:#999;", "No member pathways resolved."))
    sc <- intersect(strata_ok, names(sub))
    pc <- paste0(sc, "__padj"); pc <- pc[pc %in% names(sub)]
    sub[, sig_any := apply(.SD, 1, function(x) any(!is.na(x) & x < 0.05)), .SDcols = pc]
    cols <- list(
      mega_set_id    = colDef(show = FALSE),
      member_pathway = colDef(show = FALSE),
      sig_any        = colDef(show = FALSE),
      source = colDef(name = "src", width = 80,
                      style = function(value) {
                        col <- switch(value, Reactome="#2ca02c", KEGG="#1f77b4",
                                      Hallmark="#ff7f0e", "#666")
                        list(color = col, fontWeight = "600")
                      }),
      short_name = colDef(name = "pathway", minWidth = 280,
                          style = list(fontSize = "11px")),
      n_genes = colDef(name = "n_genes", width = 80, align = "right")
    )
    cols <- c(cols, make_nes_cols(sub, sc))
    for (s in sc) cols[[paste0(s, "__padj")]] <- colDef(show = FALSE)
    reactable(sub, columns = cols,
              rowClass = function(index) if (!sub$sig_any[index]) "nonsig-l2" else NULL,
              details = function(index) {
                pw <- sub$member_pathway[index]
                div(style = "padding:10px;background:#fafafa;",
                    tags$div(style = "font-weight:600;margin-bottom:6px;",
                             sprintf("Leading-edge genes for %s", pw)),
                    build_gene_reactable(le_lookup[[pw]]))
              },
              defaultPageSize = 20, compact = TRUE, bordered = TRUE,
              style = list(fontSize = "12px"), fullWidth = TRUE)
  }

  # L1 reactable
  l1_cols <- list(
    category    = colDef(name = "category", width = 130, style = list(fontWeight = "600")),
    subcategory = colDef(name = "subcategory", width = 150),
    mega_set_id = colDef(name = "mega-set", width = 180, style = list(fontFamily = "monospace")),
    description = colDef(name = "description", minWidth = 240,
                         style = list(fontSize = "11px", color = "#444")),
    num_genes   = colDef(name = "n_genes", width = 80, align = "right"),
    sig_any     = colDef(show = FALSE)
  )
  l1_cols <- c(l1_cols, make_nes_cols(mega_table, strata_cols))
  for (s in strata_cols) l1_cols[[paste0(s, "__padj")]] <- colDef(show = FALSE)

  l1_tbl <- reactable(mega_table, columns = l1_cols,
    groupBy = "category",
    rowClass = function(index) if (!mega_table$sig_any[index]) "nonsig-l1" else NULL,
    details = function(index) {
      mid <- mega_table$mega_set_id[index]
      div(style = "padding:10px;background:#f4f7fb;",
          tags$div(style = "font-weight:600;margin-bottom:6px;",
                   sprintf("Member pathways of %s", mid)),
          build_member_reactable(mid))
    },
    defaultPageSize = 50, searchable = TRUE, compact = TRUE,
    bordered = TRUE, striped = TRUE,
    style = list(fontSize = "13px"), fullWidth = TRUE)

  legend <- div(style = "margin:10px 0;font-size:12px;color:#555;",
    "NES color scale: ",
    span(style = "background:#2166AC;color:white;padding:2px 6px;", "neg"),
    " ", span(style = "background:white;border:1px solid #ccc;padding:2px 6px;", "0"), " ",
    span(style = "background:#B2182B;color:white;padding:2px 6px;", "pos"),
    " \u00b7 * padj<0.05, ** <0.01, *** <0.001 \u00b7 ",
    tags$span(style = "color:#2ca02c;font-weight:600;", "Reactome"), " \u00b7 ",
    tags$span(style = "color:#1f77b4;font-weight:600;", "KEGG"), " \u00b7 ",
    tags$span(style = "color:#ff7f0e;font-weight:600;", "Hallmark"))

  how_to <- div(
    style = "margin:10px 0;padding:10px;background:#f8f9fa;border-left:3px solid #2166AC;font-size:12px;",
    tags$b("How to read this:"), br(),
    "Click \u25b6 to expand a mega-set into its member pathways, then expand a pathway to see ",
    "leading-edge genes with per-stratum log2FC and padj from DESeq2 (RUVseq where available, base otherwise).")

  page <- assemble_page(contrast_id, l1_tbl, legend, how_to,
                        sprintf("%s - pathway drilldown", contrast_id),
                        default_highlight = FALSE)
  save_html(page, out_html)
  cat(sprintf("  Saved: %s (%.1f MB)\n", out_html, file.info(out_html)$size / 1024^2))
  invisible(out_html)
}

# =============================================================================
# NEGATIVE: negative-controls drilldown (LLM signatures → leading-edge genes)
# =============================================================================

cat_color <- function(x) ifelse(x == "artifact", "#D94040", "#7B3294")

generate_negative_drilldown <- function(contrast_dir, contrast_id, strata,
                                        cfg, pipeline_root, summary_dir) {
  pd     <- cfg$pathway_drilldown %||% list()
  sigs_f <- file.path(pipeline_root, pd$negative_sigs %||%
                      "resources/benchmarking/master_gene_signatures.tsv")
  out_html <- file.path(summary_dir, "negative_benchmarking.html")

  if (!file.exists(sigs_f)) {
    cat(sprintf("  [WARN] negative drilldown skipped: missing %s\n", sigs_f))
    return(invisible(NULL))
  }

  cat("  Generating negative controls drilldown HTML...\n")

  # Load per-stratum DE
  cache <- lapply(strata, function(s) load_stratum_de(contrast_dir, s))
  names(cache) <- strata
  cache <- cache[!sapply(cache, is.null)]
  strata_ok <- names(cache)
  de_by_stratum <- lapply(cache, `[[`, "de")

  # Load negative signatures
  sigs_all <- fread(sigs_f)
  sigs <- sigs_all[control_type == "negative"]

  sig_sets <- setNames(
    lapply(sigs$genes, function(g) unique(unlist(strsplit(g, "\\|")))),
    sigs$signature_id)

  class_sets <- lapply(split(sigs$signature_id, sigs$class), function(ids)
    unique(unlist(sig_sets[ids])))

  # fGSEA: class level + signature level
  class_fg <- rbindlist(lapply(strata_ok, function(s) {
    set.seed(1)
    r <- fgsea(pathways = class_sets, stats = cache[[s]]$ranks,
               minSize = 3, maxSize = 5000, eps = 0)
    r$stratum <- s; r
  }), fill = TRUE)
  class_fg[, leadingEdge := sapply(leadingEdge, paste, collapse = "|")]

  sig_fg <- rbindlist(lapply(strata_ok, function(s) {
    set.seed(1)
    r <- fgsea(pathways = sig_sets, stats = cache[[s]]$ranks,
               minSize = 3, maxSize = 5000, eps = 0)
    r$stratum <- s; r
  }), fill = TRUE)
  sig_fg[, leadingEdge := sapply(leadingEdge, paste, collapse = "|")]

  # LE lookup
  le_lookup <- build_le_lookup(sigs$signature_id, sig_fg, de_by_stratum)

  # L1 table
  class_meta <- unique(sigs[, .(category, class)])
  class_meta[, n_sigs := sapply(class, function(cl) sum(sigs$class == cl))]
  class_nes  <- dcast(class_fg, pathway ~ stratum, value.var = "NES")
  class_padj <- dcast(class_fg, pathway ~ stratum, value.var = "padj")
  setnames(class_padj, setdiff(names(class_padj), "pathway"),
           paste0(setdiff(names(class_padj), "pathway"), "__padj"))
  class_table <- merge(class_meta, class_nes,
                       by.x = "class", by.y = "pathway", all.x = TRUE)
  class_table <- merge(class_table, class_padj,
                       by.x = "class", by.y = "pathway", all.x = TRUE)
  setorder(class_table, category, class)
  strata_cols <- intersect(strata_ok, names(class_table))
  padj_cols_c <- paste0(strata_cols, "__padj")
  padj_cols_c <- padj_cols_c[padj_cols_c %in% names(class_table)]
  class_table[, sig_any := apply(.SD, 1, function(x) any(!is.na(x) & x < 0.05)),
              .SDcols = padj_cols_c]

  # L2 table
  sig_meta <- sigs[, .(signature_id, category, class, subcategory, description,
                        n_genes = num_genes)]
  sig_nes  <- dcast(sig_fg, pathway ~ stratum, value.var = "NES")
  sig_padj <- dcast(sig_fg, pathway ~ stratum, value.var = "padj")
  setnames(sig_padj, setdiff(names(sig_padj), "pathway"),
           paste0(setdiff(names(sig_padj), "pathway"), "__padj"))
  sig_table <- merge(sig_meta, sig_nes,
                     by.x = "signature_id", by.y = "pathway", all.x = TRUE)
  sig_table <- merge(sig_table, sig_padj,
                     by.x = "signature_id", by.y = "pathway", all.x = TRUE)

  # L2 builder
  build_sig_reactable <- function(cls) {
    sub <- sig_table[class == cls]
    if (nrow(sub) == 0) return(div(style = "padding:10px;color:#999;", "No signatures."))
    sc <- intersect(strata_ok, names(sub))
    pc <- paste0(sc, "__padj"); pc <- pc[pc %in% names(sub)]
    sub[, sig_any := apply(.SD, 1, function(x) any(!is.na(x) & x < 0.05)), .SDcols = pc]
    cols <- list(
      signature_id = colDef(show = FALSE), class = colDef(show = FALSE),
      sig_any = colDef(show = FALSE),
      category    = colDef(name = "type", width = 130,
                           style = function(value) list(color = cat_color(value),
                                                        fontWeight = "600")),
      subcategory = colDef(name = "subcategory", width = 150,
                           style = list(fontSize = "11px")),
      description = colDef(name = "description", minWidth = 280,
                           style = list(fontSize = "11px", color = "#444")),
      n_genes     = colDef(name = "n_genes", width = 70, align = "right")
    )
    cols <- c(cols, make_nes_cols(sub, sc))
    for (s in sc) cols[[paste0(s, "__padj")]] <- colDef(show = FALSE)
    reactable(sub, columns = cols,
              rowClass = function(index) if (!sub$sig_any[index]) "nonsig-l2" else NULL,
              details = function(index) {
                sid <- sub$signature_id[index]
                div(style = "padding:10px;background:#fafafa;",
                    tags$div(style = "font-weight:600;margin-bottom:6px;",
                             sprintf("Leading-edge genes: %s", sid)),
                    build_gene_reactable(le_lookup[[sid]]))
              },
              defaultPageSize = 20, compact = TRUE, bordered = TRUE,
              style = list(fontSize = "12px"), fullWidth = TRUE)
  }

  # L1 reactable
  l1_cols <- list(
    category = colDef(name = "category", width = 160,
                      style = function(value) list(color = cat_color(value),
                                                   fontWeight = "600")),
    class    = colDef(name = "class", minWidth = 180,
                      style = list(fontFamily = "monospace")),
    n_sigs   = colDef(name = "n_sigs", width = 70, align = "right"),
    sig_any  = colDef(show = FALSE)
  )
  l1_cols <- c(l1_cols, make_nes_cols(class_table, strata_cols))
  for (s in strata_cols) l1_cols[[paste0(s, "__padj")]] <- colDef(show = FALSE)

  l1_tbl <- reactable(class_table, columns = l1_cols,
    groupBy = "category",
    rowClass = function(index) if (!class_table$sig_any[index]) "nonsig-l1" else NULL,
    details = function(index) {
      cls <- class_table$class[index]
      div(style = "padding:10px;background:#f4f7fb;",
          tags$div(style = "font-weight:600;margin-bottom:6px;",
                   sprintf("Individual signatures: %s", cls)),
          build_sig_reactable(cls))
    },
    defaultPageSize = 50, searchable = TRUE, compact = TRUE,
    bordered = TRUE, striped = TRUE,
    style = list(fontSize = "13px"), fullWidth = TRUE)

  legend <- div(style = "margin:10px 0;font-size:12px;color:#555;",
    "NES color scale: ",
    span(style = "background:#2166AC;color:white;padding:2px 6px;", "neg"),
    " ", span(style = "background:white;border:1px solid #ccc;padding:2px 6px;", "0"), " ",
    span(style = "background:#B2182B;color:white;padding:2px 6px;", "pos"),
    " \u00b7 * padj<0.05, ** <0.01, *** <0.001 \u00b7 ",
    tags$span(style = "color:#D94040;font-weight:600;", "artifact"), " \u00b7 ",
    tags$span(style = "color:#7B3294;font-weight:600;", "subject_confounder"))

  how_to <- div(
    style = "margin:10px 0;padding:10px;background:#f8f9fa;border-left:3px solid #888;font-size:12px;",
    tags$b("How to read this:"), br(),
    "LLM-curated negative-control signatures. These should be ",
    tags$b("flat (NES \u2248 0, non-significant)"),
    " in a clean comparison. Colored/starred cells are red flags: ",
    "dissociation/ambient/ischemia = handling artifact; sex/age/BMI = donor confounder; ",
    "medication/metabolic state = clinical co-variable not in the model.")

  page <- assemble_page(contrast_id, l1_tbl, legend, how_to,
                        sprintf("%s - negative controls (should be flat)", contrast_id),
                        extra_css = "body { background: #fafbfc; } .nes-cell-nonsig { opacity: 0.55; }",
                        default_highlight = TRUE)
  save_html(page, out_html)
  cat(sprintf("  Saved: %s (%.1f MB)\n", out_html, file.info(out_html)$size / 1024^2))
  invisible(out_html)
}

# Null-coalescing operator (in case not loaded via io_utils)
`%||%` <- function(a, b) if (!is.null(a)) a else b
