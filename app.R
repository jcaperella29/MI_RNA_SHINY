# ====================== app.R ======================
suppressPackageStartupMessages({
  library(shiny); library(DT); library(plotly)
  library(readr); library(dplyr); library(ggplot2)
  library(DESeq2); library(umap); library(pwr); library(randomForest); library(pROC)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichR)
  library(multiMiR)       # targets (validated + predicted) for hs/mm
  library(gprofiler2)     # REAC/WP & fallback
  library(org.Mm.eg.db)   # mouse OrgDb
  library(AnnotationDbi)  # ID mapping
  library(digest)         # cache key
  library(R.utils)        # timeouts
  library(org.Dr.eg.db)   # zebrafish
  library(org.Dm.eg.db)   # fly
})
normalize_mirna_ids <- function(mirnas, species = c("hs","mm")) {
  species <- match.arg(species)
  prefix  <- if (species == "hs") "hsa-" else "mmu-"
  
  x <- as.character(mirnas)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  
  # Normalize separators/case
  x <- gsub("_", "-", x)
  x <- gsub("\\s+", "", x)
  x <- gsub("^mir-", "miR-", x, ignore.case = TRUE)
  x <- gsub("^mir",  "miR",  x, ignore.case = TRUE)
  
  # If IDs have no species prefix, add it
  has_prefix <- grepl("^(hsa-|mmu-|dre-|dme-)", x, ignore.case = TRUE)
  x[!has_prefix] <- paste0(prefix, x[!has_prefix])
  
  # Force canonical prefix casing
  x <- sub("^hsa-", "hsa-", x, ignore.case = TRUE)
  x <- sub("^mmu-", "mmu-", x, ignore.case = TRUE)
  
  unique(x)
}




options(timeout = 3000)  # Simple null-coalescing helper for safety
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Static CSV fallback ----
mirna_static_targets <- tryCatch(
  read.csv("data/mirna_targets.csv", stringsAsFactors = FALSE),
  error = function(e) data.frame(mirna = character(), target = character())
)



# ---- Species helpers ----
detect_species <- function(mirnas, fallback = "hs") {
  if (length(mirnas) == 0 || is.na(mirnas[1])) return(fallback)
  m <- tolower(mirnas[1])
  if (grepl("^hsa-", m)) return("hs")
  if (grepl("^mmu-", m)) return("mm")
  if (grepl("^dre-|^danrer|^dr-", m)) return("dr")
  if (grepl("^dme-|^dmel|^dm-", m)) return("dm")
  fallback
}

species_labels <- c(
  "Auto (infer from miRNA IDs)" = "auto",
  "Human (hsa-miR-…)"           = "hs",
  "Mouse (mmu-miR-…)"           = "mm",
  "Zebrafish (dre-miR-…)"       = "dr",
  "Fly (dme-miR-…)"             = "dm"
)

MM_ORG <- c(hs = "hsa", mm = "mmu")  # multiMiR org codes

# ===================== Enrichment helpers (TOP-LEVEL) =====================

# -- Safe AnnotationDbi::select wrapper (de-dup + valid keys only)
safe_select <- function(orgdb, keys, keytype, columns) {
  keys <- unique(as.character(keys))
  keys <- keys[!is.na(keys) & keys != ""]
  if (!length(keys)) return(data.frame())
  valid <- AnnotationDbi::keys(orgdb, keytype = keytype)
  keys  <- intersect(keys, valid)
  if (!length(keys)) return(data.frame())
  out <- AnnotationDbi::select(orgdb, keys = keys, keytype = keytype, columns = columns)
  dplyr::distinct(as.data.frame(out))
}

# ---- Attach per-term gene & miRNA lists ----
add_gene_mirna_cols <- function(df, genes_col, id_type = c("symbol", "entrez"),
                                target_df, symbol_col = "SYMBOL", entrez_col = "ENTREZID",
                                mirna_col = "miRNA") {
  id_type <- match.arg(id_type)
  if (is.null(df) || !nrow(df) || is.null(target_df) || !nrow(target_df)) return(df)
  if (!genes_col %in% colnames(df)) return(df)
  if (!mirna_col %in% colnames(target_df)) return(df)

  # Normalise column names in target_df
  td <- as.data.frame(target_df)
  if (!symbol_col %in% names(td))  td[[symbol_col]]  <- NA_character_
  if (!entrez_col %in% names(td))  td[[entrez_col]]  <- NA_character_
  names(td)[names(td) == mirna_col]  <- "miRNA"
  names(td)[names(td) == symbol_col] <- "SYMBOL"
  names(td)[names(td) == entrez_col] <- "ENTREZID"

  df$Genes  <- NA_character_
  df$miRNAs <- NA_character_

  for (i in seq_len(nrow(df))) {
    raw <- as.character(df[[genes_col]][i])
    if (is.na(raw) || raw == "") next

    # Enrichr uses ";", g:Profiler "/", clusterProfiler "/"
    g_ids <- unlist(strsplit(raw, "[,;/]"))
    g_ids <- g_ids[g_ids != ""]

    if (!length(g_ids)) next

    if (id_type == "symbol") {
      hit <- td[td$SYMBOL %in% g_ids | td$target %in% g_ids, , drop = FALSE]
      genes <- unique(c(hit$SYMBOL, hit$target))
    } else {
      hit <- td[td$ENTREZID %in% g_ids, , drop = FALSE]
      genes <- unique(ifelse(is.na(hit$SYMBOL) | hit$SYMBOL == "",
                             hit$ENTREZID, hit$SYMBOL))
    }

    df$Genes[i]  <- if (length(genes))  paste(genes, collapse = "; ")  else NA_character_
    df$miRNAs[i] <- if (nrow(hit))      paste(unique(hit$miRNA), collapse = "; ") else NA_character_
  }
  df
}

# -- Enrichr/g:Profiler DB menus by species
enrichr_human_choices <- c(
  "GO_Biological_Process_2021",
  "GO_Molecular_Function_2021",
  "GO_Cellular_Component_2021",
  "KEGG_2021_Human",
  "Reactome_2022",
  "WikiPathway_2021_Human"
)

gprof_sources <- c("GO:BP","KEGG","REAC","WP")

db_choices_for_species <- function(sp) {
  if (sp == "hs") enrichr_human_choices else gprof_sources
}

# -- Enrichr pull with guards
.enrichr_pull <- function(genes, db) {
  ee <- try(enrichR::enrichr(genes, databases = db), silent = TRUE)
  if (inherits(ee, "try-error") || is.null(ee)) return(data.frame())
  out <- ee[[db]]
  if (is.null(out) || !is.data.frame(out) || !nrow(out)) return(data.frame())
  out
}

# ---------------------- Mouse: targets + enrich (genes + miRNAs) ----------------------
mm_targets_df <- function(mirnas) {
  # mirnas = vector of mmu-miR-* IDs
  safe_mm <- function(org, mirs, table) {
    try(R.utils::withTimeout(
      multiMiR::get_multimir(org = org, mirna = mirs, table = table, summary = FALSE),
      timeout = 2000, onTimeout = "error"
    ), silent = TRUE)
  }

  v <- safe_mm("mmu", mirnas, "validated")
  p <- safe_mm("mmu", mirnas, "predicted")

  dat <- dplyr::bind_rows(
    if (!inherits(v,"try-error") && !is.null(v) && nrow(v@data)) tibble::as_tibble(v@data) else NULL,
    if (!inherits(p,"try-error") && !is.null(p) && nrow(p@data)) tibble::as_tibble(p@data) else NULL
  )
  if (!nrow(dat)) {
    return(tibble::tibble(miRNA   = character(),
                          SYMBOL   = character(),
                          ENTREZID = character()))
  }

  df <- dat %>%
    dplyr::transmute(
      miRNA   = as.character(mature_mirna_id),
      SYMBOL  = as.character(target_symbol),
      ENTREZID = as.character(target_entrez)
    )

  # Fill missing SYMBOL from ENTREZID
  need_sym <- which((is.na(df$SYMBOL) | df$SYMBOL == "") &
                      (!is.na(df$ENTREZID) & df$ENTREZID != ""))
  if (length(need_sym)) {
    map <- safe_select(org.Mm.eg.db,
                       unique(df$ENTREZID[need_sym]),
                       keytype = "ENTREZID",
                       columns = "SYMBOL")
    if (nrow(map)) {
      idx <- match(df$ENTREZID, map$ENTREZID)
      df$SYMBOL[need_sym] <- ifelse(
        df$SYMBOL[need_sym] == "" | is.na(df$SYMBOL[need_sym]),
        map$SYMBOL[idx][need_sym],
        df$SYMBOL[need_sym]
      )
    }
  }

  # Fill missing ENTREZID from SYMBOL
  need_ent <- which((is.na(df$ENTREZID) | df$ENTREZID == "") &
                      (!is.na(df$SYMBOL) & df$SYMBOL != ""))
  if (length(need_ent)) {
    map <- safe_select(org.Mm.eg.db,
                       unique(df$SYMBOL[need_ent]),
                       keytype = "SYMBOL",
                       columns = "ENTREZID")
    if (nrow(map)) {
      idx <- match(df$SYMBOL, map$SYMBOL)
      df$ENTREZID[need_ent] <- ifelse(
        df$ENTREZID[need_ent] == "" | is.na(df$ENTREZID[need_ent]),
        as.character(map$ENTREZID[idx][need_ent]),
        df$ENTREZID[need_ent]
      )
    }
  }

  df %>%
    dplyr::filter(!is.na(ENTREZID), ENTREZID != "") %>%
    dplyr::distinct(miRNA, SYMBOL, ENTREZID)
}

mm_enrich_from_targets <- function(target_df, db = "GO:BP", q = 0.05) {
  # target_df: miRNA, SYMBOL, ENTREZID
  if (!is.data.frame(target_df) || !nrow(target_df)) return(data.frame())
  genes <- unique(target_df$ENTREZID)
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) < 10) return(data.frame())

  if (db == "GO:BP") {
    eg <- clusterProfiler::enrichGO(
      gene          = genes,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff  = q,
      readable      = FALSE
    )
    df <- as.data.frame(eg)
    if (!nrow(df)) return(data.frame())
    df$Combined.Score <- -log10(df$p.adjust)
    df2 <- df %>%
      dplyr::transmute(
        Term             = Description,
        Adjusted.P.value = p.adjust,
        Combined.Score   = Combined.Score,
        geneID           = geneID
      )
    df2 <- add_gene_mirna_cols(
      df2,
      genes_col = "geneID",
      id_type   = "entrez",
      target_df = target_df,
      symbol_col = "SYMBOL",
      entrez_col = "ENTREZID",
      mirna_col  = "miRNA"
    )
    return(df2[order(df2$Adjusted.P.value), ])
  }

  if (db == "KEGG") {
    kk <- clusterProfiler::enrichKEGG(
      gene          = genes,
      organism      = "mmu",
      pAdjustMethod = "BH",
      qvalueCutoff  = q
    )
    df <- as.data.frame(kk)
    if (!nrow(df)) return(data.frame())
    df$Combined.Score <- -log10(df$p.adjust)
    df2 <- df %>%
      dplyr::transmute(
        Term             = Description,
        Adjusted.P.value = p.adjust,
        Combined.Score   = Combined.Score,
        geneID           = geneID
      )
    df2 <- add_gene_mirna_cols(
      df2,
      genes_col = "geneID",
      id_type   = "entrez",
      target_df = target_df,
      symbol_col = "SYMBOL",
      entrez_col = "ENTREZID",
      mirna_col  = "miRNA"
    )
    return(df2[order(df2$Adjusted.P.value), ])
  }

  # REAC/WP via g:Profiler
  gp <- try(
    gprofiler2::gost(
      query             = genes,
      organism          = "mmusculus",
      sources           = db,
      correction_method = "g_SCS"
    ),
    silent = TRUE
  )

  if (!inherits(gp, "try-error") && !is.null(gp$result) && nrow(gp$result) > 0) {
    df <- gp$result
    comb <- -log10(df$p_value) * (df$intersection_size / df$term_size)
    df2 <- data.frame(
      Term             = df$term_name,
      Adjusted.P.value = df$p_value,
      Combined.Score   = comb,
      intersection     = df$intersection,
      stringsAsFactors = FALSE
    )
    df2 <- add_gene_mirna_cols(
      df2,
      genes_col = "intersection",
      id_type   = "entrez",
      target_df = target_df,
      symbol_col = "SYMBOL",
      entrez_col = "ENTREZID",
      mirna_col  = "miRNA"
    )
    return(df2[order(df2$Adjusted.P.value), ])
  }

  data.frame()
}

# ---------------------- Zebrafish (dre): human targets → orthology → drerio ----------------------
.to_hsa_from_dre <- function(mi) sub("^dre-", "hsa-", as.character(mi), ignore.case = TRUE)

dr_targets_df <- function(mirnas) {
  # Convert dre-miR-XXX to hsa-miR-XXX for multiMiR
  hsa_mirs <- unique(stats::na.omit(.to_hsa_from_dre(mirnas)))
  if (!length(hsa_mirs)) {
    return(data.frame(miRNA = character(),
                      SYMBOL = character(),
                      ENTREZID = character(),
                      stringsAsFactors = FALSE))
  }

  pull_mm <- function(tbl) {
    try(R.utils::withTimeout(
      multiMiR::get_multimir(org = "hsa",
                             mirna = hsa_mirs,
                             table = tbl,
                             summary = FALSE),
      timeout = 2000,
      onTimeout = "error"
    ), silent = TRUE)
  }

  v <- pull_mm("validated")
  p <- pull_mm("predicted")

  hdat <- dplyr::bind_rows(
    if (!inherits(v, "try-error") && !is.null(v) && nrow(v@data)) as.data.frame(v@data) else NULL,
    if (!inherits(p, "try-error") && !is.null(p) && nrow(p@data)) as.data.frame(p@data) else NULL
  )

  if (!nrow(hdat)) {
    return(data.frame(miRNA = character(),
                      SYMBOL = character(),
                      ENTREZID = character(),
                      stringsAsFactors = FALSE))
  }

  # Keep human miRNA + human target symbol
  hdat_small <- hdat %>%
    dplyr::transmute(
      hsa_miRNA    = as.character(mature_mirna_id),
      human_symbol = as.character(target_symbol)
    ) %>%
    dplyr::filter(
      !is.na(hsa_miRNA), hsa_miRNA != "",
      !is.na(human_symbol), human_symbol != ""
    ) %>%
    dplyr::distinct()

  # Orthology: human_symbol -> zebrafish symbol
  ortho <- try(
    gprofiler2::gorth(unique(hdat_small$human_symbol), "hsapiens", "drerio"),
    silent = TRUE
  )
  if (inherits(ortho, "try-error") || is.null(ortho) || !nrow(ortho)) {
    return(data.frame(miRNA = character(),
                      SYMBOL = character(),
                      ENTREZID = character(),
                      stringsAsFactors = FALSE))
  }

  # Identify columns for human and zebrafish gene names
  human_col <- if ("name" %in% colnames(ortho)) {
    "name"
  } else if ("input" %in% colnames(ortho)) {
    "input"
  } else {
    colnames(ortho)[1]
  }

  dr_col <- if ("ortholog_name" %in% colnames(ortho)) {
    "ortholog_name"
  } else if ("target_name" %in% colnames(ortho)) {
    "target_name"
  } else if ("name" %in% colnames(ortho) && human_col != "name") {
    "name"
  } else {
    colnames(ortho)[2]
  }

  map_df <- ortho[, c(human_col, dr_col)]
  colnames(map_df) <- c("human_symbol", "dr_symbol")
  map_df <- map_df %>%
    dplyr::filter(!is.na(dr_symbol), dr_symbol != "") %>%
    dplyr::distinct()

  # Join to get dre genes per miRNA
  joined <- dplyr::inner_join(hdat_small, map_df, by = "human_symbol")
  if (!nrow(joined)) {
    return(data.frame(miRNA = character(),
                      SYMBOL = character(),
                      ENTREZID = character(),
                      stringsAsFactors = FALSE))
  }

  # Map zebrafish symbols -> ENTREZID
  ann <- safe_select(org.Dr.eg.db,
                     unique(joined$dr_symbol),
                     keytype = "SYMBOL",
                     columns = "ENTREZID")
  if (!nrow(ann)) {
    return(data.frame(miRNA = character(),
                      SYMBOL = character(),
                      ENTREZID = character(),
                      stringsAsFactors = FALSE))
  }

  out <- dplyr::inner_join(joined, ann, by = c("dr_symbol" = "SYMBOL")) %>%
    dplyr::transmute(
      miRNA   = sub("^hsa-", "dre-", hsa_miRNA, ignore.case = TRUE),
      SYMBOL  = as.character(dr_symbol),
      ENTREZID = as.character(ENTREZID)
    ) %>%
    dplyr::filter(!is.na(ENTREZID), ENTREZID != "") %>%
    dplyr::distinct()

  out
}

dr_enrich_from_targets <- function(target_df, db = "GO:BP", q = 0.05) {
  if (!is.data.frame(target_df) || !nrow(target_df)) return(data.frame())

  # Ensure columns exist
  if (!"miRNA" %in% names(target_df))  target_df$miRNA  <- NA_character_
  if (!"SYMBOL" %in% names(target_df)) target_df$SYMBOL <- NA_character_

  genes <- unique(target_df$ENTREZID)
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) < 10) return(data.frame())

  base_df <- NULL

  if (db == "GO:BP") {
    eg <- clusterProfiler::enrichGO(
      gene          = genes,
      OrgDb         = org.Dr.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff  = q,
      readable      = FALSE
    )
    base_df <- as.data.frame(eg)
  } else if (db == "KEGG") {
    kk <- clusterProfiler::enrichKEGG(
      gene          = genes,
      organism      = "dre",
      pAdjustMethod = "BH",
      qvalueCutoff  = q
    )
    base_df <- as.data.frame(kk)
  } else {
    gp <- try(
      gprofiler2::gost(
        query              = genes,
        organism           = "drerio",
        sources            = db,
        correction_method  = "g_SCS"
      ),
      silent = TRUE
    )
    if (inherits(gp, "try-error") || is.null(gp$result) || !nrow(gp$result)) {
      return(data.frame())
    }
    base_df <- gp$result
  }

  if (!nrow(base_df)) return(data.frame())

  df <- base_df

  # Normalise column names
  if (!"p.adjust" %in% names(df) && "p_value" %in% names(df))      df$p.adjust   <- df$p_value
  if (!"Description" %in% names(df) && "term_name" %in% names(df)) df$Description <- df$term_name

  # Which column carries the genes per term?
  gene_col <- if ("geneID" %in% names(df)) {
    "geneID"
  } else if ("intersection" %in% names(df)) {
    "intersection"
  } else {
    NA_character_
  }

  target_min <- target_df %>%
    dplyr::select(ENTREZID, SYMBOL, miRNA) %>%
    dplyr::distinct()

  genes_list <- character(nrow(df))
  mirna_list <- character(nrow(df))

  for (i in seq_len(nrow(df))) {
    gvec <- character(0)
    if (!is.na(gene_col)) {
      raw <- as.character(df[[gene_col]][i])
      if (!is.na(raw) && nzchar(raw)) {
        sep <- if (gene_col == "geneID") "/" else ","
        gvec <- strsplit(raw, sep, fixed = TRUE)[[1]]
        gvec <- trimws(gvec)
      }
    }

    # For clusterProfiler, geneID is ENTREZ; for gprofiler, intersection is usually symbols
    if (gene_col == "geneID") {
      subdf <- dplyr::filter(target_min, ENTREZID %in% gvec)
    } else if (gene_col == "intersection") {
      subdf <- dplyr::filter(target_min, SYMBOL %in% gvec)
    } else {
      subdf <- target_min[0, ]
    }

    genes_list[i] <- paste(sort(unique(subdf$SYMBOL)), collapse = ",")
    mirna_list[i] <- paste(sort(unique(subdf$miRNA)),  collapse = ",")
  }

  out <- data.frame(
    Term             = df$Description,
    Adjusted.P.value = df$p.adjust,
    Combined.Score   = -log10(df$p.adjust),
    Genes            = genes_list,
    miRNAs           = mirna_list,
    stringsAsFactors = FALSE
  )

  out <- out[order(out$Adjusted.P.value), ]
  out
}

# ---------------------- Fly (dme): human targets → orthology → dmel ----------------------
.to_hsa_from_dme <- function(mi) sub("^dme-", "hsa-", as.character(mi), ignore.case = TRUE)

dm_targets_df <- function(mirnas) {
  # Convert dme-miR-XXX to hsa-miR-XXX
  hsa_mirs <- unique(stats::na.omit(.to_hsa_from_dme(mirnas)))
  if (!length(hsa_mirs)) {
    return(data.frame(miRNA = character(),
                      SYMBOL = character(),
                      ENTREZID = character(),
                      stringsAsFactors = FALSE))
  }

  pull_mm <- function(tbl) {
    try(R.utils::withTimeout(
      multiMiR::get_multimir(org = "hsa",
                             mirna = hsa_mirs,
                             table = tbl,
                             summary = FALSE),
      timeout = 25,
      onTimeout = "error"
    ), silent = TRUE)
  }

  v <- pull_mm("validated")
  p <- pull_mm("predicted")

  hdat <- dplyr::bind_rows(
    if (!inherits(v, "try-error") && !is.null(v) && nrow(v@data)) as.data.frame(v@data) else NULL,
    if (!inherits(p, "try-error") && !is.null(p) && nrow(p@data)) as.data.frame(p@data) else NULL
  )

  if (!nrow(hdat)) {
    return(data.frame(miRNA = character(),
                      SYMBOL = character(),
                      ENTREZID = character(),
                      stringsAsFactors = FALSE))
  }

  hdat_small <- hdat %>%
    dplyr::transmute(
      hsa_miRNA    = as.character(mature_mirna_id),
      human_symbol = as.character(target_symbol)
    ) %>%
    dplyr::filter(
      !is.na(hsa_miRNA), hsa_miRNA != "",
      !is.na(human_symbol), human_symbol != ""
    ) %>%
    dplyr::distinct()

  ortho <- try(
    gprofiler2::gorth(unique(hdat_small$human_symbol), "hsapiens", "dmelanogaster"),
    silent = TRUE
  )
  if (inherits(ortho, "try-error") || is.null(ortho) || !nrow(ortho)) {
    return(data.frame(miRNA = character(),
                      SYMBOL = character(),
                      ENTREZID = character(),
                      stringsAsFactors = FALSE))
  }

  # Columns for human vs fly symbols
  human_col <- if ("name" %in% colnames(ortho)) {
    "name"
  } else if ("input" %in% colnames(ortho)) {
    "input"
  } else {
    colnames(ortho)[1]
  }

  dm_col <- if ("ortholog_name" %in% colnames(ortho)) {
    "ortholog_name"
  } else if ("target_name" %in% colnames(ortho)) {
    "target_name"
  } else if ("name" %in% colnames(ortho) && human_col != "name") {
    "name"
  } else {
    colnames(ortho)[2]
  }

  map_df <- ortho[, c(human_col, dm_col)]
  colnames(map_df) <- c("human_symbol", "dm_symbol")
  map_df <- map_df %>%
    dplyr::filter(!is.na(dm_symbol), dm_symbol != "") %>%
    dplyr::distinct()

  joined <- dplyr::inner_join(hdat_small, map_df, by = "human_symbol")
  if (!nrow(joined)) {
    return(data.frame(miRNA = character(),
                      SYMBOL = character(),
                      ENTREZID = character(),
                      stringsAsFactors = FALSE))
  }

  ann <- safe_select(org.Dm.eg.db,
                     unique(joined$dm_symbol),
                     keytype = "SYMBOL",
                     columns = "ENTREZID")
  if (!nrow(ann)) {
    return(data.frame(miRNA = character(),
                      SYMBOL = character(),
                      ENTREZID = character(),
                      stringsAsFactors = FALSE))
  }

  out <- dplyr::inner_join(joined, ann, by = c("dm_symbol" = "SYMBOL")) %>%
    dplyr::transmute(
      miRNA   = sub("^hsa-", "dme-", hsa_miRNA, ignore.case = TRUE),
      SYMBOL  = as.character(dm_symbol),
      ENTREZID = as.character(ENTREZID)
    ) %>%
    dplyr::filter(!is.na(ENTREZID), ENTREZID != "") %>%
    dplyr::distinct()

  out
}

dm_enrich_from_targets <- function(target_df, db = "GO:BP", q = 0.05) {
  if (!is.data.frame(target_df) || !nrow(target_df)) return(data.frame())

  if (!"miRNA" %in% names(target_df))  target_df$miRNA  <- NA_character_
  if (!"SYMBOL" %in% names(target_df)) target_df$SYMBOL <- NA_character_

  genes <- unique(target_df$ENTREZID)
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) < 10) return(data.frame())

  base_df <- NULL

  if (db == "GO:BP") {
    eg <- clusterProfiler::enrichGO(
      gene          = genes,
      OrgDb         = org.Dm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff  = q,
      readable      = FALSE
    )
    base_df <- as.data.frame(eg)
  } else if (db == "KEGG") {
    kk <- clusterProfiler::enrichKEGG(
      gene          = genes,
      organism      = "dme",
      pAdjustMethod = "BH",
      qvalueCutoff  = q
    )
    base_df <- as.data.frame(kk)
  } else {
    gp <- try(
      gprofiler2::gost(
        query              = genes,
        organism           = "dmelanogaster",
        sources            = db,
        correction_method  = "g_SCS"
      ),
      silent = TRUE
    )
    if (inherits(gp, "try-error") || is.null(gp$result) || !nrow(gp$result)) {
      return(data.frame())
    }
    base_df <- gp$result
  }

  if (!nrow(base_df)) return(data.frame())

  df <- base_df
  if (!"p.adjust" %in% names(df) && "p_value" %in% names(df))      df$p.adjust   <- df$p_value
  if (!"Description" %in% names(df) && "term_name" %in% names(df)) df$Description <- df$term_name

  gene_col <- if ("geneID" %in% names(df)) {
    "geneID"
  } else if ("intersection" %in% names(df)) {
    "intersection"
  } else {
    NA_character_
  }

  target_min <- target_df %>%
    dplyr::select(ENTREZID, SYMBOL, miRNA) %>%
    dplyr::distinct()

  genes_list <- character(nrow(df))
  mirna_list <- character(nrow(df))

  for (i in seq_len(nrow(df))) {
    gvec <- character(0)
    if (!is.na(gene_col)) {
      raw <- as.character(df[[gene_col]][i])
      if (!is.na(raw) && nzchar(raw)) {
        sep <- if (gene_col == "geneID") "/" else ","
        gvec <- strsplit(raw, sep, fixed = TRUE)[[1]]
        gvec <- trimws(gvec)
      }
    }

    if (gene_col == "geneID") {
      subdf <- dplyr::filter(target_min, ENTREZID %in% gvec)
    } else if (gene_col == "intersection") {
      subdf <- dplyr::filter(target_min, SYMBOL %in% gvec)
    } else {
      subdf <- target_min[0, ]
    }

    genes_list[i] <- paste(sort(unique(subdf$SYMBOL)), collapse = ",")
    mirna_list[i] <- paste(sort(unique(subdf$miRNA)),  collapse = ",")
  }

  out <- data.frame(
    Term             = df$Description,
    Adjusted.P.value = df$p.adjust,
    Combined.Score   = -log10(df$p.adjust),
    Genes            = genes_list,
    miRNAs           = mirna_list,
    stringsAsFactors = FALSE
  )

  out <- out[order(out$Adjusted.P.value), ]
  out
}

# ---------------------- Generic cross-species enrich (fallback) ----------------------
enrich_species <- function(genes, db, species) {
  genes <- unique(stats::na.omit(genes))
  genes <- genes[genes != ""]
  if (length(genes) < 2) return(data.frame())

  if (species == "hs" && db %in% enrichr_human_choices) {
    ee <- .enrichr_pull(genes, db)
    if (nrow(ee)) return(ee)
    db <- switch(db,
                 "GO_Biological_Process_2021" = "GO:BP",
                 "GO_Molecular_Function_2021" = "GO:MF",
                 "GO_Cellular_Component_2021" = "GO:CC",
                 "KEGG_2021_Human"            = "KEGG",
                 "Reactome_2022"              = "REAC",
                 "WikiPathway_2021_Human"     = "WP",
                 "GO:BP")
  }

  org <- switch(
    species,
    hs = "hsapiens",
    mm = "mmusculus",
    dr = "drerio",
    dm = "dmelanogaster",
    "hsapiens"
  )

  gp <- try(
    gprofiler2::gost(
      query             = genes,
      organism          = org,
      sources           = db,
      correction_method = "g_SCS"
    ),
    silent = TRUE
  )

  if (!inherits(gp, "try-error") && !is.null(gp$result) && nrow(gp$result) > 0) {
    df   <- gp$result
    comb <- -log10(df$p_value) * (df$intersection_size / df$term_size)
    out  <- data.frame(
      Term             = df$term_name,
      Adjusted.P.value = df$p_value,
      Combined.Score   = comb,
      stringsAsFactors = FALSE
    )
    return(out[order(out$Adjusted.P.value), ])
  }

  data.frame()
}

# ===================== Generic target retrieval (kept, no miRNAtap) =====================
get_targets <- function(mirnas, species) {
  mirnas <- unique(stats::na.omit(mirnas))
  if (!length(mirnas)) {
    return(data.frame(mirna = character(), target = character(),
                      evidence = character(), source = character(),
                      stringsAsFactors = FALSE))
  }
  
  sp <- tolower(as.character(species))
  mm_org <- MM_ORG[[sp]] %||% NA_character_
  
  # Only human/mouse supported by this generic multiMiR path
  if (is.na(mm_org) || !mm_org %in% c("hsa", "mmu")) {
    return(data.frame(mirna = character(), target = character(),
                      evidence = character(), source = character(),
                      stringsAsFactors = FALSE))
  }
  
  # Normalize IDs BEFORE querying multiMiR
  sp2 <- if (mm_org == "hsa") "hs" else "mm"
  mirnas_norm <- normalize_mirna_ids(mirnas, species = sp2)
  
  safe_multimir_chunked <- function(org, mirnas, table, chunk_size = 50) {
    mirnas <- unique(stats::na.omit(mirnas))
    if (!length(mirnas)) return(NULL)
    
    chunks <- split(mirnas, ceiling(seq_along(mirnas) / chunk_size))
    out <- list()
    
    for (i in seq_along(chunks)) {
      m <- chunks[[i]]
      res <- try(
        R.utils::withTimeout(
          multiMiR::get_multimir(org = org, mirna = m, table = table, summary = FALSE),
          timeout = 2000, onTimeout = "error"
        ),
        silent = TRUE
      )
      if (!inherits(res, "try-error") && !is.null(res) && nrow(res@data) > 0) {
        out[[length(out) + 1]] <- as.data.frame(res@data)
      }
    }
    
    if (!length(out)) return(NULL)
    do.call(rbind, out)
  }
  
  out <- list()
  
  mm_val_df <- safe_multimir_chunked(mm_org, mirnas_norm, "validated")
  if (!is.null(mm_val_df) && nrow(mm_val_df) > 0) {
    df <- mm_val_df[, c("mature_mirna_id", "target_symbol", "database"), drop = FALSE]
    names(df) <- c("mirna", "target", "source")
    df$evidence <- "validated"
    out$validated <- df[, c("mirna", "target", "evidence", "source")]
  }
  
  mm_pred_df <- safe_multimir_chunked(mm_org, mirnas_norm, "predicted")
  if (!is.null(mm_pred_df) && nrow(mm_pred_df) > 0) {
    df <- mm_pred_df[, c("mature_mirna_id", "target_symbol", "database"), drop = FALSE]
    names(df) <- c("mirna", "target", "source")
    df$evidence <- "predicted"
    out$predicted <- df[, c("mirna", "target", "evidence", "source")]
  }
  
  if (!length(out)) {
    return(data.frame(mirna = character(), target = character(),
                      evidence = character(), source = character(),
                      stringsAsFactors = FALSE))
  }
  
  ans <- do.call(rbind, out)
  ans <- ans[!duplicated(ans[, c("mirna", "target", "evidence")]), ]
  rownames(ans) <- NULL
  ans
}

 
  

# ---- Plot/table helpers ----
renderEnrichPlot <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- as.data.frame(df)
  names(df) <- tolower(gsub("\\.|\\s","_",names(df)))
  adj_p_col <- which(names(df) %in% c("adjusted_p_value","adjusted.p.value","adj_p","padj","p_value"))[1]
  term_col  <- which(names(df) %in% c("term","description","term_name"))[1]
  if (is.na(adj_p_col) || is.na(term_col)) return(NULL)
  df <- df %>% dplyr::rename(adj_p = !!names(df)[adj_p_col], term = !!names(df)[term_col])
  top <- df %>% arrange(adj_p) %>% slice_head(n = 10)
  plot_ly(top, x = ~reorder(term, -log10(adj_p)), y = ~-log10(adj_p), type = "bar",
          hovertemplate = "Term: %{x}<br>-log10(p): %{y:.2f}<extra></extra>") %>%
    layout(xaxis = list(title=""), yaxis = list(title="-log10 adjusted p"), margin = list(b=120))
}

renderEnrichTable <- function(df) {
  if (is.null(df) || nrow(df) == 0)
    return(data.frame(Message = "No enrichment results"))

  d  <- as.data.frame(df)
  nm <- tolower(gsub("\\.", "_", names(d)))

  # Detect key columns (robust to different tools)
  ap <- which(nm %in% c("adjusted_p_value", "adjusted.p.value", "p_value", "padj"))[1]
  tm <- which(nm %in% c("term", "description", "term_name"))[1]
  gn <- which(nm %in% c("genes", "gene", "geneid", "gene_id", "overlap_genes"))[1]

  if (!is.na(ap)) names(d)[ap] <- "Adjusted_P_value"
  if (!is.na(tm)) names(d)[tm] <- "Term"
  if (!is.na(gn)) names(d)[gn] <- "Genes"

  # Decide which columns to show
  keep <- c("Term", "Adjusted_P_value", "Genes")
  keep <- keep[keep %in% names(d)]

  d <- d[, keep, drop = FALSE]

  if ("Adjusted_P_value" %in% names(d))
    d$Adjusted_P_value <- signif(d$Adjusted_P_value, 4)

  d
}


# ===================== UI =====================
ui <- fluidPage(
  titlePanel("JCAP miRNA-SEQ (multi-species)"),
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
  sidebarLayout(
    sidebarPanel(
      fileInput("countsFile", "Upload Count Matrix CSV", accept = ".csv"),
      fileInput("metaFile", "Upload Sample Metadata CSV", accept = ".csv"),
      selectInput("species", "Species", choices = species_labels, selected = "auto"),
      uiOutput("db_picker"),
      actionButton("runAnalysis", "Run DE Analysis"),
      hr(),
      actionButton("volcanoBtn", "Volcano"), actionButton("pcaBtn", "PCA"),
      actionButton("umapBtn", "UMAP"), actionButton("heatmapBtn", "Heatmap"),
      actionButton("barplotBtn", "Top miRNAs Barplot"),
      hr(),
      actionButton("enrichAllBtn", "Enrichment: All DE (targets)"),
      actionButton("enrichUpBtn",  "Enrichment: Up (targets)"),
      actionButton("enrichDownBtn","Enrichment: Down (targets)"),
      actionButton("plotEnrichAllBtn",  "Plot All"),
      actionButton("plotEnrichUpBtn",   "Plot Up"),
      actionButton("plotEnrichDownBtn", "Plot Down"),
      hr(),
      actionButton("runRF", "Random Forest"), actionButton("runPower", "Power Analysis"),
      hr(),
      downloadButton("downloadTop", "Download Top DE miRNAs")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Top miRNAs", tableOutput("topTable")),
        tabPanel("Top miRNA Barplot", plotlyOutput("barplotPlot")),
        tabPanel("Volcano", plotlyOutput("volcanoPlot")),
        tabPanel("PCA", plotlyOutput("pcaPlot")),
        tabPanel("UMAP", plotlyOutput("umapPlot")),
        tabPanel("Heatmap", plotlyOutput("heatmapPlot")),
        tabPanel("Enrich: All DE",
                 tabsetPanel(
                   tabPanel("Table",  downloadButton("downloadEnrichAll", "Download CSV"), br(), tableOutput("enrichAllTable")),
                   tabPanel("Barplot", plotlyOutput("enrichAllPlot"))
                 )),
        tabPanel("Enrich: Upregulated",
                 tabsetPanel(
                   tabPanel("Table",  downloadButton("downloadEnrichUp", "Download CSV"), br(), tableOutput("enrichUpTable")),
                   tabPanel("Barplot", plotlyOutput("enrichUpPlot"))
                 )),
        tabPanel("Enrich: Downregulated",
                 tabsetPanel(
                   tabPanel("Table",  downloadButton("downloadEnrichDown", "Download CSV"), br(), tableOutput("enrichDownTable")),
                   tabPanel("Barplot", plotlyOutput("enrichDownPlot"))
                 )),
        tabPanel("Random Forest",
                 tabsetPanel(
                   tabPanel("Predictions", downloadButton("downloadRFpred","Download Predictions"), br(), tableOutput("rfPredictions")),
                   tabPanel("Metrics",     downloadButton("downloadRFmetrics","Download Metrics"), br(), tableOutput("rfMetrics")),
                   tabPanel("Importance",  downloadButton("downloadRFimportance","Download Importance Table"), br(), tableOutput("rfImportance"))
                 )),
        tabPanel("Power Analysis",
                 tabsetPanel(
                   tabPanel("Power Table", downloadButton("downloadPowerTable","Download Power Table (CSV)"), br(), DTOutput("power_table")),
                   tabPanel("Power Plot",  plotOutput("power_plot"))
                 )),
        tabPanel("README", verbatimTextOutput("readmeText"))
      )
    )
  )
)

# ==================== SERVER ====================
server <- function(input, output, session) {
  resultsDF   <- reactiveVal(NULL)
  ddsData     <- reactiveVal(NULL)
  vsdData     <- reactiveVal(NULL)
  power_matrix <- reactiveVal(NULL)

  enrichAll  <- reactiveVal(NULL)
  enrichUp   <- reactiveVal(NULL)
  enrichDown <- reactiveVal(NULL)
  deAllDF <- reactiveVal(NULL)   # FULL DESeq2 results (all miRNAs)
  
  rf_preds      <- reactiveVal(NULL)
  rf_metrics    <- reactiveVal(NULL)
  rf_importance <- reactiveVal(NULL)

  target_cache <- reactiveValues()
  cached_get_targets <- function(mirnas, species) {
    key <- paste0(species, "_", digest::digest(sort(mirnas)))
    if (!is.null(target_cache[[key]])) return(target_cache[[key]])
    tg <- get_targets(mirnas, species); target_cache[[key]] <- tg; tg
  }

  output$db_picker <- renderUI({
    sp <- input$species
    choices <- db_choices_for_species(if (sp %in% c("hs","mm","dr","dm")) sp else "hs")
    selectInput("selectedDB", "Pathway Database", choices = choices, selected = choices[[1]])
  })

  observeEvent(input$runAnalysis, {
  req(input$countsFile, input$metaFile)
  showNotification("🧬 DESeq2 running…", type = "message")

  
# counts
ext_counts <- tools::file_ext(input$countsFile$name)
count_data <- if (ext_counts == "csv") {
  read.csv(input$countsFile$datapath, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  read.table(input$countsFile$datapath, header = TRUE, sep = "\t",
             check.names = FALSE, stringsAsFactors = FALSE)
}

# first column is feature IDs (miRNAs)
first_col_counts <- names(count_data)[1]
count_matrix <- as.matrix(count_data[, -1, drop = FALSE])
rownames(count_matrix) <- count_data[[first_col_counts]]

# metadata
ext_meta <- tools::file_ext(input$metaFile$name)
meta_data <- if (ext_meta == "csv") {
  read.csv(input$metaFile$datapath, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  read.table(input$metaFile$datapath, header = TRUE, sep = "\t",
             check.names = FALSE, stringsAsFactors = FALSE)
}

# ---- NOW do minimal, safe column cleanup ----
names(meta_data) <- trimws(names(meta_data))
names(meta_data) <- sub("^\ufeff", "", names(meta_data))   # remove BOM if present
names(meta_data) <- tolower(names(meta_data))
names(meta_data) <- gsub("[^a-z0-9]+", "_", names(meta_data))
names(meta_data) <- gsub("^_+|_+$", "", names(meta_data))

# Recover the common failure mode you observed
names(meta_data)[names(meta_data) %in% c("ample_id","ampleid")] <- "sample_id"

# require expected columns
if (!("sample_id" %in% names(meta_data))) {
  showNotification(
    paste0("Metadata must have 'sample_id'. I see: ",
           paste(names(meta_data), collapse = ", ")),
    type = "error"
  )
  return(NULL)
}

rownames(meta_data) <- meta_data$sample_id

  # accept condition as group
  if (!("group" %in% names(meta_data)) && ("condition" %in% names(meta_data))) {
    meta_data$group <- meta_data$condition
  }

  # debug: show what the app sees
  showNotification(paste("Meta cols:", paste(names(meta_data), collapse=" | ")),
                   type="message", duration=15)

  # accept condition as group
  if (!("group" %in% names(meta_data)) && ("condition" %in% names(meta_data))) {
    meta_data$group <- meta_data$condition
  }

  if (!("sample_id" %in% names(meta_data))) {
    showNotification("Metadata must have a 'sample_id' column.", type="error")
    return(NULL)
  }
  if (!("group" %in% names(meta_data))) {
    showNotification("Metadata must have a 'group' (or 'condition') column.", type="error")
    return(NULL)
  }

  rownames(meta_data) <- meta_data$sample_id

  common <- intersect(colnames(count_matrix), rownames(meta_data))
  if (length(common) < 4) {
    showNotification("Need ≥4 overlapping samples.", type="error")
    return(NULL)
  }

  count_matrix <- count_matrix[, common, drop = FALSE]
  meta_data    <- meta_data[common, , drop = FALSE]

  # ... keep the rest of your DESeq2 code exactly as you had it ...# ... keep the rest of your DESeq2 code exactly as you had it ..


# --- Dynamically detect groups from metadata$group ---
grp_levels <- levels(factor(meta_data$group))

# Require exactly 2 groups for DESeq2 contrast
if (length(grp_levels) != 2) {
  showNotification(
    paste0(
      "Group column 'group' must have exactly 2 levels for DE analysis. Found: ",
      paste(grp_levels, collapse = ", ")
    ),
    type = "error"
  )
  return(NULL)
}

# Make 'condition' the factor DESeq2 will use
meta_data$condition <- factor(meta_data$group, levels = grp_levels)

# Optional batch covariate
use_batch <- "batch" %in% names(meta_data) && nlevels(factor(meta_data$batch)) > 1
if (use_batch) meta_data$batch <- factor(meta_data$batch)

# Check replicates per group
reps <- table(meta_data$condition)
if (any(reps < 2)) {
  showNotification(
    paste("Need ≥2 replicates per group. Found:",
          paste(names(reps), reps, collapse = " / ")),
    type = "error"
  )
  return(NULL)
}

design_formula <- if (use_batch) ~ batch + condition else ~ condition

    # ---------- DESeq2 (robust with dispersion fallbacks) ----------
    dds <- DESeqDataSetFromMatrix(
      countData = round(count_matrix),
      colData   = meta_data,
      design    = design_formula
    )
    dds <- dds[rowSums(counts(dds)) > 0, ]

    # Full-rank check (avoid “coefficients == samples”)
    X <- model.matrix(design(dds), data = as.data.frame(colData(dds)))
    if (qr(X)$rank == nrow(X)) {
      showNotification("Design is full rank (coefficients == samples). Simplify (drop batch).", type = "error")
      return(NULL)
    }

    run_deseq2_res <- function(dds) {
      # 1) Try default (parametric fit)
      res <- try({
        dds1 <- DESeq(dds, fitType = "parametric", minReplicatesForReplace = Inf)
        list(dds = dds1, method = "parametric")
      }, silent = TRUE)
      if (!inherits(res, "try-error")) return(res)

      # 2) Try local fit
      res <- try({
        dds2 <- DESeq(dds, fitType = "local", minReplicatesForReplace = Inf)
        list(dds = dds2, method = "local")
      }, silent = TRUE)
      if (!inherits(res, "try-error")) return(res)

      # 3) Try mean fit
      res <- try({
        dds3 <- DESeq(dds, fitType = "mean", minReplicatesForReplace = Inf)
        list(dds = dds3, method = "mean")
      }, silent = TRUE)
      if (!inherits(res, "try-error")) return(res)

      # 4) Final fallback: use gene-wise dispersions, then Wald test
      message("⚠️ Falling back to gene-wise dispersions -> nbinomWaldTest")
      dds4 <- estimateSizeFactors(dds)
      dds4 <- estimateDispersionsGeneEst(dds4)
      dispersions(dds4) <- mcols(dds4)$dispGeneEst
      dds4 <- nbinomWaldTest(dds4)
      list(dds = dds4, method = "gene-wise")
    }

    fit <- run_deseq2_res(dds)
    dds <- fit$dds
    showNotification(paste("DESeq2 dispersion fit:", fit$method), type = "message", duration = 4)

# ----- Build results object (robust; no hard-coded group names) -----
rn <- resultsNames(dds)

# Prefer the DESeq2-generated condition coef(s), e.g. "condition_Treated_vs_Control"
cond_coef <- rn[grepl("^condition_.+_vs_.+$", rn)]

# Determine factor levels from the dataset itself
cond_levels <- levels(colData(dds)$condition)

# Choose a sensible default contrast: (level2 vs level1)
# (DESeq2’s coef is also typically “levelX_vs_level1” where level1 is reference.)
if (length(cond_levels) >= 2) {
  ref_level <- cond_levels[1]
  test_level <- cond_levels[2]
} else {
  ref_level <- NA_character_
  test_level <- NA_character_
}

res <- tryCatch({
  if (length(cond_coef) == 1) {
    # We have a single unambiguous coefficient; use it
    if (requireNamespace("apeglm", quietly = TRUE) && fit$method != "gene-wise") {
      DESeq2::lfcShrink(dds, coef = cond_coef, type = "apeglm")
    } else {
      DESeq2::results(dds, name = cond_coef)
    }
  } else {
    # Multiple (or zero) coefficients: fall back to explicit contrast from factor levels
    if (is.na(test_level) || is.na(ref_level)) {
      stop("Could not determine two condition levels for contrast.")
    }
    DESeq2::results(dds, contrast = c("condition", test_level, ref_level))
  }
}, error = function(e) {
  # Final fallback: if coefficient matching failed, still try contrast from levels
  if (!is.na(test_level) && !is.na(ref_level)) {
    DESeq2::results(dds, contrast = c("condition", test_level, ref_level))
  } else {
    stop(e)
  }
})

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("miRNA") %>%
  dplyr::mutate(sig = !is.na(padj) & padj < 0.1)

deAllDF(res_df)  # <-- ADD THIS LINE



    # ---------- RandomForest ranking on top-20 by padj ----------
    top20 <- res_df %>%
      dplyr::filter(!is.na(padj)) %>%
      dplyr::arrange(padj) %>%
      dplyr::slice_head(n = 20)

    if (nrow(top20) < 2) {
      showNotification("Not enough DE features to rank.", type = "message"); return(NULL)
    }

    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    vst <- assay(vsd)

    # Harmonize names for safe column handling
    vst_ren <- vst
    rownames(vst_ren) <- make.names(rownames(vst))
    top20$miRNA_clean <- make.names(top20$miRNA)

    sub <- vst_ren[rownames(vst_ren) %in% top20$miRNA_clean, , drop = FALSE]
    if (nrow(sub) < 2) {
      showNotification("Top DE features not found in VST matrix.", type = "error"); return(NULL)
    }

    df_rf <- as.data.frame(t(sub))
    df_rf$condition <- colData(vsd)$condition

    set.seed(1)
    rf <- randomForest::randomForest(condition ~ ., data = df_rf, importance = TRUE)
    imp <- randomForest::importance(rf, type = 2)
    imp_df <- data.frame(
      miRNA_clean = rownames(imp),
      Importance  = imp[, 1]
    ) %>%
      dplyr::arrange(dplyr::desc(Importance)) %>%
      dplyr::slice_head(n = 11)

    map <- data.frame(miRNA = top20$miRNA, miRNA_clean = top20$miRNA_clean)
    imp_df <- merge(imp_df, map, by = "miRNA_clean", all.x = TRUE)

    final_hits <- res_df %>% dplyr::filter(miRNA %in% imp_df$miRNA)
    if (!nrow(final_hits)) {
      showNotification("No RF-ranked hits.", type = "error"); return(NULL)
    }

    resultsDF(final_hits)
    ddsData(dds)
    vsdData(vsd)
    showNotification("✅ DE + RF complete", type = "message")
  })

  # ===================== Enrichment server block =====================
  get_species_current <- reactive({
    sp <- input$species
    if (sp == "auto" && !is.null(resultsDF()))
      detect_species(resultsDF()$miRNA, "hs")
    else
      sp
  })

  do_enrich <- function(which_set = c("all","up","down")) {
    req(deAllDF())
    df <- deAllDF()
    
    which_set <- match.arg(which_set)

    mirnas <- switch(
      which_set,
      all = df %>% dplyr::filter(!is.na(padj), padj < 0.1) %>% dplyr::pull(miRNA),
      
      up   = df %>% dplyr::filter(!is.na(padj), padj < 0.1, log2FoldChange > 0) %>% dplyr::pull(miRNA),
      down = df %>% dplyr::filter(!is.na(padj), padj < 0.1, log2FoldChange < 0) %>% dplyr::pull(miRNA)
    )
    mirnas <- unique(stats::na.omit(mirnas))
    if (!length(mirnas)) return(data.frame())

    sp <- tryCatch(tolower(as.character(get_species_current())), error = function(e) "hs")
    if (length(sp) != 1 || is.na(sp) || sp == "") sp <- "hs"

    if (sp %in% c("mm","mouse","mmu")) {
      showNotification("Fetching mouse targets (multiMiR)…", type="message", duration=3)
      tg_df <- try(mm_targets_df(mirnas), silent = TRUE)
      if (inherits(tg_df,"try-error") || !is.data.frame(tg_df) || !nrow(tg_df)) {
        showNotification("No mouse targets returned.", type="error"); return(data.frame())
      }
      showNotification(paste("Targets found:", nrow(tg_df)), type="message", duration=3)
      ans <- try(mm_enrich_from_targets(tg_df, db = input$selectedDB, q = 0.05), silent = TRUE)
      return(if (inherits(ans,"try-error") || is.null(ans)) data.frame() else ans)
    }

    if (sp %in% c("dr","dre","zebrafish")) {
      showNotification("Fetching zebrafish targets (orthology)…", type="message", duration=3)
      tg_df <- try(dr_targets_df(mirnas), silent = TRUE)
      if (inherits(tg_df,"try-error") || !is.data.frame(tg_df) || !nrow(tg_df)) {
        showNotification("No zebrafish targets returned.", type="error"); return(data.frame())
      }
      showNotification(paste("Targets found:", nrow(tg_df)), type="message", duration=3)
      ans <- try(dr_enrich_from_targets(tg_df, db = input$selectedDB, q = 0.05), silent = TRUE)
      return(if (inherits(ans,"try-error") || is.null(ans)) data.frame() else ans)
    }

    if (sp %in% c("dm","dme","fly")) {
      showNotification("Fetching fly targets (orthology)…", type="message", duration=3)
      tg_df <- try(dm_targets_df(mirnas), silent = TRUE)
      if (inherits(tg_df,"try-error") || !is.data.frame(tg_df) || !nrow(tg_df)) {
        showNotification("No fly targets returned.", type="error"); return(data.frame())
      }
      showNotification(paste("Targets found:", nrow(tg_df)), type="message", duration=3)
      ans <- try(dm_enrich_from_targets(tg_df, db = input$selectedDB, q = 0.05), silent = TRUE)
      return(if (inherits(ans,"try-error") || is.null(ans)) data.frame() else ans)
    }

    # Default (human/other): use cached_get_targets → enrich_species
    tg <- try(cached_get_targets(mirnas, species = sp), silent = TRUE)
    if (inherits(tg,"try-error") || is.null(tg) || !is.data.frame(tg) || !"target" %in% names(tg))
      return(data.frame())
    genes <- unique(stats::na.omit(tg$target))
    if (length(genes) < 2) return(data.frame())

    if (sp %in% c("hs","hsa","human") && input$selectedDB %in% enrichr_human_choices) {
      ee <- .enrichr_pull(genes, input$selectedDB)
      if (nrow(ee)) return(ee)
      db_gp <- switch(
        input$selectedDB,
        "GO_Biological_Process_2021" = "GO:BP",
        "GO_Molecular_Function_2021" = "GO:MF",
        "GO_Cellular_Component_2021" = "GO:CC",
        "KEGG_2021_Human"            = "KEGG",
        "Reactome_2022"              = "REAC",
        "WikiPathway_2021_Human"     = "WP",
        "GO:BP"
      )
      return(enrich_species(genes, db_gp, "hs"))
    } else {
      return(enrich_species(genes, input$selectedDB, sp))
    }
  }

  observeEvent(input$enrichAllBtn,  { showNotification("Enrichment: all DE",  type="message"); enrichAll(do_enrich("all"))  })
  observeEvent(input$enrichUpBtn,   { showNotification("Enrichment: up",      type="message"); enrichUp(do_enrich("up"))   })
  observeEvent(input$enrichDownBtn, { showNotification("Enrichment: down",    type="message"); enrichDown(do_enrich("down")) })

  output$enrichAllTable  <- renderTable({ renderEnrichTable(enrichAll())  })
  output$enrichUpTable   <- renderTable({ renderEnrichTable(enrichUp())   })
  output$enrichDownTable <- renderTable({ renderEnrichTable(enrichDown()) })

  output$downloadEnrichAll  <- downloadHandler(
    filename = function() { "enrichment_all.csv"  },
    content  = function(f) { write.csv(enrichAll(),  f, row.names=FALSE) }
  )
  output$downloadEnrichUp   <- downloadHandler(
    filename = function() { "enrichment_up.csv"   },
    content  = function(f) { write.csv(enrichUp(),   f, row.names=FALSE) }
  )
  output$downloadEnrichDown <- downloadHandler(
    filename = function() { "enrichment_down.csv" },
    content  = function(f) { write.csv(enrichDown(), f, row.names=FALSE) }
  )

  observeEvent(input$plotEnrichAllBtn,  { req(enrichAll());  output$enrichAllPlot  <- renderPlotly(renderEnrichPlot(enrichAll()))  })
  observeEvent(input$plotEnrichUpBtn,   { req(enrichUp());   output$enrichUpPlot   <- renderPlotly(renderEnrichPlot(enrichUp()))   })
  observeEvent(input$plotEnrichDownBtn, { req(enrichDown()); output$enrichDownPlot <- renderPlotly(renderEnrichPlot(enrichDown())) })

  # ====== Plots ======
  output$topTable <- renderTable({
    req(resultsDF()); df <- resultsDF()
    df <- df[order(df$padj), ]; df <- head(df, 20)
    df$padj <- signif(df$padj, 3); df$log2FoldChange <- signif(df$log2FoldChange, 3)
    df
  })

  observeEvent(input$volcanoBtn, {
    req(resultsDF())
    df <- resultsDF() %>%
      mutate(log10padj = -log10(padj),
             sig_label = case_when(
               padj < 0.05 & log2FoldChange >  1 ~ "Up",
               padj < 0.05 & log2FoldChange < -1 ~ "Down",
               TRUE ~ "Not Sig"
             ))
    output$volcanoPlot <- renderPlotly({
      plot_ly(df, x=~log2FoldChange, y=~log10padj, color=~sig_label, type="scatter", mode="markers",
              text=~paste(miRNA, "<br>log2FC:", signif(log2FoldChange,3), "<br>padj:", signif(padj,3))) %>%
        layout(xaxis=list(title="log2FC"), yaxis=list(title="-log10(padj)"))
    })
  })

  observeEvent(input$barplotBtn, {
    req(resultsDF())
    df <- resultsDF() %>%
      arrange(padj) %>%
      filter(!is.na(padj)) %>%
      slice_head(n=20)
    df$miRNA <- factor(df$miRNA, levels = rev(df$miRNA))
    output$barplotPlot <- renderPlotly({
      plot_ly(df, x=~log2FoldChange, y=~miRNA, type="bar", orientation="h",
              text=~paste("padj:", signif(padj,3)), hoverinfo="text") %>%
        layout(xaxis=list(title="log2FC"), yaxis=list(title=""), margin=list(l=140))
    })
  })

observeEvent(input$pcaBtn, {
  req(vsdData())

  # Grab the DESeqTransform object once
  vsd_obj <- vsdData()

  # samples x features matrix
  vst <- t(assay(vsd_obj))

  # keep features with non-trivial variance
  keep <- apply(vst, 2, function(x) sd(x, na.rm = TRUE) > 1e-8)
  vst  <- vst[, keep, drop = FALSE]

  # require at least 2 variable features (NO validate/need here)
  if (ncol(vst) < 2) {
    showNotification("Too few variable features for PCA", type = "error")
    return(NULL)
  }

  # (optional but recommended) use top 500 most variable features
  vars <- apply(vst, 2, var, na.rm = TRUE)
  sel  <- order(vars, decreasing = TRUE)[seq_len(min(500, length(vars)))]
  vst  <- vst[, sel, drop = FALSE]

  # run PCA
  pca <- prcomp(vst, scale. = TRUE)

  pc <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  pc$Sample <- rownames(pc)

  # get group/condition from colData
  pc$Group <- as.factor(colData(vsd_obj)$condition)

  output$pcaPlot <- renderPlotly({
    plot_ly(
      pc,
      x     = ~PC1,
      y     = ~PC2,
      color = ~Group,
      text  = ~Sample,
      type  = "scatter",
      mode  = "markers"
    )
  })
})


  observeEvent(input$umapBtn, {
    req(vsdData())
    vst <- t(assay(vsdData()))
    um <- umap::umap(vst); ud <- as.data.frame(um$layout)
    names(ud) <- c("UMAP1","UMAP2"); ud$Sample <- rownames(vst); ud$Group <- vsdData()$condition
    output$umapPlot <- renderPlotly({
      plot_ly(ud, x=~UMAP1, y=~UMAP2, color=~Group, text=~Sample,
              type="scatter", mode="markers")
    })
  })

  observeEvent(input$heatmapBtn, {
    req(resultsDF(), vsdData())
    mat <- assay(vsdData())
    sig <- intersect(resultsDF()$miRNA, rownames(mat))
    if (length(sig) < 2) { showNotification("Not enough DE miRNAs for heatmap", type="error"); return(NULL) }
    mat <- t(scale(t(mat[sig,,drop=FALSE]))); mat[is.na(mat)] <- 0
    ro <- hclust(dist(mat))$order; co <- hclust(dist(t(mat)))$order; mat <- mat[ro,co,drop=FALSE]
    output$heatmapPlot <- renderPlotly({
      plot_ly(z=mat, x=colnames(mat), y=rownames(mat), type="heatmap", colorscale="Viridis",
              hovertemplate="miRNA: %{y}<br>Sample: %{x}<br>Z: %{z:.2f}<extra></extra>")
    })
  })

  # ====== Random Forest ======
  observeEvent(input$runRF, {
    req(resultsDF(), vsdData())
    vst <- assay(vsdData()); sig <- resultsDF()$miRNA
    vst_ren <- vst; rownames(vst_ren) <- make.names(rownames(vst))
    sub <- vst_ren[rownames(vst_ren) %in% make.names(sig), , drop=FALSE]
    df_rf <- as.data.frame(t(sub)); df_rf$condition <- factor(vsdData()$condition)
    set.seed(123); tr <- sample(seq_len(nrow(df_rf)), size = 0.7*nrow(df_rf))
    trd <- df_rf[tr,]; ted <- df_rf[-tr,]
    rf <- randomForest(condition ~ ., data = trd, importance = TRUE)
    preds <- predict(rf, ted, type="response")
    probs <- try(predict(rf, ted, type="prob"), silent=TRUE)
    rf_preds(data.frame(Sample = rownames(ted), Predicted = preds))
    cm <- table(preds, ted$condition); acc <- sum(diag(cm))/sum(cm)
    auc <- NA
    if (!inherits(probs,"try-error") && ncol(probs) == 2) {
      pos <- colnames(probs)[2]; auc <- as.numeric(pROC::auc(ted$condition, probs[,pos]))
    }
    rf_metrics(data.frame(Metric=c("Accuracy","AUC"), Value=c(acc, auc)))
    imp <- importance(rf, type=2)
    rf_importance(data.frame(
      miRNA = rownames(imp),
      Importance = imp[,1]
    ) %>% arrange(desc(Importance)))
    showNotification("✅ Random Forest complete", type="message")
  })

  output$rfPredictions <- renderTable({ req(rf_preds()); rf_preds() })
  output$rfMetrics     <- renderTable({ req(rf_metrics()); rf_metrics() })
  output$rfImportance  <- renderTable({ req(rf_importance()); rf_importance() })

  output$downloadTop   <- downloadHandler(
    filename=function(){paste0("top_de_miRNAs_",Sys.Date(),".csv")},
    content=function(f){ write.csv(resultsDF(), f, row.names=FALSE) }
  )
  output$downloadRFpred <- downloadHandler(
    filename=function(){paste0("rf_predictions_",Sys.Date(),".csv")},
    content=function(f){ write.csv(rf_preds(), f, row.names=FALSE) }
  )
  output$downloadRFmetrics <- downloadHandler(
    filename=function(){paste0("rf_metrics_",Sys.Date(),".csv")},
    content=function(f){ write.csv(rf_metrics(), f, row.names=FALSE) }
  )
  output$downloadRFimportance <- downloadHandler(
    filename=function(){paste0("rf_importance_",Sys.Date(),".csv")},
    content=function(f){ write.csv(rf_importance(), f, row.names=FALSE) }
  )

  # ====== Power ======
  observeEvent(input$runPower, {
    req(vsdData())
    vst <- assay(vsdData()); cond <- factor(vsdData()$condition)
    g1 <- which(cond == levels(cond)[1]); g2 <- which(cond == levels(cond)[2])
    n_range <- seq(5, 100, 5); genes <- rownames(vst)
    pw <- matrix(NA, nrow=length(genes), ncol=length(n_range),
                 dimnames=list(genes, paste0("n=",n_range)))
    for (i in seq_along(genes)) {
      vals <- vst[genes[i], ]; sdv <- sd(vals); if (sdv == 0) next
      d <- abs(mean(vals[g1]) - mean(vals[g2]))/sdv
      pw[i,] <- sapply(n_range, function(n) pwr.t.test(n=n, d=d, sig.level=0.05, type="two.sample")$power)
    }
    power_matrix(as.data.frame(pw))
  })

  output$power_table <- renderDT({
    req(power_matrix())
    datatable(power_matrix(), options=list(pageLength=10, scrollX=TRUE))
  })

  output$power_plot  <- renderPlot({
    req(power_matrix())
    df <- power_matrix()
    ss <- as.numeric(gsub("n=","",colnames(df)))
    mins <- apply(as.matrix(df),1,function(p){
      i <- which(p >= 0.8)[1]
      if(is.na(i)) NA else ss[i]
    })
    boxplot(na.omit(mins), main="Sample size for 80% power",
            ylab="Per group", col="skyblue")
    abline(h = median(na.omit(mins)), col="red", lty=2)
  })

  output$readmeText <- renderPrint({
    if (file.exists("readme.txt"))
      cat(readLines("readme.txt"), sep="\n")
    else
      cat("README not found.")
  })
}

shinyApp(ui, server)





