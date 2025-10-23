# ====================== app.R ======================
suppressPackageStartupMessages({
  library(shiny); library(DT); library(plotly)
  library(readr); library(dplyr); library(ggplot2)
  library(DESeq2); library(umap); library(pwr); library(randomForest); library(pROC)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichR)
  library(multiMiR)       # human/mouse validated targets
  library(miRNAtap)       # predicted targets (all 4 species)
  library(gprofiler2)     # multi-species enrichment
})

# ---- Static CSV fallback (kept from your app) ----
mirna_static_targets <- tryCatch(
  read.csv("data/mirna_targets.csv", stringsAsFactors = FALSE),
  error = function(e) data.frame(mirna=character(), target=character())
)

get_targets_from_csv <- function(mirnas, static_df = mirna_static_targets) {
  if (nrow(static_df) == 0) return(static_df)
  static_df %>% filter(tolower(mirna) %in% tolower(mirnas))
}

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
# multiMiR org codes (validated only for hs/mm)
MM_ORG <- c(hs = "hsa", mm = "mmu")
# miRNAtap taxon IDs (predicted)
TAP_TAX <- c(hs = 9606, mm = 10090, dr = 7955, dm = 7227)

# ---- Target retrieval (validated + predicted + CSV fallback) ----
get_targets <- function(mirnas, species) {
  mirnas <- unique(na.omit(mirnas))
  if (length(mirnas) < 1) return(data.frame(mirna=character(), target=character()))
  
  out <- list()
  
  # 1) Validated (multiMiR) for human/mouse
  if (species %in% c("hs","mm") && MM_ORG[[species]] %in% c("hsa","mmu")) {
    mm <- try(multiMiR::get_multimir(
      org   = MM_ORG[[species]],
      mirna = mirnas,
      table = "validated",
      summary = FALSE
    ), silent = TRUE)
    if (!inherits(mm, "try-error") && !is.null(mm) && nrow(mm@data) > 0) {
      df <- as.data.frame(mm@data)
      df <- df[, c("mature_mirna_id","target_symbol","database"), drop=FALSE]
      names(df) <- c("mirna","target","source")
      df$evidence <- "validated"
      out$validated <- df[, c("mirna","target","evidence","source")]
    }
  }
  
  # 2) Predicted (miRNAtap) for all 4 species
  if (!is.null(TAP_TAX[[species]])) {
    tx <- TAP_TAX[[species]]
    for (mi in mirnas) {
      pr <- try(miRNAtap::getPredictedTargets(
        miRNA   = mi,
        species = tx,
        method  = c("targetscan","pictar","diana","mirdb","mirmap"),
        min_src = 2
      ), silent = TRUE)
      if (!inherits(pr, "try-error") && is.data.frame(pr) && nrow(pr) > 0) {
        out[[paste0("pred_", mi)]] <- data.frame(
          mirna = mi,
          target = rownames(pr),
          evidence = "predicted",
          source = "miRNAtap",
          row.names = NULL
        )
      }
    }
  }
  
  # 3) CSV fallback (always available)
  csv_df <- get_targets_from_csv(mirnas)
  if (nrow(csv_df) > 0) {
    csv_df$evidence <- "csv"
    csv_df$source <- "csv"
    out$csv <- csv_df[, c("mirna","target","evidence","source")]
  }
  
  if (length(out) == 0) return(data.frame(mirna=character(), target=character()))
  ans <- do.call(rbind, out)
  ans <- ans[!duplicated(ans[, c("mirna","target","evidence")]), ]
  rownames(ans) <- NULL
  ans
}

# ---- Enrichment DB menus by species ----
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

# ---- Enrichment runners ----
enrich_species <- function(genes, db, species) {
  genes <- unique(na.omit(genes)); genes <- genes[genes != ""]
  if (length(genes) < 2) return(data.frame())
  
  # Human → Enrichr
  if (species == "hs" && db %in% enrichr_human_choices) {
    ee <- try(enrichR::enrichr(genes, databases = db), silent = TRUE)
    if (!inherits(ee, "try-error") && !is.null(ee[[db]]) && nrow(ee[[db]]) > 0) {
      return(ee[[db]])
    }
    # if Enrichr failed, fall through to g:Profiler human sources
    species <- "hs"  # keep org
    db <- switch(db,
                 "GO_Biological_Process_2021" = "GO:BP",
                 "GO_Molecular_Function_2021" = "GO:MF",
                 "GO_Cellular_Component_2021" = "GO:CC",
                 "KEGG_2021_Human"            = "KEGG",
                 "Reactome_2022"              = "REAC",
                 "WikiPathway_2021_Human"     = "WP",
                 "GO:BP"
    )
  }
  
  # All species → g:Profiler
  org <- switch(species,
                hs = "hsapiens",
                mm = "mmusculus",
                dr = "drerio",
                dm = "dmelanogaster",
                "hsapiens"
  )
  src <- db
  gp <- try(gprofiler2::gost(query = genes, organism = org,
                             sources = src, correction_method = "g_SCS"),
            silent = TRUE)
  if (!inherits(gp, "try-error") && !is.null(gp$result) && nrow(gp$result) > 0) {
    df <- gp$result
    # harmonize shape to look like Enrichr
    comb <- -log10(df$p_value) * (df$intersection_size / df$term_size)
    out <- data.frame(
      Term = df$term_name,
      Adjusted.P.value = df$p_value,
      Combined.Score = comb,
      stringsAsFactors = FALSE
    )
    out[order(out$Adjusted.P.value), ]
  } else {
    data.frame()
  }
}

# ---- Plot/table helpers (kept) ----
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
  if (is.null(df) || nrow(df) == 0) return(data.frame(Message="No enrichment results"))
  d <- as.data.frame(df)
  # try to standardize names
  nm <- tolower(gsub("\\.","_",names(d)))
  ap <- which(nm %in% c("adjusted_p_value","adjusted.p.value","p_value","padj"))[1]
  tm <- which(nm %in% c("term","description","term_name"))[1]
  if (!is.na(ap)) names(d)[ap] <- "Adjusted_P_value"
  if (!is.na(tm)) names(d)[tm] <- "Term"
  keep <- c("Term","Adjusted_P_value")
  keep <- keep[keep %in% names(d)]
  d <- d[, keep, drop=FALSE]
  if ("Adjusted_P_value" %in% names(d)) d$Adjusted_P_value <- signif(d$Adjusted_P_value,4)
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
  resultsDF <- reactiveVal(NULL)
  ddsData   <- reactiveVal(NULL)
  vsdData   <- reactiveVal(NULL)
  power_matrix <- reactiveVal(NULL)
  
  enrichAll <- reactiveVal(NULL); enrichUp <- reactiveVal(NULL); enrichDown <- reactiveVal(NULL)
  rf_preds <- reactiveVal(NULL); rf_metrics <- reactiveVal(NULL); rf_importance <- reactiveVal(NULL)
  
  
  
  # ---- Caching for multiMiR / miRNAtap results ----
  target_cache <- reactiveValues()
  
  cached_get_targets <- function(mirnas, species) {
    # Build reproducible key
    key <- paste0(species, "_", digest::digest(sort(mirnas)))
    # Use cache if exists
    if (!is.null(target_cache[[key]])) {
      message("✅ Using cached targets for ", species, " (", length(mirnas), " miRNAs)")
      return(target_cache[[key]])
    }
    # Otherwise run and store
    tg <- get_targets(mirnas, species)
    target_cache[[key]] <- tg
    tg
  }
  
  
  # ----- DB menu reacts to species -----
  output$db_picker <- renderUI({
    sp <- input$species
    choices <- db_choices_for_species(if (sp %in% c("hs","mm","dr","dm")) sp else "hs")
    selectInput("selectedDB", "Pathway Database", choices = choices, selected = choices[[1]])
  })
  
  # ====== Run DE + RF ranking ======
  observeEvent(input$runAnalysis, {
    req(input$countsFile, input$metaFile)
    showNotification("🧬 DESeq2 running…", type = "message")
    
    count_data <- read_csv(input$countsFile$datapath, show_col_types = FALSE)
    count_matrix <- as.matrix(count_data[,-1]); rownames(count_matrix) <- count_data[[1]]
    
    meta_data <- read_csv(input$metaFile$datapath, show_col_types = FALSE)
    meta_data <- as.data.frame(meta_data)
    rownames(meta_data) <- meta_data[[1]]; meta_data <- meta_data[,-1, drop=FALSE]
    meta_data$condition <- as.factor(meta_data[[ncol(meta_data)]])
    meta_data <- meta_data[colnames(count_matrix), , drop = FALSE]
    
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = meta_data, design = ~ condition)
    dds <- DESeq(dds); ddsData(dds)
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE); vsdData(vsd)
    
    res <- results(dds); res$miRNA <- rownames(res)
    res_df <- as.data.frame(res) %>% mutate(sig = !is.na(padj) & padj < 0.1)
    # RF ranking on top 20 by padj
    top20 <- res_df %>% arrange(padj) %>% filter(!is.na(padj)) %>% slice_head(n=20)
    vst <- assay(vsd)
    vst_ren <- vst; rownames(vst_ren) <- make.names(rownames(vst))
    top20$miRNA_clean <- make.names(top20$miRNA)
    sub <- vst_ren[rownames(vst_ren) %in% top20$miRNA_clean, , drop=FALSE]
    df_rf <- as.data.frame(t(sub)); df_rf$condition <- factor(colData(vsd)$condition)
    rf <- randomForest(condition ~ ., data = df_rf, importance = TRUE)
    imp <- importance(rf, type = 2); imp_df <- data.frame(miRNA_clean = rownames(imp), Importance = imp[,1]) %>%
      arrange(desc(Importance)) %>% slice_head(n=11)
    map <- data.frame(miRNA = top20$miRNA, miRNA_clean = top20$miRNA_clean)
    imp_df <- merge(imp_df, map, by="miRNA_clean", all.x = TRUE)
    final_hits <- res_df %>% filter(miRNA %in% imp_df$miRNA)
    if (!nrow(final_hits)) { showNotification("No RF-ranked hits.", type="error"); return(NULL) }
    resultsDF(final_hits)
    
    showNotification("✅ DE + RF complete", type="message")
    
    # If species is auto, infer now and refresh DB menu
    if (input$species == "auto") {
      sp <- detect_species(final_hits$miRNA, fallback = "hs")
      updateSelectInput(session, "species", selected = sp)
      updateSelectInput(session, "selectedDB", choices = db_choices_for_species(sp))
    }
  })
  
  # ====== Enrichment wrappers (All/Up/Down) ======
  get_species_current <- reactive({
    sp <- input$species
    if (sp == "auto" && !is.null(resultsDF())) detect_species(resultsDF()$miRNA, "hs") else sp
  })
  
  do_enrich <- function(which_set=c("all","up","down")) {
    req(resultsDF())
    df <- resultsDF()
    which_set <- match.arg(which_set)
    if (which_set == "all") {
      mirnas <- df$miRNA
    } else if (which_set == "up") {
      mirnas <- df %>% filter(!is.na(padj), padj < 0.1, log2FoldChange > 0) %>% pull(miRNA)
    } else {
      mirnas <- df %>% filter(!is.na(padj), padj < 0.1, log2FoldChange < 0) %>% pull(miRNA)
    }
    if (length(mirnas) < 1) return(data.frame())
    
    sp <- get_species_current()
    tg <- cached_get_targets(mirnas, species = sp)   # <---- cached version here!
    genes <- unique(tg$target)
    if (length(genes) < 2) return(data.frame())
    
    enrich_species(genes, input$selectedDB, sp)
  }
  
  observeEvent(input$enrichAllBtn,  { showNotification("Enrichment: all DE",  type="message"); enrichAll(do_enrich("all")) })
  observeEvent(input$enrichUpBtn,   { showNotification("Enrichment: up",      type="message"); enrichUp(do_enrich("up")) })
  observeEvent(input$enrichDownBtn, { showNotification("Enrichment: down",    type="message"); enrichDown(do_enrich("down")) })
  
  output$enrichAllTable  <- renderTable({ renderEnrichTable(enrichAll()) })
  output$enrichUpTable   <- renderTable({ renderEnrichTable(enrichUp()) })
  output$enrichDownTable <- renderTable({ renderEnrichTable(enrichDown()) })
  
  output$downloadEnrichAll  <- downloadHandler(filename=function(){ "enrichment_all.csv"  }, content=function(f){ write.csv(enrichAll(), f, row.names=FALSE) })
  output$downloadEnrichUp   <- downloadHandler(filename=function(){ "enrichment_up.csv"   }, content=function(f){ write.csv(enrichUp(), f, row.names=FALSE) })
  output$downloadEnrichDown <- downloadHandler(filename=function(){ "enrichment_down.csv" }, content=function(f){ write.csv(enrichDown(), f, row.names=FALSE) })
  
  observeEvent(input$plotEnrichAllBtn,  { req(enrichAll());  output$enrichAllPlot  <- renderPlotly(renderEnrichPlot(enrichAll())) })
  observeEvent(input$plotEnrichUpBtn,   { req(enrichUp());   output$enrichUpPlot   <- renderPlotly(renderEnrichPlot(enrichUp())) })
  observeEvent(input$plotEnrichDownBtn, { req(enrichDown()); output$enrichDownPlot <- renderPlotly(renderEnrichPlot(enrichDown())) })
  
  # ====== Heatmap / PCA / UMAP / Volcano / Barplot (unchanged) ======
  output$topTable <- renderTable({
    req(resultsDF()); df <- resultsDF()
    df <- df[order(df$padj), ]; df <- head(df, 20)
    df$padj <- signif(df$padj, 3); df$log2FoldChange <- signif(df$log2FoldChange, 3)
    df
  })
  observeEvent(input$volcanoBtn, {
    req(resultsDF()); df <- resultsDF() %>%
      mutate(log10padj = -log10(padj),
             sig_label = case_when(padj < 0.05 & log2FoldChange >  1 ~ "Up",
                                   padj < 0.05 & log2FoldChange < -1 ~ "Down",
                                   TRUE ~ "Not Sig"))
    output$volcanoPlot <- renderPlotly({
      plot_ly(df, x=~log2FoldChange, y=~log10padj, color=~sig_label, type="scatter", mode="markers",
              text=~paste(miRNA, "<br>log2FC:", signif(log2FoldChange,3), "<br>padj:", signif(padj,3))) %>%
        layout(xaxis=list(title="log2FC"), yaxis=list(title="-log10(padj)"))
    })
  })
  observeEvent(input$barplotBtn, {
    req(resultsDF()); df <- resultsDF() %>% arrange(padj) %>% filter(!is.na(padj)) %>% slice_head(n=20)
    df$miRNA <- factor(df$miRNA, levels = rev(df$miRNA))
    output$barplotPlot <- renderPlotly({
      plot_ly(df, x=~log2FoldChange, y=~miRNA, type="bar", orientation="h",
              text=~paste("padj:", signif(padj,3)), hoverinfo="text") %>%
        layout(xaxis=list(title="log2FC"), yaxis=list(title=""), margin=list(l=140))
    })
  })
  observeEvent(input$pcaBtn, {
    req(vsdData()); vst <- t(assay(vsdData()))
    vst <- vst[, apply(vst,2,function(x) sd(x,na.rm=TRUE) > 1e-8)]
    validate(need(ncol(vst) >= 2, "Too few variable features for PCA"))
    pca <- prcomp(vst, scale.=TRUE); pc <- as.data.frame(pca$x[,1:2])
    pc$Sample <- rownames(pc); pc$Group <- vsdData()$condition
    output$pcaPlot <- renderPlotly({ plot_ly(pc, x=~PC1, y=~PC2, color=~Group, text=~Sample,
                                             type="scatter", mode="markers") })
  })
  observeEvent(input$umapBtn, {
    req(vsdData()); vst <- t(assay(vsdData()))
    um <- umap::umap(vst); ud <- as.data.frame(um$layout)
    names(ud) <- c("UMAP1","UMAP2"); ud$Sample <- rownames(vst); ud$Group <- vsdData()$condition
    output$umapPlot <- renderPlotly({ plot_ly(ud, x=~UMAP1, y=~UMAP2, color=~Group, text=~Sample,
                                              type="scatter", mode="markers") })
  })
  observeEvent(input$heatmapBtn, {
    req(resultsDF(), vsdData()); mat <- assay(vsdData())
    sig <- intersect(resultsDF()$miRNA, rownames(mat))
    if (length(sig) < 2) { showNotification("Not enough DE miRNAs for heatmap", type="error"); return(NULL) }
    mat <- t(scale(t(mat[sig,,drop=FALSE]))); mat[is.na(mat)] <- 0
    ro <- hclust(dist(mat))$order; co <- hclust(dist(t(mat)))$order; mat <- mat[ro,co,drop=FALSE]
    output$heatmapPlot <- renderPlotly({
      plot_ly(z=mat, x=colnames(mat), y=rownames(mat), type="heatmap", colorscale="Viridis",
              hovertemplate="miRNA: %{y}<br>Sample: %{x}<br>Z: %{z:.2f}<extra></extra>")
    })
  })
  
  # ====== Random Forest (unchanged, simplified metrics) ======
  observeEvent(input$runRF, {
    req(resultsDF(), vsdData())
    vst <- assay(vsdData()); sig <- resultsDF()$miRNA
    vst_ren <- vst; rownames(vst_ren) <- make.names(rownames(vst))
    sub <- vst_ren[rownames(vst_ren) %in% make.names(sig), , drop=FALSE]
    df_rf <- as.data.frame(t(sub)); df_rf$condition <- factor(vsdData()$condition)
    set.seed(123); tr <- sample(seq_len(nrow(df_rf)), size = 0.7*nrow(df_rf))
    trd <- df_rf[tr,]; ted <- df_rf[-tr,]
    rf <- randomForest(condition ~ ., data = trd, importance = TRUE)
    preds <- predict(rf, ted, type="response"); probs <- try(predict(rf, ted, type="prob"), silent=TRUE)
    rf_preds(data.frame(Sample = rownames(ted), Predicted = preds))
    cm <- table(preds, ted$condition); acc <- sum(diag(cm))/sum(cm)
    auc <- NA
    if (!inherits(probs,"try-error") && ncol(probs) == 2) {
      pos <- colnames(probs)[2]; auc <- as.numeric(pROC::auc(ted$condition, probs[,pos]))
    }
    rf_metrics(data.frame(Metric=c("Accuracy","AUC"), Value=c(acc, auc)))
    imp <- importance(rf, type=2)
    rf_importance(data.frame(miRNA = rownames(imp), Importance = imp[,1]) %>% arrange(desc(Importance)))
    showNotification("✅ Random Forest complete", type="message")
  })
  output$rfPredictions <- renderTable({ req(rf_preds()); rf_preds() })
  output$rfMetrics     <- renderTable({ req(rf_metrics()); rf_metrics() })
  output$rfImportance  <- renderTable({ req(rf_importance()); rf_importance() })
  output$downloadTop   <- downloadHandler(filename=function(){paste0("top_de_miRNAs_",Sys.Date(),".csv")},
                                          content=function(f){ write.csv(resultsDF(), f, row.names=FALSE) })
  output$downloadRFpred <- downloadHandler(filename=function(){paste0("rf_predictions_",Sys.Date(),".csv")},
                                           content=function(f){ write.csv(rf_preds(), f, row.names=FALSE) })
  output$downloadRFmetrics <- downloadHandler(filename=function(){paste0("rf_metrics_",Sys.Date(),".csv")},
                                              content=function(f){ write.csv(rf_metrics(), f, row.names=FALSE) })
  output$downloadRFimportance <- downloadHandler(filename=function(){paste0("rf_importance_",Sys.Date(),".csv")},
                                                 content=function(f){ write.csv(rf_importance(), f, row.names=FALSE) })
  
  # ====== Power analysis (kept) ======
  observeEvent(input$runPower, {
    req(vsdData()); vst <- assay(vsdData()); cond <- factor(vsdData()$condition)
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
  output$power_table <- renderDT({ req(power_matrix()); datatable(power_matrix(), options=list(pageLength=10, scrollX=TRUE)) })
  output$power_plot  <- renderPlot({
    req(power_matrix()); df <- power_matrix(); ss <- as.numeric(gsub("n=","",colnames(df)))
    mins <- apply(as.matrix(df),1,function(p){ i <- which(p >= 0.8)[1]; if(is.na(i)) NA else ss[i] })
    boxplot(na.omit(mins), main="Sample size for 80% power", ylab="Per group", col="skyblue")
    abline(h = median(na.omit(mins)), col="red", lty=2)
  })
  
  # README
  output$readmeText <- renderPrint({ if (file.exists("readme.txt")) cat(readLines("readme.txt"), sep="\n") else cat("README not found.") })
}

shinyApp(ui, server)


